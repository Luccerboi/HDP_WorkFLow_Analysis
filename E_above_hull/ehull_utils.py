from collections import Counter
from itertools import combinations

from atomate2.forcefields.jobs import ForceFieldRelaxMaker
from atomate2.forcefields import MLFF
from jobflow import Flow, run_locally  # type: ignore
from pymatgen.core import Structure, Element

from os import PathLike
from pathlib import Path
import os
import json
import pandas as pd
import numpy as np
import glob


def get_chemical_subsystems(chemsys: str) -> list:
    """Return every non-empty chemical subsystem represented by ``chemsys``.

    Args:
        chemsys: Hyphen-separated element symbols, such as ``"Li-Fe-O"``.

    Returns:
        Sorted, hyphen-separated subsystem names, including the elemental
        endpoints and the complete chemical system.
    """
    subsystems = []
    elements = chemsys.split("-")
    elements_set = set(elements)
    for i in range(len(elements_set)):
        for els in combinations(elements_set, i + 1):
            subsystems.append("-".join(sorted(els)))
    return subsystems


def mlip_relax_and_get_energies(
    force_field_name: str | MLFF,
    structure_dict: pd.Series,
    calculator_kwargs: dict = {},
    relax_kwargs: dict = {},
    optimizer_kwargs: dict = {},
) -> dict:  # calc kwargs need to be adapted for nequip
    """Relax structures with a force-field workflow and collect energies.

    Args:
        force_field_name: Name of the force field or an ``MLFF`` enum value.
        structure_dict: Series mapping structure identifiers to pymatgen
            structures.
        calculator_kwargs: Keyword arguments passed to the force-field
            calculator.
        relax_kwargs: Keyword arguments controlling relaxation convergence.
        optimizer_kwargs: Keyword arguments for the relaxation optimizer.

    Returns:
        Mapping from structure identifier to relaxed energy, a convergence
        failure description, or ``None`` when the structure is unsupported.
    """
    from jobflow import SETTINGS
    store = SETTINGS.JOB_STORE
    store.connect()

    ff_maker = ForceFieldRelaxMaker(
        force_field_name=force_field_name, calculator_kwargs=calculator_kwargs,fix_symmetry=True, relax_kwargs= relax_kwargs, steps= 250000,
        optimizer_kwargs = optimizer_kwargs,
    )

    job_uuid_to_idfr = {}
    jobs = []
    for mpid, structure in structure_dict.items():
        job = ff_maker.make(structure)
        jobs.append(job)
        job_uuid_to_idfr[job.uuid] = mpid

    flow = Flow(jobs)
    resp = run_locally(flow,create_folders=True, root_dir=f"{force_field_name}", store=store)

    flow_output = {}
    for uuid, mpid in job_uuid_to_idfr.items():
        # A missing response means this force field could not process the job.
        try:
            output = resp[uuid][1].output
        except (
            KeyError
        ):  # I want to ignore elements that are not covered by a specific MLIP
            flow_output[mpid] = None
            continue

        if (
            output is not None
            and output.is_force_converged
            # and abs(output.output.energy) < energy_tol
        ):
            # Beware: differently structured dict than in direct phase diagram computation util below
            flow_output[mpid] = output.output.energy

        elif (
            output is not None
            and not output.is_force_converged
        ):
            flow_output[mpid] = f"{uuid} Failed Conv.{output.output.n_steps}, rem force: {np.max(output.output.forces)}"
        else:
            flow_output[mpid] = None
    return flow_output

def mlip_relax_2step(
    force_field_name: str | MLFF,
    structure_dict: pd.Series,
    pre_relax_kwargs: dict = {"fmax":0.03,"maxstep":0.05},
    main_relax_kwargs: dict = {"fmax":0.005, "maxstep":0.01},
    calculator_kwargs: dict = {},
    optimizer1_kwargs: dict = {"optimizer": "FIRE"},
    optimizer2_kwargs: dict = {"optimizer":"LBFGS"},
) -> dict:  # calc kwargs need to be adapted for nequip
    """Run a pre-relaxation followed by a main relaxation for each structure.
    Was attempted, but it did not help. IMPLEMENTATION IS NOT COMPLETE.

    Args:
        force_field_name: Name of the force field or an ``MLFF`` enum value.
        structure_dict: Series mapping structure identifiers to pymatgen
            structures.
        pre_relax_kwargs: Convergence settings for the first relaxation.
        main_relax_kwargs: Convergence settings for the second relaxation.
        calculator_kwargs: Keyword arguments passed to both calculators.
        optimizer1_kwargs: Keyword arguments for the pre-relaxation
            optimizer.
        optimizer2_kwargs: Keyword arguments for the main-relaxation
            optimizer.

    Returns:
        An empty dictionary. The current implementation reports each final
        relaxation to stdout but does not yet persist per-structure results.
    """
    from jobflow import SETTINGS
    store = SETTINGS.JOB_STORE
    store.connect()

    pre_rel_maker = ForceFieldRelaxMaker(
        force_field_name=force_field_name, calculator_kwargs=calculator_kwargs,fix_symmetry=True, relax_kwargs= pre_relax_kwargs, steps= 250000,
        optimizer_kwargs = optimizer1_kwargs, name="pre_relax", store_trajectory=True
    )

    main_rel_maker = ForceFieldRelaxMaker(
        force_field_name=force_field_name, calculator_kwargs=calculator_kwargs,fix_symmetry=True, relax_kwargs= main_relax_kwargs, steps= 250000,
        optimizer_kwargs = optimizer2_kwargs, name="main_relax"
    )

    job_uuid_to_idfr = {}
    jobs = []
    for mpid, structure in structure_dict.items():
        # The second job consumes the relaxed structure produced by the first.
        pre_rel_job = pre_rel_maker.make(structure=structure)
        main_rel_job = main_rel_maker.make(structure= pre_rel_job.output.structure)
        rel_flow = Flow([pre_rel_job,main_rel_job], output=main_rel_job.output, name="two-stage relax")

        response = run_locally(rel_flow, create_folders=True, root_dir=f"{force_field_name}", store=store)
        final_doc = response[main_rel_job.uuid][1].output
        print("Final max force converged:", final_doc.is_force_converged)
        print("Final energy:", final_doc.output.energy)


    return {}

def mlip_relax_batched(
    force_field_name: str | MLFF,
    structure_dict: pd.Series,
    calculator_kwargs: dict = {},
    relax_kwargs: dict = {},
    optimizer_kwargs: dict = {},
) -> dict:  # calc kwargs need to be adapted for nequip
    """Relax a batch of structures in one force-field job.

    Args:
        force_field_name: Name of the force field or an ``MLFF`` enum value.
        structure_dict: Series whose index contains structure identifiers and
            whose values are pymatgen structures.
        calculator_kwargs: Keyword arguments passed to the force-field
            calculator.
        relax_kwargs: Keyword arguments controlling relaxation convergence.
        optimizer_kwargs: Keyword arguments for the relaxation optimizer.

    Returns:
        Mapping from structure identifier to relaxed energy or a convergence
        failure description. Structures whose result cannot be read map to
        ``None``.
    """
    from jobflow import SETTINGS
    store = SETTINGS.JOB_STORE
    store.connect()

    ff_maker = ForceFieldRelaxMaker(
        force_field_name=force_field_name, calculator_kwargs=calculator_kwargs,fix_symmetry=True, relax_kwargs= relax_kwargs, steps= 250000,
        optimizer_kwargs = optimizer_kwargs,
    )

    job = ff_maker.make(structure_dict.to_list())


    flow = Flow(job)
    resp = run_locally(flow,create_folders=True, root_dir=f"{force_field_name}", store=store)

    # The output is a list of ForceFieldTaskDocuments, same order as input
    try:
        docs = resp[job.uuid][1].output  # this is a list
    except KeyError:
        print('key error in responses')
        print(resp)
        print(resp.keys())

    flow_output = {}
    for mpid, doc in zip(structure_dict.index.to_list(), docs):
        # Batched output preserves the input order, so pair it with the index.
        #print(mpid, doc.output.energy, doc.output.forces)
        try:
            force_relaxed = np.max([np.linalg.norm(x) for x in doc.output.forces]) < relax_kwargs['fmax']
        except Exception as e:
            print(e)
            flow_output[mpid] = None
            continue
        
        if force_relaxed:
            flow_output[mpid] = doc.output.energy
        else:
            print(f'{mpid} failed to converge {np.max([np.linalg.norm(x) for x in doc.output.forces])}')
            flow_output[mpid] = f"{mpid} Failed Conv.{doc.output.n_steps}, rem force: {np.max(doc.output.forces)}"
        

    return flow_output

def combine_batches_to_series(
    mlip_name: str,
    data_dir: str | PathLike = "MLIP_data",
    filename_root: str = "mpid_energy_dict",
) -> pd.Series:
    """Combine JSON energy dictionaries written by separate batch runs.
    Uses glob to find the right JSON dictionaries.

    Args:
        mlip_name: Force-field name embedded in each batch filename.
        data_dir: Directory containing the batch JSON files.
        filename_root: Common prefix used by the batch files.

    Returns:
        A Series indexed by structure identifier and containing the collected
        energies.
    """
    glob_string = f"{filename_root}_{mlip_name}_batch*"
    batch_list = glob.iglob(glob_string, root_dir=Path(data_dir))

    data = {}
    for batchfile in batch_list:
        with open(os.path.join(data_dir, batchfile), "r") as f:
            data.update(json.load(f))

    return pd.Series(data)


def collect_mlip_energies_to_df(
        structures_df: str | PathLike | pd.DataFrame = "hdp_mp_subsysandhdps.csv",
        mlip_list: list[str] = ["MACE-MPA-0", "MACE-MP-0b3", "MatPES-r2SCAN", "SevenNet"],
        data_dir: str | PathLike = "MLIP_data",
        batchfn_root: str = "mpid_energy_dict",
        output_fn: str | PathLike | None = None,
):
    """Load MLIP energies into a structures DataFrame.

    Args:
        structures_df: Input DataFrame or path to a CSV file containing the
            structures and metadata.
        mlip_list: Force-field names whose JSON energy files should be loaded.
        data_dir: Directory containing the energy JSON files.
        batchfn_root: Common prefix used by the energy files.
        output_fn: Optional path at which to write the selected output columns.

    Returns:
        A copy of the input DataFrame with one ``E_<MLIP>`` column per loaded
        force field.
    """
    if type(structures_df) == pd.DataFrame:
        struc_df = structures_df.copy()
    else:
        struc_df = pd.read_csv(structures_df, index_col=0)

    for mlip_name in mlip_list:
        # s_mlip = combine_batches_to_series(
        #     mlip_name= mlip_name,
        #     data_dir= data_dir,
        #     filename_root= batchfn_root,
        # )
        with open(os.path.join(data_dir,f'{batchfn_root}_{mlip_name}.json'),'r') as f:
            data = json.load(f)
        s_mlip = pd.Series(data=data)
        struc_df[f"E_{mlip_name}"] = s_mlip

    if output_fn is not None:
        # try:
        struc_df[['nsites','chemsys'] + [x for x in struc_df.columns if x.startswith('E_')] + ['structure_dict']].to_csv(output_fn)
        # except KeyError:
        #     struc_df.to_csv(output_fn)
            


    return struc_df

def construct_phase_diagrams(
        hdp_df: pd.DataFrame,
        subsys_MLIPenergy_df: pd.DataFrame,
        dataframe_savedir: str | PathLike | None,
        phasediagram_savedir: str | PathLike | None = None,
):
    """Construct phase diagrams and calculate formation and hull energies.

    Args:
        hdp_df: DataFrame of target HDPs with a ``subsystems`` column.
        subsys_MLIPenergy_df: DataFrame containing subsystem structures,
            chemical systems, and ``E_<MLIP>`` energy columns.
        dataframe_savedir: Optional directory for the resulting CSV files.
        phasediagram_savedir: Optional directory for serialized phase
            diagrams.

    Returns:
        A tuple containing the e-above-hull and formation-energy DataFrames.
    """
    from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry, plotly_layouts
    from monty.serialization import dumpfn


    if 'structure' not in subsys_MLIPenergy_df.columns:
        try:
            subsys_MLIPenergy_df['structure'] = subsys_MLIPenergy_df.apply(lambda row: Structure.from_dict(eval(row['structure_dict'])),axis=1)
        except:
            raise IndexError(f"Could not find/parse pymatgen.Structure objects in the MLIP energy data DataFrame.\n Only columns are: {subsys_MLIPenergy_df.columns} ")

    mlip_names = [x.removeprefix('E_') for x in subsys_MLIPenergy_df.columns if x.startswith("E_")]
    ehull_collection = []
    eform_collection = []
        
    for comp_id in hdp_df.index:
        s_ehull = pd.Series(index=mlip_names, name=comp_id)
        s_eform = pd.Series(index=mlip_names, name=comp_id)
        diagram_data = {}

        subsys_list = eval(hdp_df.loc[comp_id]['subsystems'])
        elemental_endpoints = [Element(x) for x in subsys_list if "-" not in x]


        for mlip in mlip_names:
            # A phase diagram is valid only when every required endpoint and
            # subsystem has an energy for this MLIP.
            rel_subsys = subsys_MLIPenergy_df[subsys_MLIPenergy_df['chemsys'].isin(subsys_list)]
            # print(comp_id, '\t', rel_subsys.loc[comp_id][f"E_{mlip}"])

            try:
                incomplete_data = True in rel_subsys[f"E_{mlip}"].isna().values #np.isnan(rel_subsys.loc[comp_id][f"E_{mlip}"])
                # print(comp_id, ':\t',rel_subsys['chemsys'].values,'\n\t\t', rel_subsys[f"E_{mlip}"].isna().values)
                atomic_mass_list = [x.Z for x in elemental_endpoints]
                forbidden_masses = [84, 85, 86, 87, 88, 95, 96, 97]
                too_heavy_list = [x in forbidden_masses for x in atomic_mass_list]
                if True in too_heavy_list:
                    incomplete_data = True

            except TypeError:
                incomplete_data = True



            if not incomplete_data:
                pd_entries = [
                    PDEntry(rel_subsys.loc[idx]['structure'].composition, rel_subsys.loc[idx][f"E_{mlip}"], name=idx) 
                    for idx in rel_subsys.index if rel_subsys.loc[idx][f"E_{mlip}"] is not None
                    ]

                # Need to extract the HDP entry in order to get E above hull and Formation energy
                my_entry = [x for x in pd_entries if x.name == comp_id][0]
                # pd_entries.remove(my_entry)

                phaseD = PhaseDiagram(pd_entries,elements=elemental_endpoints)
                # print(my_entry)

                s_ehull[mlip] = phaseD.get_e_above_hull(my_entry)
                s_eform[mlip] = phaseD.get_form_energy_per_atom(my_entry)
                diagram_data.update({mlip: phaseD.as_dict()})
            else:
                s_ehull[mlip] = None
                s_eform[mlip] = None
                diagram_data.update({mlip:{}})

        ehull_collection.append(s_ehull)
        eform_collection.append(s_eform)

        if phasediagram_savedir is not None:
            os.makedirs(phasediagram_savedir,exist_ok=True)
            filename = Path(phasediagram_savedir) / f"{comp_id}_PhaseDiagrams_{len(mlip_names)}MLIPS.json.gz"
            dumpfn(diagram_data,filename)

    ehull_df = pd.DataFrame(ehull_collection)
    eform_df = pd.DataFrame(eform_collection)

    if dataframe_savedir is not None:
        os.makedirs(dataframe_savedir,exist_ok=True)
        ehull_df.to_csv(os.path.join(dataframe_savedir,f"Ehull_data_{'_'.join(mlip_names)}.csv"))
        eform_df.to_csv(os.path.join(dataframe_savedir,f"Eform_data_{'_'.join(mlip_names)}.csv"))

    return ehull_df, eform_df
            

def summarize_results(
        mlip_df: pd.DataFrame,
        ehull_df: pd.DataFrame,

):
    """Print missing-data counts and stability shares for each MLIP.

    Args:
        mlip_df: DataFrame containing the raw ``E_<MLIP>`` energy columns.
        ehull_df: DataFrame containing e-above-hull values by MLIP.
    """
    print('MLIP:', '\t','NaN vals:')
    for mlip in [x for x in mlip_df.columns if x.startswith('E_')]:
        all_ens = mlip_df[mlip].dropna()
        dropped_comps = len(mlip_df) - len(all_ens)
        print(mlip, '\t', dropped_comps)

    print('MLIP:','\t','Missing HDPs:', '\t', 'share stable (<=100meV/atom):')

    for mlip in ehull_df.columns:
        sclean = ehull_df[mlip].dropna()
        dropped_cols = len(ehull_df) - len(sclean)
        stalbe = sclean[sclean<=0.1]
        print(mlip, '\t', dropped_cols, '\t', len(stalbe)/len(sclean))

    print('MLIP:','\t','Missing HDPs:', '\t', 'share stable (<=150meV/atom):')

    for mlip in ehull_df.columns:
        sclean = ehull_df[mlip].dropna()
        dropped_cols = len(ehull_df) - len(sclean)
        stalbe = sclean[sclean<=0.15]
        print(mlip, '\t', dropped_cols, '\t', len(stalbe)/len(sclean))

    print('MLIP:','\t','Missing HDPs:', '\t', 'share stable (<=200meV/atom):')

    for mlip in ehull_df.columns:
        sclean = ehull_df[mlip].dropna()
        dropped_cols = len(ehull_df) - len(sclean)
        stalbe = sclean[sclean<=0.2]
        print(mlip, '\t', dropped_cols, '\t', len(stalbe)/len(sclean))

    return

    



        




