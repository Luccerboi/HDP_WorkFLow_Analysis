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
    """Code adapted from mp_api.client.mprester.py."""
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

def mlip_relax_batched(
    force_field_name: str | MLFF,
    structure_dict: pd.Series,
    calculator_kwargs: dict = {},
    relax_kwargs: dict = {},
    optimizer_kwargs: dict = {},
) -> dict:  # calc kwargs need to be adapted for nequip
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
                pd_entries.remove(my_entry)

                phaseD = PhaseDiagram(pd_entries,elements=elemental_endpoints)

                s_ehull[mlip] = phaseD.get_e_above_hull(my_entry, allow_negative=True)
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
    print('MLIP:', '\t','NaN vals:')
    for mlip in [x for x in mlip_df.columns if x.startswith('E_')]:
        all_ens = mlip_df[mlip].dropna()
        dropped_comps = len(mlip_df) - len(all_ens)
        print(mlip, '\t', dropped_comps)

    print('MLIP:','\t','Missing HDPs:', '\t', 'share unstable (>50meV/atom):')

    for mlip in ehull_df.columns:
        sclean = ehull_df[mlip].dropna()
        dropped_cols = len(ehull_df) - len(sclean)
        non_stalbe = sclean[sclean>0.05]
        print(mlip, '\t', dropped_cols, '\t', len(non_stalbe)/len(sclean))

    print('MLIP:','\t','Missing HDPs:', '\t', 'share unstable (>100meV/atom):')

    for mlip in ehull_df.columns:
        sclean = ehull_df[mlip].dropna()
        dropped_cols = len(ehull_df) - len(sclean)
        non_stalbe = sclean[sclean>0.1]
        print(mlip, '\t', dropped_cols, '\t', len(non_stalbe)/len(sclean))

    print('MLIP:','\t','Missing HDPs:', '\t', 'share unstable (>200meV/atom):')

    for mlip in ehull_df.columns:
        sclean = ehull_df[mlip].dropna()
        dropped_cols = len(ehull_df) - len(sclean)
        non_stalbe = sclean[sclean>0.2]
        print(mlip, '\t', dropped_cols, '\t', len(non_stalbe)/len(sclean))

    return

    



        




