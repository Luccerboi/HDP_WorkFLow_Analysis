from collections import Counter
from itertools import combinations
from atomate2.forcefields.jobs import ForceFieldRelaxMaker
from atomate2.forcefields import MLFF
from jobflow import Flow, run_locally # type: ignore 
from pymatgen.core import Structure
import os
import pandas as pd

def get_chemical_subsystems(chemsys: str)->list:
    """Code adapted from mp_api.client.mprester.py."""
    subsystems = []
    elements = chemsys.split("-")
    elements_set = set(elements)
    for i in range(len(elements_set)):
        for els in combinations(elements_set, i + 1):
            subsystems.append("-".join(sorted(els)))
    return subsystems

def mlip_relax_and_get_energies(force_field_name: str | MLFF,
                                structure_dict: pd.Series[Structure],
                                calculator_kwargs: dict ={"default_dtype": "float64"}) -> dict:  # calc kwargs need to be adapted for nequip
    ff_maker = ForceFieldRelaxMaker(force_field_name=force_field_name,
                                    calculator_kwargs=calculator_kwargs)

    job_uuid_to_idfr = {}
    jobs = []
    for mpid, structure in structure_dict.items():
        job = ff_maker.make(structure)
        jobs.append(job)
        job_uuid_to_idfr[job.uuid] = mpid

    flow = Flow(jobs)
    resp = run_locally(flow)

    flow_output = {}
    for uuid, mpid in job_uuid_to_idfr.items():
        output = resp[uuid][1].output
        if (
                output is not None
                and output.is_force_converged
                # and abs(output.output.energy) < energy_tol
        ):
            # Beware: differently structured dict than in direct phase diagram computation util below
            flow_output[mpid] = output.output.energy
    return flow_output

