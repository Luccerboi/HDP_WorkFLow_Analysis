from ehull_utils import mlip_relax_and_get_energies, mlip_relax_batched
import pandas as pd
from pymatgen.core import Structure
import os
import json
#import torch
import sys
from math import ceil
from filelock import FileLock

mlip = sys.argv[1]
batchsize = int(sys.argv[2])
batchnum = int(sys.argv[3])

lock = FileLock(f"mpid_energy_dict_{mlip}.json.lock")

to_compute = ["MACE-Gabor","SevenNet","MatterSim","MatPES-r2SCAN", "MACE-MP-0b3", "MACE-MPA-0", "Nequip"]  # Nequip -> Nequip-OAM-XL
#assert mlip in to_compute, f"MLIP argument must be in {to_compute}."

## Read chemical systems to simulate and convert structure dicts to structure objects
#subsys_df = pd.read_csv('hdp_mp_subsysandhdps.csv',index_col=0)
subsys_df = pd.read_csv(f'MissingStrucs_{mlip}.csv',index_col=0)

relax_kwargs = {
    # "steps" : 1000,
    "fmax" : 0.005,
    "maxstep": 0.001,
}
optimizer_kwargs = {

    }

mlip_specifications = {
    "SevenNet-omat24" : {
        "model_name" : "SevenNet",
        "mlip_kwargs" : {
            "model" : "7net-omni-i12",
            "modal" : "omat24",
            "enable_cueq" : True,
            "device" : "cuda",
            "default_dtype" : "float64"
        }
    },
    "SevenNet-MPALOE" : {
        "model_name" : "SevenNet",
        "mlip_kwargs" : {
            "model" : "7net-omni-i12",
            "modal" : "mp_r2scan",
            "enable_cueq" : True,
            "device" : "cuda",
            "default_dtype" : "float64"
        }
    },
    "GRACE-3L-OMAT-large-ft-AM" : {
        "model_name" : {"@module": "tensorpotential.calculator.foundation_models", "@callable": "grace_fm"}, #placeholder name to forward GRACE callable
        "mlip_kwargs" : {
            "model" : "GRACE-3L-OMAT-large-ft-AM",
            "device" : "cuda",
        #    "default_dtype" : "float64", (GRACE 3L is only available in 32 bit float
        }

    },
    "PET-OAM-XL" :{
        "model_name" : "UPET",
        "mlip_kwargs" : {
            "model" : "pet-oam-xl",
            "version" : "1.0.0",
            "device" : "cuda",
            "default_dtype" : "float64"
        }
    },
    "Nequip-OAM-XL" : {
        "model_name" : "Nequip",
        "mlip_kwargs" : {
            "compile_path" : "/home/lwalterb/hdp_ehull/COMPILED_MODELS/mir-group__NequIP-OAM-XL__0.1.nequip.pt2",
            "device" : "cuda",
            "default_dtype" : "float64",
        }
    },
    "MatterSim" : {
        "model_name" : {"@module": "mattersim.forcefield.potential","@callable": "MatterSimCalculator"},
        "mlip_kwargs": {
            "potential" : "mattersim-v1.0.0-5m",
            "device" : "cuda",
            "default_dtype" : "float64"
        },
    },
    "Equiformer_v3" : {
        "model_name" : {"@module": "fairchem.core.common.relaxation.ase_utils", "@callable": "OCPCalculator"},
        "mlip_kwargs" : dict(
            checkpoint_path= "/home/lwalterb/.cache/huggingface/hub/models--mirror-physics--equiformer_v3/snapshots/ca3ce91a7109ecce14cd198b67616a9b290ef64a/checkpoint/omat24-mptrj-salex_gradient.pt",
            local_cache="pretrained_models",
            cpu=False,
        )
    }
}



data_dir = "ExpensiveMLIPs"
os.makedirs(data_dir, exist_ok=True)

n_batches = ceil(len(subsys_df)/batchsize)
try:
    batch_df = subsys_df.iloc[batchnum * batchsize : (batchnum+1)*batchsize]
except IndexError:
    batch_df = subsys_df.iloc[batchnum*batchsize:]

batch_df['structure'] = batch_df.apply(lambda row: Structure.from_dict(eval(row['structure_dict'])).to_conventional(),axis=1)

print(f"Starting {mlip} for batchnumber {batchnum}")


mpid_energy_dict = mlip_relax_and_get_energies(structure_dict= batch_df['structure'], 
                                                    force_field_name=mlip_specifications[mlip]["model_name"],
                                                    calculator_kwargs=mlip_specifications[mlip]["mlip_kwargs"],
                                                    relax_kwargs= relax_kwargs,
                                                    optimizer_kwargs = optimizer_kwargs,
                                               )

# with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}_batch_{batchnum}.json"), "w") as f:
    # json.dump(mpid_energy_dict, f)
with lock.acquire(timeout=30):
    try:
        with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}.json"), "r") as f:
            d =  json.load(f)
    except FileNotFoundError:
        print(f"starting new data json for {mlip}")
        d = {}
        pass 
    d.update(mpid_energy_dict)
    with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}.json"), "w") as f:
        json.dump(d,f)
 
     # with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}.json"), "a") as f:
         # json.dump(mpid_energy_dict, f)

print(f"Finished {mlip} flow no. {batchnum}")

