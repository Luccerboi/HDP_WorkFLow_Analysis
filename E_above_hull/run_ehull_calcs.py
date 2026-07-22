from ehull_utils import mlip_relax_and_get_energies
import pandas as pd
from pymatgen.core import Structure
import os
import json
import torch
import sys
from math import ceil


mlip = sys.argv[1]
batchsize = int(sys.argv[2])
batchnum = int(sys.argv[3])

to_compute = ["MACE-Gabor","SevenNet","MatterSim","MatPES-r2SCAN", "MACE-MP-0b3", "MACE-MPA-0", "Nequip"]  # Nequip -> Nequip-OAM-XL
#assert mlip in to_compute, f"MLIP argument must be in {to_compute}."

## Read chemical systems to simulate and convert structure dicts to structure objects
subsys_df = pd.read_csv('hdp_mp_subsysandhdps.csv',index_col=0)

relax_kwargs = {
    # "steps" : 1000,
    "fmax" : 0.05,
    "maxstep": 0.03,
}
optimizer_kwargs = {

    }

if mlip == "Nequip":
    relax_kwargs["maxstep"] = 0.001
    mlip_name = "Nequip"
    mlip_kwargs = {
            "compile_path": "/home/lwalterb/NequIP_OAM_XL/mir-group__NequIP-OAM-XL__0.1.nequip.pt2",
            "device":"cuda",
            "default_dtype":"float64",
        }
elif mlip.startswith('MatPES'):
    mlip_name = mlip
    mlip_kwargs = {
            'device':'cuda',
            'default_dtype':'float64'
                   }

    relax_kwargs['maxstep'] = 0.02 #MatPes has some issues so smaller increments
    #mlip_kwargs = {"default_dtype":"float64"}

elif mlip == "MACE-Gabor":
    mlip_name = "MACE"
    mlip_kwargs = {
        'model' : "/home/lwalterb/hdp_ehull/MACE_ALT_GABOR/MACE-MP-0b3_GaborShift.model",
        'device' : 'cuda',
        'default_dtype' : 'float64',
    }
elif mlip == "MACE-Wang":
    mlip_name = "MACE"
    mlip_kwargs = {
        'model' : "/home/lwalterb/hdp_ehull/MACE_ALT_GABOR/MACE-MP-0b3_WangShift.model",
        'device' : 'cuda',
        'default_dtype' : 'float64',
    }

else:
    mlip_name = mlip
    mlip_kwargs = {
            "default_dtype":"float64",
            'device':'cuda',
                   }


data_dir = "MLIP_race"
os.makedirs(data_dir, exist_ok=True)

n_batches = ceil(len(subsys_df)/batchsize)
try:
    batch_df = subsys_df.iloc[batchnum * batchsize : (batchnum+1)*batchsize]
except IndexError:
    batch_df = subsys_df.iloc[batchnum*batchsize:]

batch_df['structure'] = batch_df.apply(lambda row: Structure.from_dict(eval(row['structure_dict'])),axis=1)

print(f"Starting {mlip} for batchnumber {batchnum}")


mpid_energy_dict = mlip_relax_and_get_energies(structure_dict= batch_df['structure'], 
                                                    force_field_name=mlip_name,
                                                    calculator_kwargs=mlip_kwargs,
                                                    relax_kwargs= relax_kwargs,
                                                    optimizer_kwargs = optimizer_kwargs,
                                               )

with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}_batch_{batchnum}.json"), "w") as f:
    json.dump(mpid_energy_dict, f)


print(f"Finished {mlip} flow no. {batchnum}")

