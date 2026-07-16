from ehull_utils import mlip_relax_and_get_energies
import pandas as pd
from pymatgen.core import Structure
import os
import json
import sys
from math import ceil


mlip = sys.argv[1]
batchnum = int(sys.argv[2])
batchsize = 500

to_compute = ["MatterSim","MatPES-r2SCAN", "MACE-MP-0b3", "MACE-MPA-0", "Nequip"]  # Nequip -> Nequip-OAM-XL
assert mlip in to_compute, f"MLIP argument must be in {to_compute}."

## Read chemical systems to simulate and convert structure dicts to structure objects
subsys_df = pd.read_csv('hdp_mp_sybsysandhdps.csv',index_col=0)
subsys_df['structure'] = subsys_df.apply(lambda row: Structure.from_dict(eval(row['structure_dict'])),axis=1)


data_dir = "MLIP_data"
os.makedirs(data_dir, exist_ok=True)

n_batches = ceil(len(subsys_df)/batchsize)
try:
    batch_df = subsys_df.iloc[batchnum * batchsize : (batchnum+1)*batchsize]
except IndexError:
    batch_df = subsys_df.iloc[batchnum*batchsize:]
else:
    print("WEIRD ERROR IN SPLITTING DATA")


print(f"Starting {mlip} for batchnumber {batchnum}")

try:
    mpid_energy_dict = mlip_relax_and_get_energies(structure_dict= subsys_df['structure'], 
                                                    force_field_name=mlip )
except Exception as e:
    print(f"Encountered {str(e)} in {mlip}  exiting now!")
    exit(1)

with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}_batch_{batchnum}.json"), "w") as f:
    json.dump(mpid_energy_dict, f)


print(f"Finished {mlip} flow no. {batchnum}")




# # "float64" for geom. optimization
# calculator_kwargs = {"default_dtype": "float64", 
#                     #  "compile_path": "/home/htc/kueltzen/HK_hydride_project/nequip-model/NequIP-OAM-XL.nequip.pt2", 
#                      "device": "cpu"} \
#     if mlip == "Nequip" else {"default_dtype": "float64"}

