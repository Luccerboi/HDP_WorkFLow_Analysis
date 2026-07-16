from ehull_utils import mlip_relax_and_get_energies
import pandas as pd
from pymatgen.core import Structure
import os
import json
import sys



mlip = sys.argv[1]


to_compute = ["MatterSim","MatPES-r2SCAN", "MACE-MP-0b3", "MACE-MPA-0", "Nequip"]  # Nequip -> Nequip-OAM-XL
assert mlip in to_compute, f"MLIP argument must be in {to_compute}."

## Read chemical systems to simulate and convert structure dicts to structure objects
subsys_df = pd.read_csv('hdp_mp_subsys_df.csv',index_col=0)
subsys_df['structure'] = subsys_df.apply(lambda row: Structure.from_dict(eval(row['structure_dict'])),axis=1)

hdp_df = pd.read_csv('hdp_hdps_mlipenergies.csv', index_col=0)
hdp_df['structure'] = hdp_df.apply(lambda row: Structure.from_dict(eval(row['structure_dict'])))

data_dir = "MLIP_data"
os.makedirs(data_dir, exist_ok=True)


print(f"Starting {mlip} for subsystems")

try:
    mpid_energy_dict = mlip_relax_and_get_energies(structure_dict= subsys_df['structure'], 
                                                    force_field_name=mlip )
except Exception as e:
    print(f"Encountered {str(e)} in {mlip}  exiting now!")
    exit(1)

with open(os.path.join(data_dir, f"mpid_energy_dict_{mlip}.json"), "w") as f:
    json.dump(mpid_energy_dict, f)

print(f"Finished subsystems with {mlip}")

try:
    hdp_energy_dict = mlip_relax_and_get_energies(
        structure_dict=hdp_df['structure'],
        force_field_name=mlip
    )
except Exception as e:
    print(f"Encountered {str(e)} in {mlip}  exiting now!")
    exit(1)

with open(os.path.join(data_dir,f"hdp_energy_dict_{mlip}.json"),"w") as f:
    json.dump(hdp_energy_dict,f)

print(f"Finished {mlip} flow no.")




# # "float64" for geom. optimization
# calculator_kwargs = {"default_dtype": "float64", 
#                     #  "compile_path": "/home/htc/kueltzen/HK_hydride_project/nequip-model/NequIP-OAM-XL.nequip.pt2", 
#                      "device": "cpu"} \
#     if mlip == "Nequip" else {"default_dtype": "float64"}

