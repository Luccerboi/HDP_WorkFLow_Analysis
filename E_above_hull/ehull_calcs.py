import json
import os
import pandas as pd
from pymatgen.core import Structure
import sys


from utils_hydride import mlip_relax_and_get_energies


mlip = sys.argv[1]
split_id = sys.argv[2]

to_compute = ["MatPES-r2SCAN", "MACE-MP-0b3", "MACE-MPA-0", "Nequip"]  # Nequip -> Nequip-OAM-XL
assert mlip in to_compute, f"MLIP argument must be in {to_compute}."


# "float64" for geom. optimization
calculator_kwargs = {"default_dtype": "float64", 
                     "compile_path": "/home/htc/kueltzen/HK_hydride_project/nequip-model/NequIP-OAM-XL.nequip.pt2", 
                     "device": "cpu"} \
    if mlip == "Nequip" else {"default_dtype": "float64"}

data_dir = "data"
os.makedirs(mlip, exist_ok=True)

hydride_df = pd.read_json(os.path.join(data_dir,
                                       "df_new_and_existing_binary_ternary_hydride-chemical-systems_mp_full-data.json"))

flow_size = 30
hydride_df_list = [hydride_df[i:i+flow_size] for i in range(0, hydride_df.shape[0], flow_size)]
assert int(split_id) in list(range(0, len(hydride_df_list)))

subdf = hydride_df_list[int(split_id)]
input_dict = {mpid: Structure.from_dict(r["structure_dict"]) for mpid, r in subdf.iterrows()}

print(f"Starting {mlip} flow no. {split_id}"
      )
try:
    mpid_energy_dict = mlip_relax_and_get_energies(structure_dict=input_dict, 
                                                    force_field_name=mlip, 
                                                    calculator_kwargs=calculator_kwargs)
except Exception as e:
    print(f"Encountered {str(e)} in {mlip} flow no. {split_id}, exiting now!")
    exit(1)

with open(os.path.join(mlip, f"mpid_energy_dict_{mlip}_flow{split_id}.json"), "w") as f:
    json.dump(mpid_energy_dict, f)
print(f"Finished {mlip} flow no. {split_id}.")
 