import subprocess
from math import ceil
import pandas as pd
import sys

worker_id = int(sys.argv[1])
n_workers = int(sys.argv[2])
mlip = sys.argv[3]

batchsize = 10 
full_df = pd.read_csv('hdp_mp_subsysandhdps.csv',index_col=0)
n_batches = ceil(len(full_df)/batchsize)

for i in range(135 + worker_id , n_batches, n_workers):
    subprocess.run(
            ["python", "run_grace_calcs.py", mlip, str(batchsize), str(i)], check=True
            )
