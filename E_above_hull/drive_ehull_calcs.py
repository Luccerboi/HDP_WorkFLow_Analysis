import subprocess
from math import ceil
import pandas as pd
import sys

worker_id = int(sys.argv[1])
n_workers = int(sys.argv[2])
mlip = sys.argv[3]

batchsize = 25
full_df = pd.read_csv(f'MissingStrucs_{mlip}.csv',index_col=0)
n_batches = ceil(len(full_df)/batchsize)

for i in range( worker_id , n_batches, n_workers):
    subprocess.run(
            ["python", "run_ehull_calcs.py", mlip, str(batchsize), str(i)], check=True
            )
