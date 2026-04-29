from CsBBCl_ht_new import *
from job_handler import submit_and_wait, check_availability
import os

# from math import ceil, floor
# import json
# import concurrent.futures as cf
# import subprocess
# import time

import time
import pandas as pd
JobsFile = "/home/lwalterb/hdp_project/ht_programfiles/HDP_JOBS_v2507.json" #location of the jobs file
InputCompsFile = "/home/lwalterb/hdp_project/ht_programfiles/AllCompsInput" #location of the input compounds
DBlocation = "/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/AllCompsNewWF" #location of the database
restart = False #if True, the database will be deleted and reinitialized
ledgerfile = "BAMadapted_NoCs6s.csv" #name of the ledger file

p = ProcessLedger(JobsFile,StartPath=DBlocation,ledger_filename=ledgerfile)
ov = p.GetCompletionOverview()

def GetQueueLength(queue:dict):
    length = 0
    for key, val in queue.items():
        length += len(val)

    return length

def print_errors():
    p = ProcessLedger(JobsFile,StartPath=DBlocation,ledger_filename=ledgerfile)
    q = p.GetQueue()
    ov = p.GetCompletionOverview()

    ov2 = ov[['1Rel', '2Spins', '3Pre','4HSE']]

    errors = ov2[((ov2 != "0") & (ov2 != "1")).any(axis=1)]
    lob_done = ov[((ov['5LOB'] == "-1") | (ov['5LOB']== "1"))]

    print(ov.loc[errors.index])
    return



def renew_overview(queue,overview):
    q = queue
    ov = overview

    print(f"Current progress: {GetQueueLength(q)} in queue in total \n\t1Rel:\t {len(q['1Rel'])} \n\t2Spins:\t {len(q['2Spins'])} \n\t3Pre:\t {len(q['3Pre'])} \n\t4HSE:\t {len(q['4HSE'])}  \n\t5LOB:\t {len(q['5LOB'])}\n\n")

    ov2 = ov[['1Rel', '2Spins', '3Pre','4HSE']]

    errors = ov2[((ov2 != "0") & (ov2 != "1")).any(axis=1)]
    lob_done = ov[((ov['5LOB'] == "-1") | (ov['5LOB']== "1"))]

    
    print(f"Fully completed Comps:\t {len(lob_done)}")
    print(f"Erroneous Comps:\t {len(errors)}\n")
    print(f"------------------------------------------------------------------------------------------------------------------------------------------------\n")
    return ov

def get_current_jobs():
	import subprocess
	queproc = subprocess.run(["squeue","-u","lucw"],text=True,shell=True,capture_output=True)
	linearray = queproc.stdout.strip().split('\n')
	compidslist = []
	compnameslist = []
	for line in linearray:
		for word in line.strip().split():
			if '_Cs' in word:
				compidslist.append(word.split('_')[0])
				compnameslist.append(word.split('_')[-1])
	compidlist = list(map(int, compidslist))
	compnameslist
	return compidlist

def cycle_ovdf(overview: pd.DataFrame, offset: int = 0):
	return(oververview.iloc[offset])	
	

if __name__ == "__main__":
	JobsFile = "/home/lwalterb/hdp_project/ht_programfiles/HDP_JOBS_v2507.json" #location of the jobs file
	InputCompsFile = "/home/lwalterb/hdp_project/ht_programfiles/AllCompsInput" #location of the input compounds
	DBlocation = "/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/AllCompsNewWF" #location of the database
	restart = False #if True, the database will be deleted and reinitialized
	ledgerfile = "BAMadapted_NoCs6s.csv" #name of the ledger file

	
	p = ProcessLedger(JobsFile,StartPath=DBlocation,ledger_filename=ledgerfile)

	while True:
		os.system('cls' if os.name == 'nt' else 'clear')
		ov = p.GetCompletionOverview()
		q = p.GetQueue()
		runningjobs = get_current_jobs()
		runningcomps = [ov.iloc[x-1000].name for x in runningjobs]
		a = compound()
		for ix, job in enumerate(runningjobs):
			#a.CompID = runningcomps[ix]

			#p.SetSingleValue(a,'4HSE','completed','CheckAfterOSupdate')
			
			renew_overview(q,ov)
			print(f"Currently running comps:\n{runningcomps}\n")
			a = subprocess.run(["squeue","-u","lucw"],capture_output =True,text=True)
			print(a.stdout)
			print(ov.iloc[job-1005:].head(10))
			time.sleep(10)
			os.system('cls' if os.name == 'nt' else 'clear')
