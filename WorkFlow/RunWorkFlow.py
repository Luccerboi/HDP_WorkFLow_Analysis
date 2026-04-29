from CsBBCl_ht_new import *
from job_handler import submit_and_wait, check_availability
import os
import json
import concurrent.futures as cf
import subprocess
import time
import warnings
warnings.filterwarnings("ignore")

JobsFile = "/home/lwalterb/hdp_project/ht_programfiles/HDP_JOBS_v2507.json" #location of the jobs file
InputCompsFile = "/home/lwalterb/hdp_project/ht_programfiles/AllCompsInput" #location of the input compounds
DBlocation = "/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/AllCompsNewWF" #location of the database
restart = False #if True, the database will be deleted and reinitialized
ledgerfile = "BAMadapted_NoCs6s.csv" #name of the ledger file

def create_single_job(JobClass:str,ledger:ProcessLedger, comp:compound, JobName:str, server="ccp20,ccp22"):
    if JobClass == "Relaxation":
        return Relaxation(ledger, comp, JobName, server)
    elif JobClass == "SpinStateScan":
        return SpinStateScan(ledger, comp,JobName, server)
    elif JobClass == "SimpleVasp":
        return SimpleVasp(ledger, comp, JobName, server)
    elif JobClass == "SimpleLobster":
        return SimpleLobster(ledger, comp, JobName, server)
    elif JobClass == "PreLobster":
        return PreLobster(ledger, comp, JobName, server)


def AssignNextJobs(queue:dict, available:dict,list_of_compounds:list,ledger:ProcessLedger, JobSettings:dict,max_new_jobs=4):

    node_usage = {}    
    for job in queue.keys():
        if len(queue[job])>0:
            node_usage.update({job:JobSettings[job]['nnodes']})

    job_list = []
    usage_list = list(node_usage.items())
    usage_list.reverse()
    
    
    
    for server, nnodes in available.items():
        freenodes = nnodes

        while freenodes>=0:

            for ii in range(len(usage_list)):

                while freenodes >= usage_list[ii][1]:
                    try:
                        next_compID = queue[usage_list[ii][0]].pop(0)
                    except IndexError:
                        #print(f"Tried to assign comp for job {usage_list[ii][0]}, but queue was empty...")
                        break

                    #now create a job
                    #print(usage_list[ii][0])
                    newjob = create_single_job(JobSettings[usage_list[ii][0]]["JobClass"],ledger, [x for x in list_of_compounds if x.CompID == next_compID][0],usage_list[ii][0],server)
                    #print(f"Assigned {next_compID} to step {usage_list[ii][0]}")
                    
                    job_list.append(newjob)
                    max_new_jobs -= 1
                    if max_new_jobs == 0:
                        return job_list, queue
                    freenodes -= usage_list[ii][1]

            break
        
    # if len(job_list) == 0:
    #     print("Cluster is full, I'll wait for 5 minutes...")
    #     time.sleep(300)


    return job_list, queue

def GetQueueLength(queue:dict):
    length = 0
    for key, val in queue.items():
        length += len(val)

    return length

def RunJob(job:GeneralJob):
    print(f"{time.strftime('%H:%M:%S')}: starting {job.AssignedCompName} for job {job.JobName} on {job.AssignedServer}")
    out = job.Run()
    return out



ledger = ProcessLedger(JobsFile,StartPath=DBlocation,ledger_filename=ledgerfile)

with open(JobsFile,'r') as f:
    jobdata = json.load(f)

if restart:
    loc = initialize_compounds(InputCompsFile)
    os.chdir(DBlocation)
    os.system("rm -r *")
    ledger.StartNewLedger(list_of_compounds= loc)
else:
    
    loc = ledger.RestartLedger()

os.chdir(DBlocation)

ov = ledger.GetCompletionOverview()
#print(ov)



queue = ledger.GetQueue()
#print(queue)
avail = check_availability()
print(avail)
queuelength = GetQueueLength(queue)

max_procs = 4

while queuelength > 1:
    avail = check_availability()
    JobList, queue = AssignNextJobs(queue,avail, loc,ledger, jobdata,max_new_jobs=max_procs)
    print(JobList)
    for job in JobList:
        print(f"{time.strftime('%H:%M:%S')}: starting job {job}")
        RunJob(job)
    # JobList = []
    # for ii in range(max_procs):
        # next_compID = queue['5LOB'].pop(0)
        # newjob = create_single_job("SimpleLobster",ledger, [x for x in loc if x.CompID == next_compID][0],'5LOB')
        # JobList.append(newjob)

    # RunJob(JobList[1])
    with cf.ProcessPoolExecutor(max_workers=max_procs) as executor:
        future_to_comps = {executor.submit(RunJob, job): job for job in JobList}
        
        while future_to_comps:

            for future in cf.wait(future_to_comps, return_when=cf.FIRST_COMPLETED )[0]:
                comp = future_to_comps[future]
                try:
                    data = future.result()
                except Exception as exc:
                    print('%r generated an exception: %s' % (comp, exc))
                else:
                    print("%s: \t %r returned %s " % (time.strftime('%H:%M:%S'), comp, data) )
                
                # Remove the completed future
                del future_to_comps[future]

                #backup Ledger
                os.system(f"cp {ledger.LedgerFile} {ledger.LedgerBackup}")
                
                # Check if there are more jobs to submit
                if GetQueueLength(queue) > 1:
                    try:
                        next_compID = queue['5LOB'].pop(0)
                    except IndexError:
                        print('end of Queue reached, reinitializing Queue now...')
                        queue = ledger.GetQueue()
                        continue
                    
                    new_job = create_single_job("SimpleLobster",ledger, [x for x in loc if x.CompID == next_compID][0],'5LOB')
                    future_to_comps[executor.submit(RunJob, new_job)] = new_job
                    if GetQueueLength(queue) <= 2:
                        print('end of Queue reached, reinitializing Queue now...')
                        queue = ledger.GetQueue()

                    # if new_job:
                    #     for job in new_job:
                    #         future_to_comps[executor.submit(RunJob, job)] = job

        total_queue = ledger.GetQueue()
        queuelength = GetQueueLength(total_queue)
        if GetQueueLength(queue) <= 2:
            queue = ledger.GetQueue()



