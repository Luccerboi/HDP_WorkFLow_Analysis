"""
RunWorkFlow.py
Automated workflow manager for HDP calculations on HPC clusters.
Handles job queuing, submission, and monitoring across multiple clusters (ccp20, ccp22, ccp23).
"""

from HTProcessLedger import ProcessLedger, initialize_compounds
from HDPCompound import Compound
from JobDefinitions import *
from job_handler import submit_and_wait, check_availability
import os
import json
import concurrent.futures as cf
import subprocess
import time
import warnings

warnings.filterwarnings("ignore")

# ============================================================================
# Use python-dotenv to load settings from .env file.
# ============================================================================
from dotenv import load_dotenv

load_dotenv()
JobsFile = os.environ["JOB_JSON_PATH"]
DBlocation = os.environ["DATABASE_PATH"]
ledgerfile = os.environ["PROCESSLEDGER_FILENAME"]
MAX_PROCS = os.environ["MAX_NUM_JOBS"]

# Load ProcessLedger and regain the List of Compounds.
p = ProcessLedger(JobsFile, StartPath=DBlocation, ledger_filename=ledgerfile)
loc = p.RestartLedger()
_ = p.ReadFullLedger()

# ============================================================================
# CHANGE ONLY IF YOU WANT TO DELETE ALL DATA AND RESTART DATABASE
# ============================================================================
restart = False
InputCompsFile = "AllCompsInput"


def create_single_job(
    JobClass: str,
    ledger: ProcessLedger,
    comp: Compound,
    JobName: str,
    server="ccp20,ccp22",
):
    """
    Factory function to create a job instance based on the specified job class type.

    Args:
        JobClass (str): Type of job to create ('Relaxation', 'SpinStateScan', 'SimpleVasp',
                        'SimpleLobster', or 'PreLobster')
        ledger (ProcessLedger): The process ledger for tracking job status
        comp (Compound): The Compound object to associate with the job
        JobName (str): Name/identifier for the job
        server (str): Target server(s) for job submission (default: "ccp20,ccp22")

    Returns:
        GeneralJob: An instance of the appropriate job class
    """
    if JobClass == "Relaxation":
        return Relaxation(ledger, comp, JobName, server)
    elif JobClass == "SpinStateScan":
        return SpinStateScan(ledger, comp, JobName, server)
    elif JobClass == "SimpleVasp":
        return SimpleVasp(ledger, comp, JobName, server)
    elif JobClass == "SimpleLobster":
        return SimpleLobster(ledger, comp, JobName, server)
    elif JobClass == "PreLobster":
        return PreLobster(ledger, comp, JobName, server)


def AssignNextJobs(
    queue: dict,
    available: dict,
    list_of_compounds: list,
    ledger: ProcessLedger,
    JobSettings: dict,
    max_new_jobs=4,
):
    """
    Assign the next batch of jobs to available nodes on the HPC clusters.

    Intelligently distributes jobs based on:
    - Available nodes per cluster
    - Node requirements for each job type
    - Job queue priorities

    Args:
        queue (dict): Dictionary mapping job types to lists of Compound IDs waiting to be processed
        available (dict): Dictionary of cluster names to number of available nodes
        list_of_compounds (list): List of all Compound objects
        ledger (ProcessLedger): Process ledger for tracking job progress
        JobSettings (dict): Dictionary containing job configuration (node requirements, etc.)
        max_new_jobs (int): Maximum number of new jobs to assign in this call (default: 4)

    Returns:
        tuple: (job_list, updated_queue) - list of newly created jobs and updated queue dictionary
    """
    # Build a dictionary of node requirements for each job type
    node_usage = {}
    for job in queue.keys():
        if len(queue[job]) > 0:
            node_usage.update({job: JobSettings[job]["nnodes"]})

    job_list = []
    usage_list = list(node_usage.items())
    usage_list.reverse()  # Process job types in reverse order

    # Iterate through each available cluster
    for server, nnodes in available.items():
        freenodes = nnodes

        while freenodes >= 0:
            # For each job type, try to fit jobs into available nodes
            for ii in range(len(usage_list)):
                # Keep assigning jobs of this type while there are enough free nodes
                while freenodes >= usage_list[ii][1]:
                    try:
                        # Get next Compound ID from the queue for this job type
                        next_compID = queue[usage_list[ii][0]].pop(0)
                    except IndexError:
                        # Queue is empty for this job type
                        break

                    # Create a new job instance for this Compound
                    newjob = create_single_job(
                        JobSettings[usage_list[ii][0]]["JobClass"],
                        ledger,
                        [x for x in list_of_compounds if x.CompID == next_compID][0],
                        usage_list[ii][0],
                        server,
                    )

                    job_list.append(newjob)
                    max_new_jobs -= 1

                    # Stop if we've reached the maximum number of new jobs to assign
                    if max_new_jobs == 0:
                        return job_list, queue

                    # Reduce available nodes on this server
                    freenodes -= usage_list[ii][1]

            break

    return job_list, queue


def GetQueueLength(queue: dict) -> int:
    """
    Calculate the total number of jobs remaining in all queues.

    Args:
        queue (dict): Dictionary mapping job types to lists of Compound IDs

    Returns:
        int: Total number of jobs across all job type queues
    """
    length = 0
    for key, val in queue.items():
        length += len(val)

    return length


def RunJob(job: GeneralJob | None):
    """
    Execute a single job and return its completion status.

    Logs the job start time and details, then runs the job on the assigned server.

    Args:
        job (GeneralJob or None): The job to execute. If None, function returns without action.

    Returns:
        The job's return/output value upon completion
    """
    if type(job) == None:
        return

    print(
        f"{time.strftime('%H:%M:%S')}: starting {job.AssignedCompName} for job {job.JobName} on {job.AssignedServer}"
    )
    out = job.Run()
    return out


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

# Initialize the process ledger to track job completion and status
ledger = ProcessLedger(JobsFile, StartPath=DBlocation, ledger_filename=ledgerfile)

# Load job specifications from JSON file
with open(JobsFile, "r") as f:
    jobdata = json.load(f)

# Initialize or restore the database and Compound list
if restart:
    # Start fresh: delete all existing data and reinitialize
    loc = initialize_compounds(InputCompsFile)
    os.chdir(DBlocation)
    os.system("rm -r *")
    ledger.StartNewLedger(list_of_compounds=loc)
else:
    # Resume from existing database state
    loc = ledger.RestartLedger()

# Change to database directory for calculations
os.chdir(DBlocation)

# Get current completion status and job queue
ov = ledger.GetCompletionOverview()
queue = ledger.GetQueue()

# Check available nodes on each cluster
avail = check_availability()
print(avail)

queuelength = GetQueueLength(queue)


# ============================================================================
# MAIN PROCESSING LOOP
# ============================================================================
while queuelength > 0:
    # Check cluster availability before assigning new jobs
    avail = check_availability()

    # Assign next batch of jobs to available nodes
    JobList, queue = AssignNextJobs(
        queue, avail, loc, ledger, jobdata, max_new_jobs=MAX_PROCS
    )
    print(JobList)

    # Log each job being started
    for job in JobList:
        print(f"{time.strftime('%H:%M:%S')}: starting job {job}")

    # Execute jobs in parallel using process pool
    with cf.ProcessPoolExecutor(max_workers=MAX_PROCS) as executor:
        # Create a mapping of futures to job objects
        future_to_comps = {executor.submit(RunJob, job): job for job in JobList}

        # Monitor job completion and submit new jobs as resources become available
        while future_to_comps:
            # Wait for at least one job to complete
            for future in cf.wait(future_to_comps, return_when=cf.FIRST_COMPLETED)[0]:
                comp = future_to_comps[future]
                try:
                    data = future.result()
                except Exception as exc:
                    print("%r generated an exception: %s" % (comp, exc))
                else:
                    print(
                        "%s: \t %r returned %s "
                        % (time.strftime("%H:%M:%S"), comp, data)
                    )

                # Remove the completed future from tracking
                del future_to_comps[future]

                # Backup ledger after each job completion to prevent data loss
                os.system(f"cp {ledger.LedgerFile} {ledger.LedgerBackup}")

                # Check if more jobs are available and assign them to free workers
                if GetQueueLength(queue) > 1:
                    JobList, queue = AssignNextJobs(
                        queue, avail, loc, ledger, jobdata, max_new_jobs=MAX_PROCS
                    )
                    for new_job in JobList:
                        future_to_comps[executor.submit(RunJob, new_job)] = new_job
                else:
                    print("end of Queue reached, reinitializing Queue now...")
                    queue = ledger.GetQueue()

        # Update queue status after all current jobs complete
        total_queue = ledger.GetQueue()
        queuelength = GetQueueLength(total_queue)
        if GetQueueLength(queue) == 0:
            queue = ledger.GetQueue()
