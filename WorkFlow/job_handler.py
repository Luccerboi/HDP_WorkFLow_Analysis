import numpy as np
import pandas as pd
import os
import subprocess
import time

# ============================================================================
# Use python-dotenv to load settings from .env file.
# ============================================================================
from dotenv import load_dotenv

load_dotenv()
# Retrieve which partitions SLURM is allowd to use.
allowed_partitions = eval(os.environ["ALLOWED_PARTITIONS"])


def submit_and_wait(submission="vasp.sub"):
    """
    This function will submit vasp.sub to SLURM an retrieve the given job ID.
    It will then go on to check the status of this job every 30 seconds,
    once the job completes, fails, or is canceled it will return
    Args:
        submission(str): the filename of the submission scripts, defaults to "vasp.sub"

    Returns:
        (str): the function returns the final state of the calculation in SLURM, which should be COMPLETED
                but it could also be FAILED or CANCELED
    """
    assert os.path.isfile(submission)
    # os.chdir(calc.path)
    submit = subprocess.Popen(["sbatch", submission], stdout=subprocess.PIPE, text=True)
    get_jobid = subprocess.Popen(
        ["awk", "{print $NF}"], stdin=submit.stdout, stdout=subprocess.PIPE, text=True
    )
    jobid, error = get_jobid.communicate()
    done = False
    while done == False:

        if "lob".lower() in submission.lower():
            # LOBSTER cacluations take significanlty less time, so check more often.
            time.sleep(10)
        else:
            time.sleep(30)

        check_state = subprocess.Popen(
            ["sacct", "-j", jobid], stdout=subprocess.PIPE, text=False
        )
        get_status = subprocess.Popen(
            ["awk", "/%i +/ {print $6}" % (int(jobid))],
            stdin=check_state.stdout,
            stdout=subprocess.PIPE,
            text=True,
        )
        state, _ = get_status.communicate()
        state = state.strip()

        if state != "RUNNING" and state != "PENDING":
            done = True

    return state


def check_availability(
    allowed_partitions: list[str] = allowed_partitions,
) -> dict[str, int]:
    """This function will prompt SLURM to check how many nodes are available on CCP20 and CCP22
    Based on the number of nodes idle for each of these it will return the number of nodes that are free per
    cluster.

    Args:
        allowed_partitions (list[str], optional): list of partition names you want to consider for submitting jobs. Defaults to ['ccp20','ccp22'].

    Returns:
        dict[str,int]: dictionary containing with the available partitions as keys and number of idle nodes as values.
    """
    p1 = subprocess.Popen(
        ["sinfo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    partition_scan_string = "|".join(allowed_partitions)

    p2 = subprocess.Popen(
        ["awk", '/(%s)/&&/idle/ {print $1 "   " $4}'.format(partition_scan_string)],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    out, _ = p2.communicate()
    if out == "":  # all nodes busy
        return {x: 0 for x in allowed_partitions}

    avail_ls = out.strip().split("\n")  # output of call line seperated

    overview = {}
    for line in avail_ls:
        avail_ws = line.split("   ")

        if avail_ws[0] in overview.keys():
            overview[avail_ws[0]] += int(avail_ws[1])
        else:
            overview.update({avail_ws[0]: int(avail_ws[1])})

    return overview


def check_scfcompletion(logfile, erroroutput="encounteredErrors.tmp"):
    """This function will read in the logfile as it has been created by VASP. It will check for each electronic SCF cycle
    if the cycle has converged, if so it returns a 1 for that step, otherwise it will return 0. If an error is encountered
    the script will write 2, and if a relaxation is run a 3 will be printed if the required accuracy has been reached.

    Args:
        logfile (str): filename of VASP log
        errorfile (str): filename of the desired file to store encountered errors

    Returns:
        list: a list of intigers indicating how the SCF cycle(s) ended
                0   -   successfully converged
                1   -   failed to converge
                2   -   special flag for relaxation successful
                ''  -   if the run was ended prematurely it will return an empty string
               -1  -   an error occurred
    """
    p1 = subprocess.run(
        [
            "awk",
            'BEGIN{rv=0; print "\\nNew run started --->\\n\\n" >> "%s"};/self-consistency was not achieved/ {rv=1}; /F=/ {print rv; rv=0}; \
                            /E{7} * R * R * R * R * O{7}/ {print "-1"}; /reached required accuracy - stopping structural energy minimisation/ {print "2"};\
                            /E{7} * R * R * R * R * O{7}/, /---->/ {if(!/E{7} * R/ && !/---->/)print $0 >> "%s"}'
            % (erroroutput, erroroutput),
            logfile,
        ],
        stdout=subprocess.PIPE,
        text=True,
    )

    numstrings = p1.stdout.strip().split("\n")
    try:
        scf_cycles = list(map(int, numstrings))
    except ValueError as err:
        # raise Exception("Unexpected value encountered when checking SCF convergence:{err}")
        return numstrings
    return scf_cycles
