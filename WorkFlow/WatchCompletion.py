"""
WatchCompletion.py
Utility to monitor HDP workflow progress and display queue status.
"""

import os
import time
import subprocess
import pandas as pd

from CsBBCl_ht_new import ProcessLedger

# ============================================================================
# Configuration
# ============================================================================

JOBS_FILE = "/home/lwalterb/hdp_project/ht_programfiles/HDP_JOBS_v2507.json"
INPUT_COMPS_FILE = "/home/lwalterb/hdp_project/ht_programfiles/AllCompsInput"
DB_LOCATION = "/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/AllCompsNewWF"
LEDGER_FILE = "BAMadapted_NoCs6s.csv"
DEFAULT_USER = "lucw"


def GetQueueLength(queue: dict) -> int:
    """Return the total number of items in all job queues."""
    return sum(len(val) for val in queue.values())

def print_errors() -> None:
    """Print compounds with intermediate errors in the current ledger overview."""
    ledger = ProcessLedger(JOBS_FILE, StartPath=DB_LOCATION, ledger_filename=LEDGER_FILE)
    overview = ledger.GetCompletionOverview()

    status_columns = ['1Rel', '2Spins', '3Pre', '4HSE']
    errors = overview[status_columns][((overview[status_columns] != "0") & (overview[status_columns] != "1")).any(axis=1)]

    print("Error summary for intermediate steps:")
    print(overview.loc[errors.index])


def renew_overview(queue: dict, overview: pd.DataFrame) -> pd.DataFrame:
    """Display a summary of queue lengths, completed jobs, and errors."""
    print(
        f"Current progress: {GetQueueLength(queue)} in queue in total\n"
        f"\t1Rel:\t {len(queue.get('1Rel', []))}\n"
        f"\t2Spins:\t {len(queue.get('2Spins', []))}\n"
        f"\t3Pre:\t {len(queue.get('3Pre', []))}\n"
        f"\t4HSE:\t {len(queue.get('4HSE', []))}\n"
        f"\t5LOB:\t {len(queue.get('5LOB', []))}\n"
    )

    status_columns = ['1Rel', '2Spins', '3Pre', '4HSE']
    errors = overview[status_columns][((overview[status_columns] != "0") & (overview[status_columns] != "1")).any(axis=1)]
    lob_done = overview[((overview['5LOB'] == "-1") | (overview['5LOB'] == "1"))]

    print(f"Fully completed comps:\t {len(lob_done)}")
    print(f"Erroneous comps:\t {len(errors)}\n")
    print("-" * 120)

    return overview

def get_current_jobs(user: str = DEFAULT_USER) -> list[int]:
    """Return the list of currently running compound IDs from squeue output."""
    result = subprocess.run(
        ["squeue", "-u", user],
        capture_output=True,
        text=True,
        shell=False,
    )

    comp_ids = []
    for line in result.stdout.strip().splitlines():
        for word in line.split():
            if '_Cs' in word:
                try:
                    comp_ids.append(int(word.split('_')[0]))
                except ValueError:
                    continue

    return comp_ids

def cycle_ovdf(overview: pd.DataFrame, offset: int = 0) -> pd.Series:
    """Return the overview row at the requested offset."""
    return overview.iloc[offset]


def clear_screen() -> None:
    """Clear the terminal screen for a clean status display."""
    os.system('cls' if os.name == 'nt' else 'clear')


def main() -> None:
    """Run the watch loop that prints current queue and running job status."""
    ledger = ProcessLedger(JOBS_FILE, StartPath=DB_LOCATION, ledger_filename=LEDGER_FILE)

    try:
        while True:
            clear_screen()

            overview = ledger.GetCompletionOverview()
            queue = ledger.GetQueue()
            running_jobs = get_current_jobs()

            renew_overview(queue, overview)

            running_comps = []
            for job_id in running_jobs:
                row_index = job_id - 1000
                if 0 <= row_index < len(overview):
                    running_comps.append(overview.iloc[row_index].name)

            print(f"Currently running comps:\n{running_comps}\n")

            squeue_output = subprocess.run(
                ["squeue", "-u", DEFAULT_USER],
                capture_output=True,
                text=True,
                shell=False,
            )
            print(squeue_output.stdout)

            if running_jobs:
                sample_index = max(0, running_jobs[-1] - 1005)
                print(overview.iloc[sample_index:sample_index + 10])
            else:
                print("No running jobs detected.")

            time.sleep(10)
    except KeyboardInterrupt:
        print("\nWatchCompletion stopped by user.")


if __name__ == "__main__":
    main()
