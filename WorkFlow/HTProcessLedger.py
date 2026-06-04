import numpy as np
import pandas as pd  # type: ignore
import os
import json
from pathlib import Path
from HDPCompound import Compound
import ast
from filelock import FileLock, Timeout


def initialize_compounds(
    file: str = "./testinput",
) -> list[Compound]:
    """
    This function will take the input file, create a Compound object for each entry, create the directory and write the general info file
    The input file shoud be parsed as follows "ABBX, 'A', 'B1', nB1, 'B2', nB2, 'X', where '*' indicates the chemical symbol of the element,
    and nB indicates the ionic state of the B site (integer)

    Since starting a ProcessLedger adds a CompID to each composition, the input can have an additional input for CompID. 
    If the input has length 8, this function will assign the final column as the CompID to facilitate reinitialization of the ProcessLedger

    The function will return list[Compound] containig all Compounds
    """
    if not os.path.isfile(file):
        raise TypeError("intialization file could not be found/interperted as file")

    list_of_compounds = []
    with open(file, "r") as f:
        for line in f.readlines():
            inputs = line.strip().strip("\n").split(",")
            new_comp = Compound(
                inputs[1],
                inputs[2],
                int(inputs[3]),
                inputs[4],
                int(inputs[5]),
                inputs[6],
            )
            # If there is an extra column, use this as the CompID
            if len(inputs) == 8:
                new_comp.CompID = str(inputs[-1]) + "_" + str(repr(new_comp))

            list_of_compounds.append(new_comp)

    return list_of_compounds


class ProcessLedger:
    """This class takes care of tracking which calculations are done, storing additional INCAR settings and calculation paths.
    The class is essentially based on large pandas.DataFrame which is saved to a .csv at every iteration.
    A back-up copy is also created by the workflow at every iterations, since on multiple occasions the .csv file got corrupted while writing.
    It also creates in UsedInput back-up, storing the list of compositions with their assigned CompIDs

    The class uses filelock to assure only one process is working on the csv at the same time.

    There are methods for:
    - Starting new Ledger
    - Appending compositions to existing Ledger
    - Restarting Ledger from file
    - Reading the full Ledger
    - Getting info on a single job
    - Get an Overview of job Completions
    - Getting a queue based on CompletionOverview
    - Setting a Single Value in the Ledger
    - Restoring the Ledger information after critical error by scanning through the database directories.
    """
    def __init__(
        self,
        JobSettingsPath: str,
        StartPath: str = ".",
        ledger_filename: str = "JobInformationLedger.csv",
    ):
        """Initialize ProcessLedger object

        Args:
            JobSettingsPath (str): Path to json conting the standard settings for the different jobs
            StartPath (str, optional): The base directory of the database, where the ledger will be stored as well. Defaults to ".".
            ledger_filename (str, optional): Filename of the ProcessLedger csv. Defaults to "JobInformationLedger.csv".
        """
        # the ledger will be based in the parent directory of all jobfiles
        if StartPath == ".":
            self.basepath = Path(os.getcwd())
        else:
            self.basepath = Path(StartPath)

        # Load in the required jobs
        assert JobSettingsPath.endswith(".json"), "JobSettings needs to be JSON file"
        self.JobsPath = JobSettingsPath

        with open(self.JobsPath, "r") as jfile:
            # this stores the full info per step
            self.JobInfo = json.load(jfile)

        self.Jobs = list(self.JobInfo.keys())

        self.LedgerFile = self.basepath / ledger_filename
        self.inputbackup = self.basepath / str("UsedInput_" + ledger_filename)
        self.LedgerBackup = self.basepath / str("BackupLedger_" + ledger_filename)

        self.lock = FileLock(str(self.LedgerFile) + ".lock", thread_local=False)
        # this is what info the ledger stores per job
        self.info_per_job = [
            "JobPath",
            "completed",
            "JobID",
            "TimeStamp",
            "IncarExtras",
        ]

    def StartNewLedger(self, list_of_compounds: list[Compound] = []) -> None:
        """This function starts a new ledger for a given list of input Compounds.
        It will create a subdirectory based on the naming in the job-settings json file.
        For LOBSTER calculations it skips making subdirectories. It does require the LOBSTER job to have specified which vasp calculation to use.
        The function immediately also writes a backup of the inputs, where the CompID is also stored.
        The ProcessLedger DataFrame is created, and the .csv is written.

        Args:
            list_of_compounds (list[Compound], optional): list of Compound objects to include in the workflow. Defaults to [].

        Raises:
            Exception: when LOBSTER calculation misses a defined 'use_output_step'
        """
           
        os.chdir(self.basepath)

        ##Setup dataframe for the ledger info
        index = pd.MultiIndex.from_product(
            [list(self.JobInfo.keys()), self.info_per_job]
        )
        ledger = pd.DataFrame(index=index)
        ledger2 = pd.DataFrame(index=index)

        # All comps get a CompID
        # which starts at 1000 for esthetic reasons
        CompID = 1000
        input_backup_file = open(self.inputbackup, "w")

        for comp in list_of_compounds:

            comp.CompID = str(CompID) + "_" + str(repr(comp))
            comp.path = str(self.basepath) + "/" + str(CompID) + "_" + str(repr(comp))
            input_backup_file.write(
                f"{repr(comp)},{comp.A},{comp.B1},{comp.nB1},{comp.B2},{comp.nB2},{comp.X},{CompID}\n"
            )
            try:
                os.mkdir(comp.path)
                os.chdir(comp.path)
            except FileExistsError:
                os.chdir(comp.path)
            except:
                raise Exception(f"Could not make directory for {repr(comp)}")

            CompID += 1

            LedgerInfo = {}
            for JobName, settings in self.JobInfo.items():

                completion = 0

                # If it is a VASP or misc job a subdirectory will be made
                # if the subdir already exists, it will be assumed that the step was
                # completed succesfully
                # For LOBSTER jobs it will set the subdir from another step as working dir
                if settings["job_type"] != "LOBSTER":
                    calc_subdir = str(JobName) + "_" + settings["name"]
                    calcdir_full = os.path.join(comp.path, calc_subdir)

                    os.mkdir(calcdir_full)

                else:
                    assert (
                        "use_output_step" in settings.keys()
                    ), "LOBSTER job is missing required VASP directory.\n Please define 'use_output_step' in JobSettingsJson"
                    subdir = (
                        settings["use_output_step"]
                        + "_"
                        + self.JobInfo[settings["use_output_step"]]["name"]
                    )
                    calcdir_full = os.path.join(comp.path, subdir)

                # store all base information for each step
                LedgerInfo.update(
                    {
                        JobName: {
                            "JobPath": str(calcdir_full),
                            "JobID": "",
                            "completed": completion,
                            "TimeStamp": "",
                            "IncarExtras": {},
                        }
                    }
                )

            # flatten the dict of ledger info and add as new column to dataframe
            flat_dict = {
                (outer_key, inner_key): value
                for outer_key, inner_dict in LedgerInfo.items()
                for inner_key, value in inner_dict.items()
            }

            ledger2[comp.CompID] = ledger.index.map(flat_dict).copy()
        self.LedgerDF = ledger2.sort_index()
        with self.lock:
            ledger2.to_csv(self.LedgerFile)
        input_backup_file.close()
        return

    def AppendNewComps(self, list_of_compounds: list[Compound] = []):
        """Does the same as StartNewLedger, only now it first reads the previous inputs, and appends the new compositions to the list and database.
        Handy if you discorver there is another subclass you would like to include.

        Args:
            list_of_compounds (list[Compound], optional): list of Compound objects to append to ProcessLedger. Defaults to [].

        Raises:
            Exception: LOBSTER calculation does not have VASP step to use defined.
        """

        os.chdir(self.basepath)

        ##Setup dataframe for the ledger info
        index = pd.MultiIndex.from_product(
            [list(self.JobInfo.keys()), self.info_per_job]
        )
        ledger = pd.DataFrame(index=index)
        ledger2 = pd.DataFrame(index=index)

        # All comps get a CompID
        # We read the last used from UsedInput_backup
        with open(self.inputbackup, "r") as f:
            lines = f.readlines()

        last_input = lines[-1].strip().split(",")

        CompID = int(last_input[-1]) + 1
        file = open(self.inputbackup, "a")

        for comp in list_of_compounds:

            comp.CompID = str(CompID) + "_" + str(repr(comp))
            comp.path = str(self.basepath) + "/" + str(CompID) + "_" + str(repr(comp))
            file.write(
                f"{repr(comp)},{comp.A},{comp.B1},{comp.nB1},{comp.B2},{comp.nB2},{comp.X},{CompID}\n"
            )
            try:
                os.mkdir(comp.path)
                os.chdir(comp.path)
            except FileExistsError:
                os.chdir(comp.path)
            except:
                raise Exception(f"Could not make directory for {repr(comp)}")

            CompID += 1

            LedgerInfo = {}
            for JobName, settings in self.JobInfo.items():

                completion = 0

                # If it is a VASP or misc job a subdirectory will be made
                # For LOBSTER jobs it will set the subdir from another step as working dir
                if settings["job_type"] != "LOBSTER":
                    calc_subdir = str(JobName) + "_" + settings["name"]
                    calcdir_full = os.path.join(comp.path, calc_subdir)
                    try:
                        os.mkdir(calcdir_full)
                    except:
                        continue


                else:
                    assert (
                        "use_output_step" in settings.keys()
                    ), "LOBSTER job is missing required VASP directory.\n Please define 'use_output_step' in JobSettingsJson"
                    subdir = (
                        settings["use_output_step"]
                        + "_"
                        + self.JobInfo[settings["use_output_step"]]["name"]
                    )
                    calcdir_full = os.path.join(comp.path, subdir)

                # store all base information for each step
                LedgerInfo.update(
                    {
                        JobName: {
                            "JobPath": str(calcdir_full),
                            "JobID": "",
                            "completed": completion,
                            "TimeStamp": "",
                            "IncarExtras": {},
                        }
                    }
                )

            # flatten the dict of ledger info and add as new column to dataframe
            flat_dict = {
                (outer_key, inner_key): value
                for outer_key, inner_dict in LedgerInfo.items()
                for inner_key, value in inner_dict.items()
            }

            ledger2[comp.CompID] = ledger.index.map(flat_dict).copy()
        ledger = ledger2.sort_index()
        with self.lock:
            prevLedger = self.ReadFullLedger()
            newLedger = prevLedger.join(ledger, how="left", validate="1:1")
            newLedger.to_csv(self.LedgerFile)

        file.close()
        return

    def RestartLedger(self) -> list[Compound]:
        """Reinitializes the list of Compounds from the UsedInput backup.
        This can be used for recontinuing the workflow.

        Returns:
            list[Compound]: list of all Compound objects related to this ProcessLedger
        """
        loc = initialize_compounds(self.inputbackup)
        return loc

    def ReadFullLedger(self) -> pd.DataFrame:
        """Reads all information from the .csv and retuns it as a DataFrame

        Returns:
            pd.DataFrame: Info contained in the ProcessLedger csv
        """
        with self.lock:
            self.LedgerDF = pd.read_csv(
                self.LedgerFile, index_col=[0, 1], header=0, dtype=str
            ).sort_index()
        return self.LedgerDF

    def GetSingleJob(self, comp: Compound, JobName: str):
        """Returns a Dict of the information stored in the process ledger of the step specified by JobName
        for a given Compound

        Args:
            comp (Compound): which composition to get the Job for. Must have CompID assigned
            JobName (str): Which Job to access the information for. Must match one of the keys in the JobInfo json

        Returns:
            dict: a dictionary of all information stored in the ProcessLedger for a specific step for a specific Compound
        """
        comp_col = self.ReadFullLedger()  # , dtype={"IncarExtras":dict})
        jobdict = comp_col.loc[JobName][str(comp.CompID)].to_dict()
        try:
            jobdict["IncarExtras"] = ast.literal_eval(jobdict["IncarExtras"])
        except:
            jobdict["IncarExtras"] = {}

        return jobdict

    def GetCompletionOverview(self) -> pd.DataFrame | pd.Series:
        """Reads the ProcessLedger and returns only the info on the completion of jobs.
        It will return a DataFrame with compositions on the index, job steps as columns.
        Data will contain:
        - 0: Job is incomplete
        - 1: Job is complete
        - -1: Job had some error

        Returns:
            pd.DataFrame | pd.Series: Oveview of all jobs with their completion status
        """
        ledger = self.ReadFullLedger()
        summary = ledger.xs("completed", level=1).T
        return summary.astype(str)

    def GetQueue(self) -> dict:
        """Creates a dictionary with for each job step, which compositions are pending to do this job.
        Used to assign jobs

        Returns:
            dict: keys are job steps, values list of compositions that are waiting for that job step.
        """
        overview = self.GetCompletionOverview()
        AllComps = overview.index

        DoneComps = []  # used to keep track which compounds have already been treated

        stepnames = self.JobInfo.keys()
        QueueDict = {}

        ##First remove any Errors (so if completed is anything else than 1 | 0 )
        errors = overview[
            ((overview != "1") & (overview != "0")).any(axis=1)
        ].index.to_list()

        if errors:
            for error in errors:
                DoneComps.append(error)

        # DoneComps.append(AllComps.to_list()[0])

        for JobName in stepnames:
            stepcomps = []
            # Remove any Compounds contained in DoneComps

            newindex = AllComps.difference(pd.Index(DoneComps), sort=False)
            overview = overview.loc[newindex]

            # Now iterate over steps and find the comps that still need to complete that step
            stepcomps = overview[overview[JobName] == "0"].index.to_list()
            # Add these comps to the queue dict and append them to DoneComps
            QueueDict[JobName] = stepcomps
            if stepcomps:
                for comp in stepcomps:
                    DoneComps.append(comp)

        return QueueDict

    def SetSingleValue(self, comp: Compound, JobName: str, Field: str, value) -> None:
        """Set a singular value in the ProcessLedger csv

        Args:
            comp (Compound): Which compositions to set a value to. Must have CompID assigned.
            JobName (str): Which Job to set the value for. Must match one of the keys in the JobInfo json
            Field (str): Which field in the Ledger to set.
            value (_type_): The value to set.
        """
        assert (
            Field in self.info_per_job
        ), f"Trying to change a non-existent field in Ledger: {Field}"

        try:
            with self.lock.acquire(timeout=10):
                ledger = self.ReadFullLedger()
                # print(f"ledger.loc[({JobName}, {Field}),str({comp.CompID})] = {value}")
                ledger.loc[JobName, Field][str(comp.CompID)] = value
                ledger.to_csv(self.LedgerFile)
        except Timeout:
            print("Timeout occurred while trying to acquire the file lock.")
            raise
        return

    def RestoreLedgerData(self):
        """This function restores the whole Ledger from the database in case of some critical error.
        It will go to the base directory, loop through all the directories and determine if DFT SCF cylcles were completed successfully.
        It does require the InputBackup to still exist.
        """
        import glob
        from job_handler import check_scfcompletion

        os.chdir(self.basepath)

        loc = self.RestartLedger()

        index = pd.MultiIndex.from_product(
            [list(self.JobInfo.keys()), self.info_per_job]
        )

        backup = pd.read_csv(
            os.path.join(self.basepath, "AllCompsCopy.csv"), index_col=[0, 1], header=0
        ).sort_index()

        ledger = pd.DataFrame(index=index, columns=backup.columns, dtype=str)
        with self.lock:
            ledger.to_csv(self.LedgerFile)
        ledger2 = pd.DataFrame(index=index, columns=backup.columns, dtype=str)

        for comp in loc:
            print(comp)
            comp.path = str(self.basepath) + "/" + str(comp.CompID)

            LedgerInfo = {}
            jnames = list(self.JobInfo.keys())
            # jnames.reverse()

            for JobName in jnames:
                os.chdir(comp.path)
                settings = self.JobInfo[JobName]

                completion = 0

                # If it is a VASP or misc job a subdirectory will be made
                # if the subdir already exists, it will be assumed that the step was
                # completed succesfully
                # For LOBSTER jobs it will set the subdir from another step as working dir
                if settings["job_type"] != "LOBSTER":
                    calc_subdir = str(JobName) + "_" + settings["name"]
                    calcdir_full = os.path.join(comp.path, calc_subdir)

                    # Here we need to figure out whether the job was completed or not
                    if JobName == "1Rel":
                        if os.path.isfile(f"{comp.CompID}_relaxedPOSCAR"):
                            completion = 1

                    elif JobName == "2Spins":
                        res_file = glob.glob("SpinResults*")
                        if res_file:
                            completion = 1

                            spinDF = pd.read_csv(res_file[-1], index_col=0)

                            try:
                                converged_runs = spinDF.groupby(
                                    "SCF_convergence"
                                ).get_group(0)
                            except KeyError:
                                completion = 0

                            else:

                                optimal_run = converged_runs["Etot_out"].idxmin()
                                magmom_use = f"{spinDF['muB1_in'].loc[optimal_run]} {spinDF['muB2_in'].loc[optimal_run]} 2*0.0 6*0.0"

                                ##now we want to put this MAGMOM initialization into the IncarExtras for all remaining VASP jobs
                                jobindex = self.Jobs.index(JobName)

                                for futurejob in self.Jobs[jobindex + 1 :]:
                                    if (
                                        self.JobInfo[futurejob]["job_type"].upper()
                                        == "VASP"
                                    ):
                                        IncarExtras = self.GetSingleJob(comp, JobName)[
                                            "IncarExtras"
                                        ]

                                        IncarExtras.update({"MAGMOM": magmom_use})
                                        self.SetSingleValue(
                                            comp,
                                            futurejob,
                                            "IncarExtras",
                                            str(IncarExtras),
                                        )

                    elif JobName == "3Pre" or JobName == "4HSE":
                        os.chdir(calcdir_full)

                        if os.path.isfile("OUTCAR"):
                            loglist = glob.glob("log*")
                            if loglist:
                                output = check_scfcompletion(loglist[-1])
                                if output[-1] == 0:
                                    completion = 1
                                else:
                                    completion = 0
                            else:
                                completion = 0

                else:
                    assert (
                        "use_output_step" in settings.keys()
                    ), "LOBSTER job is missing required VASP directory.\n Please define 'use_output_step' in JobSettingsJson"
                    subdir = (
                        settings["use_output_step"]
                        + "_"
                        + self.JobInfo[settings["use_output_step"]]["name"]
                    )
                    calcdir_full = os.path.join(comp.path, subdir)

                # store all base information for each step
                LedgerInfo.update(
                    {
                        JobName: {
                            "JobPath": str(calcdir_full),
                            "JobID": "",
                            "completed": completion,
                            "TimeStamp": "",
                            "IncarExtras": self.GetSingleJob(comp, JobName)[
                                "IncarExtras"
                            ],
                        }
                    }
                )

            # flatten the dict of ledger info and add as new column to dataframe
            flat_dict = {
                (outer_key, inner_key): value
                for outer_key, inner_dict in LedgerInfo.items()
                for inner_key, value in inner_dict.items()
            }

            ledger2[comp.CompID] = ledger.index.map(flat_dict).copy()
        self.LedgerDF = ledger2.sort_index()
        with self.lock:
            ledger2.to_csv(self.LedgerFile)
        return

