"""_
JobDefinitions.py
Contains classes with methods for several types of VASP jobs and for a LOBSTER job.
"""

from HDPCompound import Compound
from HTProcessLedger import ProcessLedger
from job_handler import submit_and_wait, check_scfcompletion
import os
import time
import subprocess
import json

class GeneralJob:
    """Base class for workflow job execution.

    GeneralJob wraps a single workflow step for a specific Compound and ledger entry.
    It provides common behavior for job initialization, INCAR/lobsterin writing, and
    SLURM submission file generation. Subclasses implement job-specific run logic.

    Attributes:
        AssignedServer (str): SLURM partition or host assignment string.
        AssignedCompName (str): Compound identifier used in logging and paths.
        AssignedComp (Compound): Compound object for this job.
        ledger (ProcessLedger): Ledger instance tracking workflow state.
        JobName (str): Name of the current job step as defined in the job settings.
        JobBasics (dict): Standard job settings loaded from the ledger JSON file.
        JobID (int): Unique numeric identifier for this run instance.
        JobSpecifics (dict): Stored per-job override settings from the ledger.
    """

    def __init__(
        self,
        ledger: ProcessLedger,
        comp: Compound,
        JobName: str = "",
        AssignedServer: str = "ccp20,ccp22",
    ):
        assert (
            JobName in ledger.Jobs
        ), "JobName passed to (General)Job CLASS not present in the Ledger"

        self.AssignedServer = AssignedServer
        self.AssignedCompName = comp.CompID
        self.AssignedComp = comp
        self.ledger = ledger
        # retrieve the basic settings of this job from the job setting JSON
        self.JobName = JobName
        self.JobBasics = ledger.JobInfo[JobName]

        self.JobID = int(float(comp.CompID.split("_")[0]) * time.time())
        ledger.SetSingleValue(comp, JobName, Field="JobID", value=self.JobID)
        ledger.SetSingleValue(
            comp, JobName, Field="TimeStamp", value=time.strftime("%Y/%m/%d_%H:%M:%S")
        )

        # Also regtrieve the extra information stored in the ledger
        self.JobSpecifics = dict(ledger.GetSingleJob(comp, JobName))

    def __repr__(self):
        return self.AssignedCompName + "_" + self.JobName

    def WriteIncar(self):
        """
        This function will take the standard incar from JobSettings JSON file, and any extra inputs in the ProcessLedger
        put it together and write it to a file. It will automatically save to 'INCAR' for 'job_type:VASP' and lobsterin for 'job_type:LOBSTER'

        for any other job_type it will take ['name'] from JobSettings JSON and save the inputs to '{name}_INPUT'
        """
        incar = self.JobBasics["std_incar"]
        incar.update(self.JobSpecifics["IncarExtras"])

        if self.JobBasics["job_type"].upper() == "VASP":
            filename = "INCAR"
            joiner = " = "
        elif self.JobBasics["job_type"].upper() == "LOBSTER":
            filename = "lobsterin"
            joiner = " "
        else:
            filename = f"{self.JobBasics['name']}_INPUT"
            joiner = " "

        with open(filename, "w") as f:
            for key, val in incar.items():

                if type(val) == list:
                    for singleval in val:
                        newline = joiner.join([str(key), str(singleval)]) + "\n"
                        f.write(newline)
                else:
                    newline = joiner.join([str(key), str(val)]) + "\n"
                    f.write(newline)

        return incar

    def WriteSubmission(self, LogFileExtra=""):
        """Write the SLURM submission script for this job.

        The generated submission file is named based on the job type (e.g. "vasp.sub" or "lobster.sub").
        The method also records the standard output and error filenames so the job can be monitored later.

        Args:
            LogFileExtra (str, optional): Suffix to append to log and error filenames.
        """
        # Set filenames for submission, log, and error files
        filename = self.JobBasics["job_type"].lower() + ".sub"
        self.SubFile = filename
        self.LogFile = f"log_{self.JobID}_{self.AssignedCompName}{LogFileExtra}"
        self.ErrorFile = f"error_{self.JobID}_{self.AssignedCompName}{LogFileExtra}"
        
        # Extract nnodes and potentially time from JobSettings JSON
        nnodes = self.JobBasics["nnodes"]
        if "max_time" in self.JobBasics.keys():
            max_time = self.JobBasics["max_time"]
        else:
            max_time = f"0{nnodes}:00:00"

        # Set modules and slurm commands
        if self.JobBasics["job_type"].upper() == "VASP":
            module = "module load vasp/6.3.0"
            runcommand = "mpirun vasp_std"
        elif self.JobBasics["job_type"].upper() == "LOBSTER":
            module = "module load lobster"
            runcommand = "lobster"
        else:
            raise ModuleNotFoundError(
                f"Could not write submission script for unknown job_type: {self.JobBasics['job_type']}"
            )

        # Write the submission file line by line :-0
        with open(filename, "w") as file:
            file.write(f"#!/bin/bash \n#SBATCH --time={max_time}\n")
            file.write(
                f"#SBATCH -p {self.AssignedServer} \n#SBATCH -N {nnodes} --exclusive --ntasks-per-node 32\n"
            )

            file.write(
                f"#SBATCH -J {self.AssignedCompName}{self.JobName} \n#SBATCH --output={self.LogFile}\n"
            )
            file.write(
                f"#SBATCH --ntasks-per-core=1 \n#SBATCH --error={self.ErrorFile}\n\n"
            )
            file.write(f"module purge \n{module} \n\nulimit -s unlimited\n")
            file.write("startingtime=$(date)\n")
            file.write(
                f"now=$(date +'%s') \n{runcommand} \nendtime=$(date +'%s') \nduration=$(($endtime-$now)) \n"
            )
            file.write(
                "echo 'Started on: ' $startingtime 'calculation ran for: ' $duration '[s]'\n"
            )

        return

    def Run(self):
        os.chdir(self.JobSpecifics["JobPath"])

        self.WriteIncar()
        self.WriteSubmission()

        slurm_exit = submit_and_wait(self.SubFile)

        if slurm_exit == "COMPLETED":
            self.ledger.SetSingleValue(
                self.AssignedComp, self.JobName, Field="completed", value=1
            )
            return

        else:
            # some weird error, save slurm exit code to the completion field in the ledger
            self.ledger.SetSingleValue(
                self.AssignedComp,
                self.JobName,
                Field="completed",
                value=f"slurm:{slurm_exit}",
            )
            return


class SimpleVasp(GeneralJob):
    """Base class for VASP-based job steps.

    SimpleVasp extends GeneralJob with VASP-specific helpers for algorithm cycling,
    SCF completion checks, input preparation, and result extraction.
    """

    def __init__(
        self,
        ledger: ProcessLedger,
        comp: Compound,
        JobName: str = "",
        AssignedServer: str = "ccp20,ccp22",
    ):
        super().__init__(ledger, comp, JobName, AssignedServer)

        # extra attribute specifically for vasp jobs: which algos are allowed
        self.AllowedAlgos = ["N", "F", "D", "VF", "A", "S"]

    def CycleAlgos(self):
        """This function will update the chosen algorithm for a VASP calculation based on the available algorithms
        in self.AllowedAlgos.

        If the "IncarExtras" in the ProcessLedger does not have an ALGO specified for the current step it will check the previous step for that keyword.
            If it is specified it there it will update the "IncarExtras" for the current step with that ALGO
            Else it will start with AllowedAlgos[1]

        If the "IncarExtras" of the current step does have an ALGO specified, it will update it in the ProcessLedger for the next algo specified in the list
        """

        if "ALGO" in self.JobSpecifics["IncarExtras"].keys():
            prev_algo = self.JobSpecifics["IncarExtras"]["ALGO"]

            try:
                self.JobSpecifics["IncarExtras"]["ALGO"] = self.AllowedAlgos[
                    self.AllowedAlgos.index(prev_algo) + 1
                ]

            except ValueError:
                self.JobSpecifics["IncarExtras"]["ALGO"] = self.AllowedAlgos[1]
            except IndexError:
                self.ledger.SetSingleValue(
                    self.AssignedComp,
                    self.JobName,
                    "completed",
                    "AlgoError OOA: OutOfAlgos",
                )
                raise IndexError(
                    f"CycleAlgos wanted to change ALGO for Compound:{self.AssignedCompName} but went beyond the max(?) index of Job.Allowed.Algos list"
                )

        else:
            try:
                prev_job_settings = self.ledger.GetSingleJob(
                    self.AssignedComp,
                    self.ledger.Jobs[self.ledger.Jobs.index(self.JobName) - 1],
                )
            except IndexError:
                self.JobSpecifics["IncarExtras"]["ALGO"] = self.AllowedAlgos[1]
            else:
                try:
                    if prev_job_settings["IncarExtras"]["ALGO"] in self.AllowedAlgos:
                        self.JobSpecifics["IncarExtras"]["ALGO"] = prev_job_settings[
                            "IncarExtras"
                        ]["ALGO"]
                    else:
                        self.JobSpecifics["IncarExtras"]["ALGO"] = self.AllowedAlgos[1]
                except KeyError:
                    self.JobSpecifics["IncarExtras"]["ALGO"] = self.AllowedAlgos[1]

        self.ledger.SetSingleValue(
            self.AssignedComp,
            self.JobName,
            "IncarExtras",
            self.JobSpecifics["IncarExtras"],
        )

        return

    def check_scfcompletion(self, erroroutput="encounteredErrors.tmp"):
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
                self.LogFile,
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

    def PrepareVASPRun(self):
        """This function simply makes sure all required input files for VASP (except for INCAR!) are present
        in the Job subfolder. I.E. It will write KPOINTS, POTCAR and POSCAR if they are not present yet.

        For POSCAR it will check if the JobSettings.json has the keyword 'use_poscar_from',
            if it is given, it will take CONTCAR from the specified JobName and copy it to the current job as POSCAR.
            else it will use the write_poscar() method from the class Compound to write the initial poscar

        If the JobSettings.json also has the keyword 'extra_input_files' it will link the specified files from the specified job
        """
        # os.chdir(self.JobSpecifics["JobPath"])

        if not os.path.exists("KPOINTS"):
            self.AssignedComp.write_kpoints()

        self.AssignedComp.write_potcar()

        if "extra_input_files" in self.JobBasics.keys():
            for job, files in self.JobBasics["extra_input_files"].items():
                sourcedir = self.ledger.GetSingleJob(self.AssignedComp, job)["JobPath"]

                if type(files) == list:
                    for file in files:
                        sourcefile = os.path.join(sourcedir, file)
                        if file == "CONTCAR":

                            os.system(f"cp -f {sourcefile} ./POSCAR")
                        else:

                            if not os.path.exists(file):
                                os.system(f"ln -f {sourcefile} .")

                else:
                    sourcefile = os.path.join(sourcedir, files)
                    if files == "CONTCAR":

                        os.system(f"cp -f {sourcefile} ./POSCAR")
                    else:

                        if not os.path.exists(files):
                            os.system(f"ln -f {sourcefile} .")

        if not os.path.exists("POSCAR"):
            self.AssignedComp.write_poscar()

        return

    def ExtractEtot(self, file="OSZICAR"):
        """Extract the total energy from an OSZICAR file.

        Args:
            file (str, optional): OSZICAR-like file to parse. Defaults to "OSZICAR".

        Returns:
            float | list[float]: Total energy value(s) extracted from the file.
        """
        # os.chdir(self.JobSpecifics['JobPath'])
        p = subprocess.run(
            ["awk", "/F=/ {print $5}", file], stdout=subprocess.PIPE, text=True
        )
        output = p.stdout.strip().split("\n")
        if len(output) > 1:
            Eout = list(map(float, output))
        else:
            Eout = float(output[0])

        return Eout

    def ExtractMAGtot(self, file="OSZICAR"):
        """Extract the total magnetization from an OSZICAR file.

        Args:
            file (str, optional): OSZICAR-like file to parse. Defaults to "OSZICAR".

        Returns:
            float | list[float]: Total magnetic moment value(s) extracted from the file.
        """
        # os.chdir(self.JobSpecifics['JobPath'])
        p = subprocess.run(
            ["awk", "/F=/ {print $NF}", file], stdout=subprocess.PIPE, text=True
        )
        output = p.stdout.strip().split("\n")
        if len(output) > 1:
            MAGMOM = list(map(float, output))
        else:
            MAGMOM = float(output[0])

        return MAGMOM

    def UsePOTCARENMAX(self, MultFactor: float = 1.2, MinValue: int = 350):
        """Set ENCUT in the INCAR extras based on POTCAR ENMAX.

        Reads the maximum ENMAX value from the POTCAR and scales it by MultFactor.
        A minimum ENCUT value is enforced to avoid too-low cutoffs.

        Args:
            MultFactor (float, optional): Multiplicative scaling factor for ENMAX. Defaults to 1.2.
            MinValue (int, optional): Minimum allowed ENCUT. Defaults to 350.

        Returns:
            float: The selected ENCUT value.
        """

        # os.chdir(self.JobSpecifics["JobPath"])
        get_enmax = subprocess.run(
            ["awk", "BEGIN{a = 0} /ENMAX/ {if($3>0+a) a=$3 fi} END{print a}", "POTCAR"],
            stdout=subprocess.PIPE,
            text=True,
        )
        enout = float(get_enmax.stdout.strip().strip(";"))

        ENCUT = MultFactor * enout

        if ENCUT < MinValue:
            ENCUT = MinValue

        self.JobSpecifics["IncarExtras"].update({"ENCUT": ENCUT})

        return ENCUT

    def Run(self):
        """Run a VASP job, check SCF convergence, and update ledger state.

        This method writes the INCAR file, submits the VASP job, and uses the SCF
        completion helper to determine whether the run succeeded, failed, or needs retries.

        Returns:
            str: Status indicator such as "succesfull", "convergence_failed", "error", or "run_failed".
        """
        os.chdir(self.JobSpecifics["JobPath"])

        self.WriteIncar()
        self.WriteSubmission()

        slurm_out = submit_and_wait(self.SubFile)
        self.ledger.SetSingleValue(
            self.AssignedComp,
            self.JobName,
            "IncarExtras",
            self.JobSpecifics["IncarExtras"],
        )
        if slurm_out != "COMPLETED":
            self.ledger.SetSingleValue(
                self.AssignedComp,
                self.JobName,
                Field="completed",
                value=f"slurm:{slurm_out}",
            )
            return "slurm"

        scf_out = self.check_scfcompletion()

        if scf_out[-1] == 0:
            # job succesfull, update this in ledger and return
            self.ledger.SetSingleValue(
                self.AssignedComp, self.JobName, Field="completed", value=1
            )
            return "succesfull"

        elif scf_out[-1] == 1:
            # convergence failed, so make useful changes
            self.CycleAlgos()
            return "convergence_failed"

        elif scf_out[-1] == -1:
            # error occured, deal with it appropriate for this step
            self.ledger.SetSingleValue(
                self.AssignedComp, self.JobName, Field="completed", value=-1
            )
            return "error"
        else:
            # something else went wrong in converging
            # convergence failed, so make useful changes
            self.CycleAlgos()
            return "run_failed"


class Relaxation(SimpleVasp):
    """VASP relaxation job step with adaptive correction logic.

    Relaxation extends SimpleVasp by using frozen-core POTCARs, updating ENCUT,
    checking for structural changes, and handling restart conditions.
    """

    def __init__(
        self,
        ledger: ProcessLedger,
        comp: Compound,
        JobName: str = "",
        AssignedServer: str = "ccp20,ccp22",
    ):
        super().__init__(ledger, comp, JobName, AssignedServer)
        self.AllowedAlgos = ["N", "A", "D", "F", "VF", "S"]

    def check_changepos(self, outcar="OUTCAR"):
        """This function is specifically for the relaxation step!
        It will take all position updates out of OUTCAR and compare the first set with the last set op coordinates
        If the relaxation reverts back to the starting position and the end and beginning position are the same
        this function will return False, and adaptations need to be made to get the correct structure.

        Args:
            outcar (str, optional): A string containing the filename of the OUTCAR file. Defaults to "OUTCAR".

        Returns:
            Bool: a boolean indicating whether the final structure is different from the initial positions
                True means the structure has changed, indicating succesfull relaxation.
        """
        p1 = subprocess.run(
            ["awk", '/TOTAL-FORCE/ {getline;getline; for (i=0; i<=9; i++)\
            print $1"   "$2"   "$3 getline; print"\\n"}', outcar],
            text=True,
            stdout=subprocess.PIPE,
        )

        all_pos = p1.stdout.strip().split("\n")
        i_empties = [i for i, x in enumerate(all_pos) if x == ""]
        try:
            begin_pos = all_pos[0 : i_empties[0]]
        except IndexError:
            if len(i_empties) == 0:
                # only one set of positions was found so no success
                return False
            else:
                return False
        except Exception as exc:
            # something else went wrong
            print("unexpected exception in check_changepos: %s" % (exc))
            return False
        else:

            end_pos = all_pos[i_empties[-1] :]
            end_pos.pop(0)
            return end_pos != begin_pos

    def check_zbrent_message(self):
        """This function checks the OUTCAR of the current job for the literal string
        'rerun with smaller EDIFF', and returns a boolean indicating whether it occured that specific
        VASP error message in the OUTCAR. This can be used to

        Args:
            outcar (str, optional): String that needs to point to a VASP OUTCAR file. Defaults to "OUTCAR".

        Returns:
            bool: True if the error in ZBRENT restart with smaller EDIFF and restart from contcar
                    False if this error message is not encountered
        """
        pzbrent = subprocess.run(
            [
                "awk",
                "BEGIN{enc=0}; /rerun with smaller EDIFF/ {enc++}; END{print enc}",
                self.LogFile,
            ],
            stdout=subprocess.PIPE,
            text=True,
        )
        zbrent_encounters = int(pzbrent.stdout.strip())
        return bool(zbrent_encounters)

    def Run(self):
        """Run the relaxation job and apply adaptive recovery strategies.

        This method submits the relaxation calculation, checks convergence,
        saves the relaxed CONTCAR on success, and updates INCAR settings if the
        relaxation fails or requires restart handling.

        Returns:
            str: Run status such as "relaxed", "slurm", or "relaxation_failed".
        """

        os.chdir(self.JobSpecifics["JobPath"])

        ##Preparation steps
        self.PrepareVASPRun()  # assures KPOINTS, POSCAR, and POTCAR are there
        self.AssignedComp.write_potcar(use_frozen=True)
        self.UsePOTCARENMAX(
            MultFactor=1.5, MinValue=self.JobBasics["std_incar"]["ENCUT"]
        )
        incar = self.WriteIncar()

        self.WriteSubmission()
        slurm_out = submit_and_wait(self.SubFile)
        if slurm_out != "COMPLETED":
            self.ledger.SetSingleValue(
                self.AssignedComp,
                self.JobName,
                Field="completed",
                value=f"slurm:{slurm_out}",
            )
            return "slurm"

        scf_out = self.check_scfcompletion()

        if scf_out[-1] == 2:
            # this is a special flag recognizing VASP successful relaxation message
            self.ledger.SetSingleValue(
                self.AssignedComp, self.JobName, Field="completed", value=1
            )
            os.system(f"cp CONTCAR ../{self.AssignedCompName}_relaxedPOSCAR")
            return "relaxed"

        poschange = self.check_changepos()
        contcarrestart = self.check_zbrent_message()

        if not poschange:
            # Relaxation returned to initial position, this can be solved by initializing with larger a_lat
            self.AssignedComp.write_poscar(a=5.7)

        if -1 in scf_out:
            # check for errors in the SCF cycles
            # many errors are solved by using ALGO = All and high PREC
            self.CycleAlgos()
            self.JobSpecifics["IncarExtras"].update({"PREC": "A", "NELM": 150})
            self.ledger.SetSingleValue(
                self.AssignedComp,
                self.JobName,
                "IncarExtras",
                self.JobSpecifics["IncarExtras"],
            )

        if contcarrestart:
            # special flag from VASP telling you to restart from CONTCAR with smaller EDIFF
            os.system("cp CONTCAR POSCAR")
            newEDIFF = incar["EDIFF"] / 10
            self.JobSpecifics["IncarExtras"].update({"EDIFF": newEDIFF})
            self.ledger.SetSingleValue(
                self.AssignedComp,
                self.JobName,
                "IncarExtras",
                self.JobSpecifics["IncarExtras"],
            )
        elif np.sum(scf_out) / len(scf_out) > 0.3:
            # more than 30% of SCF cycles did not converge, maybe cyclin algos will help
            self.CycleAlgos()

        return f"relaxation_failed, Poschange: {poschange}, Restart CONTCAR: {contcarrestart}"


class SpinStateScan(SimpleVasp):
    """VASP job step that scans magnetic spin initializations.

    SpinStateScan generates multiple MAGMOM initializations for B-site spin states,
    runs each combination, and records the lowest-energy converged result.
    """

    def __init__(
        self,
        ledger: ProcessLedger,
        comp: Compound,
        JobName: str = "",
        AssignedServer: str = "ccp20,ccp22",
    ):
        super().__init__(ledger, comp, JobName, AssignedServer)
        self.AllowedAlgos = ["N", "F", "VF", "A", "S"]

    def GetSpinCombos(self):
        """Build a list of spin moment combinations to test for the B-site cations.

        Returns:
            list[list[float | str]]: Requested MAGMOM initializations for B1 and B2.
        """
        b1HS, b1LS = self.AssignedComp.get_spinstate("B1")
        if b1HS == b1LS:
            B1s = [b1HS]

        else:
            B1s = [b1HS, b1LS]

        if self.AssignedComp.B2 == "Vac":
            B2s = [""]
        else:

            b2HS, b2LS = self.AssignedComp.get_spinstate("B2")
            ##Filter out possible duplicate spinstates
            if b2HS == b2LS:
                B2s = [b2HS]

            else:
                B2s = [b2HS, b2LS]

        combos = []
        for b1 in B1s:
            for b2 in B2s:
                combos.append([b1, b2])
                if b2 != 0 and b2 != "":
                    combos.append([b1, -b2])

        return combos

    def Run(self):
        """Run multiple spin-initialized VASP jobs and select the best converged result.

        Returns:
            str: "good_run" if at least one converged spin scan succeeds, "convergence_failed" if not,
                 or "slurm" if the submission itself fails.
        """

        os.chdir(self.JobSpecifics["JobPath"])

        self.PrepareVASPRun()
        self.UsePOTCARENMAX(
            MultFactor=1.0, MinValue=self.JobBasics["std_incar"]["ENCUT"]
        )

        spin_combos = self.GetSpinCombos()

        spinDF = pd.DataFrame(
            columns=["muB1_in", "muB2_in", "SCF_convergence", "Etot_out", "MAGtot_out"],
            dtype=str,
        )
        subrunIndex = 0
        for combo in spin_combos:
            subrunIndex += 1
            ##Set the required MAGMOM initialization
            MAGMOM_line = f"{combo[0]} {combo[1]} 2*0.0 6*0.0"
            self.JobSpecifics["IncarExtras"].update({"MAGMOM": MAGMOM_line})
            # self.ledger.SetSingleValue(self.AssignedComp,self.JobName,'IncarExtras', self.JobSpecifics['IncarExtras'])

            # write INCAR
            self.WriteIncar()
            self.WriteSubmission(LogFileExtra=f"_MAGMOM_{combo[0]}_{combo[1]}")

            slurm_out = submit_and_wait(self.SubFile)

            if slurm_out != "COMPLETED":
                # some weird error, save slurm exit code to the completion field in the ledger
                self.ledger.SetSingleValue(
                    self.AssignedComp,
                    self.JobName,
                    Field="completed",
                    value=f"slurm:{slurm_out}",
                )
                return "slurm"

            scf_out = self.check_scfcompletion()
            if scf_out == "" or scf_out == [""] or scf_out == [-1]:
                scf_out = ["1"]
                Eout = 100.0
                MagOut = ""
            else:
                Eout = self.ExtractEtot()
                MagOut = self.ExtractMAGtot()

            outputDict = {
                "muB1_in": combo[0],
                "muB2_in": combo[1],
                "SCF_convergence": scf_out[-1],
                "Etot_out": Eout,
                "MAGtot_out": MagOut,
            }

            spinDF.loc[subrunIndex] = spinDF.columns.map(outputDict)

        # safe the results of this analysis to a csv file
        spinDF.to_csv(f"../SpinResults_{self.AssignedCompName}_job{self.JobID}.csv")

        try:
            converged_runs = spinDF.groupby("SCF_convergence").get_group(0)
        except KeyError:
            # not a single run was converged, so we need to retry with a different algo
            self.CycleAlgos()
            return "convergence_failed"

        else:
            self.ledger.SetSingleValue(self.AssignedComp, self.JobName, "completed", 1)

        spinDF.fillna("", inplace=True)
        optimal_run = converged_runs["Etot_out"].idxmin()
        magmom_use = f"{spinDF['muB1_in'].loc[optimal_run]} {spinDF['muB2_in'].loc[optimal_run]} 2*0.0 6*0.0"

        ##now we want to put this MAGMOM initialization into the IncarExtras for all remaining VASP jobs
        jobindex = self.ledger.Jobs.index(self.JobName)

        for futurejob in self.ledger.Jobs[jobindex + 1 :]:
            if self.ledger.JobInfo[futurejob]["job_type"].upper() == "VASP":
                IncarExtras = self.ledger.GetSingleJob(self.AssignedComp, futurejob)[
                    "IncarExtras"
                ]

                IncarExtras.update({"MAGMOM": magmom_use})
                self.ledger.SetSingleValue(
                    self.AssignedComp, futurejob, "IncarExtras", IncarExtras
                )

        return "good_run"


class PreLobster(SimpleVasp):
    """Prepares and executes the final VASP run before LOBSTER projection.

    PreLobster ensures the VASP calculation uses enough bands for the planned
    LOBSTER projection and updates INCAR settings accordingly.
    """

    def __init__(
        self,
        ledger: ProcessLedger,
        comp: Compound,
        JobName: str = "",
        AssignedServer: str = "ccp20,ccp22",
    ):
        super().__init__(ledger, comp, JobName, AssignedServer)
        if JobName == "4HSE":
            self.AllowedAlgos = ["D", "A", "N", "S"]
        else:
            self.AllowedAlgos = ["N", "VF", "F", "D", "A", "S"]

    def Run(self):
        """Run the final pre-LOBSTER VASP step and ensure NBANDS is set correctly.

        Returns:
            str: "good_run", "convergence_failed", "error", or "run_failed".
        """
        os.chdir(self.JobSpecifics["JobPath"])
        ##Preparation steps
        self.PrepareVASPRun()  # assures KPOINTS, POSCAR, and POTCAR are there
        self.UsePOTCARENMAX(
            MultFactor=1.0, MinValue=self.JobBasics["std_incar"]["ENCUT"]
        )

        NBANDS = self.AssignedComp.get_maxNBANDS()
        self.JobSpecifics["IncarExtras"].update({"NBANDS": NBANDS})
        self.ledger.SetSingleValue(
            self.AssignedComp,
            self.JobName,
            "IncarExtras",
            self.JobSpecifics["IncarExtras"],
        )

        self.WriteIncar()
        self.WriteSubmission()

        slurm_out = submit_and_wait(self.SubFile)
        if slurm_out != "COMPLETED":
            self.ledger.SetSingleValue(
                self.AssignedComp,
                self.JobName,
                Field="completed",
                value=f"slurm:{slurm_out}",
            )
            return "slurm"

        scf_out = self.check_scfcompletion()

        if scf_out[-1] == 0:
            # job succesfull, update this in ledger and return
            self.ledger.SetSingleValue(
                self.AssignedComp, self.JobName, Field="completed", value=1
            )
            return "good_run"

        elif scf_out[-1] == 1:
            # convergence failed, so make useful changes
            self.CycleAlgos()
            return "convergence_failed"

        elif scf_out[-1] == -1:
            # error occured, deal with it appropriate for this step
            self.CycleAlgos()
            return "error"

        else:
            # something else went wrong in converging
            # convergence failed, so make useful changes
            return "run_failed"


class SimpleLobster(GeneralJob):
    """LOBSTER projection job step for chemical bonding analysis.

    SimpleLobster prepares lobsterin settings, manages basis-function combinations,
    and runs LOBSTER in separate directories while checking spillings and overlaps.
    """

    def __init__(
        self,
        ledger: ProcessLedger,
        comp: Compound,
        JobName: str = "",
        AssignedServer: str = "ccp20,ccp22",
    ):
        super().__init__(ledger, comp, JobName, AssignedServer)

    def ChangeCOHPErange(self, use_file="DOSCAR"):
        """Adjust LOBSTER energy window settings based on the DOSCAR energy range.
        Also adds COHPgenerator keyword twice for each B-X combination to the IncarExtras

        Args:
            use_file (str, optional): DOSCAR-like file to parse. Defaults to "DOSCAR".

        Returns:
            tuple[float, float]: COHP start and end energies.
        """


        with open(use_file, "r") as file:
            lines = file.readlines()

        energyline = lines[5]
        values = list(map(float, energyline.strip().split()))

        self.JobSpecifics["IncarExtras"].update(
            {
                "COHPstartEnergy": values[1] - values[3],
                "COHPendEnergy": values[0] - values[3],
            }
        )

        if self.AssignedComp.B2 == "Vac":
            self.JobSpecifics["IncarExtras"].update(
                {
                    "cohpGenerator from 0.5 to 4.": f"type {self.AssignedComp.B1} type {self.AssignedComp.X} orbitalwise "
                }
            )
        else:
            self.JobSpecifics["IncarExtras"].update(
                {
                    "cohpGenerator from 0.5 to 4.": [
                        f"type {self.AssignedComp.B1} type {self.AssignedComp.X} orbitalwise",
                        f"type {self.AssignedComp.B2} type {self.AssignedComp.X} orbitalwise ",
                    ]
                }
            )

        return values[1], values[0]

    def add_basisfunctions(self, which="basis0"):
        """Append basis function definitions for a named LOBSTER basis set to lobsterin.

        Args:
            which (str, optional): Basis combination name from get_all_basis_combos().
                Defaults to "basis0".
        """
        # os.chdir(self.JobSpecifics['JobPath'])
        all_basis = self.AssignedComp.get_all_basis_combos()
        with open("lobsterin", "a") as file:
            for el, funcs in all_basis[which].items():
                file.write(f"basisFunctions {str(el)} {str(funcs)}\n")

        return

    def check_charge_spilling(self, max_val=3):
        """Check whether LOBSTER charge spilling exceeds a threshold.

        Args:
            max_val (int, optional): Maximum allowed absolute charge spilling percent.
                Defaults to 3.

        Returns:
            bool: True if charge spilling is above the threshold.
        """

        awk_spil = subprocess.run(
            ["awk", "/abs. charge spilling/ {print $NF}", "lobsterout"],
            stdout=subprocess.PIPE,
            text=True,
        )
        spilling = awk_spil.stdout.strip().split("\n")

        bad_spilling = False
        for value in spilling:
            if float(value.strip("%")) > max_val:
                bad_spilling = True

        return bad_spilling

    def check_overlaps(self, max_val=0.05):
        """Check LOBSTER band overlap deviations against a tolerance.

        Args:
            max_val (float, optional): Maximum allowed overlap deviation value. Defaults to 0.05.

        Returns:
            bool: True if any overlap deviation exceeds the tolerance.
        """

        if os.path.exists("bandOverlaps.lobster"):
            awk_maxDev = subprocess.run(
                ["awk", "/maxDeviation/ {print $NF}", "bandOverlaps.lobster"],
                stdout=subprocess.PIPE,
                text=True,
            )

            dev_vals = list(map(float, awk_maxDev.stdout.strip().split("\n")))

        else:
            dev_vals = [0]

        for val in dev_vals:
            if val > max_val:
                return True

        return False

    def Run(self):
        """Run the LOBSTER projection across available basis-function combinations.

        Returns:
            str: "good_run" if one basis combination succeeds, "slurm" if there was an issue with the submission, otherwise "no_good_run".
        """
        os.chdir(self.JobSpecifics["JobPath"])
        self.ChangeCOHPErange()
        os.system("export OMP_NUM_THREADS=2")
        # self.WriteIncar()
        # self.WriteSubmission()
        # print(self.ledger)
        self.ledger.SetSingleValue(
            self.AssignedComp,
            self.JobName,
            "IncarExtras",
            self.JobSpecifics["IncarExtras"],
        )
        pos_bases = self.AssignedComp.get_all_basis_combos()

        with open(f"{self.AssignedCompName}_basisfunctions.json", "w") as file:
            json.dump(pos_bases, file, indent=4)

        for set_number in pos_bases.keys():
            try:
                os.mkdir(set_number)
            except FileExistsError:
                os.system(f"rm -rf {set_number}")
                os.mkdir(set_number)

            os.chdir(set_number)
            os.system("ln -f ../CHGCAR .")
            os.system("ln -f ../POSCAR .")
            os.system("ln -f ../WAVECAR .")
            os.system("ln -f ../INCAR .")
            os.system("ln -f ../POTCAR .")
            os.system("ln -f ../DOSCAR")
            os.system("ln -f ../KPOINTS .")
            os.system("ln -f ../IBZKPT .")
            os.system("ln -f ../EIGENVAL")
            os.system("ln -f ../vasprun.xml .")
            os.system("ln -f ../CONTCAR .")
            os.system("ln -f ../OUTCAR .")
            os.system("ln -f ../OSZICAR .")
            self.WriteIncar()
            self.add_basisfunctions(set_number)
            self.WriteSubmission()

            print(f"starting LOB for {self.AssignedCompName}")
            slurm_exit = submit_and_wait(self.SubFile)
            
            if slurm_exit != "COMPLETED":
                #some weird error, save slurm exit code to the completion field in the ledger
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_exit}")
                return 'slurm'

            # Check projection quality parameters
            spil_restart = self.check_charge_spilling()
            over_restart = self.check_overlaps()

            if spil_restart or over_restart:
                os.chdir(self.JobSpecifics["JobPath"])
                continue
            else:
                self.ledger.SetSingleValue(
                    self.AssignedComp, self.JobName, Field="completed", value=1
                )
                return "good_run"

        self.ledger.SetSingleValue(
            self.AssignedComp, self.JobName, Field="completed", value=-1
        )
        return "no_good_run"