import numpy as np
import pandas as pd # type: ignore
import os
import subprocess
import math
import sys
import json
import time
from pathlib import Path
from job_handler_bam import submit_and_wait
import ast
from filelock import FileLock, Timeout

chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                    'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

atomic_numbers = {'A':1}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z

atomic_spinstates = pd.read_csv('~/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/ht_programfiles/UnpairedSpins.csv', index_col=[0,1], skipinitialspace=True)


class compound:
    """This class contains all compound specific methods. In this case this means that this class contains all basic setups for treating Halide Double Perovskites (HDPs).
    For the initialization we need to know the present ions. This iteration assumes both A and X to be monovalent. 
     
        __init__() for this class requires A(str: element),B1(str: element),nB1(int:charge), B2(str: element),B2(str: element),nB2(int:charge),X(str: element)
    
    This class contains methods for 
        - getting the (electronic) spin states of the B-sites
        - writing POTCAR
        - writing KPOINTS grid
        - writing POSCAR
        - Getting basissets for LOBSTER projection
        
    """
    
    def __init__(self, A:str='Cs', B1:str='Pt', nB1:int=2, B2:str='Am', nB2:int=2, X:str='Cl'):
        self.A= A
        self.B1 = B1
        self.nB1 = nB1
        self.B2 = B2
        self.nB2 = nB2
        self.X = X

    def __repr__(self):
        return self.A + self.B1 + self.B2 + self.X

    def get_spinstate(self, which='B1'):
        """
        Helper function that returns the atomic spin bohr-magnetron number (Ignoring L and J contributions to the magnetic moment).
        The magnetron numbers are retreived from the file 'UnpairedSpins.csv', and are only defined for the transition metals and
        Lanthanides and Actinides, any other atom will return a magnetic moment of 0
        Inputs: 
              which:    Selection of which B-site to return the spinstates for
                        must be string 'B*' with *==1 or 2  
        Outputs:
                muHS, muLS  float of mu_s in Bohr magnetrons for the High and Low spinstate 
        """
        if which.upper() == 'B1':
            Zatom = atomic_numbers[self.B1]
            natom = self.nB1
        elif which.upper() == 'B2':
            Zatom = atomic_numbers[self.B2]
            natom = self.nB2
        else:
            raise ValueError("Cannot return spinstates for given B-site, fill in either 'B1' or 'B2'")
        
        try:
            muHS = atomic_spinstates.loc[Zatom, natom]["mus_HS"]
        except KeyError:
            muHS = 0.0
        
        try:    
            muLS = atomic_spinstates.loc[Zatom, natom]["mus_LS"]
        except KeyError:
            muLS = 0.0
        
        
        return muHS, muLS
    
    def write_geninfo(self):
        """writes a short textfile with an overview of the compound information known before any calculations
        """
        filename = self.name+"_geninfo.txt"
        b1HS, b1LS = self.get_spinstate('B1')
        #b2HS, b2LS = self.get_spinstate('B2')
        with open(filename, 'w') as file:
            file.write("This compound is %s \n" % (self.name))
            file.write("A-site is: %s \t with atomic number: %s \n" % (self.A, atomic_numbers[self.A]))
            file.write("B1-site is: %s \t with atomic number: %s \n" % (self.B1, atomic_numbers[self.B1]))
            file.write("\t has charge %s, and spinstates: %s [mu_B] and %s [mu_B] \n" % (self.nB1, b1HS, b1LS))
            if self.B2 == 'Vac':
                file.write("B2-site is a vacancy \n")
            else:
                b2HS, b2LS = self.get_spinstate('B2')
                file.write("B2-site is: %s \t with atomic number: %s \n" % (self.B2, atomic_numbers[self.B2]))
                file.write("\t has charge %s, and spinstates: %s [mu_B] and %s [mu_B] \n" % (self.nB2, b2HS, b2LS))
            file.write("X-site is %s \t with atomic number: %s \n" % (self.X, atomic_numbers[self.X]))
        return               
 
    def write_poscar(self, a = 5.5):
        """
        This function writes a standard primitive Fm3m unit cell POSCAR into the current directory. 
        It will check if B2-site is a vacancy and adjust the POSCAR accordingly
        Args:
            a (float): the lattice paramater of the primitive unit cell
        """
           

        with open('POSCAR', 'w') as f:
            f.write("%s \n1.00000000000000 \n%.10f    %.10f    %.10f\n" %(repr(self),0.0,a,a ))
            f.write("%.10f    %.10f    %.10f\n%.10f    %.10f    %.10f\n" %(a, 0.0, a, a, a, 0.0))
            
            if self.B2=='Vac':
                f.write("%s    %s    %s\n1    2    6\n" %(self.B1, self.A, self.X))
                f.write("Direct\n  %.16f  %.16f  %.16f\n" % (0.0, 0.0, 0.0)) 
            else:
                f.write("%s    %s    %s    %s\n1    1    2    6\n" %(self.B1, self.B2, self.A, self.X))
                f.write("Direct\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n" % (0.0, 0.0, 0.0, 0.5, 0.5, 0.5)) 
                
            f.write("  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n" %(0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75, 0.25, 0.25, 0.25, 0.75, 0.75))
            f.write("  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n" %(0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75, 0.75, 0.75, 0.25))
        return
    
    def write_kpoints(self, npoints=4):
        """Creates the file KPOINTS in the current directory. 
        This function only supports cubic equidistant gamma-centered kpoint-grids (as this is the relevant grid for this studycase) 

        Args:
            npoints (int, optional): number of k-grid points in each direction. Defaults to 4.
        """
        with open('KPOINTS', 'w') as file:
            file.write('k-points\n')
            file.write(' 0\nGamma \n %i  %i  %i \n 0  0  0' %(npoints,npoints,npoints))
        return
    
    def write_potcar(self, use_frozen: bool = False, write_file:bool = True):
        """This function concatonates the POTCARS for each species defined by the object compound. Which specific POTCAR is used for each
        species is determined by potLUT (potential Look Up Table). 
        If any species is given that is not present in the potLUT this function will fail

        Arguments:
            use_frozen(bool):   This specifies whether we want to use potential with extra electrons in the frozen core. 
                                This option is specifically used for the relaxation step to improve stability
        """
        if use_frozen:
            LUTpath = "/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/ht_programfiles/potLUT"
        else:
            LUTpath = "/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/ht_programfiles/potLUT_nofrozen"
        
        PPpath = "/home/lucw/VASP-PAW-pseudo/potpaw_PBE.64/"
        if os.path.isfile("POTCAR") and write_file:
            #print('POTCAR encountered, will be overwritten...')
            os.system("rm POTCAR")
        
        potcar_list = []
        for el in [ self.B1, self.B2, self.A, self.X]:
            #print(el)
            if el == "Vac":
                continue
            #grep_proc = subprocess.Popen(["grep", el, LUTpath], stdout=subprocess.PIPE, text=True)
            awk_proc = subprocess.Popen(["awk", "/%s\t/{print $2}"%(el), LUTpath], stdout=subprocess.PIPE, text=True)
            potname, _ = awk_proc.communicate()
            potcar_list.append(potname.strip())
            if write_file:
                potdir = PPpath + potname.strip() +  "/POTCAR"
                os.system("cat %s >> POTCAR" %(potdir))
        


        return potcar_list
    
    def get_basis_functions(self,basis="min"):
        """This function will look at BASISPBE.yaml file (adapted from the Lobsterpy package) and return the set of basisfunctions 
        that should be used by LOBSTER. There are two sets to chose from, the minimal or maximal basis set.

        Args:
            basis (str, optional): Should be either "min" or "max" and determines which basisset-file to look in. Defaults to "min".

        Returns:
            dict: a Dictonary of 'element':'basis functions' that are to be used in LOBSTER projection
        """

        basisfile = f"/home/lwalterb/hdp_project/ht_programfiles/BASIS_PBE_64_{basis}.yaml"
        LUTpath = "/home/lwalterb/hdp_project/ht_programfiles/potLUT_nofrozen"

        with open(LUTpath, 'rt') as f:
            LUTlines = f.readlines()

        with open(basisfile, 'rt') as bases:
            basislines = bases.readlines()

        basis_dict = {}

        for el in [ self.B1, self.B2, self.A, self.X]:
            if el == "Vac":
                continue

            potLUTline = [line for line in LUTlines if line.startswith(f"{el}\t")]
            assert len(potLUTline) == 1, f"Looking up the required POTCAR fo r{el} returned multilple options: {potLUTline}"
            potname = potLUTline[0].strip('\n').strip().split('\t')[1]
            # print(potname)


            basisline = [line for line in basislines if f"{potname}:" in line]
            assert len(basisline) == 1, f"Looking up the required LOBSTER basis for {potname} returned multilple options: {basisline}"
            basisfuncs = basisline[0].strip('\n').strip().split(':')[1]
            # print(basisfuncs)

            # awk_bases = subprocess.run(['awk', '-F:', '/%s:/ {print $2}'%(potname), basisfile], stdout=subprocess.PIPE, text = True)
            # print(awk_bases.stdout.strip())
            # basisfuncs = awk_bases.stdout.strip()
            # if '\n' in basisfuncs:
            #     basisfuncs = basisfuncs.split('\n')[0]

            basis_dict.update({el: basisfuncs})

        return basis_dict

    def get_all_basis_combos(self):
        """This function will return all possible basisset combinations for the LOBSTER projection. This is based on the basis sets 
        available in BASIS_PBE_64_{min/max}.yaml which is ultimately derived from the listed orbitals in the POTCAR files.

        Returns:
            dict: A nested Dictonary of every basisset combination (given iteratively, i.e., basis0:..., basis1:..., etc...) and for each basis 'element':'basis;
        """

        min_basis = self.get_basis_functions(basis="min")
        max_basis = self.get_basis_functions(basis="max")

        all_base_dict = {}
        base_index = 0
        base_name = "basisset" + str(base_index)

        all_base_dict[base_name] = min_basis
        
        extraB1 = False
        extraB2 = False

        if min_basis[self.B1] != max_basis[self.B1]:
            dummy = min_basis.copy()
            dummy.update({self.B1: max_basis[self.B1]})
            
            base_index += 1
            base_name = "basisset" + str(base_index)

            all_base_dict[base_name] = dummy
            extraB1 = True

        if self.B2 != "Vac":
            if min_basis[self.B2] != max_basis[self.B2]:
                dummy = min_basis.copy()
                dummy.update({self.B2: max_basis[self.B2]})

                base_index += 1
                base_name = "basisset" + str(base_index)   

                all_base_dict[base_name] = dummy    

                extraB2 = True

        if extraB1 and extraB2:
            base_index += 1
            base_name = "basisset" + str(base_index)
            all_base_dict[base_name] = max_basis


    
        return all_base_dict
        
    def get_maxNBANDS(self):
        """This function will take the max basis functions for a specific HDP, and calculate the number of bands
        required in the DFT calculation for the LOBSTER projection.

        Returns:
            int: The number of bands required for the maximal basis set for this HDP composition
                    This is to be set in the INCAR file as NBANDS
            
            """
        max_basis = self.get_basis_functions("max")
        
        nbands = 0
        
        for el, basis in max_basis.items():
            mult = 0  #multiplier to account for number of occurences in unit cell
            if el.split('_')[0] == self.A:
                mult += 2
            if el.split('_')[0] == self.X:
                mult += 6
            if el.split('_')[0] == self.B1:
                #necessary in case of duplicate entries (e.g. Cs2TlTlCl6 or Cs2CsVBr)
                mult+= 1
            if el.split('_')[0] == self.B2:
                mult += 1
            

            for band in basis.split():

                if band.endswith('s'):
                    nbands += 1*mult

                elif band.endswith('p'):
                    nbands += 3 * mult

                elif band.endswith('d'):
                    nbands += 5 * mult

                elif band.endswith('f'):
                    nbands += 7 * mult

        return nbands








def initialize_compounds(file:str ="/home/lucw/ht_programfiles/testinput_larger"):
    """
    This function will take the input file, create an compound object for each entry, create the directory and write the general info file
    The input file shoud be parsed as follows "ABBX, 'A', 'B1', nB1, 'B2', nB2, 'X', where '*' indicates the chemical symbol of the element,
    and nB indicates the ionic state of the B site (integer) 
    
    The function will return a list of objects of all compounds
    """
    if not os.path.isfile(file):
        raise TypeError('intialization file could not be found/interperted as file')

    
    list_of_compounds = []
    with open(file,'r') as f:
        for line in f.readlines():
            inputs = line.strip().strip('\n').split(',')
            new_comp = compound(inputs[1],inputs[2],int(inputs[3]),inputs[4],int(inputs[5]),inputs[6])
            #print(len(inputs))
            if len(inputs) == 8:
                new_comp.CompID = str(inputs[-1]) + '_' + str(repr(new_comp))
                #print(f"CompID is {new_comp.CompID}")

            
            list_of_compounds.append(new_comp)

    return list_of_compounds




class ProcessLedger:
    def __init__(self, JobSettingsPath:str,StartPath:str = '.', ledger_filename: str="JobInformationLedger.csv"):
        #the ledger will be based in the parent directory of all jobfiles
        if StartPath == '.':
            self.basepath = Path(os.getcwd())
        else:
            self.basepath = Path(StartPath)
        

        #Load in the required jobs
        assert JobSettingsPath.endswith('.json'), "JobSettings needs to be JSON file"
        self.JobsPath = JobSettingsPath



        with open(self.JobsPath,'r') as jfile:
            #this stores the full info per step
            self.JobInfo = json.load(jfile)

        self.Jobs = list(self.JobInfo.keys())

        self.LedgerFile = self.basepath / ledger_filename
        self.inputbackup = self.basepath / str("UsedInput_" + ledger_filename)
        self.LedgerBackup = self.basepath / str("BackupLedger_" + ledger_filename)

        self.lock = FileLock(str(self.LedgerFile) + ".lock",thread_local=False)
        #this is what info the ledger stores per job
        self.info_per_job = ["JobPath","completed","JobID","TimeStamp", "IncarExtras"]  
        # try:
        #     os.mkdir(self.path)
        # except:
        #     pass 
    



    def StartNewLedger(self,  list_of_compounds:list = []):    
        os.chdir(self.basepath)

        ##Setup dataframe for the ledger info
        index = pd.MultiIndex.from_product([list(self.JobInfo.keys()),self.info_per_job])    
        ledger = pd.DataFrame(index=index)
        ledger2 = pd.DataFrame(index = index)

        #All comps get a CompID
        #which starts at 1000 for esthetic reasons
        CompID = 1000
        file = open(self.inputbackup,'w')

        for comp in list_of_compounds:

            comp.CompID = str(CompID) + '_' + str(repr(comp))
            comp.path = str(self.basepath) +'/'   + str(CompID)+ "_" + str(repr(comp))
            file.write(f"{repr(comp)},{comp.A},{comp.B1},{comp.nB1},{comp.B2},{comp.nB2},{comp.X},{CompID}\n")
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

                ##If it is a VASP or misc job a subdirectory will be made
                    #if the subdir already exists, it will be assumed that the step was
                    #completed succesfully
                ##For LOBSTER jobs it will set the subdir from another step as working dir
                if settings['job_type'] != "LOBSTER":
                    calc_subdir = str(JobName) + "_" + settings['name']
                    calcdir_full = os.path.join(comp.path, calc_subdir)
                    try:
                        os.mkdir(calcdir_full)
                    except FileExistsError:
                        ##If Vasp SubDir exists job is considered completed
                        completion = 1

                    #If VASP job we already write KPOINTS and POTCAR
                    # else:
                    #     if settings['job_type'] == 'VASP':
                    #         os.chdir(calcdir_full)
                    #         comp.write_kpoints()
                    #         comp.write_potcar()


                else:
                    assert "use_output_step" in settings.keys(), "LOBSTER job is missing required VASP directory.\n Please define 'use_output_step' in JobSettingsJson"
                    subdir = settings['use_output_step'] + "_" + self.JobInfo[settings['use_output_step']]['name']
                    calcdir_full = os.path.join(comp.path, subdir)
                    
                    

                #store all base information for each step
                LedgerInfo.update({JobName: {
                    "JobPath" : str(calcdir_full),
                    "JobID" : "",
                    "completed" : completion,
                    "TimeStamp" : "",
                    "IncarExtras" : {}
                }})


            
            #flatten the dict of ledger info and add as new column to dataframe
            flat_dict = {(outer_key, inner_key): value 
                for outer_key, inner_dict in LedgerInfo.items() 
                for inner_key, value in inner_dict.items()}    

            ledger2[comp.CompID] = ledger.index.map(flat_dict).copy()
        self.LedgerDF = ledger2.sort_index()
        with self.lock:
            ledger2.to_csv(self.LedgerFile)
        file.close()
        return

    def AppendNewComps(self,  list_of_compounds:list = []):    
        os.chdir(self.basepath)

        ##Setup dataframe for the ledger info
        index = pd.MultiIndex.from_product([list(self.JobInfo.keys()),self.info_per_job])    
        ledger = pd.DataFrame(index=index)
        ledger2 = pd.DataFrame(index = index)

        #All comps get a CompID
        #We read the last used from UsedInput_backup 
        with open(self.inputbackup,'r') as f:
            lines = f.readlines()
        
        last_input = lines[-1].strip().split(',')


        CompID = int(last_input[-1]) + 1
        file = open(self.inputbackup,'a')

        for comp in list_of_compounds:

            comp.CompID = str(CompID) + '_' + str(repr(comp))
            comp.path = str(self.basepath) +'/'   + str(CompID)+ "_" + str(repr(comp))
            file.write(f"{repr(comp)},{comp.A},{comp.B1},{comp.nB1},{comp.B2},{comp.nB2},{comp.X},{CompID}\n")
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

                ##If it is a VASP or misc job a subdirectory will be made
                    #if the subdir already exists, it will be assumed that the step was
                    #completed succesfully
                ##For LOBSTER jobs it will set the subdir from another step as working dir
                if settings['job_type'] != "LOBSTER":
                    calc_subdir = str(JobName) + "_" + settings['name']
                    calcdir_full = os.path.join(comp.path, calc_subdir)
                    try:
                        os.mkdir(calcdir_full)
                    except FileExistsError:
                        ##If Vasp SubDir exists job is considered completed
                        completion = 1

                    #If VASP job we already write KPOINTS and POTCAR
                    # else:
                    #     if settings['job_type'] == 'VASP':
                    #         os.chdir(calcdir_full)
                    #         comp.write_kpoints()
                    #         comp.write_potcar()


                else:
                    assert "use_output_step" in settings.keys(), "LOBSTER job is missing required VASP directory.\n Please define 'use_output_step' in JobSettingsJson"
                    subdir = settings['use_output_step'] + "_" + self.JobInfo[settings['use_output_step']]['name']
                    calcdir_full = os.path.join(comp.path, subdir)
                    
                    

                #store all base information for each step
                LedgerInfo.update({JobName: {
                    "JobPath" : str(calcdir_full),
                    "JobID" : "",
                    "completed" : completion,
                    "TimeStamp" : "",
                    "IncarExtras" : {}
                }})


            
            #flatten the dict of ledger info and add as new column to dataframe
            flat_dict = {(outer_key, inner_key): value 
                for outer_key, inner_dict in LedgerInfo.items() 
                for inner_key, value in inner_dict.items()}    

            ledger2[comp.CompID] = ledger.index.map(flat_dict).copy()
        ledger = ledger2.sort_index()
        with self.lock:
            prevLedger = self.ReadFullLedger()
            newLedger = prevLedger.join(ledger,how='left', validate ='1:1')
            newLedger.to_csv(self.LedgerFile)



        file.close()
        return

    def RestartLedger(self):
        loc = initialize_compounds(self.inputbackup)
        return loc

    def ReadFullLedger(self):
        with self.lock: 
            self.LedgerDF = pd.read_csv(self.LedgerFile, index_col=[0,1],header=0, dtype=str).sort_index()
        return self.LedgerDF


    def GetSingleJob(self,comp:compound, JobName:str):
        """Returns a Dict of the information stored in the process ledger of the step specified by JobName
        for a given compound

        Args:
            comp (compound): _description_
            JobName (str): _description_

        Returns:
            dict: a dictionary of all information stored in the ProcessLedger for a specific step for a specific compound
        """
        comp_col = self.ReadFullLedger()#, dtype={"IncarExtras":dict})
        jobdict = comp_col.loc[JobName][str(comp.CompID)].to_dict()
        try:
            jobdict["IncarExtras"] = ast.literal_eval(jobdict["IncarExtras"])
        except:
            jobdict["IncarExtras"] = {}

        return jobdict

    def GetCompletionOverview(self):
        ledger = self.ReadFullLedger()
        summary =  ledger.xs('completed',level=1).T
        return summary.astype(str)
    
    def GetQueue(self):
        overview = self.GetCompletionOverview()
        AllComps = overview.index

        DoneComps = [] #used to keep track which compounds have already been treated
        
        stepnames = self.JobInfo.keys()
        QueueDict = {} 

        ##First remove any Errors (so if completed is anything else than 1 | 0 )
        errors = overview[((overview != '1') & (overview != '0' )).any(axis=1)].index.to_list()

        if errors:
            for error in errors: 
                DoneComps.append(error)

        #DoneComps.append(AllComps.to_list()[0])
        
        
        for JobName in stepnames:
            stepcomps = []
            #Remove any compounds contained in DoneComps
      
            newindex = AllComps.difference(pd.Index(DoneComps), sort=False)
            overview = overview.loc[newindex]

            #Now iterate over steps and find the comps that still need to complete that step
            stepcomps = overview[overview[JobName]=='0'].index.to_list()
            #Add these comps to the queue dict and append them to DoneComps
            QueueDict[JobName] = stepcomps
            if stepcomps:
                for comp in stepcomps:
                    DoneComps.append(comp)

        return QueueDict



    def SetSingleValue(self,comp:compound,JobName:str,Field:str,value):
        assert Field in self.info_per_job, f"Trying to change a non-existent field in Ledger: {Field}"
        
        
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
        import glob
        from job_handler import check_scfcompletion
        os.chdir(self.basepath)

        loc = self.RestartLedger()

        index = pd.MultiIndex.from_product([list(self.JobInfo.keys()),self.info_per_job])    
        
        backup = pd.read_csv(os.path.join(self.basepath,"AllCompsCopy.csv"), index_col=[0,1],header=0).sort_index()


        ledger = pd.DataFrame(index=index, columns = backup.columns, dtype=str)
        with self.lock:
            ledger.to_csv(self.LedgerFile)
        ledger2 = pd.DataFrame(index = index, columns = backup.columns, dtype = str)

        for comp in loc:
            print(comp)
            comp.path = str(self.basepath) +'/'   + str(comp.CompID)
            


            LedgerInfo = {}
            jnames = list(self.JobInfo.keys())
            #jnames.reverse()
            

            for JobName in jnames:
                os.chdir(comp.path)
                settings = self.JobInfo[JobName]

                completion = 0

                ##If it is a VASP or misc job a subdirectory will be made
                    #if the subdir already exists, it will be assumed that the step was
                    #completed succesfully
                ##For LOBSTER jobs it will set the subdir from another step as working dir
                if settings['job_type'] != "LOBSTER":
                    calc_subdir = str(JobName) + "_" + settings['name']
                    calcdir_full = os.path.join(comp.path, calc_subdir)
                    
                    #Here we need to figure out whether the job was completed or not
                    if JobName == "1Rel":
                        if os.path.isfile(f"{comp.CompID}_relaxedPOSCAR"):
                            completion = 1

                    elif JobName == "2Spins":
                        res_file = glob.glob("SpinResults*")
                        if res_file:
                            completion = 1

                            spinDF = pd.read_csv(res_file[-1], index_col=0)

                            try:
                                converged_runs = spinDF.groupby('SCF_convergence').get_group(0)
                            except KeyError:
                                completion = 0
                            
                            else:

                                optimal_run = converged_runs['Etot_out'].idxmin()
                                magmom_use = f"{spinDF['muB1_in'].loc[optimal_run]} {spinDF['muB2_in'].loc[optimal_run]} 2*0.0 6*0.0"


                                ##now we want to put this MAGMOM initialization into the IncarExtras for all remaining VASP jobs
                                jobindex = self.Jobs.index(JobName)

                                for futurejob in self.Jobs[jobindex+1:]:
                                    if self.JobInfo[futurejob]['job_type'].upper() == "VASP":
                                        IncarExtras = self.GetSingleJob(comp,JobName)["IncarExtras"]

                                        IncarExtras.update({"MAGMOM": magmom_use})
                                        self.SetSingleValue(comp,futurejob,"IncarExtras", str(IncarExtras))

                    
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
                    assert "use_output_step" in settings.keys(), "LOBSTER job is missing required VASP directory.\n Please define 'use_output_step' in JobSettingsJson"
                    subdir = settings['use_output_step'] + "_" + self.JobInfo[settings['use_output_step']]['name']
                    calcdir_full = os.path.join(comp.path, subdir)

                    

                #store all base information for each step
                LedgerInfo.update({JobName: {
                    "JobPath" : str(calcdir_full),
                    "JobID" : "",
                    "completed" : completion,
                    "TimeStamp" : "",
                    "IncarExtras" : self.GetSingleJob(comp, JobName)["IncarExtras"]
                }})


            
            #flatten the dict of ledger info and add as new column to dataframe
            flat_dict = {(outer_key, inner_key): value 
                for outer_key, inner_dict in LedgerInfo.items() 
                for inner_key, value in inner_dict.items()}    

            ledger2[comp.CompID] = ledger.index.map(flat_dict).copy()
        self.LedgerDF = ledger2.sort_index()
        with self.lock:
            ledger2.to_csv(self.LedgerFile)
        return


class GeneralJob:

    def __init__(self, ledger:ProcessLedger, comp:compound, JobName:str = "", AssignedServer:str = "ccp20,ccp22"):
        assert JobName in ledger.Jobs, "JobName passed to (General)Job CLASS not present in the Ledger"
        
        self.AssignedServer = AssignedServer
        self.AssignedCompName = comp.CompID
        self.AssignedComp = comp
        self.ledger = ledger
        #retrieve the basic settings of this job from the job setting JSON
        self.JobName = JobName
        self.JobBasics = ledger.JobInfo[JobName]
        
        self.JobID = int(float(comp.CompID.split('_')[0])* time.time())
        ledger.SetSingleValue(comp, JobName, Field='JobID', value=self.JobID)
        ledger.SetSingleValue(comp, JobName, Field='TimeStamp', value=time.strftime("%Y/%m/%d_%H:%M:%S"))



        #Also regtrieve the extra information stored in the ledger
        self.JobSpecifics = dict(ledger.GetSingleJob(comp, JobName))


    def __repr__(self):
        return self.AssignedCompName +'_' + self.JobName


    def WriteIncar(self):
        """
        This function will take the standard incar from JobSettings JSON file, and any extra inputs in the ProcessLedger
        put it together and write it to a file. It will automatically save to 'INCAR' for 'job_type:VASP' and lobsterin for 'job_type:LOBSTER'

        for any other job_type it will take ['name'] from JobSettings JSON and output it to '{name}_INPUT'
        """
        incar = self.JobBasics['std_incar']
        incar.update(self.JobSpecifics["IncarExtras"])

        if self.JobBasics['job_type'].upper() == "VASP":
            filename = "INCAR"
            joiner = " = "
        elif self.JobBasics['job_type'].upper() == "LOBSTER":
            filename = "lobsterin"
            joiner = " "
        else:
            filename = f"{self.JobBasics['name']}_INPUT"
            joiner = " "

        with open(filename,'w') as f :
            for key, val in incar.items():
                # try:
                #     literal_val = eval(val)
                # except Exception as exc:
                #     print(f"tried to eval: {val} as literal val but got exception {exc}")
                #     newline = joiner.join([str(key),str(val)]) + '\n'
                #     f.write(newline)
                # else:
                if type(val) == list:
                    for singleval in val:
                        newline = joiner.join([str(key),str(singleval)]) + '\n'
                        f.write(newline)
                else:
                    newline = joiner.join([str(key),str(val)]) + '\n'
                    f.write(newline)

        return incar

    def WriteSubmission(self, LogFileExtra = ''):
        
        filename = self.JobBasics['job_type'].lower() + '.sub'
        self.SubFile = filename
        self.LogFile = f"log_{self.JobID}_{self.AssignedCompName}{LogFileExtra}"
        self.ErrorFile = f"error_{self.JobID}_{self.AssignedCompName}{LogFileExtra}"
        ##extract nnodes and potentially time from JobSettings JSON
        nnodes = self.JobBasics['nnodes']
        if 'max_time' in self.JobBasics.keys():
            max_time = self.JobBasics['max_time']
        else:
            max_time = f'0{nnodes}:00:00'



        if self.JobBasics['job_type'].upper() == "VASP":
            module = "module load vasp/6.3.0"
            runcommand = "mpirun vasp_std"
        elif self.JobBasics['job_type'].upper() == "LOBSTER":
            module = "module load lobster"
            runcommand = "lobster"
        else:
            raise ModuleNotFoundError(f"Could not write submission script for unknown job_type: {self.JobBasics['job_type']}")

        with open(filename,'w') as file:
            file.write(f'#!/bin/bash \n#SBATCH --time={max_time}\n')
            file.write(f'#SBATCH -p {self.AssignedServer} \n#SBATCH -N {nnodes} --exclusive --ntasks-per-node 32\n')
        
            file.write(f'#SBATCH -J {self.AssignedCompName}{self.JobName} \n#SBATCH --output={self.LogFile}\n')
            file.write(f'#SBATCH --ntasks-per-core=1 \n#SBATCH --error={self.ErrorFile}\n\n')
            file.write(f'module purge \n{module} \n\nulimit -s unlimited\n')
            file.write('startingtime=$(date)\n')
            file.write(f"now=$(date +'%s') \n{runcommand} \nendtime=$(date +'%s') \nduration=$(($endtime-$now)) \n")
            file.write("echo 'Started on: ' $startingtime 'calculation ran for: ' $duration '[s]'\n")


        return

    def Run(self):
        os.chdir(self.JobSpecifics['JobPath'])

        self.WriteIncar()
        self.WriteSubmission()

        slurm_exit = submit_and_wait(self.SubFile)

        if slurm_exit == "COMPLETED":
            self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=1)
            return

        else:
            #some weird error, save slurm exit code to the completion field in the ledger
            self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_exit}")
            return



class SimpleVasp(GeneralJob):
    def __init__(self, ledger:ProcessLedger, comp:compound, JobName:str = "", AssignedServer:str = "ccp20,ccp22"):
        super().__init__( ledger, comp, JobName, AssignedServer)

        #extra attribute specifically for vasp jobs: which algos are allowed
        self.AllowedAlgos = [ "N","F", "D", "VF", "A", "S"]

    def CycleAlgos(self):
        """This function will update the chosen algorithm for a VASP calculation based on the available algorithms
        in GeneralJob.AllowedAlgos. 
        
        If the "IncarExtras" in the ProcessLedger does not have an ALGO specified for the current step it will check the previous step for that keyword. 
            If it is specified it there it will update the "IncarExtras" for the current step with that ALGO
            Else it will start with AllowedAlgos[1]
        
        If the "IncarExtras" of the current step does have an ALGO specified, it will update it in the ProcessLedger for the next algo specified in the list
        """

        if "ALGO" in self.JobSpecifics['IncarExtras'].keys():
            prev_algo = self.JobSpecifics['IncarExtras']["ALGO"]
            
            try:
                self.JobSpecifics['IncarExtras']["ALGO"] = self.AllowedAlgos[self.AllowedAlgos.index(prev_algo) + 1]

            except ValueError:
                self.JobSpecifics['IncarExtras']["ALGO"] = self.AllowedAlgos[1]
            except IndexError:
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName,"completed", "AlgoError OOA: OutOfAlgos")
                raise IndexError(f"CycleAlgos wanted to change ALGO for compound:{self.AssignedCompName} but went beyond the max(?) index of Job.Allowed.Algos list")

    

        else:
            try:
                prev_job_settings = self.ledger.GetSingleJob(self.AssignedComp, self.ledger.Jobs[self.ledger.Jobs.index(self.JobName) - 1])
            except IndexError:
                self.JobSpecifics['IncarExtras']["ALGO"] = self.AllowedAlgos[1]
            else:
                try:
                    if prev_job_settings["IncarExtras"]["ALGO"] in self.AllowedAlgos:
                        self.JobSpecifics['IncarExtras']["ALGO"] = prev_job_settings["IncarExtras"]["ALGO"]
                    else:
                        self.JobSpecifics['IncarExtras']["ALGO"] = self.AllowedAlgos[1]
                except KeyError:
                    self.JobSpecifics['IncarExtras']["ALGO"] = self.AllowedAlgos[1]

        self.ledger.SetSingleValue(self.AssignedComp, self.JobName, "IncarExtras", self.JobSpecifics["IncarExtras"] )

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
        p1 = subprocess.run(["awk", 'BEGIN{rv=0; print "\\nNew run started --->\\n\\n" >> "%s"};/self-consistency was not achieved/ {rv=1}; /F=/ {print rv; rv=0}; \
                            /E{7} * R * R * R * R * O{7}/ {print "-1"}; /reached required accuracy - stopping structural energy minimisation/ {print "2"};\
                            /E{7} * R * R * R * R * O{7}/, /---->/ {if(!/E{7} * R/ && !/---->/)print $0 >> "%s"}'%(erroroutput, erroroutput),\
                            self.LogFile], stdout=subprocess.PIPE, text=True)

        numstrings = p1.stdout.strip().split('\n')
        try:
           scf_cycles = list(map(int,numstrings))
        except ValueError as err:
            #raise Exception("Unexpected value encountered when checking SCF convergence:{err}")
            return numstrings
        return scf_cycles

    def PrepareVASPRun(self):
        """This function simply makes sure all required input files for VASP (except for INCAR!) are present
        in the Job subfolder. I.E. It will write KPOINTS, POTCAR and POSCAR if they are not present yet.

        For POSCAR it will check if the JobSettings.json has the keyword 'use_poscar_from', 
            if it is given, it will take CONTCAR from the specified JobName and copy it to the current job as POSCAR.
            else it will use the write_poscar() method from the class compound to write the initial poscar

        If the JobSettings.json also has the keyword 'extra_input_files' it will link the specified files from the specified job
        """
        #os.chdir(self.JobSpecifics["JobPath"])
        
        if not os.path.exists('KPOINTS'):
            self.AssignedComp.write_kpoints()

        
        self.AssignedComp.write_potcar()

        if "extra_input_files" in self.JobBasics.keys():
            for job, files in self.JobBasics['extra_input_files'].items():
                sourcedir = self.ledger.GetSingleJob(self.AssignedComp,job)["JobPath"]

                if type(files) == list:
                    for file in files:
                        sourcefile = os.path.join(sourcedir,file)
                        if file == "CONTCAR":
                            
                            os.system(f"cp -f {sourcefile} ./POSCAR")
                        else:
                            
                            if not os.path.exists(file):   
                                os.system(f"ln -f {sourcefile} .")

                else:
                    sourcefile = os.path.join(sourcedir,files)
                    if files == "CONTCAR":
                        
                        os.system(f"cp -f {sourcefile} ./POSCAR")
                    else:
                        
                        if not os.path.exists(files):
                            os.system(f"ln -f {sourcefile} .")



        if not os.path.exists('POSCAR'):
            self.AssignedComp.write_poscar()

        return
    
    def ExtractEtot(self, file="OSZICAR"):
        #os.chdir(self.JobSpecifics['JobPath'])
        p = subprocess.run(['awk', "/F=/ {print $5}", file],stdout = subprocess.PIPE, text=True)
        output = p.stdout.strip().split('\n')
        if len(output) > 1:
            Eout = list(map(float, output))
        else:
            Eout = float(output[0])

        return Eout
    
    def ExtractMAGtot(self, file="OSZICAR"):
        #os.chdir(self.JobSpecifics['JobPath'])
        p = subprocess.run(['awk', "/F=/ {print $NF}", file],stdout = subprocess.PIPE, text=True)
        output = p.stdout.strip().split('\n')
        if len(output) > 1:
            MAGMOM = list(map(float, output))
        else:
            MAGMOM = float(output[0])

        return MAGMOM

    def UsePOTCARENMAX(self, MultFactor:float=1.2 , MinValue:int= 350):

        #os.chdir(self.JobSpecifics["JobPath"])
        get_enmax = subprocess.run(["awk", "BEGIN{a = 0} /ENMAX/ {if($3>0+a) a=$3 fi} END{print a}", "POTCAR"], stdout=subprocess.PIPE, text=True)
        enout = float(get_enmax.stdout.strip().strip(';'))

        ENCUT = MultFactor * enout

        if ENCUT < MinValue:
            ENCUT = MinValue

        self.JobSpecifics["IncarExtras"].update({"ENCUT": ENCUT})
        

        return ENCUT

    
    def Run(self):
            os.chdir(self.JobSpecifics['JobPath'])

            self.WriteIncar()
            self.WriteSubmission()

            slurm_out = submit_and_wait(self.SubFile)
            self.ledger.SetSingleValue(self.AssignedComp,self.JobName,"IncarExtras", self.JobSpecifics["IncarExtras"])
            if slurm_out != "COMPLETED":
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_out}")
                return 'slurm'


            scf_out = self.check_scfcompletion()

            if scf_out[-1] == 0:
                #job succesfull, update this in ledger and return
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=1)
                return 'succesfull'

            elif scf_out[-1] == 1:
                #convergence failed, so make useful changes
                self.CycleAlgos()
                return 'convergence_failed'

            elif scf_out[-1] == -1:
                #error occured, deal with it appropriate for this step
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=-1)
                return 'error'
            else:
                #something else went wrong in converging
                #convergence failed, so make useful changes
                self.CycleAlgos()
                return 'run_failed'




class Relaxation(SimpleVasp):
    def __init__(self, ledger:ProcessLedger, comp:compound, JobName:str = "", AssignedServer:str = "ccp20,ccp22"):
        super().__init__(ledger, comp, JobName, AssignedServer)
        self.AllowedAlgos = ["N","A", "D","F", "VF",  "S"]



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
        p1 = subprocess.run(['awk',  '/TOTAL-FORCE/ {getline;getline; for (i=0; i<=9; i++)\
            print $1"   "$2"   "$3 getline; print"\\n"}',outcar],text=True, stdout=subprocess.PIPE)

        all_pos = p1.stdout.strip().split('\n')
        i_empties = [i for i, x in enumerate(all_pos) if x=='']
        try:
            begin_pos = all_pos[0:i_empties[0]]
        except IndexError:
            if len(i_empties)==0:
                #only one set of positions was found so no success
                return False
            else:
                return False
        except Exception as exc:
            #something else went wrong
            print('unexpected exception in check_changepos: %s'%(exc))
            return False
        else:

            end_pos = all_pos[i_empties[-1]:]
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
        pzbrent = subprocess.run(['awk', 'BEGIN{enc=0}; /rerun with smaller EDIFF/ {enc++}; END{print enc}',self.LogFile],stdout=subprocess.PIPE, text=True)
        zbrent_encounters = int(pzbrent.stdout.strip())
        return bool(zbrent_encounters)
    
    def Run(self):
        
        
        
        os.chdir(self.JobSpecifics['JobPath'])

        ##Preparation steps
        self.PrepareVASPRun() #assures KPOINTS, POSCAR, and POTCAR are there
        self.AssignedComp.write_potcar(use_frozen=True)
        self.UsePOTCARENMAX( MultFactor=1.5, MinValue = self.JobBasics["std_incar"]["ENCUT"])
        incar = self.WriteIncar()

        self.WriteSubmission()
        slurm_out = submit_and_wait(self.SubFile)
        if slurm_out != "COMPLETED":
            self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_out}")
            return 'slurm'
        
        scf_out = self.check_scfcompletion()

        if scf_out[-1] == 2:
            #this is a special flag recognizing VASP successful relaxation message
            self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=1)
            os.system(f"cp CONTCAR ../{self.AssignedCompName}_relaxedPOSCAR")
            return 'relaxed'
      
        poschange = self.check_changepos()
        contcarrestart = self.check_zbrent_message()

        if not poschange:
            #Relaxation returned to initial position, this can be solved by initializing with larger a_lat
            self.AssignedComp.write_poscar(a = 5.7)

        if -1 in scf_out:
            #check for errors in the SCF cycles
            #many errors are solved by using ALGO = All and high PREC
            self.CycleAlgos()
            self.JobSpecifics["IncarExtras"].update({ "PREC": "A", "NELM": 150})
            self.ledger.SetSingleValue(self.AssignedComp,self.JobName, "IncarExtras", self.JobSpecifics['IncarExtras'])


        
        if contcarrestart:
            #special flag for VASP telling you to restart from CONTCAR with smaller EDIFF
            os.system("cp CONTCAR POSCAR")
            newEDIFF = incar['EDIFF']/10
            self.JobSpecifics['IncarExtras'].update({"EDIFF": newEDIFF})
            self.ledger.SetSingleValue(self.AssignedComp,self.JobName, "IncarExtras", self.JobSpecifics['IncarExtras'])
        elif np.sum(scf_out)/len(scf_out) > 0.3:
            #more than 30% of SCF cycles did not converge, maybe cyclin algos will help
            self.CycleAlgos()

        return f'relaxation_failed, Poschange: {poschange}, Restart CONTCAR: {contcarrestart}'

class SpinStateScan(SimpleVasp):
    def __init__(self, ledger:ProcessLedger, comp:compound, JobName:str = "", AssignedServer:str = "ccp20,ccp22"):
        super().__init__(ledger, comp, JobName, AssignedServer)
        self.AllowedAlgos = [ "N","F", "VF", "A", "S"]

    def GetSpinCombos(self):
            b1HS, b1LS = self.AssignedComp.get_spinstate('B1')
            if b1HS == b1LS:
                B1s = [b1HS]


            else:
                B1s = [b1HS, b1LS]


            if self.AssignedComp.B2 == 'Vac':
                B2s = ['']
            else:

                b2HS, b2LS = self.AssignedComp.get_spinstate('B2')
                ##Filter out possible duplicate spinstates
                if b2HS == b2LS:
                    B2s = [b2HS]

                else:
                    B2s = [b2HS, b2LS]

            combos = []
            for b1 in B1s:
                for b2 in B2s:
                    combos.append([b1,b2])
                    if b2 != 0 and b2 != '':
                        combos.append([b1,-b2])

            return combos
    
    def Run(self):

        os.chdir(self.JobSpecifics['JobPath'])

        self.PrepareVASPRun()
        self.UsePOTCARENMAX(MultFactor =1.0, MinValue=self.JobBasics['std_incar']['ENCUT'])

        spin_combos = self.GetSpinCombos()

        spinDF = pd.DataFrame(columns=["muB1_in", "muB2_in", "SCF_convergence", "Etot_out", "MAGtot_out"],dtype=str)
        subrunIndex = 0 
        for combo in spin_combos:
            subrunIndex += 1
            ##Set the required MAGMOM initialization
            MAGMOM_line = f"{combo[0]} {combo[1]} 2*0.0 6*0.0"
            self.JobSpecifics['IncarExtras'].update({"MAGMOM": MAGMOM_line})
            #self.ledger.SetSingleValue(self.AssignedComp,self.JobName,'IncarExtras', self.JobSpecifics['IncarExtras'])

            #write INCAR
            self.WriteIncar()
            self.WriteSubmission(LogFileExtra = f'_MAGMOM_{combo[0]}_{combo[1]}')

            slurm_out = submit_and_wait(self.SubFile)

            if slurm_out != "COMPLETED":
                #some weird error, save slurm exit code to the completion field in the ledger
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_out}")
                return 'slurm'

            scf_out = self.check_scfcompletion()
            if scf_out == '' or scf_out == [''] or scf_out ==[-1]:
                scf_out = ['1']
                Eout = 100.
                MagOut = ''
            else:
                Eout = self.ExtractEtot()
                MagOut = self.ExtractMAGtot()

            outputDict = {"muB1_in": combo[0],
                            "muB2_in": combo[1], 
                            "SCF_convergence": scf_out[-1], 
                            "Etot_out": Eout, 
                            "MAGtot_out": MagOut
                            }
            
            spinDF.loc[subrunIndex] = spinDF.columns.map(outputDict)

        #safe the results of this analysis to a csv file
        spinDF.to_csv(f'../SpinResults_{self.AssignedCompName}_job{self.JobID}.csv')

        try:
            converged_runs = spinDF.groupby('SCF_convergence').get_group(0)
        except KeyError:
            #not a single run was converged, so we need to retry with a different algo
            self.CycleAlgos()
            return 'convergence_failed'
        
        else:
            self.ledger.SetSingleValue(self.AssignedComp,self.JobName,"completed", 1)

        spinDF.fillna('', inplace=True)
        optimal_run = converged_runs['Etot_out'].idxmin()
        magmom_use = f"{spinDF['muB1_in'].loc[optimal_run]} {spinDF['muB2_in'].loc[optimal_run]} 2*0.0 6*0.0"


        ##now we want to put this MAGMOM initialization into the IncarExtras for all remaining VASP jobs
        jobindex = self.ledger.Jobs.index(self.JobName)

        for futurejob in self.ledger.Jobs[jobindex+1:]:
            if self.ledger.JobInfo[futurejob]['job_type'].upper() == "VASP":
                IncarExtras = self.ledger.GetSingleJob(self.AssignedComp,futurejob)['IncarExtras']

                IncarExtras.update({"MAGMOM": magmom_use})
                self.ledger.SetSingleValue(self.AssignedComp,futurejob,"IncarExtras", IncarExtras)
        
        return 'good_run'

class PreLobster(SimpleVasp):
    def __init__(self, ledger:ProcessLedger, comp:compound, JobName:str = "", AssignedServer:str = "ccp20,ccp22"):
        super().__init__(ledger, comp, JobName, AssignedServer)
        if JobName == "4HSE":
            self.AllowedAlgos = ["D","A","N","S"]
        else:
            self.AllowedAlgos = ["N","VF","F","D","A", "S"]


    def Run(self):
            os.chdir(self.JobSpecifics['JobPath'])
            ##Preparation steps
            self.PrepareVASPRun() #assures KPOINTS, POSCAR, and POTCAR are there
            self.UsePOTCARENMAX(MultFactor=1.0, MinValue = self.JobBasics["std_incar"]["ENCUT"])
            
            NBANDS = self.AssignedComp.get_maxNBANDS()
            self.JobSpecifics['IncarExtras'].update({"NBANDS":NBANDS})
            self.ledger.SetSingleValue(self.AssignedComp, self.JobName, "IncarExtras", self.JobSpecifics['IncarExtras'])


            self.WriteIncar()
            self.WriteSubmission()

            slurm_out = submit_and_wait(self.SubFile)
            if slurm_out != "COMPLETED":
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_out}")
                return 'slurm'


            scf_out = self.check_scfcompletion()

            if scf_out[-1] == 0:
                #job succesfull, update this in ledger and return
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=1)
                return 'good_run'

            elif scf_out[-1] == 1:
                #convergence failed, so make useful changes
                self.CycleAlgos()
                return 'convergence_failed'
            

            elif scf_out[-1] == -1:
                #error occured, deal with it appropriate for this step
                self.CycleAlgos()
                return 'error'
            
            else:
                #something else went wrong in converging
                #convergence failed, so make useful changes
                self.CycleAlgos()
                return 'run_failed'



class SimpleLobster(GeneralJob):
    def __init__(self, ledger:ProcessLedger, comp:compound, JobName:str = "", AssignedServer:str = "ccp20,ccp22"):
        super().__init__(ledger, comp, JobName, AssignedServer)


    def ChangeCOHPErange(self,use_file="DOSCAR"):

        #os.chdir(self.JobSpecifics["JobPath"])

        with open(use_file,'r') as file:
            lines = file.readlines()

        energyline = lines[5]
        values = list(map(float,energyline.strip().split()))
        
        self.JobSpecifics['IncarExtras'].update({"COHPstartEnergy": values[1] - values[3], "COHPendEnergy":values[0] - values[3]})
        
        if self.AssignedComp.B2 == "Vac":
            self.JobSpecifics['IncarExtras'].update({"cohpGenerator from 0.5 to 4.": f"type {self.AssignedComp.B1} type {self.AssignedComp.X} orbitalwise "})
        else:
            self.JobSpecifics['IncarExtras'].update({"cohpGenerator from 0.5 to 4.": [f"type {self.AssignedComp.B1} type {self.AssignedComp.X} orbitalwise", f"type {self.AssignedComp.B2} type {self.AssignedComp.X} orbitalwise "]})

        return values[1], values[0]


    def add_basisfunctions(self,which = 'basis0'):
        #os.chdir(self.JobSpecifics['JobPath'])
        all_basis = self.AssignedComp.get_all_basis_combos()
        with open('lobsterin','a') as file:
            for el, funcs in all_basis[which].items():
                file.write(f'basisFunctions {str(el)} {str(funcs)}\n')

        return

    def check_charge_spilling(self, max_val = 3):

        awk_spil = subprocess.run(["awk","/abs. charge spilling/ {print $NF}", 'lobsterout'],stdout=subprocess.PIPE, text=True)
        spilling = awk_spil.stdout.strip().split('\n')
        
        bad_spilling = False
        for value in spilling:
            if float(value.strip('%')) > max_val:
                bad_spilling = True

        return bad_spilling

    def check_overlaps(self, max_val = 0.05):

        if os.path.exists('bandOverlaps.lobster'):
            awk_maxDev = subprocess.run(["awk",'/maxDeviation/ {print $NF}','bandOverlaps.lobster'], stdout = subprocess.PIPE, text=True)

            dev_vals = list(map(float,awk_maxDev.stdout.strip().split('\n')))

        else:
            dev_vals = [0]

        for val in dev_vals:
            if val > max_val:
                return True

        return False



    def Run(self):
        os.chdir(self.JobSpecifics['JobPath'])
        self.ChangeCOHPErange()
        os.system('export OMP_NUM_THREADS=2')
        #self.WriteIncar()
        #self.WriteSubmission()
        # print(self.ledger)
        self.ledger.SetSingleValue(self.AssignedComp,self.JobName, "IncarExtras", self.JobSpecifics["IncarExtras"])
        pos_bases = self.AssignedComp.get_all_basis_combos()

        with open(f'{self.AssignedCompName}_basisfunctions.json','w') as file:
            json.dump(pos_bases,file, indent=4)

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
            # self.WriteSubmission()

            # slurm_exit = submit_and_wait(self.SubFile)
# 
            # if slurm_exit != "COMPLETED":
                # some weird error, save slurm exit code to the completion field in the ledger
                # self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=f"slurm:{slurm_exit}")
                # return 'slurm'
            print(f"starting LOB for {self.AssignedCompName}")
            with open('extralob.log','w') as f:
                subprocess.run(['/home/lwalterb/Downloads/lobster-5.1.1/lobster-5.1.1'],shell=True,stderr=subprocess.STDOUT, stdout=f, check=True) 
            # time.sleep(40)

            spil_restart = self.check_charge_spilling()
            over_restart = self.check_overlaps()

            if spil_restart or over_restart:
                os.chdir(self.JobSpecifics['JobPath'])
                continue
            else:
                self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=1)
                return 'good_run'

        self.ledger.SetSingleValue(self.AssignedComp, self.JobName, Field="completed", value=-1)
        return 'no_good_run'








            



