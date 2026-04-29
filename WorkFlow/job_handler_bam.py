#!/software/compilers/intel/AIToolKit/intelpython/latest/bin/python3

import numpy as np
import pandas as pd
import os
import subprocess
import time
#import queue
from math import floor, ceil 
import sys
import CsBBCl_HT as ht
from collections import deque

def submit_and_wait(submission="vasp.sub"):
    pass 

def check_availability():
    pass

def write_submission(comp,time='02:00:00',nnodes=4, submission='vasp.sub'):
    
    with open(submission,'w') as file:
        file.write(f'#!/bin/bash \n#SBATCH --time={time}\n')
        file.write(f'#SBATCH -p {comp.assigned_server} \n#SBATCH -N {nnodes} --exclusive --ntasks-per-node 32\n')
        
        file.write(f'#SBATCH -J {comp.calc_step}{comp.name} \n#SBATCH --output={comp.cur_log()}\n')
        file.write(f'#SBATCH --ntasks-per-core=1 \n#SBATCH --error={comp.cur_log()}\n\n')
        file.write('module purge \nmodule load vasp/6.3.0 \n\nulimit -s unlimited\n')
        file.write('startingtime=$(date)\n')
        file.write("now=$(date +'%s') \nmpirun vasp_std \nendtime=$(date +'%s') \nduration=$(($endtime-$now)) \n")
        file.write("echo 'Started on: ' $startingtime 'calculation ran for: ' $duration '[s]'\n")

    return       


  
def write_lobsub(comp,time='01:00:00',nnodes=1, submission='lobster.sub'):
    
    with open(submission,'w') as file:
        file.write(f'#!/bin/bash \n#SBATCH --time={time}\n')
        file.write(f'#SBATCH -p {comp.assigned_server} \n#SBATCH -N {nnodes} --exclusive --ntasks-per-node 32\n')
        
        file.write(f'#SBATCH -J {comp.calc_step}{comp.name} \n#SBATCH --output={comp.cur_log()}\n')
        file.write(f'#SBATCH --ntasks-per-core=1 \n#SBATCH --error={comp.cur_log()}\n\n')
        file.write('module purge \nmodule load lobster \n\nulimit -s unlimited\n')
        file.write('startingtime=$(date)\n')
        file.write("now=$(date +'%s') \n/software/codes/lobster/5.1.1/bin/lobster \nendtime=$(date +'%s') \nduration=$(($endtime-$now)) \n")
        file.write("echo 'Started on: ' $startingtime 'calculation ran for: ' $duration '[s]'\n")

    return

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
        p1 = subprocess.run(["awk", 'BEGIN{rv=0; print "\\nNew run started --->\\n\\n" >> "%s"};/self-consistency was not achieved/ {rv=1}; /F=/ {print rv; rv=0}; \
                            /E{7} * R * R * R * R * O{7}/ {print "-1"}; /reached required accuracy - stopping structural energy minimisation/ {print "2"};\
                            /E{7} * R * R * R * R * O{7}/, /---->/ {if(!/E{7} * R/ && !/---->/)print $0 >> "%s"}'%(erroroutput, erroroutput),\
                            logfile], stdout=subprocess.PIPE, text=True)

        numstrings = p1.stdout.strip().split('\n')
        try:
           scf_cycles = list(map(int,numstrings))
        except ValueError as err:
            #raise Exception("Unexpected value encountered when checking SCF convergence:{err}")
            return numstrings
        return scf_cycles   

class CompQueue(deque):
    """This class creates a list like object with an extra function to rewrite the queue
    according to the way i like it"""
    
    def __init__(self,iterable, nsteps):
        super().__init__(iterable)
        self.nsteps = nsteps
    
    def rewrite(self, full_list, banlist):
        self.clear()
        for n in np.arange(self.nsteps,0,-1):
            partiallist = ht.filter_calcstep(full_list, n)
            #print(f'Rewriting queue with {len(partiallist)} compounds at step {n}\n')
            for el in partiallist:
                if el in banlist:
                    continue
                else:
                    self.append(el)
        return
    
    def from_file(self,filename):
        pass
    
    
