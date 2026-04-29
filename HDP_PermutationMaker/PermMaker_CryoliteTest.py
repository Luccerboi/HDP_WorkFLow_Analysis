# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 17:09:00 2024

@author: lucwa
"""

from mendeleev.fetch import fetch_ionic_radii
import numpy as np
import pandas as pd
import csv
from mendeleev import element


def calculate_tau(nA, r_A, r_B, r_X):
    tau = (r_X/r_B) - nA*(nA - (r_A/r_B)/np.log(r_A/r_B))
    return tau


def compare_ens(en1,en2):
    rat = en1/en2
    if rat < 0.9:
        return 3,1
    elif rat > 1.1:
        return 1,3
    else:
        return 2, 2

#def append_compounds(filename, an_a, nA, an_b1, nB1, an_b2, nB2, halide_irs):
#    halide_names = ('F', 'Cl', 'Br', 'I')
    
    
    

tau_cutoff = 4.18

irs_full = fetch_ionic_radii(radius='ionic_radius')[['I', 'II', 'III', 'IV',  'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIV']].interpolate(axis=1, limit_direction='both').dropna(axis=0)
irs_full = irs_full.drop([1,1])
##Seperate periodic table into relevant charges
chg1 = irs_full.groupby(['charge']).get_group(1)
#chg1 = chg1.iloc[:3,:]

#chg1.drop()
#chg1 = chg1.drop([1,1])
chg2 = irs_full.groupby(['charge']).get_group(2)
chg3 = irs_full.groupby(['charge']).get_group(3)

bsites = pd.concat([chg1,chg2,chg3]).sort_values(by='atomic_number')
#+1, +2, +3 are the groups that can form multiple combinations in B B'

chg4 = irs_full.groupby(['charge']).get_group(4)
halides = irs_full.groupby(['charge']).get_group(-1)

##safe halide information for easier use
halide_names = ('F', 'Cl', 'Br', 'I')
halide_irs = halides['VI'].values
nX = -1

##get column indices of relevant CNs
index6 = irs_full.columns.get_loc('VI')
index12 = irs_full.columns.get_loc('XII')

###If 2 charge combinations are possible we need a electronegativity (EN) scale
en_scale = 'mulliken'


###we write all new entries straight to one file

filename = 'HDPPermutations_IncludingCryolites.csv'
with open(filename, 'w') as file:
    file.write('A2BBX6,A,B1,B2,X,nA,nB1,nB2,nX,rA,rB1,rB2,rBavg,rX,tau\n')

i_asite = 0
tot_asites = len(chg1.index.get_level_values('atomic_number').unique()) + len(chg2.index.get_level_values('atomic_number').unique())

for iter_a1 in chg1.index.get_level_values('atomic_number').unique():
    print(f'treating A-site {i_asite}/{tot_asites}')
    i_asite += 1
    
    ir_a = irs_full.loc[iter_a1,1]['XII']
    nA = 1
    
    jj = 0
    bsits = len(bsites.index.get_level_values('atomic_number').unique())
    for iter_b1 in bsites.index.get_level_values('atomic_number').unique():
        jj+=1
        if jj%10==0:
            print(f'B1-site:{jj}/{bsits}')
        
        for iter_b2 in bsites.loc[iter_b1:,:].index.get_level_values('atomic_number').unique():
            ##loop over all B B' combinations
            
            
            # if iter_a1 == iter_b1  or iter_a1 == iter_b2:
            #     #Same element cannot occupy A site and a B site
            #     continue
            
            perov_bool = False
            elpas_bool = False
            
            #check available charges
            chgs_b1 = bsites.loc[iter_b1,:].index
            chgs_b2 = bsites.loc[iter_b2,:].index
            
            
            #check if perovskite or elpasolite is possible
            if 2 in chgs_b1 and 2 in chgs_b2 and iter_b1 != iter_b2:
                perov_bool = True
                nB1 = 2
                nB2 = 2
                
            if (1 in chgs_b1 and 3 in chgs_b2):
                elpas_bool = True
                nB1 = 1
                nB2 = 3
            
            elif (3 in chgs_b1 and 1 in chgs_b2):
                elpas_bool = True
                nB1 = 3
                nB2 = 1
                
            ##If both charge buildups are possible check ENs
            if perov_bool and elpas_bool:
                #print('checking ENs')
                en_b1 = element(iter_b1).electronegativity(en_scale)
                en_b2 = element(iter_b2).electronegativity(en_scale)
                
                nB1, nB2 = compare_ens(en_b1, en_b2)
            
            #if no charge combination possible continue
            if not perov_bool and not elpas_bool:
                continue
            
            ###With charges determined calculate tau based on those irs
            ir_b1 = irs_full.loc[iter_b1,nB1]['VI']
            ir_b2 = irs_full.loc[iter_b2,nB2]['VI']
            
            ir_bavg = (ir_b1+ir_b2)/2
            if ir_bavg >= ir_a:
                continue
            
            taus = calculate_tau(nA,ir_a, ir_bavg,halide_irs)
            
            nameA=element(iter_a1).symbol
            nameB1=element(iter_b1).symbol
            nameB2=element(iter_b2).symbol
            
            with open(filename, 'a') as file:
                for ii in range(len(taus)):
                    cname = nameA + '2' + nameB1 + nameB2 + halide_names[ii] +'6'
                    file.write(f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n")
            
    for iter_b1 in chg4.index.get_level_values('atomic_number').unique():
        #calculate vacancy ordered possibilities
        nB1 = 4
        nB2 = 0
        
        ir_bavg = irs_full.loc[iter_b1,nB1]['VI']
        if ir_bavg >= ir_a:
            continue
        
        taus = calculate_tau(nA,ir_a, ir_bavg,halide_irs)
        
        nameA=element(iter_a1).symbol
        nameB1=element(iter_b1).symbol
        nameB2='Vac'
        
        with open(filename, 'a') as file:
            for ii in range(len(taus)):
                cname = nameA + '2' + nameB1 + halide_names[ii] +'6'
                file.write(f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n")        


for iter_a1 in chg2.index.get_level_values('atomic_number').unique():
    ir_a = irs_full.loc[iter_a1,2]['XII']
    nA = 2
    print(f'treating A-site {i_asite}/{tot_asites}')
    i_asite += 1
    
    jj = 0
    for iter_b1 in chg1.index.get_level_values('atomic_number').unique():
        nB1 = 1
        jj+=1
        if jj%10==0:
            print(f'B1-site:{jj}/{bsits}')
        
        
        for iter_b2 in chg1.loc[iter_b1:,:].index.get_level_values('atomic_number').unique():
                nB2 = 1
                ##Making A(2)B(1)B(1)X(-1) combinations
                if iter_b1 == iter_b2:
                    continue
                
                ir_b1 = irs_full.loc[iter_b1,nB1]['VI']
                ir_b2 = irs_full.loc[iter_b2,nB2]['VI']
                
                ir_bavg = (ir_b1+ir_b2)/2
                if ir_bavg >= ir_a:
                    continue
                
                taus = calculate_tau(nA,ir_a, ir_bavg,halide_irs)
                
                nameA=element(iter_a1).symbol
                nameB1=element(iter_b1).symbol
                nameB2=element(iter_b2).symbol
                
                with open(filename, 'a') as file:
                    for ii in range(len(taus)):
                        cname = nameA + '2' + nameB1 + nameB2 + halide_names[ii] +'6'
                        file.write(f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n")
         
    for iter_b1 in chg2.index.get_level_values('atomic_number').unique():
        #calculate vacancy ordered possibilities
        nB1 = 2
        nB2 = 0
        
        ir_bavg = irs_full.loc[iter_b1,nB1]['VI']
        if ir_bavg >= ir_a:
            continue
        
        taus = calculate_tau(nA,ir_a, ir_bavg,halide_irs)
        
        nameA=element(iter_a1).symbol
        nameB1=element(iter_b1).symbol
        nameB2='Vac'
        
        with open(filename, 'a') as file:
            for ii in range(len(taus)):
                cname = nameA + '2' + nameB1 + halide_names[ii] +'6'
                file.write(f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n") 
            
#file.close()            
    

