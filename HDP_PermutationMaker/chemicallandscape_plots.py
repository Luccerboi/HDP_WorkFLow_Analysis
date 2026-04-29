# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 21:59:46 2024

@author: lucwa
"""

import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
import numpy as np
from mendeleev import element
from mendeleev.fetch import fetch_table


dposs = pd.read_csv('HDPcandidates_final_allordered.csv')
#dposs = pd.read_csv('HDPpossible_final_ordered.csv')


###create bar plot of different A sites, seperated in each X-site contribution
As = dposs['A'].unique()
Xs = dposs['X'].unique()


#find total counts
totals = []
for ii in range(len(As)):
    asite = As[ii]
    tot = len(dposs.query(f'A==\"{asite}\"'))
    totals.append({"A": asite, "tot":tot})
    
def myFunc(e):
    return e['tot']

totals.sort(key=myFunc, reverse=False)


countdict = {}
As_ordered = []
for xsite in Xs:
    counts = np.zeros(len(As))

    for ii in range(len(As)):
        asite = totals[ii]['A']
        
        if xsite == 'F':
            nA = '+' + str(dposs.query(f'A==\"{asite}\"')['nA'].iloc[0])
            aformatted = asite + '$^{%s}$'%(nA)
            
            
            As_ordered.append(aformatted)
        counts[ii] = len(dposs.query(f'A==\"{asite}\" & X==\"{xsite}\"'))
    
    countdict.update({xsite: counts})
    
#print(countdict)

fig, ax = plt.subplots(figsize=(13,9))
bottom = np.zeros(len(As))
ax.xaxis.grid()
for label, comp_count in countdict.items():
    p = ax.barh(As_ordered, comp_count, label=label, left=bottom)
    bottom += comp_count

ax.set_title(f'Considered HDP compositions per A-site and X-site ion, total: {len(dposs)}',fontsize=17)
ax.bar_label(p,fontsize=15)
ax.legend(['F$^-$','Cl$^-$','Br$^-$','I$^-$'],title='X-site choice:', loc='lower right',fontsize=15)
ax.set_xlabel('Number of considered HDP compositions',fontsize=16)
ax.set_ylabel('A-site cation',fontsize=16)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
fig.tight_layout()
plt.savefig('compsperaxallcomps', dpi=200)
plt.show()
    

#def extract_bpairing(x):
#    return element(x).block


#dposs['B1'].apply(extract_bpairing)

#dposs.query('A=="Cs"')['B1'].map(lambda x: element(x).block)

#['Fluorides','Chlorides','Bromides','Iodides']

dext = pd.read_csv('HDPpossible_extended.csv')

#define color scheme

# #allpairings = dext.groupby('Bpairing').count()['A'].keys()
# allpairings = ['Vac-d', 'Vac-f', 'Vac-p', 'Vac-s', 'd-d', 'd-f', 'd-p', 'd-s', 'f-f', 'f-p', 'f-s', 'f-dummy','p-p', 'p-s', 'p-dumm1y', 'p-dummy2',  's-s']


# def filter_pairings(pairings,starting_letter:str):
#     ##helper function to filter all possible B-B' pairings on the starting letter
#     filtered_pairings = []
#     for pair in pairings:
#         if pair.startswith(starting_letter):
#             filtered_pairings.append(pair)
#     return filtered_pairings


# dpairs = filter_pairings(allpairings,'d')
# spairs = filter_pairings(allpairings,'s')
# ppairs = filter_pairings(allpairings,'p')
# fpairs = filter_pairings(allpairings,'f')
# Vacpairs = filter_pairings(allpairings, 'V')

# colors = { label:color for label, color in zip(allpairings, colormaps['tab20c'].colors)}
# #colors.update({ label:color for label, color in zip(fpairs, colormaps['Reds'].colors)})
# #colors.update({ label:color for label, color in zip(ppairs, colormaps['Greens'].colors)})
# #colors.update({ label:color for label, color in zip(spairs, colormaps['Oranges'].colors)})
# #colors.update({ label:color for label, color in zip(Vacpairs, colormaps['Blues'].colors)})


# ###Make figure of bpairings for CsBBX for all 4 Xs
# plt.close('all')
# fig, axes = plt.subplots(2,2, figsize=(10,10))


# ax = axes[0,0]
# paircounts = dext.query('A=="Cs" & X=="F"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Cs_2 BB\' F_6$')

# ax = axes[0,1]
# paircounts = dext.query('A=="Cs" & X=="Cl"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Cs_2 BB\' Cl_6$')

# ax = axes[1,0]
# paircounts = dext.query('A=="Cs" & X=="Br"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Cs_2 BB\' Br_6$')

# ax = axes[1,1]
# paircounts = dext.query('A=="Cs" & X=="I"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Cs_2 BB\' I_6$')

# fig.suptitle('Count of different B-B\' pairings for $Cs_2BB\'X_6$ ')
# fig.tight_layout()
# plt.savefig('bpairingpie_Cs__X.png', dpi=200)
# #plt.show()
# plt.close()




# ####Same for Rb
# fig, axes = plt.subplots(2,2, figsize=(10,10))


# ax = axes[0,0]
# paircounts = dext.query('A=="Rb" & X=="F"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(),autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Rb_2 BB\' F_6$')

# ax = axes[0,1]
# paircounts = dext.query('A=="Rb" & X=="Cl"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Rb_2 BB\' Cl_6$')

# ax = axes[1,0]
# paircounts = dext.query('A=="Rb" & X=="Br"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Rb_2 BB\' Br_6$')

# ax = axes[1,1]
# paircounts = dext.query('A=="Rb"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Rb_2 BB\' X_6$')

# fig.suptitle('Count of different B-B\' pairings for $Rb_2BB\'X_6$ ')
# fig.tight_layout()
# plt.savefig('bpairingpie_Rb__X.png', dpi=200)
# #plt.show()
# plt.close()





# ####same but different entries
# fig, axes = plt.subplots(2,2, figsize=(10,10))


# ax = axes[0,0]
# paircounts = dext.query('A=="Na" ').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Na_2 BB\' X_6$')

# ax = axes[0,1]
# paircounts = dext.query('A=="Au" & X=="F"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Au_2 BB\' F_6$')

# ax = axes[1,0]
# paircounts = dext.query('A=="Ag" & X=="F"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Ag_2 BB\' F_6$')

# ax = axes[1,1]
# paircounts = dext.query('A=="Hg" & X=="F"').groupby('Bpairing').count()['A']

# ax.pie(paircounts, labels=paircounts.keys(), autopct='%1.1f%%', pctdistance=0.9, labeldistance=1.15, colors=[colors[k] for k in paircounts.keys()])
# ax.set_title('$Hg_2 BB\' F_6$')

# fig.suptitle('Count of different B-B\' pairings for $A_2BB\'F_6$ ')
# fig.tight_layout()
# plt.savefig('bpairingpie_altAsites.png',dpi=200)
# #plt.show()
# plt.close()


#fig, ax = plt.subplots()
#dbar = dext.query('A=="Cs"')[['X', 'Bpairing']].groupby('X').value_counts().sort_index(axis=0)
#dbar.plot.bar(subplots=True, color=[colors[k] for k in dbar.keys().get_level_values(1)])



def bpairing_barplot(bar_df, *, title='Counts of different B-B\' block pairings', fn=''):
    allpairings = ['Vac-d', 'Vac-f', 'Vac-p', 'Vac-s', 'd-d', 'd-f', 'd-p', 'd-s', 'f-f', 'f-p', 'f-s', 'f-dummy','p-p', 'p-s', 'p-dumm1y', 'p-dummy2',  's-s']
    colors = { label:color for label, color in zip(allpairings, colormaps['tab20c'].colors)}
    
    species = list(bar_df.keys())
    #paircounts = bar_dict.values()
    
    def return_paircounts(bar_df,pair):
        list_of_counts = []
        for key in bar_df.keys():
            list_of_counts.append(bar_df[key][pair])
        return list_of_counts
    
    fig,ax = plt.subplots(figsize=(10,6))
    width = 0.2
    offset = 0
    
    
    maxpairs = bar_df[species[0]].keys()
    x = 3*np.arange(len(species))
    multiplier = 0
        
    for pair in maxpairs:
        offset = width* multiplier
        p = ax.bar(x+offset, return_paircounts(bar_df,pair),width=width , label = pair, color=colors[pair])
        ax.bar_label(p,fontsize=12)
        multiplier += 1
        
    ax.set_xticks(x + len(maxpairs)/2*width,species, fontsize=14)
    ax.legend(ncols=2, loc='upper right')
    ax.yaxis.grid()
    
    ax.set_title(title,fontsize=16)
    ax.set_xlabel('Compound group',fontsize=14)
    ax.set_ylabel('Count',fontsize=14)
    fig.tight_layout()
    #plt.show()
    
    if fn!='':
        plt.savefig(fn, dpi=200)
        plt.close()
    else:
        plt.show()
   
    
   
    
###select data to plot
CsF = dext.query('A=="Cs" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
CsCl = dext.query('A=="Cs" & X == "Cl"')['Bpairing'].value_counts().sort_index(axis=0)
CsBr = dext.query('A=="Cs" & X == "Br"')['Bpairing'].value_counts().sort_index(axis=0)
CsI = dext.query('A=="Cs" & X == "I"')['Bpairing'].value_counts().sort_index(axis=0)

bar_df = pd.DataFrame({f'$Cs_2[BB\']F_6$\ntot:{CsF.sum()}' : CsF, f'$Cs_2[BB\']Cl_6$\ntot:{CsCl.sum()}': CsCl, f'$Cs_2[BB\']Br_6$\ntot:{CsBr.sum()}': CsBr, f'$Cs_2[BB\']I_6$\ntot:{CsI.sum()}': CsI}, index=CsF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for $Cs_2BB\'X_6$ \n total:{CsF.sum() + CsCl.sum() + CsBr.sum() + CsI.sum()}',fn='BBpair_bar_CsX.png')




###select data to plot
RbF = dext.query('A=="Rb" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
RbCl = dext.query('A=="Rb" & X == "Cl"')['Bpairing'].value_counts().sort_index(axis=0)
RbBr = dext.query('A=="Rb" & X == "Br"')['Bpairing'].value_counts().sort_index(axis=0)

bar_df = pd.DataFrame({f'$Rb_2[BB\']F_6$\ntot:{RbF.sum()}' : RbF, f'$Rb_2[BB\']Cl_6$\ntot:{RbCl.sum()}': RbCl, f'$Rb_2[BB\']Br_6$\ntot:{RbBr.sum()}': RbBr})#, f'$Rb_2[BB\']I_6$\ntot:{CsI.sum()}': CsI}, index=CsF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for $Rb_2BB\'X_6$ \n total:{RbF.sum() + RbCl.sum() + RbBr.sum()}',fn='BBpair_bar_RBX.png')



###select data to plot
TlF = dext.query('A=="Tl" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
TlCl = dext.query('A=="Tl" & X == "Cl"')['Bpairing'].value_counts().sort_index(axis=0)
TlBr = dext.query('A=="Tl" & X == "Br"')['Bpairing'].value_counts().sort_index(axis=0)

bar_df = pd.DataFrame({f'$Tl_2[BB\']F_6$\ntot:{TlF.sum()}' : TlF, f'$Tl_2[BB\']Cl_6$\ntot:{TlCl.sum()}': TlCl, f'$Tl_2[BB\']Br_6$\ntot:{TlBr.sum()}': TlBr})#, f'$Cs_2[BB\']I_6$\ntot:{CsI.sum()}': CsI}, index=CsF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for $Tl_2BB\'X_6$ \n total:{TlF.sum() + TlCl.sum() + TlBr.sum()}',fn='BBpair_bar_TlX.png')



###select data to plot
FrF = dext.query('A=="Fr" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
FrCl = dext.query('A=="Fr" & X == "Cl"')['Bpairing'].value_counts().sort_index(axis=0)
FrBr = dext.query('A=="Fr" & X == "Br"')['Bpairing'].value_counts().sort_index(axis=0)

bar_df = pd.DataFrame({f'$Fr_2[BB\']F_6$\ntot:{FrF.sum()}' : FrF, f'$Fr_2[BB\']Cl_6$\ntot:{FrCl.sum()}': FrCl, f'$Fr_2[BB\']Br_6$\ntot:{FrBr.sum()}': FrBr})#, f'$Fr_2[BB\']I_6$\ntot:{FrI.sum()}': FrI}, index=FrF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for $Fr_2BB\'X_6$ \n total:{FrF.sum() + FrCl.sum() + FrBr.sum()}',fn='BBpair_bar_FrX.png')



    
###select data to plot
NaF = dext.query('A=="Na" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
AuF = dext.query('A=="Au" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
AgF = dext.query('A=="Ag" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
HgF = dext.query('A=="Hg" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)

bar_df = pd.DataFrame({f'$Na_2[BB\']F_6$\ntot:{NaF.sum()}' : NaF, f'$Au_2[BB\']Cl_6$\ntot:{AuF.sum()}': AuF, f'$Ag_2[BB\']F_6$\ntot:{AgF.sum()}': AgF, f'$Hg_2[BB\']F_6$\ntot:{HgF.sum()}': HgF}, index=CsF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for different $A_2BB\'F_6$ \nWith A(+1) ', fn='BBpair_bar_AltA1s.png')


    
###select data to plot
BaF = dext.query('A=="Ba" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
BaCl = dext.query('A=="Ba" & X == "Cl"')['Bpairing'].value_counts().sort_index(axis=0)
RaF = dext.query('A=="Ra" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
RaCl = dext.query('A=="Ra" & X == "Cl"')['Bpairing'].value_counts().sort_index(axis=0)
RaBr = dext.query('A=="Ra" & X == "Br"')['Bpairing'].value_counts().sort_index(axis=0)

bar_df = pd.DataFrame({f'$Ba_2[BB\']F_6$\ntot:{BaF.sum()}' : BaF, f'$Ba_2[BB\']Cl_6$\ntot:{BaCl.sum()}': BaCl, f'$Ra_2[BB\']F_6$\ntot:{RaF.sum()}' : RaF, f'$Ra_2[BB\']Cl_6$\ntot:{RaCl.sum()}': RaCl, f'$Ra_2[BB\']Br_6$\ntot:{RaBr.sum()}': RaBr}, index=RaF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for $A_2BB\'X_6$ \nWith A(+2) ',fn='BBpair_bar_earthalkali.png')

###select data to plot
PbF = dext.query('A=="Pb" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
SrF = dext.query('A=="Sr" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
NdF = dext.query('A=="Nd" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)
EuF = dext.query('A=="Eu" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0)

#f'$Cs_2[BB\']F_6$\ntot:{CsF.sum()}' : CsF
#CsF = dext.query('A=="Cs" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0) 
CaF = dext.query('A=="Ca" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0) 
CdF = dext.query('A=="Cd" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0) 
AmF = dext.query('A=="Am" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0) 
SmF = dext.query('A=="Cd" & X == "F"')['Bpairing'].value_counts().sort_index(axis=0) 

bar_df = pd.DataFrame({f'$Pb_2[BB\']F_6$\ntot:{PbF.sum()}' : PbF, f'$Sr_2[BB\']F_6$\ntot:{SrF.sum()}': SrF, f'$Nd_2[BB\']F_6$\ntot:{NdF.sum()}': NdF, f'$Eu_2[BB\']F_6$\ntot:{EuF.sum()}': EuF,
                       f'$Ca_2[BB\']F_6$\ntot:{CaF.sum()}' : CaF,f'$Cd_2[BB\']F_6$\ntot:{CdF.sum()}' : CdF, f'$Am_2[BB\']F_6$\ntot:{AmF.sum()}' : AmF, f'$Sm_2[BB\']F_6$\ntot:{SmF.sum()}' : SmF}, index=PbF.keys()).fillna(0)

bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for different $A_2BB\'F_6$ \nWith A(+2) ', fn='BBpair_bar_AltA2s.png')

    
bar_df = pd.DataFrame({f'$Ra_2[BB\']F_6$\ntot:{RaF.sum()}' : RaF, f'$Ba_2[BB\']F_6$\ntot:{BaF.sum()}' : BaF, f'$Pb_2[BB\']F_6$\ntot:{PbF.sum()}' : PbF,f'$Cd_2[BB\']F_6$\ntot:{CdF.sum()}' : CdF,f'$Sm_2[BB\']F_6$\ntot:{SmF.sum()}' : SmF})
bpairing_barplot(bar_df, title=f'Counts of different B-B\' block pairings for different $A_2BB\'F_6$ \nWith A(+2) ', fn='BBpair_bar_A2report.png')