#!/software/compilers/intel/AIToolKit/intelpython/latest/bin/python3
import numpy as np
from pymatgen.electronic_structure.dos import Dos, CompleteDos, _lobster_orb_labs, LobsterCompleteDos
from pymatgen.core.periodic_table import Element, Species

from pymatgen.core import Structure, get_el_sp
from pymatgen.core.spectrum import Spectrum
from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.util.coord import get_linear_interpolated_value
from pymatgen.util.typing import SpeciesLike
    
def plot_pdos(d:dict, comp_info:dict, efermi: float,filename,*,xmax=10,ymax=15, magmom:dict={}):
    import matplotlib.pyplot as plt # type: ignore
    #import matplotlib # type: ignore
    #matplotlib.use('agg')
    import scienceplots  # type: ignore
    
    plt.style.use(['science','no-latex'])
    x = d['energy']
    

    
    colors = ('b','r','g','m')  
    if comp_info['B2']=='Vac':
        sites = ['B1','A','X']
        nplots = 4
    else:
        sites = ['B1','B2','A','X']
        nplots = 5
    
    
    fig, axes = plt.subplots( nplots, 1, figsize=(10,20), layout='constrained')
    ax = axes[0]
    
    
    ax.fill_between(x,0, d['tdos'][Spin.up],color='grey',alpha = 0.3, label='t-DOS')
    ax.fill_between(x,0, d['tdos'][Spin.down],color='grey',alpha = 0.3)
    #B1
    ax.plot(x,d['B1'][Spin.up],color=colors[0],alpha=0.8,linewidth=0.5)
    ax.plot(x,d['B1'][Spin.down],color=colors[0],alpha=0.8,linewidth=0.5)
    ax.fill_between(x,0,d['B1'][Spin.up],color=colors[0],alpha=0.1,label=comp_info['B1'])
    ax.fill_between(x,0,d['B1'][Spin.down],color=colors[0],alpha=0.1)
    #B2
    if comp_info['B2'] == 'Vac':
        pass
    else:
        ax.plot(x,d['B2'][Spin.up],color=colors[1],alpha=0.8,linewidth=0.5)
        ax.plot(x,d['B2'][Spin.down],color=colors[1],alpha=0.8,linewidth=0.5)
        ax.fill_between(x,0,d['B2'][Spin.up],color=colors[1],alpha=0.1,label=comp_info['B2'])
        ax.fill_between(x,0,d['B2'][Spin.down],color=colors[1],alpha=0.1)
    
    #A-site
    ax.plot(x,d['A'][Spin.up],color=colors[2],alpha=0.8,linewidth=0.5)
    ax.plot(x,d['A'][Spin.down],color=colors[2],alpha=0.8,linewidth=0.5)
    ax.fill_between(x,0,d['A'][Spin.up],color=colors[2],alpha=0.1,label=comp_info['A'])
    ax.fill_between(x,0,d['A'][Spin.down],color=colors[2],alpha=0.1)
    
    #x-SITE
    ax.plot(x,d['X'][Spin.up],color=colors[3],alpha=0.8,linewidth=0.5)
    ax.plot(x,d['X'][Spin.down],color=colors[3],alpha=0.5,linewidth=0.5)
    ax.fill_between(x,0,d['X'][Spin.up],color=colors[3],alpha=0.1,label=comp_info['X'])
    ax.fill_between(x,0,d['X'][Spin.down],color=colors[3],alpha=0.1)
    
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([-xmax,min(xmax,x[-1])])
    ax.vlines(0,-ymax,ymax,ls='--',color='black',alpha=0.7)
    ax.tick_params(axis='x',which='major',labelsize=13)
    ax.set_yticks([ymax/2,-ymax/2],['spin-up','spin-down'],rotation='vertical',fontsize=13,verticalalignment='center')
    #ax.set_yticklabels(['spin-up','spin-down'],rotation='vertical')
    ax.set_xlabel('$E-E_F$ [eV]',fontsize=14)
    ax.set_ylabel('Density of States [a.u.]',fontsize=14)
    ax.set_title('Total DOS',fontsize=16)
    ax.set_title('A)',fontsize=16,loc='left')
    ax.legend()
    
    
    total_states = np.sum(d['tdos'][Spin.up])+ np.sum(np.abs(d['tdos'][Spin.down]))
    ymax = ymax/2
    plot_labels= ('B','C','D', 'E')
    lo = 1 #label offset, to make sure the subplots get the rigth label
    for sp in range(1,len(axes)):
        ax = axes[sp]
        s = sites[sp-1]
        #if comp_info[s]=='Vac':
        #    lo += 1
        #    continue
        
        ax.fill_between(x,0, d['tdos'][Spin.up],color='grey',alpha = 0.3)
        ax.fill_between(x,0, d['tdos'][Spin.up],color='grey',alpha = 0.3)
        
        bands = ('s','p','d','f')
        marks = ('//////','\\\\','......', 'oo' )
        dup = d[s][Spin.up]
        ddown = d[s][Spin.up]
        for ii in range(len(bands)):
            
            band_contribution = np.sum(dup[ii,:])+ np.sum(np.abs(ddown[ii,:]))
            if band_contribution/total_states <= 0.01:
                continue
            
            ax.plot(x, dup[ii,:],color=colors[sp-1],alpha=0.8,linewidth=0.5)
            ax.plot(x, ddown[ii,:],color=colors[sp-1],alpha=0.8,linewidth=0.5)
            ax.fill_between(x, dup[ii,:],color=colors[sp-1],alpha=0.4,hatch=marks[ii],label=(comp_info[s]+f"-({bands[ii]})"))
            ax.fill_between(x, ddown[ii,:],color=colors[sp-1],alpha=0.4,hatch=marks[ii])
        
        ax.set_ylim([-ymax,ymax])
        ax.set_xlim([-xmax,min(xmax,x[-1])])
        ax.vlines(0,-ymax,ymax,ls='--',color='black',alpha=0.7)
        ax.tick_params(axis='x',which='major',labelsize=13)
        ax.set_yticks([ymax/2,-ymax/2],['spin-up','spin-down'],fontsize=14, rotation='vertical',verticalalignment='center')
        #ax.set_yticklabels(['spin-up','spin-down'],rotation='vertical')
        ax.set_xlabel('$E-E_F$ [eV]',fontsize=14)
        ax.set_ylabel('Density of States [a.u.]',fontsize=14)
        if len(magmom) == 0:
            title_string = f"{s}-site ({comp_info[s]}) pDOS"
        else:
            try:
                title_string = f"{s}-site ({comp_info[s]}) pDOS, $\\mu({s})={magmom[s]:.2f} [\\mu_B$]"
            except:
                title_string = f"{s}-site ({comp_info[s]}) pDOS"
        ax.set_title(title_string,fontsize=16)
        ax.set_title(f'{plot_labels[sp-lo]})',fontsize=16,loc='left')
        ax.legend()
        
    fig.suptitle(f'{comp_info["A"]+ "2" +comp_info["B1"]+comp_info["B2"]+comp_info["X"] + "6"} projected DOS\n',fontsize=20, y=1.035)
    if len(magmom) == 0:
        fig.text(0.5,1.01, f'$E_F$={efermi:.2f} [eV]',fontsize=15,ha='center',va='baseline')
    else:
        fig.text(0.5,1.01, f'$E_F$={efermi:.2f} [eV],'+ " $\\mu_{tot}$" + f'={magmom["total"]:.2f} [$\\mu_B$]',fontsize=15,ha='center',va='baseline')
    
    
    plt.savefig(filename,dpi=250)
    #fig.save('trial.png')
    #plt.show()

def plot_tdos(d:dict, comp_info:dict, efermi: float, filename,*,xmax=10,ymax=15, magmom:dict={}):
    import matplotlib.pyplot as plt # type: ignore
    
    #matplotlib.use('agg')
    import scienceplots  # type: ignore
    
    plt.style.use(['science','grid'])
    x = d['energies'] - efermi 
    element_map = comp_info['input_species']
    #symmetric limits to the axes

    
    colors = ('b','r','g','m')
    fig, axes = plt.subplots( 1, 1, figsize=(10,5), layout='constrained')
    ax = axes
    
    
    ax.fill_between(x,0, d['tdos'][Spin.up],color='grey',alpha = 0.3, label='t-DOS')
    ax.fill_between(x,0, d['tdos'][Spin.down]*(-1),color='grey',alpha = 0.3)
    #B1
    ax.plot(x,d['B1']['densities'][Spin.up],color=colors[0],alpha=0.8,linewidth=0.5)
    ax.plot(x,d['B1']['densities'][Spin.down]*(-1),color=colors[0],alpha=0.8,linewidth=0.5)
    ax.fill_between(x,0,d['B1']['densities'][Spin.up],color=colors[0],alpha=0.1,label=str(element_map['B1'].element))
    ax.fill_between(x,0,d['B1']['densities'][Spin.down]*(-1),color=colors[0],alpha=0.1)
    #B2
    if element_map['B2'] == None:
        pass
    else:
        ax.plot(x,d['B2']['densities'][Spin.up],color=colors[1],alpha=0.8,linewidth=0.5)
        ax.plot(x,d['B2']['densities'][Spin.down]*(-1),color=colors[1],alpha=0.8,linewidth=0.5)
        ax.fill_between(x,0,d['B2']['densities'][Spin.up],color=colors[1],alpha=0.1,label=str(element_map['B2'].element))
        ax.fill_between(x,0,d['B2']['densities'][Spin.down]*(-1),color=colors[1],alpha=0.1)
    
    #A-site
    ax.plot(x,d['A']['densities'][Spin.up],color=colors[2],alpha=0.8,linewidth=0.5)
    ax.plot(x,d['A']['densities'][Spin.down]*(-1),color=colors[2],alpha=0.8,linewidth=0.5)
    ax.fill_between(x,0,d['A']['densities'][Spin.up],color=colors[2],alpha=0.1,label=str(element_map['A'].element))
    ax.fill_between(x,0,d['A']['densities'][Spin.down]*(-1),color=colors[2],alpha=0.1)
    
    #x-SITE
    ax.plot(x,d['X']['densities'][Spin.up],color=colors[3],alpha=0.8,linewidth=0.5)
    ax.plot(x,d['X']['densities'][Spin.down]*(-1),color=colors[3],alpha=0.5,linewidth=0.5)
    ax.fill_between(x,0,d['X']['densities'][Spin.up],color=colors[3],alpha=0.1,label=str(element_map['X'].element))
    ax.fill_between(x,0,d['X']['densities'][Spin.down]*(-1),color=colors[3],alpha=0.1)
    
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([-xmax,min(xmax,x[-1])])
    ax.vlines(0,-ymax,ymax,ls='--',color='black',alpha=0.7)
    ax.tick_params(axis='x',which='major',labelsize=12)
    ax.set_yscale("linear")
    ax.set_yticks([ymax/2,-ymax/2],['spin-up','spin-down'],rotation='vertical',fontsize=13,verticalalignment='center')
    #ax.set_yticklabels(['spin-up','spin-down'],rotation='vertical')
    ax.set_xlabel('$E-E_F$ [eV]',fontsize=14)
    ax.set_ylabel('Density of States [a.u.]',fontsize=14)
    if element_map['B2'] != None:
        
        ax.set_title(f'{str(element_map["A"].element)+"2"+str(element_map["B1"].element)+str(element_map["B2"].element)+str(element_map["X"].element)+ "6"} projected DOS',fontsize=18,va='baseline', y=1.05)
    else:
        ax.set_title(f'{str(element_map["A"].element)+"2"+ str(element_map["B1"].element)+ "Vac" +str(element_map["X"].element)+ "6"} projected DOS',fontsize=18,va='baseline', y=1.05)
        
    if len(magmom) == 0:
        fig.text(0.5,0.92, f'$E_F$={efermi:.2f} [eV]',fontsize=14,ha='center',va='baseline')
    else:
        fig.text(0.5,0.92,  f'$E_F=${efermi:.2f} [eV],'+ " $\\mu_{tot}$" + f'={magmom["total"]:.2f} [$\\mu_B$]',fontsize=14,ha='center',va='baseline')
    
    #ax.set_title('A)',fontsize=16,loc='left')
    ax.legend()  
    #print('tdos done')
    plt.savefig(filename,dpi=200)  
