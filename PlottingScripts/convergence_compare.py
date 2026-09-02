import plotly.graph_objects as go
import pandas as pd
import numpy as np
from HDP_PlotslyPlots import plot_dos
from pathlib import Path
from os import makedirs
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin
from monty.serialization import dumpfn
from monty.io import zopen
from monty.json import MontyEncoder, jsanitize
import json
from plotly.subplots import make_subplots
import subprocess

hm_list = ['2347_CsYbAmF',
 '1198_CsNiSmCl',
 '3551_CsPtAmI',
 '3443_CsFeTlI',
 '2828_CsMnAuBr']

mag_list = [
    '2314_CsNaFeF',
 '2160_CsCrCdF',
 '3066_CsYbAmBr',
 '1426_CsNaAgCl',
 '1606_CsLiErCl']

compare_list = mag_list + hm_list

info_df_path = Path("/home/lwalterb/hdp_project/HDP_WorkFlow_Analysis/AnalysisResults/HDP_CombinedInfo_260510.csv")
dcomb = pd.read_csv(info_df_path, index_col=0)
# data_dump_path = Path("/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/ReviewReport/ConvCheck/DOSJson/")
data_dump_path = Path("./convergence_data/")
image_path = data_dump_path / 'images/'
makedirs(data_dump_path,exist_ok=True)
makedirs(image_path, exist_ok=True)

# encoder = MontyEncoder()

for comp_id in compare_list:

    # /home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/AllCompsNewWF/2314_CsNaFeF/4HSE_HSEDOS
    compare_path = Path(f"/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/ReviewReport/ConvCheck/{comp_id}_4HSE/")
    base_data_path = Path(f"/home/lwalterb/RZ-Dienste/hpc-user/lwalterb/HDP_project/UTwente_backup/home/lucw/AllCompsNewWF/{comp_id}/4HSE_HSEDOS/")

    vrun1 = Vasprun(filename=base_data_path/"vasprun.xml")
    
    dos1 = vrun1.complete_dos
    parse_dict = {}
    parse_dict['tdos'] = {}
    parse_dict['tdos_per_site'] = {}
    parse_dict['tdos']['energies'] = dos1.energies - dos1.efermi
    parse_dict['tdos']['densities'] = {'1': dos1.densities[Spin.up], '-1': dos1.densities[Spin.down]}
    for ii, site in enumerate(dos1.structure.sites):
        sitedos = dos1.get_site_dos(site).as_dict()
        parse_dict['tdos_per_site'][ii] = {'densities': {
            '1': sitedos['densities']['1'],
            '-1': sitedos['densities']['-1']
        }}

    # dumpfn(encoder.default(parse_dict),f"{comp_id}_DOSediff4.json.gz")
    with zopen(data_dump_path/f"{comp_id}_DOSediff4.json.gz",'wt') as f:
        json.dump(jsanitize(parse_dict),f)


    vrun2 = Vasprun(filename=compare_path/"vasprun.xml")
    dos2 = vrun2.complete_dos
    parse_dict = {}
    parse_dict['tdos'] = {}
    parse_dict['tdos_per_site'] = {}
    parse_dict['tdos']['energies'] = dos2.energies - dos2.efermi

    tdos = dos2.get_densities(spin=None)
    below_fermi = [
                    i
                    for i in range(len(dos2.energies))
                    if dos2.energies[i] < dos2.efermi and tdos[i] > 1e-4
                ]

    vbm_start = max(below_fermi)
    vbm_energy = dos2.energies[vbm_start]
    parse_dict['tdos']['energies'] = dos2.energies  + dcomb.loc[comp_id]['VBM'] - vbm_energy



    parse_dict['tdos']['densities'] = {'1': dos2.densities[Spin.up], '-1': dos2.densities[Spin.down]}
    for ii, site in enumerate(dos2.structure.sites):
        sitedos = dos2.get_site_dos(site).as_dict()
        parse_dict['tdos_per_site'][ii] = {'densities': {
            '1': sitedos['densities']['1'],
            '-1': sitedos['densities']['-1']
        }}
    
    # dumpfn(encoder.default(parse_dict),f"{comp_id}_DOSediff6.json.gz")
    with zopen(data_dump_path/f"{comp_id}_DOSediff6.json.gz",'wt') as f:
        json.dump(jsanitize(parse_dict),f)


    dosfig1 = plot_dos(
        comp_id= comp_id,
        info_df_path= info_df_path,
        dos_path= data_dump_path,
        dos_extension= 'DOSediff4.json.gz'
    )

    dosfig2 = plot_dos(
        comp_id= comp_id,
        info_df_path= info_df_path,
        dos_path= data_dump_path,
        dos_extension= 'DOSediff6.json.gz'
    )

    # after dosfig1 and dosfig2 are created

    mag_basis = subprocess.run(['awk', '/mag=/ {print $NF}',base_data_path/"OSZICAR"],text=True, capture_output=True).stdout
    mag_comp = subprocess.run(['awk', '/mag=/ {print $NF}',compare_path/"OSZICAR"],text=True,capture_output=True).stdout
    
    fig = make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True,
        shared_xaxes=True,
        subplot_titles=(
            f"{comp_id}: EDIFF=1e-4 mag={mag_basis}",
            f"{comp_id}: EDIFF=1e-6 mag={mag_comp}",
        ),
    )
    fig.update_annotations(font=dict(size=22))

    for tr in dosfig1.data:
        fig.add_trace(tr, row=1, col=1)

    for tr in dosfig2.data:
        fig.add_trace(tr, row=1, col=2)

    fig.update_xaxes(title_text="Dens. of States (#/eV)", row=1, col=1)
    fig.update_xaxes(title_text="Dens. of States (#/eV)", row=1, col=2)
    fig.update_yaxes(title_text="Energy (eV)", row=1, col=1)
    fig.update_yaxes(title_text="Energy (eV)", row=1, col=2)
    fig.update_yaxes(range=[-10, 10], row=1, col=1)
    fig.update_yaxes(range=[-10, 10], row=1, col=2)
    fig.update_xaxes(range=[-30, 30], row=1, col=1)
    fig.update_xaxes(range=[-30, 30], row=1, col=2)
    fig.update_xaxes(
        tickfont={'size': 20},
        title={'font':{'size':22}},
        row=1, col=1
        )
    fig.update_xaxes(
        tickfont={'size': 20},
        title={'font':{'size':22}},
        row=1, col=2
        )
    fig.update_yaxes(
        tickfont={'size': 20},
        title={'font':{'size':22}},
        row=1, col=1
        )
    fig.update_yaxes(
        tickfont={'size': 20},
        title={'font':{'size':22}},
        row=1, col=2
        )


    fig.update_layout(
        # title={"text":f"{comp_id} DOS comparison",
        #        "font":{"size":24}},
        template="plotly_white",
        legend=dict(x=1.02, y=1),
        width=1400,
        height=600,
    )

    fig.show()
    fig.write_image(image_path/f"{comp_id}_DOScomparison.png",width=1400,height=600,scale=2.2)

    if comp_id in hm_list:
        # fig.update_layout(
        #     title={"text":f"{comp_id} DOS comparison zoomed",
        #            "font":{"size":24}},

        # )
        fig.update_yaxes(range=[-1, 1], row=1, col=1)
        fig.update_yaxes(range=[-1, 1], row=1, col=2)
        fig.update_xaxes(range=[-5, 5], row=1, col=1)
        fig.update_xaxes(range=[-5, 5], row=1, col=2)
        fig.show()
        fig.write_image(image_path/f"{comp_id}_DOScomparison_zoomed.png",width=1400,height=600,scale=2.2)


