""" 
HDP_PlotslyPlots.py
Contains several plotting functions based on Plotly. These plots were mainly used for the interactive UMAP explorer.
Plotting functions for:
- Density of States
- COHP/COBI
- Periodic Table Heatplot with several options to select the colorscheme on.
"""

import plotly.express as px
import plotly.graph_objects as go
from typing import Literal
from plotly.subplots import make_subplots
import pymatviz as pmv
import pandas as pd
import numpy as np
from monty.io import zopen
import json
from pathlib import Path
import os


def plot_dos(
    comp_id: str,
    info_df_path: Path,
    dos_path: Path = Path("./lsodos_smeared"),
    dos_extension: str = "lsosmeareddos_persite.json.gz",
):
    """Plots the DOS using plotly. Total DOS is plotted in grey area, pDOS is plotted in color.
    Assures the VBM lies at 0 eV, and that the majority spin channel is plotted as spin-up.

    Args:
        comp_id (str): the CompID of the composition to plot, used to index CombinedInfo DataFrame and locate DOS file.
        info_df_path (Path): Path to saved CombinedInfo DataFrame
        dos_path (Path, optional): path to where parsed DOSs are saved. Defaults to Path("./lsodos_orignals").
        dos_extension (str, optional): extension of DOS json files, will open files according to {CompID}_{extension}. Defaults to "lsodos_persite.json.gz".

    Returns:
        go.Figure
    """

    info_df = pd.read_csv(info_df_path, index_col=0)

    dos_file = dos_path / f"{comp_id}_{dos_extension}"
    with zopen(dos_file, "rt") as f:
        dos_data = json.load(f)

    if info_df.loc[comp_id, "popdiff.total"] >= 0:
        up = "1"
        down = "-1"
    else:
        up = "-1"
        down = "1"

    # make filled traces for up and down total dos
    energie = dos_data["tdos"]["energies"] - info_df.loc[comp_id, "VBM"]
    data = [
        go.Scatter(
            x=np.array(dos_data["tdos"]["densities"][up]),
            y=energie,
            fill="tozerox",
            fillcolor="lightgrey",
            line=dict(color="grey"),
            name="TDOS-up",
            showlegend=True,
        ),
        go.Scatter(
            x=-1 * np.array(dos_data["tdos"]["densities"][down]),
            y=energie,
            fill="tozerox",
            fillcolor="lightgrey",
            line=dict(color="grey"),
            name="TDOS-down",
            showlegend=True,
        ),
    ]

    info_df["element.B2"] = info_df["element.B2"].fillna("Vac")
    if info_df.loc[comp_id, "element.B2"] == "Vac":
        site_map = {
            "A": ["1", "2"],
            "X": ["3", "4", "5", "6", "7", "8"],
            "B1": ["0"],
        }
    else:
        site_map = {
            "A": ["2", "3"],
            "X": ["4", "5", "6", "7", "8", "9"],
            "B1": ["0"],
            "B2": ["1"],
        }

    color_map = {
        "B1": ("255", "0", "0"),
        "B2": ("0", "0", "255"),
        "A": ("20", "145", "20"),
        "X": ("230", "90", "230"),
    }

    fillalpha = 0.1

    for site, site_indices in site_map.items():
        site_dos_up = np.zeros_like(dos_data["tdos"]["energies"])
        site_dos_down = np.zeros_like(dos_data["tdos"]["energies"])
        for idx in site_indices:
            site_dos_up += np.array(dos_data["tdos_per_site"][idx]["densities"][up])
            site_dos_down += np.array(dos_data["tdos_per_site"][idx]["densities"][down])
        data.append(
            go.Scatter(
                x=site_dos_up,
                y=energie,
                fill="tozerox",
                fillcolor=f"rgba({','.join(color_map[site])},{fillalpha})",
                opacity=0.0,
                line=dict(width=2, color=f"rgb({','.join(color_map[site])})"),
                name=f'{info_df.loc[comp_id,f"element.{site}"]}-up',
                showlegend=True,
            )
        )
        data.append(
            go.Scatter(
                x=-1 * site_dos_down,
                y=energie,
                fill="tozerox",
                fillcolor=f"rgba({','.join(color_map[site])},{fillalpha})",
                opacity=0.0,
                line=dict(width=2, color=f"rgb({','.join(color_map[site])})"),
                name=f'{info_df.loc[comp_id,f"element.{site}"]}-down',
                showlegend=True,
            )
        )

    ymax = info_df.loc[comp_id, "bandgap"] + 5.0

    fig = go.Figure(data=data)
    fig.update_layout(
        title=f"Total DOS",
        xaxis_title="Dens. of States (#/eV)",
        yaxis_title="Energy (eV)",
        yaxis=dict(
            range=[-5, ymax], zeroline=True, zerolinewidth=2, zerolinecolor="LightPink"
        ),
        xaxis=dict(
            range=[-10, 10], zeroline=True, zerolinewidth=2, zerolinecolor="Red"
        ),
        showlegend=True,
    )

    return fig


def plot_coxx(
    comp_id: str,
    info_df_path: Path,
    coxx_type: Literal["cohp", "cobi"] = "cohp",
    coxx_path: Path = Path("./cohps_smeared"),
    coxx_extenstion: str = "smearedCOHP.json.gz",
):
    """Plots the COHP or COBI for a specific composition.
    Assures the VBM for each composition is set to 0 eV and that the majority spin channel is plotted as spin-up

    Args:
        comp_id (str): CompID of the composition to plot. Used to index CombinedInfo DataFrame and locate the COHP/COBI file
        info_df_path (Path): Path to saved CombinedInfo DataFrame
        coxx_type (Literal[&quot;cohp&quot;, &quot;cobi&quot;], optional): whether to plot COHP/COBI. Defaults to "cohp".
        coxx_path (Path, optional): Path where COHP/COBI files are stored. Defaults to Path("./cohps_original").
        coxx_extenstion (str, optional): extension used in opening standardized filenames. Will open files according to {CompID}_{extension}. Defaults to "smearedCOHP.json.gz".

    Returns:
        go.Figure
    """

    info_df = pd.read_csv(info_df_path, index_col=0)
    info_df["element.B2"] = info_df["element.B2"].fillna("Vac")

    if coxx_type == "cohp":
        datafile = coxx_path / f"{comp_id}_{coxx_extenstion}"
        prefactor = -1
        xlim = 2
    else:
        datafile = coxx_path / f"{comp_id}_{coxx_extenstion}"
        prefactor = 1
        xlim = 1

    with zopen(datafile, "rt") as f:
        coxx_data = json.load(f)

    energie = (
        np.array(coxx_data["energies"])
        - coxx_data["efermi"]
        - float(info_df.loc[comp_id, "VBM"])
    )

    if info_df.loc[comp_id, "popdiff.total"] >= 0:
        up = "1"
        down = "-1"
    else:
        up = "-1"
        down = "1"

    b1_traces = [
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1avg"][up]),
            name=f'{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, "element.X"]}-avg-up',
            line=dict(color="red"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1avg"][down]),
            name=f"{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, 'element.X']}-avg-down",
            line=dict(color="red", dash="dash"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1x"][up]),
            name=f'{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, "element.X"]}-x-up',
            line=dict(color="blue"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1x"][down]),
            name=f"{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, 'element.X']}-x-down",
            line=dict(color="blue", dash="dash"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1y"][up]),
            name=f'{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, "element.X"]}-y-up',
            line=dict(color="green"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1y"][down]),
            name=f"{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, 'element.X']}-y-down",
            line=dict(color="green", dash="dash"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1z"][up]),
            name=f'{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, "element.X"]}-z-up',
            line=dict(color="purple"),
            showlegend=True,
        ),
        go.Scatter(
            y=energie,
            x=prefactor * np.array(coxx_data["coxx_b1z"][down]),
            name=f"{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, 'element.X']}-z-down",
            line=dict(color="purple", dash="dash"),
            showlegend=True,
        ),
    ]
    if info_df.loc[comp_id, "element.B2"] == "Vac":
        b2_traces = []
    else:
        b2_traces = [
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2avg"][up]),
                name=f'{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, "element.X"]}-avg-up',
                line=dict(color="red"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2avg"][down]),
                name=f"{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, 'element.X']}-avg-down",
                line=dict(color="red", dash="dash"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2x"][up]),
                name=f'{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, "element.X"]}-x-up',
                line=dict(color="blue"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2x"][down]),
                name=f"{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, 'element.X']}-x-down",
                line=dict(color="blue", dash="dash"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2y"][up]),
                name=f'{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, "element.X"]}-y-up',
                line=dict(color="green"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2y"][down]),
                name=f"{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, 'element.X']}-y-down",
                line=dict(color="green", dash="dash"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2z"][up]),
                name=f'{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, "element.X"]}-z-up',
                line=dict(color="purple"),
                showlegend=True,
            ),
            go.Scatter(
                y=energie,
                x=prefactor * np.array(coxx_data["coxx_b2z"][down]),
                name=f"{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, 'element.X']}-z-down",
                line=dict(color="purple", dash="dash"),
                showlegend=True,
            ),
        ]

    # fig = make_subplots(rows=1, cols=2, subplot_titles=(f"{str(prefactor).strip('1')}{coxx_type.upper()}-{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, 'element.X']}", f"{str(prefactor).strip('1')}{coxx_type.upper()}-{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, 'element.X']}"),shared_yaxes=True)
    fig1 = go.Figure(data=b1_traces)
    fig2 = go.Figure(data=b2_traces)

    # for trace in b1_traces:
    #     fig.add_trace(trace, row=1, col=1)
    # for trace in b2_traces:
    #     fig.add_trace(trace, row=1, col=2)
    fig1.update_layout(
        title=f"{str(prefactor).strip('1')}{coxx_type.upper()}-{info_df.loc[comp_id, 'element.B1']}-{info_df.loc[comp_id, 'element.X']}",
        xaxis_title=f"{str(prefactor).strip('1')}{coxx_type.upper()}",
        xaxis=dict(range=[-xlim, xlim]),
        yaxis_title="Energy (eV)",
        yaxis=dict(
            zeroline=True,
            zerolinewidth=2,
            range=[-5, float(info_df.loc[comp_id, "bandgap"]) + 5],
            zerolinecolor="LightPink",
        ),
        showlegend=True,
    )

    fig2.update_layout(
        title=f"{str(prefactor).strip('1')}{coxx_type.upper()}-{info_df.loc[comp_id, 'element.B2']}-{info_df.loc[comp_id, 'element.X']}",
        xaxis_title=f"{str(prefactor).strip('1')}{coxx_type.upper()}",
        xaxis=dict(range=[-xlim, xlim]),
        yaxis_title="Energy (eV)",
        yaxis=dict(
            zeroline=True,
            zerolinewidth=2,
            range=[-5, float(info_df.loc[comp_id, "bandgap"]) + 5],
            zerolinecolor="LightPink",
        ),
        showlegend=True,
    )
    return fig1, fig2


def plot_ptable(
    info_df_orig: pd.DataFrame,
    color_value: str = "Element Counts",
    use_log: bool = False,
):
    """Creates a periodic table heatplot showing average values or count of elements.

    Args:
        info_df_orig (pd.DataFrame): the CombinedInfo DataFrame to base the heatmaps on.
        color_value (str, optional): what value to color the table by. Must be from:
        [Element Count, Average -ICOHP, Average ICOBI, Average Dir.Asym.Index (ICOHP/ICOBI), Average popdiff, Average Band Gap (eV), Average B-Charge, Average X-Charge]. Defaults to "Element Counts".
        use_log (bool, optional): Whether to apply logarithmic colorscaling. Defaults to False.

    Raises:
        ValueError: color_value is not implimented

    Returns:
        go.Figure
    """

    info_df = info_df_orig.copy()
    info_df["element.B2"] = info_df["element.B2"].fillna("Vac")
    b1counts = info_df.groupby("element.B1").count()
    b2counts = info_df.groupby("element.B2").count()
    xcounts = info_df.groupby("element.X").count()
    dcounts = b1counts.add(b2counts, fill_value=0)
    dcounts = dcounts.add(xcounts, fill_value=0)["compID_num"]
    dcounts.name = "Element Counts"
    dcounts["H"] = dcounts["Vac"]

    col_dict = {
        "Element Counts": {"format": ".3g", "exclude_element": ["F", "Cl", "Br", "I"]},
        "Average -ICOHP": {
            "colB1": "Icohp.B1.avg",
            "colB2": "Icohp.B2.avg",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average ICOBI": {
            "colB1": "Icobi.B1.avg",
            "colB2": "Icobi.B2.avg",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Dir.Asym.Index (ICOHP)": {
            "colB1": "Icohp.B1.directional_asym_index",
            "colB2": "Icohp.B2.directional_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Axial.Asym.Index (ICOHP)": {
            "colB1": "Icohp.B1.axial_asym_index",
            "colB2": "Icohp.B2.axial_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Dir.Asym.Index (ICOBI)": {
            "colB1": "Icobi.B1.directional_asym_index",
            "colB2": "Icobi.B2.directional_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Axial.Asym.Index (ICOBI)": {
            "colB1": "Icobi.B1.axial_asym_index",
            "colB2": "Icobi.B2.axial_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average popdiff": {
            "colB1": "popdiff.B1",
            "colB2": "popdiff.B2",
            "format": ".2f",
            "exclude_element": ["F", "Cl", "Br", "I"],
        },
        "Average Band Gap (eV)": {
            "colB1": "bandgap",
            "colB2": "bandgap",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average B-charge": {
            "colB1": "charge.B1",
            "colB2": "charge.B2",
            "format": ".2g",
            "exclude_element": [],
        },
        "Average X-charge": {
            "colB1": "charge.X",
            "colB2": "charge.X",
            "format": ".2g",
            "exclude_element": [],
        },
    }
    if color_value == "Element Counts":
        dcustom = dcounts
        # exclude_elements=['F','Cl','Br','I']

    elif color_value in col_dict.keys():
        info_df["b1col"] = np.abs(info_df[col_dict[color_value]["colB1"]].round(3))
        info_df["b2col"] = np.abs(info_df[col_dict[color_value]["colB2"]].round(3))
        b1vals = info_df.groupby("element.B1")["b1col"].sum()
        b2vals = info_df.groupby("element.B2")["b2col"].sum()
        x1vals = info_df.groupby("element.X")["b1col"].sum()
        x2vals = info_df.groupby("element.X")["b2col"].sum()
        xvals = x1vals.add(x2vals, fill_value=0)
        # print(xvals.count())
        xvals = xvals / 2
        vals = b1vals.add(b2vals, fill_value=0)
        vals = vals.add(xvals, fill_value=0)
        dcustom = vals.div(dcounts, fill_value=0)
        dcustom["H"] = dcustom["Vac"]
        dcustom.name = color_value
        exclude_elements = []
    else:
        raise ValueError(
            f"Could not parse requested color value for P-table: {color_value}\nPlease use 'Element Counts' or one of {list(col_dict.keys())}"
        )

    def format(value):
        return f"{value:.2f}"

    # print(dcounts.loc['H'])
    fig = pmv.ptable_heatmap_plotly(
        dcustom,
        fmt=col_dict[color_value]["format"],
        colorbar={"tickformat": col_dict[color_value]["format"]},
        element_symbol_map={"H": "Vac"},
        log=use_log,
        exclude_elements=col_dict[color_value]["exclude_element"],
        opacity=0.7,
        customdata=dcustom,
    )
    fig.update_layout(clickmode="event+select", font=dict(size=16))
    # Increase colorbar title and tick font sizes for better readability
    cbar_title_size = 24
    cbar_tick_size = 20
    for tr in fig.data:
        if hasattr(tr, "colorbar") and tr.colorbar is not None:
            try:
                title_text = (
                    tr.colorbar.title.text
                    if (tr.colorbar.title and "text" in tr.colorbar.title)
                    else color_value
                )
            except Exception:
                title_text = color_value
            # print(type(tr.colorbar))
            tr.colorbar.title = dict(
                text=title_text,
                font=dict(size=cbar_title_size, weight="bold"),
                side="top",
            )
            tr.colorbar.tickfont = dict(size=cbar_tick_size)
            tr.colorbar.exponentformat = "none"
    return fig


if __name__ == "__main__":

    info_df_path = Path("../AnalysisResults/HDP_CombinedInfo_260510.csv")
    dcomb = pd.read_csv(info_df_path, index_col=0)

    # fig.show()
    saving_dir = Path("../AnalysisResults/hdp_plotly_plots")
    saving_dir.mkdir(exist_ok=True)



    col_dict = {
        "Element Counts": {"format": ".3g", "exclude_element": ["F", "Cl", "Br", "I"]},
        "Average -ICOHP": {
            "colB1": "Icohp.B1.avg",
            "colB2": "Icohp.B2.avg",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average ICOBI": {
            "colB1": "Icobi.B1.avg",
            "colB2": "Icobi.B2.avg",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Dir.Asym.Index (ICOHP)": {
            "colB1": "Icohp.B1.directional_asym_index",
            "colB2": "Icohp.B2.directional_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Axial.Asym.Index (ICOHP)": {
            "colB1": "Icohp.B1.axial_asym_index",
            "colB2": "Icohp.B2.axial_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Dir.Asym.Index (ICOBI)": {
            "colB1": "Icobi.B1.directional_asym_index",
            "colB2": "Icobi.B2.directional_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average Axial.Asym.Index (ICOBI)": {
            "colB1": "Icobi.B1.axial_asym_index",
            "colB2": "Icobi.B2.axial_asym_index",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average popdiff": {
            "colB1": "popdiff.B1",
            "colB2": "popdiff.B2",
            "format": ".2f",
            "exclude_element": ["F", "Cl", "Br", "I"],
        },
        "Average Band Gap (eV)": {
            "colB1": "bandgap",
            "colB2": "bandgap",
            "format": ".2f",
            "exclude_element": [],
        },
        "Average B-charge": {
            "colB1": "charge.B1",
            "colB2": "charge.B2",
            "format": ".2g",
            "exclude_element": [],
        },
        "Average X-charge": {
            "colB1": "charge.X",
            "colB2": "charge.X",
            "format": ".2g",
            "exclude_element": [],
        },
    }

    for plot_type in col_dict.keys():
        ptable = plot_ptable(dcomb, color_value="Average Band Gap (eV)")
        ptable.write_image(file=saving_dir/f"PTable_{plot_type.replace(" ","_")}.svg")

    # fig = plot_coxx(test_comp1, info_df_path, coxx_type='cohp')
    # fig[0].write_html(saving_dir/f'{test_comp1}_B1COHPplot.html')
    # fig[1].write_html(saving_dir/f'{test_comp1}_B2COHPplot.html')
#
# fig = plot_coxx(test_comp2, info_df_path, coxx_type='cohp')
# fig[0].write_html(saving_dir/f'{test_comp2}_B1COHPplot.html')
# fig[1].write_html(saving_dir/f'{test_comp2}_B2COHPplot.html')
