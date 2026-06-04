""" 
CsBBX_Plotting.py
Contains several plotting functions for the Halide Double Perovskite DataBase. 
"""

import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
import numpy as np
from pathlib import Path
import time
from os import PathLike

def plot_bgmagmom(
    combined_info: pd.DataFrame,
    bandgap_df: pd.DataFrame,
    suptitle: str | None = None,
    filename: str | PathLike | None = None,
    use_complex: bool = False,
    use_transition: bool = True,
):
    """Creates a scatter plot of net-spin versus bandgap for the dataset. The datapoints can either be colored by main character of VBM and CBM (p-d, indicates a transition from a p-band to a d-band),
    or by the block-pairing (p-d in this instance indicates a pairing between a B-site from the p-block and a B-site from the d-block).
    Also a more complex glyph formatting was implemented, which changes the symbol based on which sites have the main contribution and half-filled color based on the character,
    but this was not further used.


    Args:
        combined_info (pd.DataFrame): DataFrame containing CombinedInfo on the HDP dataset
        bandgap_df (pd.DataFrame): the DataFrame describing bandedges, used to access the 'secondary gaps' for half-metals
        suptitle (str | None, optional): A supertitle to add to the plot. Defaults to None.
        filename (str | None, optional): Filename to save the plot. If None is given plot will be shown instead. Defaults to None.
        use_complex (bool, optional): Whether to use the more complex glyph formatting, results in very busy and complex plot. Defaults to False.
        use_transition (bool, optional): Whether to base datapoint coloring on description of VBM-CBM. If False will use B-site block-pairing instead. Defaults to True.
    """
    plt.rcParams.update(
        {
            "text.usetex": True,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "legend.handletextpad": 0.4,
            "legend.handlelength": 0.6,
            "legend.columnspacing": 0.6,
        }
    )
    # combine all info into one DF for plotting

    # we also want to plot the minority spin gap in half-metals
    halfm_index = combined_info[combined_info["cond_type"] == "half-metal"].index
    hm_channel_gaps = bandgap_df.loc[halfm_index].drop("combined", level=1)
    secondary_gaps = hm_channel_gaps[hm_channel_gaps["bandgap"] > 0.1]

    metallic_index = combined_info[combined_info["cond_type"] == "metallic"].index

    ##Start plotting
    fig, axes = plt.subplots(
        2,
        2,
        sharex="col",
        sharey="row",
        figsize=(10, 8),
        width_ratios=(5, 1),
        height_ratios=(1, 4),
        layout="constrained",
    )

    ##plot histogram of bandgaps
    bg_binning = [
        0.0,
        0.1,
        0.25,
        0.5,
        0.75,
        1.0,
        1.25,
        1.5,
        1.75,
        2.0,
        2.25,
        2.5,
        2.75,
        3.0,
        3.25,
        3.5,
        3.75,
        4.0,
        4.25,
        4.5,
        4.75,
        5.0,
        5.25,
        5.5,
        5.75,
        6.0,
        6.25,
        6.5,
        6.75,
        7.0,
        7.25,
        7.5,
        7.75,
        8.0,
        8.25,
        8.5,
        8.75,
        9.0,
        9.25,
        9.5,
        9.75,
        10.0,
        10.25,
    ]
    axes[0][0].hist(
        "bandgap",
        color="#198B15E2",
        bins=bg_binning,
        data=combined_info.drop(metallic_index).drop(halfm_index),
        alpha=0.7,
        label="Bandgap values",
    )
    axes[0][0].hist(
        "bandgap",
        color="#8D800BDF",
        bins=bg_binning,
        data=combined_info[combined_info["cond_type"] == "metallic"],
        alpha=1,
        label="Metallic",
    )
    axes[0][0].hist(
        "bandgap",
        color="#4C1D72FB",
        bins=bg_binning,
        data=secondary_gaps,
        alpha=0.7,
        label="Half-metal gaps",
    )
    axes[0][0].set_ylabel("Counts", fontsize=24)

    axes[0][0].set_ylim([0, 55])
    # axes[0][0].set_xlabel('Bandgap (eV)',fontsize=14)

    # axes[0][0].set_title('Bandgap counts', fontsize=16)
    axes[0][0].legend(fontsize=20, loc="upper right")

    axes[0][1].axis("off")

    if use_transition:
        legendtitle = "Character of \n VBM-CBM:"
    else:
        legendtitle = "Block pairing"

    # plot histogram of magmom0
    magmom_binning = np.arange(0, 13, 1)
    axes[1][1].hist(
        abs(combined_info["popdiff.total"]),
        color="#001D7C9B",
        bins=magmom_binning,
        alpha=0.7,
        align="left",
        label="Net polarization",
        orientation="horizontal",
    )
    # axes[1][1].set_ylabel('Mag. moment ($\\mu_B$)',fontsize=14)
    axes[1][1].set_xlabel("Counts", fontsize=24)
    axes[1][1].set_xlim([0, 350])
    # axes[1][1].set_title('Magmom counts', fontsize=16)
    # axes[1][1].legend(fontsize=13)

    # Plot scatter of magmom vs. bg
    # determe whether to use vasp or lob magmom

    bandcolors = {
        "s": "tab:orange",
        "p": "tab:purple",
        "d": "tab:blue",
        "f": "tab:green",
    }
    #
    transshapes = {
        "B-B": "s",
        "B-X": "P",
        "X-B": "X",
        "X-X": "o",
        "B-A": "^",
        "A-B": "v",
        "A-A": "d",
        "X-A": "<",
        "A-X": ">",
    }

    trans_bands = [
        "d-d",
        "d-f",
        "d-p",
        "d-s",
        "s-d",
        "s-f",
        "s-p",
        "s-s",
        "f-d",
        "f-f",
        "f-p",
        "f-s",
        "p-d",
        "p-f",
        "p-p",
        "p-s",
    ]
    transband_colormap = {
        label: color for label, color in zip(trans_bands, colormaps["tab20c"].colors)
    }

    allpairings = [
        "d-d",
        "d-f",
        "d-p",
        "d-s",
        "s-s",
        "s-dum",
        "s-dum1",
        "s-dum2",
        "f-f",
        "f-p",
        "f-s",
        "f-dummy",
        "p-p",
        "p-s",
        "p-dummy1",
        "p-dummy2",
        "Vac-d",
        "Vac-f",
        "Vac-p",
        "Vac-dummy",
    ]
    blockpair_colormap = {
        label: color for label, color in zip(allpairings, colormaps["tab20c"].colors)
    }

    had_transbands = []
    # combined_info = determineTransitionBands(combined_info,include_x=True)

    # combined_info = combined_info.loc[combined_info.index]

    ii = 0
    for idx in combined_info.index:
        complex_marker_style = {
            "marker": transshapes[combined_info.loc[idx, "transition_sites"]],
            # 'fillstyle':'left',
            "markersize": 6
            + combined_info.loc[idx, "VBMtotcontr.X"] * 2
            + combined_info.loc[idx, "CBMtotcontr.X"] * 2,
            "markerfacecolor": bandcolors[combined_info.loc[idx, "vbmorbchar"]],
            "markerfacecoloralt": bandcolors[combined_info.loc[idx, "cbmorbchar"]],
            "markeredgecolor": "black",
            "alpha": (
                combined_info.loc[idx, "vbmorbcontr"]
                + combined_info.loc[idx, "cbmorbcontr"]
            )
            / 2
            - 0.2,
        }

        simplified_marker_style = {
            "marker": "o",  # transshapes[combined_info.loc[idx,'transition_sites']],
            # 'fillstyle' : 'full',
            "markersize": 7,
            "markerfacecolor": (
                transband_colormap[combined_info.loc[idx, "transition_bands"]]
                if use_transition
                else blockpair_colormap[combined_info.loc[idx, "block_pairing"]]
            ),
            "markeredgecolor": "black",
            "markeredgewidth": 0.1,
            "linewidth": 0,
            "alpha": 0.8,
        }

        bandgap = (
            secondary_gaps.loc[idx, "bandgap"].min()
            if combined_info.loc[idx, "cond_type"] == "half-metal"
            else combined_info.loc[idx, "bandgap"]
        )
        # print(idx, bandgap)
        if (
            combined_info.loc[idx, "transition_bands"] not in had_transbands
            and use_transition
        ):
            # ii+=1
            # if ii >= len(trans_bands):
            #     ii = len(trans_bands) -1
            had_transbands.append(combined_info.loc[idx, "transition_bands"])
            axes[1][0].plot(
                bandgap,
                abs(combined_info.loc[idx, "popdiff.total"]),
                fillstyle="full",
                **simplified_marker_style,
                label=combined_info.loc[idx, "transition_bands"],
            )

        elif (
            not use_transition
            and combined_info.loc[idx, "block_pairing"] not in had_transbands
            and not use_complex
        ):
            had_transbands.append(combined_info.loc[idx, "block_pairing"])
            axes[1][0].plot(
                bandgap,
                abs(combined_info.loc[idx, "popdiff.total"]),
                fillstyle="full",
                **simplified_marker_style,
                label=combined_info.loc[idx, "block_pairing"],
            )
        elif use_complex:
            axes[1][0].plot(
                bandgap,
                abs(combined_info.loc[idx, "popdiff.total"]),
                fillstyle="left",
                **complex_marker_style,
            )
        else:
            axes[1][0].plot(
                bandgap,
                abs(combined_info.loc[idx, "popdiff.total"]),
                **simplified_marker_style,
            )

    handles, labels = axes[1][0].get_legend_handles_labels()
    order = []

    bandorder = [
        "s-d",
        "s-f",
        "s-p",
        "s-s",
        "p-d",
        "p-f",
        "p-p",
        "p-s",
        "d-d",
        "d-f",
        "d-p",
        "d-s",
        "f-d",
        "f-f",
        "f-p",
        "f-s",
    ]
    pairorder = [
        "s-s",
        "s-dum",
        "s-dum1",
        "s-dum2",
        "p-p",
        "p-s",
        "p-dummy1",
        "p-dummy2",
        "d-d",
        "d-f",
        "d-p",
        "d-s",
        "f-f",
        "f-p",
        "f-s",
        "f-dummy",
        "Vac-d",
        "Vac-f",
        "Vac-p",
        "Vac-dummy",
    ]

    orderedlist = bandorder if use_transition else pairorder

    for bands in orderedlist:
        try:
            order.append(labels.index(bands))
        except:
            continue
    # print(order)
    axes[1][0].legend(
        [handles[idx] for idx in order],
        [labels[idx] for idx in order],
        ncols=2,
        title=legendtitle,
        title_fontsize=22,
        fontsize=22,
        frameon=True,
        framealpha=0.9,
        loc="upper right",
    )

    axes[1][0].set_xlabel("Band Gap [eV]", fontsize=24)
    axes[1][0].set_ylabel("Spin Magnetic Moment ($\\mu_s$)", fontsize=24)
    axes[1][0].set_ylim([-0.2, 14.2])
    axes[1][0].set_xlim([-0.1, 9.5])

    if suptitle:
        fig.suptitle(suptitle, fontsize=28)

    if filename:
        plt.savefig(filename, dpi=200)
    else:
        plt.show()

    return


def plot_icohpbg(
    combined_info: pd.DataFrame,
    use_transition: bool = False,
    includex_transition: bool = True,
    use_icobi: bool = False,
    filename: str | None = None,
):

    trans_bands = [
        "d-d",
        "d-f",
        "d-p",
        "d-s",
        "s-d",
        "s-f",
        "s-p",
        "s-s",
        "f-d",
        "f-f",
        "f-p",
        "f-s",
        "p-d",
        "p-f",
        "p-p",
        "p-s",
    ]
    transband_colormap = {
        label: color for label, color in zip(trans_bands, colormaps["tab20c"].colors)
    }
    cbm_colormap = {
        "s": "tab:orange",
        "p": "tab:purple",
        "d": "tab:blue",
        "f": "tab:green",
    }

    allpairings = [
        "d-d",
        "d-f",
        "d-p",
        "d-s",
        "s-s",
        "s-dum",
        "s-dum1",
        "s-dum2",
        "f-f",
        "f-p",
        "f-s",
        "f-dummy",
        "p-p",
        "p-s",
        "p-dummy1",
        "p-dummy2",
        "Vac-d",
        "Vac-f",
        "Vac-p",
        "Vac-dummy",
    ]
    blockpair_colormap = {
        label: color for label, color in zip(allpairings, colormaps["tab20c"].colors)
    }

    if use_transition:
        assert "vbmorbchar" in combined_info.columns, "Tried to use transition, but combined_info DataFrame does not contain right information"
        groupcolumn = "vbmorbchar"
        colorcolumn = "cbmorbchar"
        labelcolumn = "transition_bands"
        cmap = cbm_colormap

    else:
        groupcolumn = "block_pairing"
        colorcolumn = "block_pairing"
        labelcolumn = "block_pairing"
        cmap = cbm_colormap

    if use_icobi:
        combined_info["Icobi.total"] = combined_info["Icobi.B1.sum"].fillna(
            0
        ) + combined_info["Icobi.B2.sum"].fillna(0)
        plotcol = "Icobi.total"
        ylabel = "ICOBI"
        ymax = 7
    else:
        combined_info["Icohp.total"] = -combined_info["Icohp.B1.sum"].fillna(
            0
        ) - combined_info["Icohp.B2.sum"].fillna(0)
        combined_info["Icohp.diff"] = abs(
            -combined_info["Icohp.B1.sum"].fillna(0)
            + combined_info["Icohp.B2.sum"].fillna(0)
        )
        plotcol = "Icohp.diff"
        ylabel = "-ICOHP difference (eV)"
        ymax = combined_info[plotcol].max() * 1.1

    if not use_transition:
        starts = ["Vac", "s", "p", "d", "f"]
        groups = {
            "Vac": combined_info[combined_info["block_pairing"].str.startswith("Vac")],
            "s": combined_info[combined_info["block_pairing"].str.startswith("s")],
            "p": combined_info[combined_info["block_pairing"].str.startswith("p")],
            "d": combined_info[combined_info["block_pairing"].str.startswith("d")],
            "f": combined_info[combined_info["block_pairing"].str.startswith("f")],
        }
    else:
        groups = {
            "s": combined_info[combined_info[groupcolumn] == "s"],
            "p": combined_info[combined_info[groupcolumn] == "p"],
            "d": combined_info[combined_info[groupcolumn] == "d"],
            "f": combined_info[combined_info[groupcolumn] == "f"],
        }

    for name, group in groups.items():
        fig, axes = plt.subplots(
            2, 2, layout="constrained", sharex=True, sharey=True, figsize=(8, 8)
        )
        axmap = {"F": axes[0][0], "Cl": axes[0][1], "Br": axes[1][0], "I": axes[1][1]}
        for xatom, bla in group.groupby("element.X"):
            data = group.groupby("element.X").get_group(xatom)
            hadlabels = []
            ax = axmap[xatom]
            for idx in data.index:
                style = {
                    "marker": "o",
                    "markersize": 7,
                    "markerfacecolor": (
                        cmap[data.loc[idx, colorcolumn]]
                        if use_transition
                        else cmap[data.loc[idx, colorcolumn].split("-")[-1]]
                    ),
                    "markeredgecolor": "black",
                    "markeredgewidth": 0.1,
                    "linewidth": 0,
                    "alpha": 0.75,
                }
                if data.loc[idx, colorcolumn] not in hadlabels:
                    hadlabels.append(data.loc[idx, colorcolumn])
                    ax.plot(
                        data.loc[idx, "bandgap"],
                        data.loc[idx, plotcol],
                        label=data.loc[idx, labelcolumn],
                        **style,
                    )
                else:

                    ax.plot(data.loc[idx, "bandgap"], data.loc[idx, plotcol], **style)
            ax.set_title(f"$X$ = {xatom}", fontsize=16)
            ax.set_xlabel("Bandgap [eV]", fontsize=14)
            ax.set_ylabel(ylabel, fontsize=14)
            ax.legend(
                fontsize=12, title=labelcolumn.replace("_", "\n"), title_fontsize=11
            )
            ax.set_xlim([-0.1, 10.3])
            ax.set_ylim([0.0, ymax])
            ax.grid(True, linestyle="--", alpha=0.4)
            fig.suptitle(f"{name} - {ylabel} vs. Bandgap", fontsize=20)
        if filename:
            plt.savefig(f"{filename}_{name}.png", dpi=200)
        else:
            fig.show()

    return


# def plot_icohpmagmom(combined_info:pd.DataFrame):
def plot_rowicohp(
    combined_info: pd.DataFrame, use_icobi: bool = False, filename: str | None = None
):
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "Helvetica",
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.handletextpad": 0.4,
            "legend.handlelength": 0.6,
            "legend.columnspacing": 0.8,
        }
    )
    # combine all info into one DF for plotting
    fig, axes = plt.subplots(2, 2, figsize=(8, 8), layout="constrained")
    combined_info["Icobi.total"] = combined_info["Icobi.B1.sum"].fillna(
        0
    ) + combined_info["Icobi.B2.sum"].fillna(0)
    combined_info["Icohp.total"] = -combined_info["Icohp.B1.sum"].fillna(
        0
    ) - combined_info["Icohp.B2.sum"].fillna(0)

    if use_icobi:
        plotcol = ["Icobi", "avg"]
        ylabel = "ICOBI$_{(B.X)avg}$"
        ylim = [0, combined_info[f"{plotcol[0]}.B1.{plotcol[1]}"].max() * 1.1]
        prefactor = 1
    else:
        plotcol = ["Icohp", "avg"]
        ylabel = "-ICOHP$_{(B.X)avg}$ (eV)"
        ylim = [
            -combined_info[f"{plotcol[0]}.B1.{plotcol[1]}"].max() * 1.1,
            -combined_info[f"{plotcol[0]}.B1.{plotcol[1]}"].min() * 1.1,
        ]
        prefactor = -1

    block_colormap = {
        "s": "tab:orange",
        "p": "tab:purple",
        "d": "tab:blue",
        "f": "tab:green",
        "Vac": "tab:gray",
    }
    xatom_symbolmap = {"F": "o", "Cl": "s", "Br": "^", "I": "D"}
    axmap = {"F": axes[0][0], "Cl": axes[0][1], "Br": axes[1][0], "I": axes[1][1]}
    offsetmap = {"s": -0.2, "p": -0.1, "d": 0.0, "f": 0.1, "Vac": 0.15}

    for xatom, data in combined_info.groupby("element.X"):
        hadlabels = []
        ax = axmap[xatom]
        # style = {
        #     'marker' : xatom_symbolmap[xatom],
        #     'markersize' : 5,
        #     'markeredgecolor' : 'black',
        #     'markeredgewidth' : 0.1,
        #     'linewidth' : 0,
        #     'alpha' : 0.65
        # }
        # print(data['row.B1']+ [offsetmap[x] for x in data['block.B1']])
        for idx in data.index:
            if data.loc[idx, "block.B1"] not in hadlabels:
                hadlabels.append(data.loc[idx, "block.B1"])
                ax.scatter(
                    data.loc[idx, "row.B1"] + offsetmap[data.loc[idx, "block.B1"]],
                    prefactor * data.loc[idx, f"{plotcol[0]}.B1.{plotcol[1]}"],
                    marker="o",
                    s=5,
                    alpha=0.6,
                    c=block_colormap[data.loc[idx, "block.B1"]],
                    label=data.loc[idx, "block.B1"],
                )
                ax.scatter(
                    data.loc[idx, "row.B2"] + offsetmap[data.loc[idx, "block.B2"]],
                    prefactor * data.loc[idx, f"{plotcol[0]}.B2.{plotcol[1]}"],
                    marker="o",
                    s=5,
                    alpha=0.6,
                    c=block_colormap[data.loc[idx, "block.B2"]],
                )
            else:
                ax.scatter(
                    data.loc[idx, "row.B1"] + offsetmap[data.loc[idx, "block.B1"]],
                    prefactor * data.loc[idx, f"{plotcol[0]}.B1.{plotcol[1]}"],
                    marker="o",
                    s=5,
                    alpha=0.6,
                    c=block_colormap[data.loc[idx, "block.B1"]],
                )
                ax.scatter(
                    data.loc[idx, "row.B2"] + offsetmap[data.loc[idx, "block.B2"]],
                    prefactor * data.loc[idx, f"{plotcol[0]}.B2.{plotcol[1]}"],
                    marker="o",
                    s=5,
                    alpha=0.6,
                    c=block_colormap[data.loc[idx, "block.B2"]],
                )
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_xlabel("row of B", fontsize=14)
        ax.set_title(f"X={xatom}", fontsize=16)
        ax.grid("minor", linestyle="--", alpha=0.4)
        ax.legend(fontsize=11, ncols=2)
        ax.set_xlim([1, combined_info["row.B1"].max() * 1.1])
        ax.set_ylim(ylim)

    if filename:
        plt.savefig(filename, dpi=200)
    else:
        fig.show()
    return


def plot_icohpviolin_perspecies(
    combined_info: pd.DataFrame, use_icobi: bool = False, filedir: Path | None = None
):
    plt.rcParams.update(
        {
            "text.usetex": True,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.handletextpad": 0.4,
            "legend.handlelength": 0.6,
            "legend.columnspacing": 0.8,
        }
    )

    unique_specie_list = combined_info["element.B1"].unique().tolist() + (
        list(
            set(combined_info["element.B2"].dropna().unique())
            - set(combined_info["element.B1"].dropna().unique())
        )
    )
    print(unique_specie_list)
    for specie in unique_specie_list:

        specie_data1 = combined_info[combined_info["element.B1"] == specie]
        specie_data2 = combined_info[combined_info["element.B2"] == specie]

        fig, axs = plt.subplots(
            ncols=4, nrows=2, sharey="row", figsize=(15, 9), layout="constrained"
        )
        icohpaxmap = {"F": axs[0, 0], "Cl": axs[0, 1], "Br": axs[0, 2], "I": axs[0, 3]}
        icobiaxmap = {"F": axs[1, 0], "Cl": axs[1, 1], "Br": axs[1, 2], "I": axs[1, 3]}

        if not specie_data1.empty and not specie_data2.empty:
            icohpymax = (
                min(
                    specie_data1["Icohp.B1.avg"].dropna().min(),
                    specie_data2["Icohp.B2.avg"].dropna().min(),
                )
                * -1.1
            )
            icohpymin = (
                max(
                    specie_data1["Icohp.B1.avg"].dropna().max(),
                    specie_data2["Icohp.B2.avg"].dropna().max(),
                )
                * -1
                - 0.05
            )

            icobiymin = (
                min(
                    specie_data1["Icobi.B1.avg"].dropna().min(),
                    specie_data2["Icobi.B2.avg"].dropna().min(),
                )
                * 0.9
            )
            icobiymax = (
                max(
                    specie_data1["Icobi.B1.avg"].dropna().max(),
                    specie_data2["Icobi.B2.avg"].dropna().max(),
                )
                * 1.1
            )
        elif specie_data1.empty and not specie_data2.empty:
            icohpymax = specie_data2["Icohp.B2.avg"].dropna().min() * -1.1
            icohpymin = specie_data2["Icohp.B2.avg"].dropna().max() * -1 - 0.05

            icobiymin = specie_data2["Icobi.B2.avg"].dropna().min() * 0.9
            icobiymax = specie_data2["Icobi.B2.avg"].dropna().max() * 1.1
        elif not specie_data1.empty and specie_data2.empty:
            icohpymax = specie_data1["Icohp.B1.avg"].dropna().min() * -1.1
            icohpymin = specie_data1["Icohp.B1.avg"].dropna().max() * -1 - 0.05

            icobiymin = specie_data1["Icobi.B1.avg"].dropna().min() * 0.9
            icobiymax = specie_data1["Icobi.B1.avg"].dropna().max() * 1.1
        else:
            print(f"No data for specie {specie}, skipping...")
            continue


        for xatom in combined_info["element.X"].unique():
            ax1 = icohpaxmap[xatom]
            ax2 = icobiaxmap[xatom]
            icohpdata_per_block = {
                "all": [],
                "s": [],
                "p": [],
                "d": [],
                "f": [],
                "Vac": [],
            }

            icobidata_per_block = {
                "all": [],
                "s": [],
                "p": [],
                "d": [],
                "f": [],
                "Vac": [],
            }

            if xatom in specie_data1["element.X"].values:
                for idx in specie_data1.groupby("element.X").get_group(xatom).index:
                    icohpdata_per_block["all"].append(
                        -specie_data1.loc[idx, "Icohp.B1.avg"]
                    )
                    icohpdata_per_block[specie_data1.loc[idx, "block.B2"]].append(
                        -specie_data1.loc[idx, "Icohp.B1.avg"]
                    )

                    icobidata_per_block["all"].append(
                        specie_data1.loc[idx, "Icobi.B1.avg"]
                    )
                    icobidata_per_block[specie_data1.loc[idx, "block.B2"]].append(
                        specie_data1.loc[idx, "Icobi.B1.avg"]
                    )

            if xatom in specie_data2["element.X"].values:
                for idx in specie_data2.groupby("element.X").get_group(xatom).index:
                    icohpdata_per_block["all"].append(
                        -specie_data2.loc[idx, "Icohp.B2.avg"]
                    )
                    icohpdata_per_block[specie_data2.loc[idx, "block.B1"]].append(
                        -specie_data2.loc[idx, "Icohp.B2.avg"]
                    )

                    icobidata_per_block["all"].append(
                        specie_data2.loc[idx, "Icobi.B2.avg"]
                    )
                    icobidata_per_block[specie_data2.loc[idx, "block.B1"]].append(
                        specie_data2.loc[idx, "Icobi.B2.avg"]
                    )

            icohplabels = [
                (
                    f"{key}\n{len(val)}\n{round(np.mean(val),2)}\n{round(np.median(val),2)}\n{round(np.std(val),2)}"
                    if key != "all"
                    else f"{key}\nn:     {len(val)}\nmean: {round(np.mean(val),2)}\nmedian: {round(np.median(val),2)}\nstd: {round(np.std(val),2)}"
                )
                for key, val in icohpdata_per_block.items()
            ]
            icobilabels = [
                (
                    f"{key}\n{len(val)}\n{round(np.mean(val),2)}\n{round(np.median(val),2)}\n{round(np.std(val),2)}"
                    if key != "all"
                    else f"{key}\nn:     {len(val)}\nmean: {round(np.mean(val),2)}\nmedian: {round(np.median(val),2)}\nstd: {round(np.std(val),2)}"
                )
                for key, val in icobidata_per_block.items()
            ]

            for key, val in icohpdata_per_block.items():
                if len(val) == 0:
                    icohpdata_per_block[key] = list(np.zeros(10))
                    icobidata_per_block[key] = list(np.zeros(10))

            ax1.violinplot(
                icohpdata_per_block.values(), showmeans=False, showmedians=True
            )
            ax1.set_title(f"-ICOHP for $X$={xatom}", fontsize=16)
            ax1.set_xticks(
                np.arange(1, len(icohpdata_per_block.keys()) + 1), labels=icohplabels
            )
            for ii, label in enumerate(ax1.get_xticklabels()):
                # label.set_rotation(0)
                if ii == 0:
                    label.set_ha("right")
            ax1.set_xlabel("Block of other B atom\n  ", fontsize=14)
            ax1.set_ylabel("-ICOHP (eV)", fontsize=14)
            ax1.yaxis.grid(True)
            try:
                ax1.set_ylim([icohpymin, icohpymax])
            except Exception as e:
                print(
                    f"Could not set ylim for -ICOHP plot of {specie} with X={xatom}, due to {e}"
                )

            ax2.violinplot(
                icobidata_per_block.values(), showmeans=False, showmedians=True
            )
            ax2.set_title(f"ICOBI for $X$={xatom}", fontsize=16)
            ax2.set_xticks(
                np.arange(1, len(icobidata_per_block.keys()) + 1), labels=icobilabels
            )
            for ii, label in enumerate(ax2.get_xticklabels()):
                # label.set_rotation(0)
                if ii == 0:
                    label.set_ha("right")
            ax2.set_xlabel("Block of other B atom", fontsize=14)
            ax2.set_ylabel("ICOBI", fontsize=14)
            ax2.yaxis.grid(True)
            ax2.set_ylim([icobiymin, icobiymax])

        fig.suptitle(
            f"Bonding indicator distribution for bonds involving {specie}", fontsize=22
        )
        if filedir:
            plt.savefig(filedir / f"IcohpDist_{specie}.png", dpi=200)
        else:
            fig.show()
            input("press enter to continue...")

    return


def compare_alattice(
    combined_info: pd.DataFrame, icsd_df: pd.DataFrame, make_fig: bool = False
):
    plt.rcParams.update(
        {
            "text.usetex": False,
            "font.family": "sans-serif",
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
        }
    )
    icsd_df["cleanedformula"] = icsd_df.apply(
        lambda row: row["StructuredFormula"].replace("(", "").replace(")", ""), axis=1
    )
    print("Total ICSD entries:", icsd_df.shape[0])
    icsd_df = icsd_df[icsd_df["Pressure"] < 1]
    print("ICSD entries after dropping high P measurements", icsd_df.shape[0])

    ###extract exact B-site and X-site makeup of ICSD data due to changing ordering and duplicate entries
    bsites = []
    xsites = []
    simplenames = []
    for val in icsd_df["cleanedformula"].astype(str).values:
        lb = [x for x in val.split(" ") if "2" not in x and "6" not in x]
        lx = [x for x in val.split(" ") if "6" in x]
        xsites.append(lx[0].replace("6", ""))
        lb.sort()
        simplenames.append("Cs" + "".join(lb) + lx[0].replace("6", ""))
        if len(lb) == 1:
            lb.append(np.nan)
        bsites.append(set(lb))

    # print(bsites)
    icsd_df["bsites"] = bsites
    icsd_df["xsite"] = xsites
    icsd_df["simplename"] = simplenames
    print(f"Total hdp compositions {len(icsd_df['simplename'].unique())}")

    ##filter out all compounds for which no Fm3m entries exist
    noncubic_df = icsd_df[icsd_df["HMS"] != "F m -3 m"]
    cubic_df = icsd_df[icsd_df["HMS"] == "F m -3 m"]
    noncubic_df["unique"] = noncubic_df.apply(
        lambda row: row["simplename"] not in list(cubic_df["simplename"].unique()),
        axis=1,
    )
    noncubic_df = noncubic_df[noncubic_df["unique"] == True]
    noncubic_df.to_csv(f'ICSD_NonCubicHDPdata{time.strftime("%d%m%y")}.csv')
    print(
        f"Dropped {len(noncubic_df['simplename'].unique())} compounds for having no Fm-3m ICSD entry"
    )
    print(f"{len(cubic_df['simplename'].unique())} compounds left for a_lat comparison")

    # print(cubic_df['simplename'].unique())
    # make set of bsites so its unordered
    combined_info["bsites"] = combined_info.apply(
        lambda row: set([row["element.B1"], row["element.B2"]]), axis=1
    )

    # extract a-lat from ICSD data
    cubic_df["a_lat"] = cubic_df.apply(
        lambda row: row["StandardisedCellParameter"].split(" ")[0], axis=1
    )

    # combine ICSD data to avg, min, max and add simulated values
    comparison_df = pd.DataFrame(index=cubic_df["simplename"].unique())
    for hdpname, idxs in cubic_df.groupby("simplename").groups.items():
        comparison_df.loc[hdpname, "bsites"] = str(cubic_df.loc[idxs[0], "bsites"])
        comparison_df.loc[hdpname, "xsite"] = cubic_df.loc[idxs[0], "xsite"]

        comparison_df.loc[hdpname, "icsd_a_lat_avg"] = (
            cubic_df.loc[idxs, "a_lat"].astype(np.float64).mean()
        )
        comparison_df.loc[hdpname, "icsd_a_lat_min"] = (
            cubic_df.loc[idxs, "a_lat"].astype(np.float64).min()
        )
        comparison_df.loc[hdpname, "icsd_a_lat_max"] = (
            cubic_df.loc[idxs, "a_lat"].astype(np.float64).max()
        )

        comparison_df.loc[hdpname, "a_lat_diff_min"] = (
            comparison_df.loc[hdpname, "icsd_a_lat_avg"]
            - comparison_df.loc[hdpname, "icsd_a_lat_min"]
        )
        comparison_df.loc[hdpname, "a_lat_diff_max"] = (
            -comparison_df.loc[hdpname, "icsd_a_lat_avg"]
            + comparison_df.loc[hdpname, "icsd_a_lat_max"]
        )

        try:
            comparison_df.loc[hdpname, "sim_a_lat"] = combined_info[
                (
                    (combined_info["bsites"] == cubic_df.loc[idxs[0], "bsites"])
                    & (combined_info["element.X"] == cubic_df.loc[idxs[0], "xsite"])
                )
            ].iloc[0]["lattice_a_conventional"]
        except IndexError:
            comparison_df.loc[hdpname, "sim_a_lat"] = np.nan
    print(
        f"Could not find simulated data for {len(comparison_df.index) - len(comparison_df.dropna().index)} compounds"
    )

    comparison_df.to_csv(f'ICSD_CubicHDPdataComparison{time.strftime("%d%m%y")}.csv')
    # print(icsd_df['bsites'].unique())
    comparison_df = comparison_df.dropna()
    # Calculate mean absolute error
    abs_err_sum = np.abs(
        comparison_df["sim_a_lat"] - comparison_df["icsd_a_lat_avg"]
    ).sum()
    sq_err_sum = np.sum(
        np.array(comparison_df["sim_a_lat"] - comparison_df["icsd_a_lat_avg"]) ** 2
    )
    mae = abs_err_sum / len(comparison_df["icsd_a_lat_avg"])
    rmse = (sq_err_sum / len(comparison_df["sim_a_lat"])) ** (1 / 2)
    print("Comparison complete, error calculation...")
    print(f"{abs_err_sum} MAE is :{mae}")
    print(f"{sq_err_sum} RMSE is : {rmse}")

    if make_fig:
        fig, ax = plt.subplots(figsize=(6, 6), layout="constrained")
        compgrouped = comparison_df.groupby("xsite")

        for xsite in compgrouped.groups.keys():
            data = compgrouped.get_group(xsite)
            # print([list(data['icsd_a_lat_min'].values), list(data['icsd_a_lat_max'].values)])
            ax.errorbar(
                x=data["sim_a_lat"],
                y=data["icsd_a_lat_avg"],
                yerr=[
                    list(data["a_lat_diff_min"].values),
                    list(data["a_lat_diff_max"].values),
                ],
                label=xsite,
                fmt="o",
            )
        ax.plot([1, 20], [1, 20], "k--")
        ax.set_xlabel("Simulated lattice constant [Ang]", fontsize=16)
        ax.set_ylabel("Experimental lattice constant [Ang]", fontsize=16)
        # ax.set_title('Comparison experimental and Simulated lattice constants',fontsize=16)
        ax.set_ylim(8.5, 12.5)
        ax.set_xlim(8.5, 12.5)
        ax.legend(title="X-atom", fontsize=16, title_fontsize=18)
        plt.savefig(f"ICSD_ComparisonPlot.png", dpi=250)
    return cubic_df


anal_dir = Path("../AnalysisResults/")
basicinfo_csv = anal_dir / "HDP_BasicInfo_260510.csv"
edgeinfo_csv = anal_dir / "HDP_bandedgeInfo_lsodos_260510.csv"
lobinfo_csv = anal_dir / "HDP_LobsterInfo_260510.csv"
combinfo_csv = anal_dir / "HDP_CombinedInfo_260510.csv"
strucutralinfo_csv = anal_dir  / "HDP_StructuralInfo_260510.csv"
dbasic = pd.read_csv(basicinfo_csv, index_col=0)
dedge = pd.read_csv(edgeinfo_csv, index_col=[0, 1])
dlob = pd.read_csv(lobinfo_csv, index_col=0)
dstruc = pd.read_csv(strucutralinfo_csv, index_col=0)
dcomb = pd.read_csv(combinfo_csv, index_col=0)


dicsd = pd.read_csv(
    "DATA NOT INCLUDED DUE TO LICENSE. SEE AnalysisResults/FoundICSD_HDPs.txt for ICSD ENTRIES",
    index_col=0,
    sep="\t",
)

if __name__ == "__main__":


    dcomps_perx = dcomb.groupby("element.X")
    for xatom in dcomps_perx.groups.keys():
        data = dcomps_perx.get_group(xatom)
        plot_bgmagmom(
            data,
            dedge.loc[data.index],
            suptitle=f"$Cs_2[B,B']{xatom}_6$",
            filename=anal_dir / f"BG_Mag_{xatom}transition.png",
        )
        plot_bgmagmom(
            data,
            dedge.loc[data.index],
            suptitle=f"$Cs_2[B,B']{xatom}_6$",
            use_transition=False,
            filename=anal_dir / f"BG_Mag_{xatom}blockpairs.png",
        )
    #
    # Make BG-MAGMOM plot for paper
    data = dcomb.groupby("element.X").get_group("Cl")
    plot_bgmagmom(
        data,
        dedge.loc[data.index],
        filename=anal_dir / "BG_Mag_Cl_Transition_Print.png",
    )