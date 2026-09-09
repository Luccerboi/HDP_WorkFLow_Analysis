import pandas as pd
import numpy as np
from pymatgen.core import Species


def comparative_histograms(
    hdp_df: pd.DataFrame,
    base_indices: list | pd.Index | pd.Series,
    compare_indices: list | pd.Index | pd.Series | dict,
    base_label: str = "full set",
    compare_label: str = "unstable comps",
    normalize: bool = False,
    **kwargs,
):
    import plotly.graph_objects as go
    from plotly.express import colors
    from pymatgen.core import Element

    col = hdp_df.columns.to_list()[0]
    hdp_df.fillna({"element.B2": "Vac"}, inplace=True)
    base_count_b1 = hdp_df.loc[base_indices].groupby("element.B1").count()[col]
    base_count_b2 = hdp_df.loc[base_indices].groupby("element.B2").count()[col]

    elem_list = list(set(base_count_b2.index.to_list() + base_count_b1.index.to_list()))
    try:
        elem_map = {elem: Element(elem).Z for elem in elem_list if elem != "Vac"}
    except ValueError:
        elem_map = {elem: Species(elem).Z*2 + Species(elem).oxi_state for elem in elem_list if elem != "Vac"}
        # print(elem_map)
    elem_map.update({"Vac": 1})

    elem_idx = pd.Index(dict(sorted(elem_map.items(), key=lambda item: item[1])).keys())
    # print(elem_idx)
    elem_idx = pd.Index([str(x) for x in elem_idx])

    base_count_ser = base_count_b1.reindex(elem_idx).fillna(0) + base_count_b2.reindex(
        elem_idx
    ).fillna(0)

    if isinstance(compare_indices, dict):
        df_dict = {}
        for label, idx in compare_indices.items():
            compare_count_b1 = hdp_df.loc[idx].groupby("element.B1").count()[col]
            # print(compare_count_b1)
            # print(compare_count_b1.reindex(elem_idx).fillna(0))
            compare_count_b2 = hdp_df.loc[idx].groupby("element.B2").count()[col]
            compare_count_ser = pd.Series(data=np.zeros(len(elem_idx)),index=elem_idx).add(compare_count_b1,fill_value=0).add(compare_count_b2,fill_value=0)
            # print(compare_count_ser.head(10))
            # print(compare_count_ser.index)
            df_dict.update({label: compare_count_ser})
        # print(df_dict)
        counter_df = pd.DataFrame(df_dict)
        # Source - https://stackoverflow.com/a/22650162
        # Posted by 8one6, modified by community. See post 'Timeline' for change history
        # Retrieved 2026-08-31, License - CC BY-SA 4.0

        counter_df = counter_df.loc[~(counter_df == 0).all(axis=1)]

        new_idx = pd.Index(dict(sorted({elem: elem_map[elem] for elem in counter_df.index}.items(), key= lambda item: item[1])).keys()),
        if '5LOB' not in counter_df.columns:
            total_counts = counter_df.loc[new_idx].fillna(0).sum(axis=1)
        else:
            total_counts = counter_df.loc[new_idx]['5LOB']

        counter_df = counter_df.loc[new_idx]
        # print(counter_df.head())
        base_count_ser = base_count_ser.loc[new_idx]
        elem_idx = new_idx
        # print(elem_idx)
        # print(counter_df.loc[new_idx]['1Rel'])
        # print(base_count_ser.loc[new_idx])
        # print(counter_df.head(10))

    else:
        compare_count_b1 = (
            hdp_df.loc[compare_indices].groupby("element.B1").count()[col]
        )
        compare_count_b2 = (
            hdp_df.loc[compare_indices].groupby("element.B2").count()[col]
        )
        compare_count_ser = compare_count_b1.reindex(elem_idx).fillna(
            0
        ) + compare_count_b2.reindex(elem_idx).fillna(0)

        counter_df = pd.DataFrame(
            {compare_label: compare_count_ser}
        )


    if normalize:

        fig = go.Figure(
            data=[
                # go.Bar(
                #     x=elem_idx,
                #     y=counter_df[base_label] / np.sum(counter_df[base_label]) * 100,
                #     name=base_label,
                #     marker_color="#575554",
                # ),
                go.Bar(
                    x = [f"{elem} tot: {int(count)}" for elem, count in total_counts.items()],
                    y=counter_df[key] / base_count_ser * 100,
                    name=key,
                    # marker_color="#E55107",
                )
                for key in counter_df.columns
            ]
        )
        fig.update_layout(
            xaxis={
                "title": {"text": "Element", "font": {"size": 28}},
                "tickfont": {"size": 18},
            },
            yaxis={
                "title": {"text": "Occurrence Share (%)", "font": {"size": 28}},
                "tickfont": {"size": 22},
            },
            barmode="overlay",
        )

    else:
        fig = go.Figure(
            data=[
                go.Bar(
                    x=counter_df.index,
                    y=base_count_ser,
                    name=base_label,
                    marker_color="#8D8D8D",
                ),
            ]
        )
        fig.add_traces(
            data=[
                go.Bar(
                    x = counter_df.index,
                    y = counter_df[key],
                    name = key,
                )
                    for key in counter_df.columns
            ]
        )
        fig.update_layout(
            barmode="overlay",
            xaxis={
                "title": {"text": "Element", "font": {"size": 28}},
                "tickfont": {
                    "size": 17,
                    #  'family': "Arial Black"
                },
            },
            yaxis={
                "title": {"text": "Occurrence Count", "font": {"size": 28}},
                "tickfont": {"size": 22},
            },
        )

    fig.update_layout(
        legend={
            "font": {"size": 22},
            "x": 0.98,
            "y": 0.98,
            "xanchor": "right",
            "yanchor": "top",
            "bgcolor": "rgba(255,255,255,0.7)",
            "bordercolor": "black",
            "borderwidth": 1,
        },
        colorway=colors.qualitative.Plotly
    )
    fig.update_layout(
        margin={'t':0,'l':0,'b':0,'r':0}
    )


    fig.update_traces(opacity=1)
    return fig


def stable_counts(
    ehull_df: pd.DataFrame, cutoff_values: list[float] = [100, 150, 200]
) -> pd.DataFrame:
    counts_df = pd.DataFrame(index=ehull_df.index)
    counts_df["MLIP_entries"] = ehull_df.count(axis=1)
    for cutoff in cutoff_values:
        val = cutoff / 1000
        counts_df[f"Ehull <= {cutoff}meV"] = ehull_df.where(ehull_df < val).count(
            axis=1
        )

    return counts_df


if __name__ == "__main__":

    hdp_df = pd.read_csv(
        "../E_above_hull/HDP_CombinedInfo_WithStructures.csv",
        index_col=0,
    )
    hdp_df.set_index("comp", inplace=True)

    ehull_df = pd.read_csv(
        "../AnalysisResults/HDP_Ehull_overview.csv",
        index_col=0,
    )

    fig1 = comparative_histograms(
        hdp_df=hdp_df,
        base_indices=ehull_df["Ehull_avg"].dropna().index,
        # compare_indices=ehull_df.dropna()[ehull_df["Ehull_avg"] > 0.15].index,
        compare_indices=ehull_df.dropna()[ehull_df["Ehull <= 150meV"] <= 2].index,
        base_label="Full Dataset",
        compare_label="Unstable HDPs",
        normalize=False,
    )
    fig1.update_layout(height=500)
    fig1.write_image("ehull_images/MLIPStab_Ehull150.png", width=1800, height=500, scale=2.2)
    fig1.write_html("ehull_images/MLIPStab_Ehull150.html")
    fig1.show()

    input_data = pd.read_csv(
        "../WorkFlow/ImplementedWorkFlow/UsedInput_HDPLedger_BkCfFmMd_removed.csv",
        header=None,
    )
    input_data["compID"] = input_data.apply(lambda row: f"{row[7]}_{row[0]}", axis=1)
    input_data.set_index("compID", inplace=True)
    # input_data["element.B1"] = Species(f"{input_data[2]}{input_data[3]}+")

    input_data["element.B1"] = input_data.apply(lambda row: f"{row[2]}{row[3]}+" if row[2]!="Vac" else "Vac",axis=1)
    input_data["element.B2"] = input_data.apply(lambda row: f"{row[4]}{row[5]}+" if row[4]!="Vac" else "Vac",axis=1)
    input_data["element.X"] = input_data[6]
    # print(input_data)
    d1 = pd.read_csv(
        "../HDP_WorkFlow_Analysis/WorkFlow/ImplementedWorkFlow/HDPLedger_NoCs6s.csv",
        index_col=[0, 1],
    ).T
    d1 = d1.loc[input_data.index]
    # converged_idx = d1[d1[("3Pre", "completed")].str.contains("1")].index
    # nonconverged_idx = pd.Index(set(d1.index) - set(converged_idx))
    conv_overview = {
        step: d1[d1[(step, "completed")].str.contains("1")].index
        for step in sorted(list(set(d1.columns.get_level_values(0))), reverse=True)
    }
    nonconv_overview = {
        key: pd.Index(set(d1.index) - set(val)) for key, val in conv_overview.items()
    }
    # print(nonconv_overview)

    fig2 = comparative_histograms(
        input_data,
        d1.index,
        nonconv_overview,
        base_label="All input HDPs",
        compare_label="Non-converged HDPs",
        normalize=False,
    )
    fig2.update_layout(legend={"title": {"text": "Calc. Step:", "font": {"size": 28}}})
    fig2.update_xaxes(tickfont={"size": 24})
    fig2.update_layout(
        xaxis={"title": {"font": {"size": 32}}},
        yaxis={"title": { "font": {"size": 32}}},
        legend={"font":{"size":26}},
    )

    fig2.write_image(
        "ehull_images/NonConvOverview_CalcStep.png", width=1800, height=900, scale=2.2
    )
    fig2.write_html("ehull_images/NonConvOverview_CalcStep.html")
    fig2.show()

    xcount = input_data.loc[nonconv_overview['5LOB']].groupby('element.X')
    xidx = ['F','Cl','Br','I']
    # xcoded_index = {xel : xcount.get_group(xel).index for xel in xidx}
    xcum_index = {}
    for ii in range(len(xidx)):
        if ii == 0:
            xcum_index[xidx[ii]] = pd.Index(list(xcount.get_group(xidx[ii]).index))
        else:
            xcum_index[xidx[ii]] = pd.Index(list(xcount.get_group(xidx[ii]).index) + list(xcum_index[xidx[ii-1]]))
    print(xcum_index)
    xcum_reversed = dict(sorted(xcum_index.items(), key=lambda item: len(item[1]),reverse=True))

    fig4 = comparative_histograms(
        input_data,
        d1.index,
        xcum_reversed,
        base_label="All Input HDPs",
        compare_label="Non-converged HDPs",
        normalize=False,
    )
    fig4.update_layout(legend={"title": {"text": "X-Element:", "font": {"size": 26}}})
    fig4.update_xaxes(tickfont={"size": 24})
    fig4.update_layout(
        xaxis={"title": {"font": {"size": 32}}},
        yaxis={"title": { "font": {"size": 32}}},
        legend={"font":{"size":26}},
    )
    # fig4.update_layout(barmode="stack")
    fig4.show()
    fig4.write_image("ehull_images/NonConvOverview_Xel.png",width=1800,height=900,scale=2.2)
    fig4.write_html("ehull_images/NonConvOverview_Xel.html")

    import plotly.express as px


    fig3 = px.scatter(
        x=hdp_df.loc[ehull_df.dropna().index]["tau_factor"],
        y=ehull_df.dropna()["Ehull_avg"],
        color=hdp_df.loc[ehull_df.dropna().index]["element.X"],
        hover_name= ehull_df.dropna().index.to_list(),
    )
    fig3.update_layout(
        # xaxis={
        # "title": {"text": "$\\tau-\\text{factor}$", "font": {"size": 22}},
        # "tickfont": {"size": 18},
        # },
        yaxis={
            "title": {"text": "Energy above Hull (eV/atom)", "font": {"size": 28}},
            "tickfont": {"size": 22},
            "maxallowed": 1.0,
            "minallowed": -0.01,
        },
        legend={
            "title": {"text": "  X-site:", "font": {"size": 23}},
            "font": {"size": 22},
            "x": 0.98,
            "y": 0.98,
            "xanchor": "right",
            "yanchor": "top",
            "bgcolor": "rgba(255,255,255,0.7)",
            "bordercolor": "black",
            "borderwidth": 1,
        },
        height=500,
    )
    fig3.update_xaxes(
        # title_text=r"$\tau-\text{factor}$",
        # title_font={"size": 60},
        tickfont={"size": 22},
    )
    fig3.update_xaxes(
        title={
            "text": "<i>τ</i> - factor",
            "font": {"size": 28},
        }
    )
    fig3.write_image("ehull_images/Ehull_vs_Tau.png", width=1800, height=500, scale=2.2)
    fig3.write_html("ehull_images/Ehull_vs_Tau.html")
    fig3.show()
