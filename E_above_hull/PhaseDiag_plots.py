import plotly.graph_objects as go
import pandas as pd
import numpy as np
import os
from os import PathLike
from pathlib import Path
import glob
from copy import deepcopy
from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PhaseDiagram


def plot_phasediagram(comp_id: str, data_dir: str | PathLike):

    data_file = glob.glob(pathname=f'{comp_id}*',root_dir=data_dir)[0]

    pd_data = loadfn(os.path.join(data_dir,data_file))

    button_list = []
    for mlip, diagram in pd_data.items():
        button_list.append(
            dict(
                label = mlip,
                method = "update",
                args = [{"data" : diagram.get_plot().to_dict()['data'], "layout": diagram.get_plot().to_dict()['layout']}],
            )
        )

    # fig = go.Figure()
    fig = pd_data['SevenNet-MPALOE'].get_plot()

    fig.update_layout(
    updatemenus=[
        dict(
            type="buttons",
            direction="right",
            active=0,
            x=0.57,
            y=1.2,
            buttons= button_list,
        )
    ])

    fig.show()
    return

# def plot_phasediagram(comp_id: str, data_dir: str | PathLike):

#     data_file = glob.glob(pathname=f'{comp_id}*',root_dir=data_dir)[0]

#     pd_data = loadfn(os.path.join(data_dir,data_file))

#     fig = pd_data['SevenNet-MPALOE'].get_plot()

#     false_list = [False] * (len(pd_data.keys()) + 1)
#     ii = 1

#     button_list = []

#     for mlip, diagram in pd_data.items():
#         fig.add_trace(go.Figure(
#             data=diagram.get_plot().to_dict()['data'],
#             layout=diagram.get_plot().to_dict()['layout'],
#         ))

#         vis_list = deepcopy(false_list)
#         vis_list[ii] = True
#         button_list.append(
#             dict(
#                 label = mlip,
#                 method = "update",
#                 args = [{"visible": vis_list, "title": mlip}],
#             )
#         )

#     # fig = go.Figure()
    

#     fig.update_layout(
#     updatemenus=[
#         dict(
#             type="buttons",
#             direction="right",
#             active=0,
#             x=0.57,
#             y=1.2,
#             buttons= button_list,
#         )
#     ])

#     fig.show()
#     return