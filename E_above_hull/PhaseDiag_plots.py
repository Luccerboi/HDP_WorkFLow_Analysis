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


def plot_phasediagram(comp_id: str, data_dir: str | PathLike, mlip: str | int = 0):

    data_file = glob.glob(pathname=f'{comp_id}*',root_dir=data_dir)[0]

    pd_data = loadfn(os.path.join(data_dir,data_file))



    # fig = go.Figure()
    if isinstance(mlip, str):
        try:
            mlip_name = mlip
            fig = pd_data[mlip_name].get_plot()
        except Exception as e:
            raise KeyError(f"failed to load {mlip} phasediagram for {comp_id}.\n Got Exception {e}")
    elif isinstance(mlip,int):
        mlip_list = list(pd_data.keys())
        mlip_name = mlip_list[mlip]
        fig = pd_data[mlip_name].get_plot()

    fig.update_layout({'title':f"{comp_id} PhaseDiagram {mlip_name}"})



    fig.show()
    return pd_data

