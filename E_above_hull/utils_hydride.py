def mlip_relax_and_get_energies(force_field_name: str | MLFF,
                                structure_dict: dict[str, Structure],
                                calculator_kwargs: dict ={"default_dtype": "float64"}) -> dict:  # calc kwargs need to be adapted for nequip
    ff_maker = ForceFieldRelaxMaker(force_field_name=force_field_name,
                                    calculator_kwargs=calculator_kwargs)

    job_uuid_to_idfr = {}
    jobs = []
    for mpid, structure in structure_dict.items():
        job = ff_maker.make(structure)
        jobs.append(job)
        job_uuid_to_idfr[job.uuid] = mpid

    flow = Flow(jobs)
    resp = run_locally(flow)

    flow_output = {}
    for uuid, mpid in job_uuid_to_idfr.items():
        output = resp[uuid][1].output
        if (
                output is not None
                and output.is_force_converged
                # and abs(output.output.energy) < energy_tol
        ):
            # Beware: differently structured dict than in direct phase diagram computation util below
            flow_output[mpid] = output.output.energy
    return flow_output


def make_comp_phase_diagram(
        elemental_endpoints: list[str] | tuple[str] | set[str],
        pd_entries: list[PDEntry],
        force_field_name: str | MLFF,
        structures: dict[str, Structure] | None = None,
        energy_tol: float | None = 1e6,
):
    """
    Modified version of Aaron Kaplan's tutorial
    https://github.com/Matgenix/atomate2_cecam_school/blob/main/notebooks/hands_on_examples/force_field_basics.ipynb
    without pulling structures via API.
    Also, avoid unnecessary relaxation calls by storing PD entries per chemical system.
    """
    energy_tol = energy_tol or float("inf")

    if structures:
        ff_maker = ForceFieldRelaxMaker(force_field_name=force_field_name,
                                        calculator_kwargs={"default_dtype": "float64"})

        job_uuid_to_idfr = {}
        jobs = []
        for idfr, structure in structures.items():
            job = ff_maker.make(structure)
            jobs.append(job)
            job_uuid_to_idfr[job.uuid] = idfr

        flow = Flow(jobs)
        resp = run_locally(flow)

        flow_output = {}
        for uuid, idfr in job_uuid_to_idfr.items():
            output = resp[uuid][1].output
            if (
                    output is not None
                    and output.is_force_converged
                    and abs(output.output.energy) < energy_tol
            ):
                flow_output[idfr] = {
                    "energy": output.output.energy,
                    "composition": output.structure.composition,
                }

        pd_entries.extend([PDEntry(entry["composition"], entry["energy"], name=mpid)for mpid, entry in flow_output.items()])

    pd = PhaseDiagram(pd_entries, [Element(e) for e in elemental_endpoints])

    pd_entries = [PDEntry.from_dict(pde) for pde in pd.as_dict()["all_entries"]]
    for i, pde in enumerate(pd_entries):
        pd_entries[i].entry_id = pde.name

    pd = PhaseDiagram(pd_entries, [Element(e) for e in elemental_endpoints])

    hull_dict = {e.entry_id: pd.get_decomp_and_e_above_hull(entry=e, check_stable=False)[1] for e in pd.all_entries}

    fig = pd.get_plot(ternary_style="3d")
    fig.update_layout(
        height=600,
        width=1000,
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(color="black")
    )
        
    return fig, pd, hull_dict