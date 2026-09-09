# E_above_hull

This directory contains scripts, input tables, jobflow configuration, and saved results for calculating energies with machine-learned interatomic potentials (MLIPs), then using those energies to construct phase diagrams and calculate formation energies and energy-above-hull for HDP structures.

The main workflow is:

1. Read structures and chemical-system metadata from CSV files.
2. Relax structures with an MLIP through `atomate2`, `jobflow`, and `pymatgen`.
3. Save relaxed energies in MLIP-specific JSON dictionaries.
4. Combine the energies with structure metadata in `PD_Construction.ipynb`.
5. Construct phase diagrams and save energy-above-hull and formation-energy tables.
6. Optionally save serialized phase diagrams for plotting or later inspection.

## Python utilities

### `ehull_utils.py`

Shared functions for the calculation and analysis workflow:

- `get_chemical_subsystems` generates all non-empty chemical subsystems for a hyphen-separated chemical system.
- `mlip_relax_and_get_energies` runs separate relaxation jobs and returns an energy or failure value for each structure.
- `mlip_relax_batched` runs a batch of structures in one force-field job and returns energies keyed by structure identifier.
- `mlip_relax_2step` is an experimental two-stage relaxation helper. It is currently incomplete and returns an empty dictionary after printing final relaxation information.
- `combine_batches_to_series` combines JSON files produced by separate batch runs.
- `collect_mlip_energies_to_df` adds `E_<MLIP>` columns to a structures DataFrame from MLIP JSON files.
- `construct_phase_diagrams` builds `pymatgen` phase diagrams and calculates formation energies and e-above-hull values.
- `summarize_results` prints missing-energy counts and the fractions of structures below 100, 150, and 200 meV/atom above hull.

## Relaxation and batch drivers

The calculation scripts are organized into two general roles:

### `run_*.py` scripts

The `run_*.py` scripts execute a particular batch of MLIP calculations. They generally:

1. Read either the full subsystem table (`hdp_mp_subsysandhdps.csv`) or an MLIP-specific missing-structure table (`MissingStrucs_<MLIP>.csv`).
2. Convert serialized structure dictionaries into `pymatgen.Structure` objects.
3. Select the MLIP model and its calculator settings from an `mlip_specifications` mapping.
4. Relax the structures using either the batched or individual relaxation helper in `ehull_utils.py`.
5. Print progress and, for the production batch paths, merge energies into an MLIP-specific JSON dictionary while holding a file lock.

The scripts differ by calculation variant rather than by overall workflow:

- Some use `mlip_relax_batched` to process a list of structures as one job.
- Some use `mlip_relax_and_get_energies` to create separate jobs for each structure.
- Some convert structures to conventional cells before relaxation.
- Some target a particular MLIP or resume a particular portion of the full subsystem calculation.
- The batch-oriented path writes results to `ExpensiveMLIPs/mpid_energy_dict_<MLIP>.json`; other scripts may print results or have output code disabled for testing.

The usual command-line interface for a batch script is:

```text
python run_<variant>.py <MLIP> <batch size> <batch number>
```

The exact variant name and supported MLIP names must be checked in the selected script's `mlip_specifications` mapping. The available `run_*.py` files include the standard e-hull batch path, a conventional-cell/convergence path, a full-subsystem batch path, and the GRACE-specific path.

### `drive_*.py` scripts

The `drive_*.py` scripts distribute batches across workers. Each driver reads the relevant input table, calculates the number of batches, assigns every worker a subset of batch numbers, and launches a corresponding `run_*.py` script through `subprocess.run`.

The usual command-line interface is:

```text
python drive_<variant>.py <worker id> <number of workers> <MLIP>
```

The Slurm submission scripts in this directory use these drivers to start multiple workers on GPU nodes. `check_completion.sh` can be used to inspect the resulting jobs and output files.

## Analysis notebooks and plotting

### `PD_Construction.ipynb`

Primary analysis notebook. It loads the HDP and subsystem tables, reads the four MLIP energy JSON files in `MLIP_Results/`, creates the combined energy table, constructs phase diagrams, writes e-above-hull and formation-energy CSVs to `../AnalysisResults`, saves phase-diagram data under `MLIP_Results/PhaseDiagramData`, and generates MLIP-specific missing-structure CSVs.

The MLIPs used in the notebook are:

- `GRACE-3L-OMAT-large-ft-AM`
- `SevenNet-MPALOE`
- `SevenNet-omat24`
- `PET-OAM-XL`

### `query_strucs.ipynb`

Exploratory notebook for querying and preparing structure tables. It reads the combined HDP data and subsystem-related CSV data, including `HDP_CombinedInfo_WithStructures.csv`.

### `PhaseDiag_plots.py`

Provides `plot_phasediagram`, which loads a serialized phase diagram from a file in a supplied directory, selects an MLIP by name or by index, displays the Plotly figure, and returns the loaded data.

## Input CSV and JSON data

### Structure and HDP tables

- `HDP_CombinedInfo_WithStructures.csv`: combined HDP information, including composition, structural/electronic descriptors, bonding-related columns, and a serialized `structure_dict` column.
- `hdp_mp_subsysandhdps.csv`: subsystem structure table used by the relaxation and phase-diagram workflows. Its key columns are `structure_dict`, `chemsys`, and `nsites`.
- `hdp_subsystemstosimulate_mp.json`: JSON input describing the subsystem systems selected for simulation.

### Missing-structure tables

- `MissingStrucs_GRACE-3L-OMAT-large-ft-AM.csv`
- `MissingStrucs_PET-OAM-XL.csv`
- `MissingStrucs_SevenNet-MPALOE.csv`
- `MissingStrucs_SevenNet-omat24.csv`

Each contains `chemsys`, `nsites`, and `structure_dict` for structures without a usable energy for that MLIP. These files are consumed by `run_ehull_calcs.py` and `run_conv_ehull.py`.

### Summary tables

- `HDP_Ehull_MissingStrucs.csv`: MLIP-by-MLIP missing-entry information, with MLIP columns and a `composition` column.
- `HDP_Ehull_overview.csv`: summary statistics for formation energy and energy-above-hull results, including means, standard deviations, entry counts, and stability thresholds at 100, 150, and 200 meV/atom.

### `MLIP_Results/`

Saved energy and phase-diagram results:

- `hdp_mlipenergies_4MLIPS.csv`: combined subsystem table with energy columns for the four MLIPs used in `PD_Construction.ipynb`.
- `mpid_energy_dict_<MLIP>.json`: structure-ID-to-energy mappings for each MLIP.
- `PhaseDiagramData/`: serialized phase-diagram data written for individual HDP compositions.

The relaxation scripts may also create an `ExpensiveMLIPs/` directory at runtime. It is not included in the attached directory listing but is referenced as the default output location by those scripts.

## Batch submission and monitoring

- `submit_grace.sh`: Slurm setup and launch for GRACE.
- `submit_pet.sh`: Slurm setup and launch for PET-OAM-XL with four concurrent workers.
- `submit_svnnet.sh`: Slurm setup and launch for SevenNet-MPALOE with four concurrent workers.
- `submit_svnomat.sh`: Slurm setup and launch for SevenNet-omat24 with four concurrent workers.
- `check_completion.sh`: reports running Slurm jobs, completed batches, failed entries, runtime, and recent log output.

The submission scripts activate MLIP-specific micromamba environments, set CUDA/CPU environment variables, select a jobflow configuration, and launch the corresponding driver. Their `JOBFLOW_CONFIG_FILE` values point to `/home/lwalterb/hdp_ehull/jobflow_settings/...`; verify those absolute paths before using the scripts from another checkout or directory layout.

## Environments and jobflow configuration

### `env_requirements/`

Pinned package requirements for the MLIP environments:

- `reqs_grace_p11.txt`
- `reqs_pet.txt`
- `reqs_sevennet.txt`

### `jobflow_settings/`

Jobflow store configurations for the four MLIP variants:

- `grace_jobflow.yaml`
- `pet_jobflow.yaml`
- `7net_mpaloe_jobflow.yaml`
- `7netomat_jobflow.yaml`

Each file configures a `MontyStore` document store and blob/additional store. The configured database paths are outside this repository under `/home/lwalterb/hdp_ehull/jobstores/`.

## Notes on running

- Run scripts from `E_above_hull/` unless you update their relative paths.
- The relaxation scripts expect the requested MLIP name to match one of the keys in their `mlip_specifications` dictionaries.
- MLIP relaxation requires the corresponding environment, model support, and usually a CUDA-capable GPU.
- Failed convergence is represented in intermediate dictionaries by strings containing `Failed`; the analysis notebook masks those entries before phase-diagram construction.
- Phase diagrams are skipped when required energies are missing or when an elemental endpoint is in the explicitly excluded atomic-number list in `construct_phase_diagrams`.
