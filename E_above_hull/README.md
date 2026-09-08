# E_above_hull

This directory contains scripts, input tables, jobflow configuration, and saved results for calculating energies with machine-learned interatomic potentials (MLIPs), then using those energies to construct phase diagrams and calculate formation energies and energies above the convex hull for HDP structures.

The main workflow is:

1. Read structures and chemical-system metadata from CSV files.
2. Relax structures with an MLIP through `atomate2`, `jobflow`, and `pymatgen`.
3. Save relaxed energies in MLIP-specific JSON dictionaries.
4. Combine the energies with structure metadata in `PD_Construction.ipynb`.
5. Construct phase diagrams and save e-above-hull and formation-energy tables.
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

### `run_ehull_calcs.py`

Runs one batch of MLIP relaxations for structures listed in `MissingStrucs_<MLIP>.csv`. It converts serialized structure dictionaries into `pymatgen.Structure` objects, calls `mlip_relax_batched`, and merges the returned energies into `ExpensiveMLIPs/mpid_energy_dict_<MLIP>.json` while holding a file lock.

Command-line arguments:

```text
python run_ehull_calcs.py <MLIP> <batch size> <batch number>
```

### `drive_ehull_calcs.py`

Splits the missing-structure input for an MLIP across workers. Each worker calls `run_ehull_calcs.py` for every batch assigned to it.

```text
python drive_ehull_calcs.py <worker id> <number of workers> <MLIP>
```

### `run_batch_ehull.py`

A related batch driver that reads the full `hdp_mp_subsysandhdps.csv` input rather than the MLIP-specific missing-structure CSV. It uses batched relaxation and prints the returned energy dictionary; its JSON-writing code is currently commented out.

### `run_conv_ehull.py`

Runs a batch using `mlip_relax_and_get_energies` rather than the batched helper. It reads `MissingStrucs_<MLIP>.csv`, converts structures to conventional cells, and merges results into an MLIP JSON dictionary under a file lock.

### `drive_eqv3_calcs.py`

Worker launcher for an Equiformer-v3-style calculation path. It calls `run_eqv3_ehull.py`; that target script is not present in this directory, so this launcher is not self-contained here.

### `drive_grace_calcs.py`

Worker launcher for GRACE calculations. It reads the full subsystem table to determine the number of batches, starts at batch offset 135, and calls `run_grace_calcs.py`.

### `run_grace_calcs.py`

Runs GRACE relaxation batches for the full subsystem table and merges results into `ExpensiveMLIPs/mpid_energy_dict_<MLIP>.json`.

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
- `HDP_Ehull_overview.csv`: summary statistics for formation energy and e-above-hull results, including means, standard deviations, entry counts, and stability thresholds at 100, 150, and 200 meV/atom.

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
- `compile_nequip.sh`: Slurm script that compiles the NequIP-OAM-XL model for CUDA/ASE use.

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
