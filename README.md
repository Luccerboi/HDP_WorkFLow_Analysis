# HDP WorkFlow Analysis

This repository contains tools used for the creation of our database published as [Spin-Polarized Electronic Structure and Chemical Bonding Data for 2,500+ Halide Double Perovskites](LINK TO OUR PAPER). It contains functionality for creating input files, monitoring completion, and analyzing results. The main automation logic lives under the `WorkFlow/` folder. The full dataset is available over [NOMAD](LINK TO NOMAD DOI)

## Repository structure

- `AnalysisResults/` - Contains output data: Summarized data in .csv format. Scatter plots for spin magnetic moment vs. bandgap. Periodic Table heatmaps with several quantities averaged per element.
- `HDP_PermutationMaker/` - Scripts used to generate all possible Halide Double Perovskite permutations, based on the tolerance factor from [Bartel et al.](https://www.science.org/doi/10.1126/sciadv.aav0693).
- `PlottingScripts/` - plotting functions and notebooks for visualization.
- `WorkFlow/` - main workflow automation engine for job creation, submission, queue management, and monitoring.

## `WorkFlow/` contents

- `.env` - environment configuration used by the workflow scripts.
- `AllCompsInput` - input compound list used to initialize the workflow ledger.
- `CsBBX_Analyzer/` - additional analyzer scripts or tools.
- `HDPCompound.py` - compound-level logic, including POTCAR generation and compound metadata.
- `HTProcessLedger.py` - process ledger implementation for job state tracking and queue management.
- `JobDefinitions.py` - job class definitions and submission template generation for VASP/LOBSTER jobs.
- `RunWorkFlow.py` - main workflow driver script that assigns jobs, submits them, and monitors progress.
- `WatchCompletion.py` - monitoring utility that displays current queue status and running jobs.
- `job_handler.py` - SLURM helper functions for job submission and availability checks.
- `testinput` - sample or test input data for the workflow.
- `ImplementedWorkFlow/` - workflow database directory containing job JSON, ledger files, and backups.
- `assets/` - auxiliary assets: LOBSTER basisfunction set YAML files, POTCAR look-up-tables linking element to specific VASP PAWs, and high and low spin-only magnetic moments for all d- and f- elements.

### How the workflow works

The workflow is centered in `WorkFlow/RunWorkFlow.py`:

1. `RunWorkFlow.py` loads environment variables from `WorkFlow/.env` using `python-dotenv`.
2. It reads `JOB_JSON_PATH`, `DATABASE_PATH`, and `PROCESSLEDGER_FILENAME` to locate the job definitions and ledger.
3. `HTProcessLedger.ProcessLedger` is used to restore or initialize the workflow state and retrieve the compound list.
4. The workflow periodically checks cluster availability using `job_handler.check_availability()`.
5. `AssignNextJobs()` chooses jobs from the queue based on available resources and node requirements.
6. Jobs are created from classes in `JobDefinitions.py`, which wrap VASP or LOBSTER execution steps.
7. Assigned jobs are executed in parallel using `concurrent.futures.ProcessPoolExecutor` up to `MAX_NUM_JOBS`.
8. After each job completes, the ledger is backed up and the queue is refreshed.

### Job definitions

`JobDefinitions.py` defines the common `GeneralJob` class and specialized subclasses for different workflow stages such as relaxations, spin scans, HSE, and LOBSTER.

- It reads `VASP_COMMAND` and `LOBSTER_COMMAND` from `.env`.
- It writes submission scripts, INCAR/LOBSTERIN files, and manages job-specific behavior.

### Compound and ledger management

- `HDPCompound.py` handles compound data and file generation.
  - It reads `PATH_TO_PSEUDO` from `.env` to locate the VASP PAW pseudopotentials.
- `HTProcessLedger.py` stores workflow state, job progress, and queue assignments.
- `AllCompsInput` is used to initialize the ledger when starting from scratch. **NOTE:** Due to some errors in previous iterations of the HDP permutation determination, there are some duplicate entries in this input liste (e.g., Cs2PdHgCl6 and Cs2HgPdCl6), or compositions that were assumed stable based on a different oxidation state assignment before using the difference in ElectronNegativies in the permutation maker. These duplicate/unstable compositions have been removed from the results, but are still included in the used input and process ledger.

### How `.env` works

The `.env` file in `WorkFlow/` contains runtime settings for the workflow scripts. Each script loads this file at startup and reads environment variables through `os.environ`.

Current `.env` variables:

- `ALLOWED_PARTITIONS` - Python list of SLURM partitions to consider for job scheduling.
- `DEFAULT_USERNAME` - username used by `WatchCompletion.py` to query `squeue`.
- `MAX_NUM_JOBS` - maximum number of concurrent jobs the workflow will attempt to run.
- `PATH_TO_PSEUDO` - path to the VASP PAW pseudopotential directory.
- `LOBSTER_COMMAND` - command or path to the LOBSTER executable.
- `VASP_COMMAND` - command used to run VASP (e.g. `mpirun vasp_std`).
- `JOB_JSON_PATH` - path to the JSON job settings file used by the workflow.
- `DATABASE_PATH` - path to the workflow database directory where ledger and outputs are stored.
- `PROCESSLEDGER_FILENAME` - filename of the CSV ledger used to track jobs.

### Important notes

- The `.env` file is required by the workflow scripts. If variables are missing, the scripts will fail when trying to read the configuration.
- The workflow should normally be launched from the `WorkFlow/` directory or with the working directory set so that `python-dotenv` can locate `WorkFlow/.env`.
- Keep `.env` private and system-specific, since it contains local paths and cluster settings.

## Running the workflow

From the `WorkFlow/` directory:

```bash
cd /home/lwalterb/hdp_project/HDP_WorkFlow_Analysis/WorkFlow
python RunWorkFlow.py
```

To monitor progress in a separate terminal:

```bash
python WatchCompletion.py
```

## Quick start

1. Open `WorkFlow/.env` and update the values for your environment:
   - `JOB_JSON_PATH`
   - `DATABASE_PATH`
   - `PROCESSLEDGER_FILENAME`
   - `VASP_COMMAND`
   - `LOBSTER_COMMAND`
   - `PATH_TO_PSEUDO`
   - `DEFAULT_USERNAME`
   - `MAX_NUM_JOBS`

2. Confirm the workflow data folder exists and contains the files referenced by `DATABASE_PATH`.
3. If you are starting from scratch, make sure `AllCompsInput` is populated with the compounds you want to run.
4. Launch the workflow from `WorkFlow/`:

```bash
cd /home/lwalterb/hdp_project/HDP_WorkFlow_Analysis/WorkFlow
python RunWorkFlow.py
```

5. Use `WatchCompletion.py` in another terminal to watch queue progress and running jobs.

6. To restart the workflow and reload the ledger, stop the current run, change the `Restart` parameter to `True` in `RunWorkFlow.py` or ledger state as needed, and rerun `python RunWorkFlow.py`. **Note**: Setting Restart to True will delete all data from the DataBase path.

### Notes

- `RunWorkFlow.py` uses `job_handler.py` to submit jobs to SLURM and to query job states.
- `WatchCompletion.py` is intended to show current jobs and queue progress, but it depends on the same `.env` settings.

## Plotting functions in `PlottingScripts/`

The Python scripts in `PlottingScripts/` provide plotting utilities for both interactive and publication-style figures.

- `HDP_PlotslyPlots.py`
  - `plot_dos()` — builds a plotly DOS figure for a single compound, showing total DOS and projected site contributions with the VBM aligned at 0 eV.
  - `plot_coxx()` — creates COHP/COBI plots for a compound, returning separate figures for B1 and B2 bonding contributions.
  - `plot_ptable()` — generates a periodic table heatmap using [pymatviz](https://github.com/janosh/pymatviz/tree/main) from `CombinedInfo`, supporting element counts, average ICOHP/ICOBI, band gap, charge, and other aggregate metrics. 
  - These functions were mainly used for the [Interactive UMAP Explorer](https://github.com/Luccerboi/HDP_UMAP_Explorer). The DOS/COHP/COBI plots are omitted in this repository, but can easily be reproduced using the functionalities contained here. All periodic table heatmap plots can be found in `AnalysisResults/hdp_plotly_plots/`.

- `CsBBX_Plotting.py`
  - `plot_bgmagmom()` — plots net spin magnetic moment versus bandgap, including bandgap and magnetization histograms, with optional coloring by transition character or block pairing. All bandgap-magnetic moment scatter plots for each halide anion are contained in `AnalysisResults/BandGap_Magmom_ScatterPlots/`.
  - `compare_alattice()` — compares simulated HDP lattice constants against ICSD experimental values, filters cubic entries, and produces comparison summary output. Was used to create the structural comparison figure in the publication.
  ### Implemented plots that were not used in the end:
  - `plot_icohpbg()` — plots ICOHP or ICOBI versus bandgap, optionally grouped by transition character or block pairing and split by halide species. 
  - `plot_rowicohp()` — plots average ICOHP/ICOBI by B-site row across different halides, useful for comparing bonding trends across the periodic table.
  - `plot_icohpviolin_perspecies()` — makes violin plots of ICOHP/ICOBI distributions for each atomic species and halide type.

- `PtablePlot.ipynb`
  - Was used to create the periodic table plot showing which elements are featured on which sites in all possible permutations.
  - Also used to create the bar plot showing the B-B block pairing counts for the stable Cs-HDPs.
  - Both these figures are featured in Figure 1 of [our publication](LINK TO OUR PUBLICATION)

- `ChargeSpill_FigMaker.ipynib`
  - Was used to make the histogram showing the charge spilling from the LOBSTER projection for all compositions.

## Data in `AnalysisResults`

When the workflow was run initially, there were some duplicate entries (e.g., Cs2HgPdCl6 and Cs2PdHgCl6) and we later added the method of using differences in ElectroNegativities to assign oxidation states (same methodology as Bartel et al.) to B-sites that had multiple options to satisfy charge neutrality. This later addition made that some compositions were later predicted to be unstable based on updated ionic radii. These duplicate and unstable compositions have been removed from all results, but are still featured in the UsedInput and ProcessLedger. This means we have:
- `WorkFlow/AllCompsInput` and `/WorkFlow/ImplementedWorkFlow/UsedInput_HDPLedger_NoCs6s.csv` contain the duplicates and later unstable predicted compositions.
- `Completed_HDPComps_NoDuplicates.txt` is the list of compositions that completed the workflow after removing duplicate and (unstable) unstable compositions.
- `DuplicatedOrUnstable_HDPComp_list.txt` the list of compositions that completed the workflow but were actually duplicate or (predicted) unstable compositions.
- `Duplicated_Unstable_HDPs_CombinedInfo.csv` the info contained in the CombinedInfo.csv for the compositions that were later filtered out based on the updates to the permutation maker.
- `NonConverged_HDPComp_list.txt` list of all compositions that did not complete the full workflow. All ran out of the maximum number of tries in some step of the workflow. The full data can be extracted by retrieving the CompletionOverview form the ProcessLedger. This list may still contain some duplicates or (predicted) unstable compositions.

- Several .csv files containing data extracted using `CsBBX_Analyzer.py` and outlined in the [Data Descriptor](LINK TO OUR PAPER)
  - `HDP_BasicInfo_260510.csv` Contains info on input species, used LOBSTER basisset, and some elemental data.
  - `HDP_StructuralInfo_260510.csv` Contains some interatomic distances from the relaxed *Fm3m* structure.
  - `HDP_bandedgeInfo_lsodos_260510.csv` Contains VBM, CBM, bandgaps, and pDOS contributions to 0.5eV around the bandedges for each element. Analyzed from DOSCAR.lso.lobster file. Data is given for spin-up, spin-down, and combined spin-channels.
  - `HDP_LobsterInfo_260510.csv` Contains data extracted from LOBSTER projection (ICOHP/ICOBI values), and a bonding descriptor describing asymmetry in ICOHP/ICOBI along x,y,z,-axes.

- `FoundICSD_HDPs.txt` gives an overview of all ICSD entries we found matching the Cesium-HDP chemical formula. It features CollectionCode, Chemical formula, Title, Authors, and Reference for all found matches.