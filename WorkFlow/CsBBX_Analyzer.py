"""
CsBBX_Analyzer.py
Contains the SingleHDPanalysis and GroupedAnalysis classes. SingleHDPanalysis is used to collect and analyse outputs from a single HDP compound.
GroupedAnalysis is used to analyse all compounds from the DataBase.
The main loop contains several calls to run analysis functions.
"""

from HDPCompound import Compound
from HTProcessLedger import ProcessLedger, initialize_compounds
import pandas as pd
import numpy as np
import os
from mendeleev.fetch import fetch_ionic_radii
from monty.io import zopen
from os import PathLike
from tqdm.autonotebook import tqdm
from pathlib import Path
import json
import glob
from typing import Literal
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import multiprocessing as mp
from monty.json import MontyDecoder, MontyEncoder, jsanitize

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from scipy.ndimage import gaussian_filter1d
import warnings
import time

from lobsterpy.plotting import (
    InteractiveCohpPlotter,
    IcohpDistancePlotter,
    PlainCohpPlotter,
    PlainDosPlotter,
    get_style_list,
)
from lobsterpy.featurize.core import FeaturizeIcoxxlist
from pymatgen.electronic_structure.dos import (
    Dos,
    CompleteDos,
    _lobster_orb_labs,
    LobsterCompleteDos,
)
from pymatgen.core.periodic_table import Element, Species

from pymatgen.core import Structure, get_el_sp
from pymatgen.core.spectrum import Spectrum
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.io.lobster import (
    Charge,
    Grosspop,
    Icohplist,
    Doscar,
    MadelungEnergies,
    SitePotential,
)
from pymatgen.util.coord import get_linear_interpolated_value
from pymatgen.util.typing import SpeciesLike

warnings.filterwarnings("ignore")

# ============================================================================
# Use python-dotenv to load settings from .env file.
# ============================================================================
from dotenv import load_dotenv

load_dotenv()
JobsFile = os.environ["JOB_JSON_PATH"]
DBlocation = os.environ["DATABASE_PATH"]
ledgerfile = os.environ["PROCESSLEDGER_FILENAME"]

# Load ProcessLedger and regain the List of Compounds.
p = ProcessLedger(JobsFile, StartPath=DBlocation, ledger_filename=ledgerfile)
loc = p.RestartLedger()
_ = p.ReadFullLedger()


# Code for caluclating tau-factor and generalized t-factor
irs_full = (
    fetch_ionic_radii(radius="ionic_radius")[
        [
            "I",
            "II",
            "III",
            "IV",
            "V",
            "VI",
            "VII",
            "VIII",
            "IX",
            "X",
            "XI",
            "XII",
            "XIV",
        ]
    ]
    .interpolate(axis=1, limit_direction="both")
    .dropna(axis=0)
)
irs_full = irs_full.drop([1, 1])


class SingleHDPanalysis:
    """This class contains functions used to collect and analyze data for a spicific composition
    from the database.
    """

    def __init__(self, comp: Compound, LedgerInfo: pd.Series):
        """Initializes the class aimed at analyzing data for a specific compostion from the HDP database.
        Paths to relevant output files will be set depending on the job completion.
        Input data and orbitals will be set as attributes.

        Args:
            comp (Compound): Instance of Compound class containing the elemental and input info.
            LedgerInfo (pd.Series): Row of the ProcessLedger corresponding to *comp* containing the required computational info.
        """
        self.comp = comp
        self.ledgerinfo = LedgerInfo
        self.joblist = LedgerInfo.index.get_level_values(0).unique().to_list()
        self.jobcompletion = LedgerInfo.xs("completed", level=1).to_dict()
        self.jobpaths = LedgerInfo.xs("JobPath", level=1).to_dict()
        self._suborbitals = {}

        for orbtype in OrbitalType:
            self._suborbitals[orbtype] = [
                x for x in Orbital if str(x).startswith(str(orbtype))
            ]

        if self.jobcompletion["1Rel"] == "1":
            path_to_poscar = Path(self.jobpaths["1Rel"]) / "CONTCAR"
            self.structure = Structure.from_file(path_to_poscar)

        if self.jobcompletion["4HSE"] == "1":
            self.path_to_doscar = Path(self.jobpaths["4HSE"]) / "DOSCAR"

        if bool(int(self.jobcompletion["5LOB"])):
            self.used_basis = "basisset0"
            self.path_to_lsodos = (
                Path(self.jobpaths["5LOB"]) / self.used_basis / "DOSCAR.LSO.lobster"
            )
            self.path_to_lobdos = (
                Path(self.jobpaths["5LOB"]) / self.used_basis / "DOSCAR.lobster"
            )
        self.input_dict = self._get_WFinput_dict()

    def update_lobbasis(self, basis: str):
        """This function will be used by GroupedHDPanalysis in order to change the basisset of the LOBSTER projection if required.
        Used to change the paths to LOBSTER output files for analysis functions.

        Args:
            basis (str): The name of the basis set to be used (basisset*).
                         Needs to match the subdirectory name for downstream functionality!
        """
        self.used_basis = basis
        self.path_to_lsodos = (
            Path(self.jobpaths["5LOB"]) / self.used_basis / "DOSCAR.LSO.lobster"
        )
        self.path_to_lobdos = (
            Path(self.jobpaths["5LOB"]) / self.used_basis / "DOSCAR.lobster"
        )

    def _get_WFinput_dict(self) -> dict:
        """Will return a dict of the input information used to start the workflow.

        Returns:
            dict: Dictionary of inputs, specifically: input species, used POTCARS, high/low spinstates, LOBSTER basisset overview.
        """

        spinsB1 = list(self.comp.get_spinstate(which="B1"))
        if self.comp.B2 != "Vac":
            spinsB2 = list(self.comp.get_spinstate(which="B2"))
        else:
            spinsB2 = [None, None]

        WFinput: dict = {
            "input_species": {
                "A": Species(self.comp.A, oxidation_state=1),
                "B1": Species(self.comp.B1, oxidation_state=self.comp.nB1),
                "B2": (
                    Species(self.comp.B2, oxidation_state=self.comp.nB2)
                    if self.comp.B2 != "Vac"
                    else None
                ),
                "X": Species(self.comp.X, oxidation_state=-1),
            },
            "potcar_list_frozen": self.comp.write_potcar(
                use_frozen=True, write_file=False
            ),
            "potcar_list_nofrozen": self.comp.write_potcar(
                use_frozen=False, write_file=False
            ),
            "expected_spinstates": {
                "B1": {"HighSpin": spinsB1[0], "LowSpin": spinsB1[1]},
                "B2": {"HighSpin": spinsB2[0], "LowSpin": spinsB2[1]},
            },
            "all_lob_basis": self.comp.get_all_basis_combos(),
            "used_basis_num": self.used_basis,
            "used_basis_funcs": self.comp.get_all_basis_combos()[self.used_basis],
        }

        return WFinput

    def get_ionicradii_based_data(self) -> dict:
        """Calculates Bartel's tau-factor, and M.R. Filip's generalized t-factor with geometric limits for this composition.

        Returns:
            dict: Dictionary of used ionic radii, tau-factor, t-factor, octahedral factor, octahedral mismatch,
                    whether geometric limits are broken, and if so, which.
        """

        def calculate_tau(nA, r_A, r_B, r_X):
            tau = (r_X / r_B) - nA * (nA - (r_A / r_B) / np.log(r_A / r_B))
            return tau

        def calculate_gen_t(r_A, mu_avg, delta_mu, r_X):
            t = (r_A / r_X + 1) / (2 * (mu_avg + 1) ** 2 + delta_mu**2) ** (1 / 2)
            return t

        def calculate_geom_stability(t, mu_avg, delta_mu):
            SL = bool(t <= 1)
            OL = bool(mu_avg >= np.sqrt(2) - 1 + delta_mu)
            TL1 = bool(t >= (0.44 * mu_avg + 1.37) / (np.sqrt(2) * (mu_avg + 1)))
            TL2 = bool(t >= (0.73 * mu_avg + 1.13) / (np.sqrt(2) * (mu_avg + 1)))
            CL1 = bool(t <= 2.46 / (2 * (mu_avg + 1) ** 1 + delta_mu**2) ** (1 / 2))
            CL2 = bool(mu_avg <= 1.14)
            if False in [SL, OL, TL1, TL2, CL1, CL2]:
                false_index = [SL, OL, TL1, TL2, CL1, CL2].index(False)
                broken_cond = ["SL", "OL", "TL1", "TL2", "CL1", "CL2"][false_index]
                return False, broken_cond
            else:
                return True, None

        species = self.input_dict["input_species"]
        r_A = irs_full.loc[(species["A"].Z, species["A"].oxi_state), "XII"].astype(
            np.float64
        )
        r_B1 = irs_full.loc[(species["B1"].Z, species["B1"].oxi_state), "VI"].astype(
            np.float64
        )
        r_B2 = (
            irs_full.loc[(species["B2"].Z, species["B2"].oxi_state), "VI"].astype(
                np.float64
            )
            if species["B2"] is not None
            else r_B1
        )
        r_X = irs_full.loc[(species["X"].Z, species["X"].oxi_state), "VI"].astype(
            np.float64
        )

        tau = calculate_tau(1, r_A, (r_B1 + r_B2) / 2, r_X)

        mu_avg = (r_B1 + r_B2) / (2 * r_X)
        delta_mu = np.abs(r_B1 - r_B2) / (2 * r_X)
        gen_t = calculate_gen_t(r_A, mu_avg, delta_mu, r_X)

        geom_stabil, condition = calculate_geom_stability(gen_t, mu_avg, delta_mu)

        return {
            "r_A": r_A,
            "r_B1": r_B1,
            "r_B2": r_B2 if species["B2"] is not None else None,
            "r_X": r_X,
            "oct_factor": mu_avg,
            "oct_mismatch": delta_mu,
            "gen_t_factor": gen_t,
            "geom_stable": geom_stabil,
            "broken_condition": condition,
            "tau_factor": tau,
        }

    ###DOS related methods

    def _get_bandedges_eigenval(
        self, cald_dir: PathLike, filename: str = "EIGENVAL"
    ) -> dict:
        """Determines VBM and CBM edges from the EIGENVAL file in stead of DOSCAR.
        Was used to investigate some discrepancies, but was not used for the final analysis.

        Args:
            cald_dir (PathLike): Path to the calculation containing EIGENVAL file
            filename (str, optional): File to read and determine VBM/CBM energies from (must be parsed like VASP EIGENVAL). Defaults to "EIGENVAL".

        Raises:
            ValueError: If filename could not be resolved

        Returns:
            dict: Dictionary of VBM/CBM energy and corresponding k-point
        """
        if cald_dir.is_dir():
            file = cald_dir / filename
        elif cald_dir.is_file():
            file = cald_dir
        else:
            raise ValueError(
                f"Could Not read read EIGENVALUE file {filename} from {cald_dir} for {self.comp}"
            )
        # print(file)
        with zopen(file, "rt") as f:
            lines = f.readlines()
        vbm_up = -100.0
        vbm_down = -100.0
        kp_vbm_up = []
        kp_vbm_down = []
        cbm_up = 100.0
        cbm_down = 100.0
        kp_cbm_up = []
        kp_cbm_down = []
        for line in lines[7:]:
            newline = line.strip().split()
            values = list(map(float, newline))
            # print(values)
            # print(len(newline))
            if len(newline) == 4:
                cur_kpoint = values[:3]

            if len(newline) == 5:
                if values[3] > 0.1 and values[1] >= vbm_up:
                    # Spin up is occupied
                    vbm_up = values[1]
                    kp_vbm_up = cur_kpoint
                if values[4] > 0.1 and values[2] >= vbm_down:
                    # Spin down is occupied
                    vbm_down = values[2]
                    kp_vbm_down = cur_kpoint
                if values[3] < 0.1 and values[1] <= cbm_up:
                    # Spin up is unoccupied
                    cbm_up = values[1]
                    kp_cbm_up = cur_kpoint
                if values[4] < 0.1 and values[2] <= cbm_down:
                    # Spin down is unoccupied
                    cbm_down = values[2]
                    kp_cbm_down = cur_kpoint
        return {
            Spin.up: {
                "VBM": vbm_up,
                "kp_vbm": kp_vbm_up,
                "CBM": cbm_up,
                "kp_cbm": kp_cbm_up,
            },
            Spin.down: {
                "VBM": vbm_down,
                "kp_vbm": kp_vbm_down,
                "CBM": cbm_down,
                "kp_cbm": kp_cbm_down,
            },
        }

    def _get_bandedges_dos(
        self, dos_type: Literal["vasp", "lso", "lobster"], tol: float = 1e-3
    ) -> dict:
        """Determine VBM CBM edges for spin-up and spin-down separately based on the total DOS.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): Type of DOS to analyze (from VASP or LOBSTER (standard or LSO basis))
            tol (float, optional): Tolerance to use when determining whether states are occupied. Defaults to 1e-3.

        Raises:
            ValueError: DOSCAR type was not understood, could not open.

        Returns:
            dict: Dictionary containing VBM CBM energy and the index of edges for each spin channel separately.
        """
        tdos = self.get_tdos(dos_type)
        if dos_type.lower() == "vasp":

            efermi = tdos["efermi"]
        elif dos_type.lower() == "lso":
            efermi = 0
        elif dos_type.lower() == "lobster":
            efermi = 0
        else:
            raise ValueError(f"Error processsing DOSCAR dos_type {dos_type}")

        edge_dict = {}
        for spin, tdensity in tdos["densities"].items():
            below_fermi = [
                i
                for i in range(len(tdos["energies"]))
                if tdos["energies"][i] < efermi and tdensity[i] > tol
            ]
            above_fermi = [
                i
                for i in range(len(tdos["energies"]))
                if tdos["energies"][i] > efermi and tdensity[i] > tol
            ]
            # print(self.comp.CompID, below_fermi,tdensity[5000:5300])
            vbm_start = max(below_fermi)
            vbm_energy = tdos["energies"][vbm_start]
            cbm_start = min(above_fermi)
            cbm_energy = tdos["energies"][cbm_start]

            edge_dict.update(
                {
                    spin: {
                        "VBM": vbm_energy,
                        "i_vbm": vbm_start,
                        "CBM": cbm_energy,
                        "i_cbm": cbm_start,
                    }
                }
            )

        return edge_dict

    def get_rel_dos_contributions(
        self,
        spin: Spin | None,
        emin: float,
        emax: float,
        dos_type: Literal["vasp", "lso", "lobster"],
        tol: float = 1e-3,
    ) -> dict:
        """Analyze the relative contribution from each species to the total DOS.
        If a species contributes more than 1% of the total states it will also determine contributions from individual bands.
        Function analyzes one spin-channel at a time, or sums boths spins together before analyzing

        Args:
            spin (Spin | None): Which spin channel to analyze, Spin.up or Spin.down. If None then the DOS for spin-up and down will be summed.
            emin (float): Lowest energy to include for analysis
            emax (float): Highest energy to include for analysis
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): Which DOSCAR file to use for analysis
            tol (float, optional): Tolerance when determining which states are occupied. Defaults to 1e-3.

        Returns:
            dict: Dictionary containing the relative contribution of each species, and their respective orbitals
        """

        if spin == None:
            summed_spins = True
            spin = Spin.up
        else:
            summed_spins = False

        assert emin < emax

        tdos = self.get_tdos(dos_type, summed_spins)
        tdos_species = self.get_tdos_per_species(dos_type, summed_spins)
        spddos_species = self.get_spddos_per_species(dos_type, summed_spins)

        E_full = tdos["energies"]

        include_index = [
            i
            for i in range(len(E_full))
            if (E_full[i] >= emin and E_full[i] <= emax)
            and tdos["densities"][spin][i] > tol
        ]

        total_states = np.sum(tdos["densities"][spin][include_index])
        # print(f'total_states: {total_states} for {spin}')
        rel_contr_species = {}
        tdos_species.pop("efermi")
        tdos_species.pop("energies")
        for specie, tdos_spec in tdos_species.items():

            specie_tot_contr = (
                np.sum(tdos_spec["densities"][spin][include_index]) / total_states
            )
            # print(f'specie {specie} contributes {specie_tot_contr}')
            # print(f"{specie} has {specie_tot_contr} states")
            if specie_tot_contr > 0.01:
                rel_contr_species.update({specie: {"tot": specie_tot_contr}})
                for band, pdos in spddos_species[specie].items():
                    if band in ["energies", "efermi"]:
                        continue
                    band_contr = (
                        np.sum(pdos["densities"][spin][include_index]) / total_states
                    )
                    rel_contr_species[specie].update({band: band_contr})

        return rel_contr_species

    def analyze_bandedges(
        self,
        dos_type: Literal["vasp", "lso", "lobster"],
        use_eigenval: bool = False,
        integration_width=0.5,
        tol: float = 1e-3,
    ) -> dict:
        """Calls self._get_bandedges_[dos/eigenval] and self.get_rel_dos_contributions() and collects the relevant info in a dict.
        The dict contains the VBM, CBM, bandgap, and relative contributions for each spin-channel separately and a combined
        channel where the overall bandgap is determined and the corresponding orbital contributions are saved.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): Which DOSCAR type is used for analysis
            use_eigenval (bool, optional): Whether to use EIGENVAL for VBM and CBM determination. Defaults to False.
            integration_width (float, optional): The width of analysis window from band edgeswith which to determine relative contributions. Defaults to 0.5.
            tol (float, optional): Tolerance to use when determining which states are occupied. Defaults to 1e-3.

        Returns:
            dict: Dictionary containing VBM, CBM and bandgap energies, and valence and conduction states.
                  high level keys are ['Spin.up','Spin.down','combined'] where combined contains the info from the spin channel with
                  highest VBM and lowest CBM.
        """

        if use_eigenval:
            print(f"Note {self.comp.CompID} is using EIGENVAL for band edge detection")
            edgedict = self._get_bandedges_eigenval(Path(self.jobpaths["4HSE"]))
        else:
            edgedict = self._get_bandedges_dos(dos_type)

        bandedge_info = {Spin.up: {}, Spin.down: {}, "combined": {}}

        # dummy values to determine which spin has the lowest CBM and highest VBM
        Evbm = -100
        Ecbm = 100

        for spin in edgedict.keys():
            # print(f'doing valence for {spin}')
            valence_states = self.get_rel_dos_contributions(
                spin,
                edgedict[spin]["VBM"] - integration_width,
                edgedict[spin]["VBM"],
                dos_type,
                tol,
            )
            # print(f'doing cond for {spin}')
            conduction_states = self.get_rel_dos_contributions(
                spin,
                edgedict[spin]["CBM"],
                edgedict[spin]["CBM"] + integration_width,
                dos_type,
                tol,
            )

            bandedge_info[spin].update({"VBM": edgedict[spin]["VBM"]})
            bandedge_info[spin].update({"CBM": edgedict[spin]["CBM"]})
            gap = edgedict[spin]["CBM"] - edgedict[spin]["VBM"]
            bandedge_info[spin].update({"bandgap": gap})

            if use_eigenval:  # we have KPOINT information in this case, so lets safe it
                bandedge_info[spin].update({"kp_vbm": edgedict[spin]["kp_vbm"]})
                bandedge_info[spin].update({"kp_cbm": edgedict[spin]["kp_cbm"]})
                bandedge_info[spin]["direct"] = (
                    bandedge_info[spin]["kp_vbm"] == bandedge_info[spin]["kp_cbm"]
                )

            bandedge_info[spin].update(
                {
                    "valence_states": valence_states,
                    "conduction_states": conduction_states,
                }
            )

            if edgedict[spin]["VBM"] > Evbm:
                Evbm = edgedict[spin]["VBM"]
                bandedge_info["combined"].update(
                    {
                        "spin_vbm": spin,
                        "VBM": edgedict[spin]["VBM"],
                        "valence_states": valence_states,
                    }
                )

            if edgedict[spin]["CBM"] < Ecbm:

                Ecbm = edgedict[spin]["CBM"]
                bandedge_info["combined"].update(
                    {
                        "spin_cbm": spin,
                        "CBM": edgedict[spin]["CBM"],
                        "conduction_states": conduction_states,
                    }
                )

            combined_gap = (
                bandedge_info["combined"]["CBM"] - bandedge_info["combined"]["VBM"]
            )
            bandedge_info["combined"]["bandgap"] = combined_gap
            bandedge_info["combined"]["spin_forbidden"] = (
                bandedge_info["combined"]["spin_vbm"]
                != bandedge_info["combined"]["spin_cbm"]
            ) and (
                (
                    abs(bandedge_info[Spin.up]["VBM"] - bandedge_info[Spin.down]["VBM"])
                    > 0.05
                )
                or (
                    abs(bandedge_info[Spin.up]["CBM"] - bandedge_info[Spin.down]["CBM"])
                    > 0.05
                )
            )

            if use_eigenval:

                bandedge_info["combined"].update(
                    {
                        "kp_vbm": edgedict[bandedge_info["combined"]["spin_vbm"]][
                            "kp_vbm"
                        ]
                    }
                )
                bandedge_info["combined"].update(
                    {
                        "kp_cbm": edgedict[bandedge_info["combined"]["spin_cbm"]][
                            "kp_cbm"
                        ]
                    }
                )
                bandedge_info["combined"]["direct"] = (
                    bandedge_info[spin]["kp_vbm"] == bandedge_info[spin]["kp_cbm"]
                )

        return bandedge_info

    def analyze_edgecontr(
        self,
        dos_type: Literal["vasp", "lso", "lobster"],
        use_eigenval: bool = False,
        integration_width=0.5,
        tol: float = 1e-3,
    ) -> dict:
        """Similar to self.analyze_bandedges() only saves the relative contributions of each species in a way thats easier to combine with pd.Dataframes.
        Instead of nested dictionairy containing orbital contributions, it saves the total contribution of a species and the largest contributing orbital.

        Calls self._get_bandedges_[dos/eigenval] and self.get_rel_dos_contributions() and collects the relevant info in a dict.
        The dict contains the VBM, CBM, bandgap, and relative contributions for each spin-channel separately and a combined
        channel where the overall bandgap is determined and the corresponding orbital contributions are saved.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): Which DOSCAR type is used for analysis
            use_eigenval (bool, optional): Whether to use EIGENVAL for VBM and CBM determination. Defaults to False.
            integration_width (float, optional): The width of analysis window from band edgeswith which to determine relative contributions. Defaults to 0.5.
            tol (float, optional): Tolerance to use when determining which states are occupied. Defaults to 1e-3.

        Returns:
            dict: Dictionary containing VBM, CBM and bandgap energies, and valence and conduction states.
                  high level keys are ['Spin.up','Spin.down','combined'] where combined contains the info from the spin channel with
                  highest VBM and lowest CBM.
        """
        edge_dict = self.analyze_bandedges(
            dos_type, use_eigenval, integration_width, tol
        )
        input_species = self._get_WFinput_dict()["input_species"]

        output_dict = {}

        for channel, edges in edge_dict.items():
            vbm_states = edges["valence_states"]
            cbm_states = edges["conduction_states"]
            spin_dict = {}
            spin_dict.update(
                {"VBM": edges["VBM"], "CBM": edges["CBM"], "bandgap": edges["bandgap"]}
            )

            if channel == "combined":
                spin_dict.update(
                    {
                        "spin_vbm": edges["spin_vbm"],
                        "spin_cbm": edges["spin_cbm"],
                        "spin_forbidden": edges["spin_forbidden"],
                    }
                )

            for site, species in input_species.items():
                try:
                    if species == None:
                        raise KeyError
                    site_vbmstates = vbm_states[site]

                except KeyError:
                    spin_dict.update(
                        {
                            f"VBMtotcontr.{site}": 0,
                            f"VBMorbital.{site}": None,
                            f"VBMorbchar.{site}": None,
                            f"VBMorborder.{site}": None,
                            f"VBMorbcontr.{site}": 0,
                        }
                    )

                else:
                    tmpdict = site_vbmstates.copy()
                    tmpdict.__delitem__("tot")
                    orb = ""
                    contr = 0.0

                    for orbit, contri in tmpdict.items():
                        # find max contributing orbital of this species
                        if contri > contr:
                            contr = contri
                            orb = orbit

                    spin_dict.update(
                        {
                            f"VBMtotcontr.{site}": vbm_states[site]["tot"],
                            f"VBMorbital.{site}": orb,
                            f"VBMorbchar.{site}": orb[-1],
                            f"VBMorborder.{site}": int(orb[0]),
                            f"VBMorbcontr.{site}": contr,
                        }
                    )

                try:
                    if species == None:
                        raise KeyError
                    site_cbmstates = cbm_states[site]

                except KeyError:
                    spin_dict.update(
                        {
                            f"CBMtotcontr.{site}": 0,
                            f"CBMorbital.{site}": None,
                            f"CBMorbchar.{site}": None,
                            f"CBMorborder.{site}": None,
                            f"CBMorbcontr.{site}": 0,
                        }
                    )

                else:
                    tmpdict = site_cbmstates.copy()
                    tmpdict.__delitem__("tot")
                    orb = ""
                    contr = 0.0

                    for orbit, contri in tmpdict.items():
                        # find max contributing orbital of this species
                        if contri > contr:
                            contr = contri
                            orb = orbit

                    spin_dict.update(
                        {
                            f"CBMtotcontr.{site}": cbm_states[site]["tot"],
                            f"CBMorbital.{site}": orb,
                            f"CBMorbchar.{site}": orb[-1],
                            f"CBMorborder.{site}": int(orb[0]),
                            f"CBMorbcontr.{site}": contr,
                        }
                    )
            output_dict.update({channel: spin_dict})

        return output_dict

    def _parse_doscar(self, path_to_doscar: PathLike) -> dict:
        """Parses DOSCAR into dictionary containing total DOS and pDOS per orbital for each site.
        Sub-orbitals (i.e., 2p_x, 2p_y, and 2p_z) are summed together in to 2p

        Args:
            path_to_doscar (PathLike): Path to the DOSCAR type to parse

        Raises:
            Exception: Something went wrong while parsing the DOSCAR file

        Returns:
            dict: Collection of energies, total DOS, pDOS per orbital per site.
        """
        with zopen(path_to_doscar, "rt") as f:
            lines = f.readlines()
        headline = lines[5].strip().split()
        nedos = int(headline[2])
        e_fermi_vasp = float(headline[3])

        if "LOBSTER" in lines[4]:  # LOBSTER DOS EFERMI HS BEEN SHIFTED
            EFERMI = 0
        else:
            EFERMI = e_fermi_vasp

        energies = np.zeros(nedos)
        tdos_up = np.zeros(nedos)
        tdos_down = np.zeros(nedos)
        for nval, line in enumerate(lines[6 : 6 + nedos]):
            words = line.strip().split()
            energies[nval] = words[0]
            tdos_up[nval] = words[1]
            tdos_down[nval] = words[2]
        tdos = {
            "efermi": EFERMI,
            "energies": energies,
            "densities": {Spin.up: tdos_up, Spin.down: tdos_down},
        }
        site_dos_dict = {}
        for natom, site in enumerate(self.structure.sites):
            site_dos_dict.update({natom: {"species": site.species}})
            line_offset = 6 + (nedos + 1) * (
                natom + 1
            )  # offset of lines to next projected DOS
            if "LOBSTER" in lines[4]:
                orbline = lines[line_offset - 1].strip().split(";")[-1].strip().split()
                norbs = len(orbline)
                mainbands = []
                suborbdict = {}
                for orb in orbline:
                    band = orb.split("_")[0]
                    if band not in mainbands:
                        mainbands.append(band)
                        # print(band)
                        suborbs = [x for x in Orbital if str(x).startswith(band[-1])]
                        # print(suborbs)
                        suborbdict.update({band: suborbs})

            else:
                suborbdict = self._suborbitals

            projdos_lines = lines[line_offset : line_offset + nedos]
            nbands = int((len(projdos_lines[2].strip().split()) - 1) / 2)
            projdos_energies = np.zeros(nedos)
            projdos_spinup = np.zeros((nedos, nbands))
            projdos_spindown = np.zeros((nedos, nbands))
            for nval, nextline in enumerate(projdos_lines):
                try:
                    line = nextline.strip().split()
                    ncols = len(line)
                    spin_up_cols = np.arange(1, ncols, 2)
                    spin_down_cols = np.arange(2, ncols, 2)
                    projdos_energies[nval] = line[0]
                    projdos_spinup[nval, :] = [line[x] for x in spin_up_cols.tolist()]
                    projdos_spindown[nval, :] = [
                        line[x] for x in spin_down_cols.tolist()
                    ]
                except Exception as exc:
                    raise Exception(
                        f"{self.comp.CompID} raised an exception when parsing {path_to_doscar} at line {line_offset + nval}: {exc}"
                    )

            maxorbs = len(spin_up_cols.tolist())
            iband = 0
            # print(f"Parsing Doscar {str(path_to_doscar).split('/')[-1]}, for {site.species_string} if have {ncols} columns, {len(spin_up_cols)} spinup cols, and bands is {suborbdict}")
            for mainband, suborb_list in suborbdict.items():
                # maxorbs -= len(suborb_list)
                # if maxorbs < 0:
                #     break
                site_dos_dict[natom].update({mainband: {}})
                for suborb in suborb_list:

                    site_dos_dict[natom][mainband].update(
                        {
                            suborb: {
                                "efermi": EFERMI,
                                "energies": projdos_energies,
                                "densities": {
                                    Spin.up: projdos_spinup[:, iband],
                                    Spin.down: projdos_spindown[:, iband],
                                },
                            }
                        }
                    )
                    iband += 1

        full_dos_dict = {
            "e_fermi_proper": EFERMI,
            "e_fermi_vasp": e_fermi_vasp,
            "energies": energies,
            "tdos": tdos,
            "dos_per_orbital": site_dos_dict,
        }
        return full_dos_dict

    def get_tdos(
        self, dos_type: Literal["vasp", "lso", "lobster"], summed_spins: bool = False
    ) -> dict:
        """Parses the DOSCAR and returns the total DOS.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): DOSCAR from which basis should be analyzed.
            summed_spins (bool, optional): Whether to sum the spin-up and down channels. Defaults to False.

        Raises:
            ValueError: DOSCAR type was not understood

        Returns:
            dict: Containing the total DOS for spin up and down separately or summed.
        """
        if dos_type.lower() == "vasp":
            full_dos_dict = self._parse_doscar(self.path_to_doscar)
        elif dos_type.lower() == "lso":
            full_dos_dict = self._parse_doscar(self.path_to_lsodos)
        elif dos_type.lower() == "lobster":
            full_dos_dict = self._parse_doscar(self.path_to_lobdos)
        else:
            raise ValueError(f"Error processsing DOSCAR dos_type {dos_type}")

        if summed_spins:
            dens_sep = full_dos_dict["tdos"]["densities"]

            tdos_dict = full_dos_dict["tdos"]["densities"] = {
                Spin.up: dens_sep[Spin.up] + dens_sep[Spin.down]
            }
            return tdos_dict

        else:
            return full_dos_dict["tdos"]

    def get_spddos_per_site(
        self, dos_type: Literal["vasp", "lso", "lobster"], summed_spins: bool = False
    ) -> dict:
        """Parse the DOSCAR and return the pDOS per orbital for each site in the composition.
        If VASP DOSCAR is analyzed orbitals are resolved on the level of s,p,d,f.
        If LOBSTER or LOBSTER.LSO DOSCAR is analyzed also orbital order is resolved.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): The DOSCAR of which basis should be analyzed.
            summed_spins (bool, optional): Whether spin-up and down states should be summed. Defaults to False.

        Raises:
            ValueError: Given DOSCAR type was not understood.

        Returns:
            dict[site_index, orbital, Dos]: Dictionary of site resolved orbital pDOS.
        """
        if dos_type.lower() == "vasp":
            full_dos_dict = self._parse_doscar(self.path_to_doscar)
        elif dos_type.lower() == "lso":
            full_dos_dict = self._parse_doscar(self.path_to_lsodos)
        elif dos_type.lower() == "lobster":
            full_dos_dict = self._parse_doscar(self.path_to_lobdos)
        else:
            raise ValueError(f"Error processsing DOSCAR dos_type {dos_type}")

        energies = full_dos_dict["energies"]
        dos_per_orbit = full_dos_dict["dos_per_orbital"]
        efermi = full_dos_dict["e_fermi_proper"]
        site_spd_dos = {"efermi": efermi, "energies": energies}
        for natom, nested_dict in dos_per_orbit.items():
            site_spd_dos[natom] = {"efermi": efermi, "energies": energies}
            for nested_key, nested_val in nested_dict.items():
                if nested_key == "species":
                    continue
                # site_spd_dos.update({natom: {'nested_key})
                site_spd_dos[natom].update({nested_key: {}})
                sum_up = np.zeros(len(energies))
                sum_down = np.zeros(len(energies))

                for suborbital in nested_val.keys():
                    density = nested_val[suborbital]["densities"]

                    if summed_spins:
                        sum_up += density[Spin.up] + density[Spin.down]

                    else:
                        sum_up += density[Spin.up]
                        sum_down += density[Spin.down]

                if summed_spins:
                    densities = {Spin.up: sum_up}
                else:
                    densities = {Spin.up: sum_up, Spin.down: sum_down}
                site_spd_dos[natom][nested_key].update(
                    {
                        "efermi": full_dos_dict["e_fermi_proper"],
                        "energies": energies,
                        "densities": densities,
                    }
                )

        return site_spd_dos

    def get_tdos_per_site(
        self, dos_type: Literal["vasp", "lso", "lobster"], summed_spins: bool = False
    ) -> dict:
        """Retrieve the 'total' pDOS per site, determined by summing all orbital pDOS.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): The DOSCAR of which basis should be analyzed.
            summed_spins (bool, optional): Whether spin-up and down states should be summed. Defaults to False.

        Returns:
            dict[site_index,Dos]: Dictionary containing total pDOS per site.
        """

        site_spd_dos = self.get_spddos_per_site(dos_type, summed_spins=False)
        energies = site_spd_dos["energies"]
        efermi = site_spd_dos["efermi"]

        site_tdos = {"efermi": efermi, "energies": energies}

        site_spd_dos.pop("efermi")
        site_spd_dos.pop("energies")
        for natom, dos_per_band in site_spd_dos.items():
            site_tdos.update({natom: {}})
            dos_per_band.pop("efermi")
            dos_per_band.pop("energies")
            sum_up = 0
            sum_down = 0
            for band, dos in dos_per_band.items():
                densitie = dos["densities"]

                sum_up += densitie[Spin.up]
                sum_down += densitie[Spin.down]

            if summed_spins:
                new_dens = {Spin.up: sum_up + sum_down}
            else:
                new_dens = {Spin.up: sum_up, Spin.down: sum_down}

            site_tdos[natom].update(
                {"efermi": efermi, "energies": energies, "densities": new_dens}
            )

        return site_tdos

    def get_spddos_per_species(
        self, dos_type: Literal["vasp", "lso", "lobster"], summed_spins: bool = False
    ) -> dict:
        """Retrieve orbital resolved pDOS per species (A/B1/B2/X). Retrieved by summing pDOS over symmetrically equivalent sites.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): The DOSCAR of which basis should be analyzed.
            summed_spins (bool, optional): Whether spin-up and down states should be summed. Defaults to False.

        Returns:
            dict[(A/B1/B2/X), orbital, Dos]: Dictionary containing the orbital resolved pDOS per species of this composition.
        """

        dos_per_site = self.get_spddos_per_site(dos_type, summed_spins=False)
        energies = dos_per_site["energies"]
        efermi = dos_per_site["efermi"]

        specie_spd_dos = {"efermi": efermi, "energies": energies}
        dos_per_site.pop("energies")
        dos_per_site.pop("efermi")

        if self.comp.B2 == "Vac":
            assert len(dos_per_site.keys()) == 9
            specie_sites = {"B1": [0], "A": [1, 2], "X": [3, 4, 5, 6, 7, 8]}
        else:
            assert len(dos_per_site.keys()) == 10
            specie_sites = {"B1": [0], "B2": [1], "A": [2, 3], "X": [4, 5, 6, 7, 8, 9]}

        for specie, sites in specie_sites.items():
            specie_spd_dos.update({specie: {}})

            for site in sites:
                sitedos = dos_per_site[site]
                sitedos.pop("efermi")
                sitedos.pop("energies")

                for band in sitedos.keys():
                    if band not in specie_spd_dos[specie].keys():
                        specie_spd_dos[specie].update({band: {}})
                        new_dens = sitedos[band]["densities"]

                        if summed_spins:
                            specie_spd_dos[specie][band].update(
                                {
                                    "efermi": efermi,
                                    "energies": energies,
                                    "densities": {
                                        Spin.up: new_dens[Spin.up] + new_dens[Spin.down]
                                    },
                                }
                            )
                        else:
                            specie_spd_dos[specie][band].update(
                                {
                                    "efermi": efermi,
                                    "energies": energies,
                                    "densities": new_dens,
                                }
                            )

                    else:
                        old_dens = specie_spd_dos[specie][band]["densities"]
                        new_dens = sitedos[band]["densities"]
                        if summed_spins:
                            density = (
                                old_dens[Spin.up]
                                + new_dens[Spin.up]
                                + new_dens[Spin.down]
                            )
                            specie_spd_dos[specie][band].update(
                                {"densities": {Spin.up: density}}
                            )
                        else:
                            density = {
                                Spin.up: old_dens[Spin.up] + new_dens[Spin.up],
                                Spin.down: old_dens[Spin.down] + new_dens[Spin.down],
                            }
                            specie_spd_dos[specie][band].update({"densities": density})

        return specie_spd_dos

    def get_tdos_per_species(
        self, dos_type: Literal["vasp", "lso", "lobster"], summed_spins: bool = False
    ) -> dict:
        """Retrieve the 'total' pDOS per species (A/B1/B2/X). Retrieved by summing orbital resolved pDOS per species.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): The DOSCAR of which basis should be analyzed.
            summed_spins (bool, optional): Whether spin-up and down states should be summed. Defaults to False.

        Returns:
            dict[(A/B1/B2/X),Dos]: Dictionary containing the total pDOS per symmetrically inequivalent species.
        """
        specie_spd_dos = self.get_spddos_per_species(dos_type, summed_spins=False)
        energies = specie_spd_dos["energies"]
        efermi = specie_spd_dos["efermi"]

        specie_tdos = {"efermi": efermi, "energies": energies}

        specie_spd_dos.pop("efermi")
        specie_spd_dos.pop("energies")
        for specie, dos_per_band in specie_spd_dos.items():
            specie_tdos.update({specie: {}})
            sum_up = 0
            sum_down = 0
            for band, dos in dos_per_band.items():
                densitie = dos["densities"]

                sum_up += densitie[Spin.up]
                sum_down += densitie[Spin.down]

            if summed_spins:
                new_dens = {Spin.up: sum_up + sum_down}
            else:
                new_dens = {Spin.up: sum_up, Spin.down: sum_down}

            specie_tdos[specie].update(
                {"efermi": efermi, "energies": energies, "densities": new_dens}
            )

        return specie_tdos

    def get_CompleteDosObject(
        self, dos_type: Literal["vasp", "lso", "lobster"]
    ) -> CompleteDos:
        """Retrieve pymatgen CompleteDos object based on specified DOS type.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): The DOSCAR of which basis should be analyzed.

        Raises:
            ValueError: Given DOSCAR type was not understood.

        Returns:
            CompleteDos: Pymatgen CompleteDos object corresponding to this composition and DOSCAR type.
        """
        if dos_type == "vasp":
            from pymatgen.io.vasp.outputs import Vasprun

            ##Pymatgen doesnt have a built-in method for parsing VASP Doscars...
            vrunfile = Path(self.jobpaths["5LOB"]) / self.used_basis / "vasprun.xml"
            compdos = Vasprun(
                filename=vrunfile,
                parse_eigen=False,
                parse_projected_eigen=False,
                parse_potcar_file=False,
                separate_spins=True,
            ).complete_dos
            return compdos
        elif dos_type == "lso":
            dosfile = "DOSCAR.LSO.lobster"
        elif dos_type == "lobster":
            dosfile = "DOSCAR.lobster"
        else:
            raise ValueError(
                f"Please give valid doscar type [vasp/lso/lobster], got {dos_type} instead"
            )

        dospath = Path(self.jobpaths["5LOB"]) / self.used_basis / dosfile
        strucpath = (
            Path(self.jobpaths["5LOB"]) / self.used_basis / "POSCAR.lobster.vasp"
        )

        compdos = Doscar(doscar=dospath, structure_file=strucpath).completedos

        return compdos

    def parse_dos_withsmearing(
        self, dos_type: Literal["vasp", "lso", "lobster"], smearing_width: float = 0.1
    ) -> dict:
        """Parse the DOS with additional Gaussian smearing. Returned DOS contains total DOS and 'total' pDOS per site.

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): The DOSCAR in which basis should be analyzed.
            smearing_width (float, optional): How much smearing should be applied (synonomous to VASP's SIGMA). Defaults to 0.1.

        Returns:
            dict: Dictionary containing total DOS and pDOS per site.
        """

        cdos = self.get_CompleteDosObject(dos_type)
        dos_dict = {}
        dos_dict.update(
            {"tdos": {"densities": cdos.get_smeared_densities(smearing_width)}}
        )
        dos_dict["tdos"].update({"energies": cdos.energies})
        dos_dict.update({"tdos_per_site": {}})
        dos_dict["tdos_per_site"].update({"energies": cdos.energies})

        for ii, site in enumerate(cdos.structure.sites):
            dos_dict["tdos_per_site"].update(
                {
                    ii: {
                        "densities": cdos.get_site_dos(site).get_smeared_densities(
                            smearing_width
                        )
                    }
                }
            )

        return dos_dict

    ###LOBSTER projection related methods

    def _get_qualanalysis_object(self, basis_dir: PathLike | str) -> dict:
        """Returns the Lobsterpy Analysis quality summary for specified lobster basis.

        Args:
            basis_dir (PathLike | str): Which subdirectory (corresponding to a lobster basis set) should be analyzed.

        Raises:
            Exception: If something went wrong in Analysis object construction

        Returns:
            dict: Dictionary containing charge spillings, and bandOverlaps (if existent)
        """
        if type(basis_dir) == str:
            if basis_dir.startswith("basis"):
                basisdir = Path(self.jobpaths["5LOB"] + "/" + basis_dir)
            else:
                basisdir = Path(basis_dir)

        if os.path.exists(basisdir / "bandOverlaps.lobster"):
            fbandoverlap = basisdir / "bandOverlaps.lobster"
        else:
            fbandoverlap = None
        try:
            quality_obj = Analysis.get_lobster_calc_quality_summary(
                path_to_poscar=basisdir / "POSCAR.lobster.vasp",
                path_to_charge=basisdir / "CHARGE.lobster",
                path_to_potcar=basisdir / "POTCAR",
                path_to_bandoverlaps=fbandoverlap,
                path_to_doscar=basisdir / "DOSCAR.LSO.lobster",
                path_to_lobsterin=basisdir / "lobsterin",
                path_to_lobsterout=basisdir / "lobsterout",
                path_to_vasprun=Path(self.jobpaths["4HSE"]) / "vasprun.xml",
                e_range=[-10, 0],
                dos_comparison=True,
                n_bins=1000,
                bva_comp=True,
            )
        except Exception as exc:
            raise Exception(
                f"{self.comp.CompID} quality analysis object raised and unexpected error {exc}"
            )
        return quality_obj

    def _get_COXXanalysis_object(
        self,
        basis_dir: PathLike | str,
        coxx_type: Literal["COOP", "COHP", "COBI"],
        charge_type: Literal["Mulliken", "Loewdin"],
    ):
        """Retrieve Lobsterpy Analysis object for this composition and specified basisset.

        Args:
            basis_dir (PathLike | str): Which subdirectory (corresponding to a lobster basis set) should be analyzed.
            coxx_type (Literal[&quot;COOP&quot;, &quot;COHP&quot;, &quot;COBI&quot;]): Whether COOP, COHP, or COBI should be analyzed.
            charge_type (Literal[&quot;Mulliken&quot;, &quot;Loewdin&quot;]): Whether Mulliken or Loewdin charges should be used.

        Returns:
            Analysis: Lobsterpy Analysis object.
        """

        if type(basis_dir) == str:
            if basis_dir.startswith("basisset"):
                basis_dir = Path(self.jobpaths["5LOB"] + "/" + basis_dir)
            else:
                basis_dir = Path(basis_dir)

        analyse_coxx = Analysis(
            path_to_poscar=basis_dir / "POSCAR",
            path_to_icohplist=basis_dir / f"I{coxx_type}LIST.lobster",
            path_to_cohpcar=basis_dir / f"{coxx_type}CAR.lobster",
            path_to_charge=basis_dir / "CHARGE.lobster",
            orbital_resolved=True,
            are_cobis=coxx_type == "COBI",
            are_coops=coxx_type == "COOP",
            summed_spins=False,
            noise_cutoff=0.000001 if coxx_type == "COBI" else 0.00001,
            orbital_cutoff=0.01,
            cutoff_icohp=0.00001,
            which_bonds="all",
            type_charge=charge_type,
        )

        return analyse_coxx

    def _get_quality_all_basis(self) -> dict:
        """Collect projection quality data for all LOBSTER basissets for which calculations are performed.
        Detecting performed calculations is based on presence of subdirectories with fixed naming scheme (basisset*)

        Returns:
            dict[basisset#,quality_overview]: Dictionary containing charge spillings, and band overlaps for each basisset.
        """
        basesdirs = glob.glob(self.jobpaths["5LOB"] + "/basisset*")
        basesdirs.sort()

        quality_all = {}
        for ii, subdir in enumerate(basesdirs):
            # print(ii)
            quality_analysis = self._get_qualanalysis_object(subdir)

            if quality_analysis["band_overlaps_analysis"]["file_exists"]:
                bo_maxdev = quality_analysis["band_overlaps_analysis"]["max_deviation"]
                perc_dev = quality_analysis["band_overlaps_analysis"][
                    "percent_kpoints_abv_limit"
                ]
            else:
                bo_maxdev = "NaN"
                perc_dev = "NaN"

            quality_all.update(
                {
                    subdir.split("/")[-1]: {
                        "abs_charge_spilling": quality_analysis["charge_spilling"][
                            "abs_charge_spilling"
                        ],
                        "bandOverlap_exists": quality_analysis[
                            "band_overlaps_analysis"
                        ]["file_exists"],
                        "bandOverlap_maxDev": bo_maxdev,
                        "bandOverlap_perc_kpts": perc_dev,
                    }
                }
            )

        # self.lobquality_all = quality_all

        return quality_all

    def analyze_lobster_outputs(
        self,
        charge_type: Literal["Mulliken", "Loewdin"],
        basis_dir: PathLike | str | None = None,
    ):
        """Analyze the LOBSTER outputs; namely the ICOHPLIST, ICOBILIST, CHARGES, SitePotentials, and GROSSPOP.
        Mostly summarizes the results on a species (A/B1/B2/X) level. Also calculates axial and directional asymmetry indices.
        Also reports the projection quality for the basisset.

        Args:
            charge_type (Literal[&quot;Mulliken&quot;, &quot;Loewdin&quot;]): Whether Mulliken or Loewdin charges should be analyzed.
            basis_dir (PathLike | str | None, optional): Which basisset (and also subdirectory) the analysis should be based on
                                                        If None is given it will resort to self.used_basis. Defaults to None.

        Returns:
            dict: Data extracted from the various LOBSTER outputs.
        """

        if basis_dir == None:
            basis_dir = self.used_basis

        lobdata_dict = self._get_quality_all_basis()[self.used_basis]

        lobsterdir = Path(self.jobpaths["5LOB"]) / basis_dir

        gpop_list = Grosspop(lobsterdir / "GROSSPOP.lobster").list_dict_grosspop
        netpopdiff_list = [
            gpop_list[x][charge_type + " GP"]["total"][Spin.up]
            - gpop_list[x][charge_type + " GP"]["total"][Spin.down]
            for x in range(len(gpop_list))
        ]

        lobdata_dict.update(
            {
                "popdiff.total": round(np.sum(netpopdiff_list), 2),
                "popdiff.B1": round(netpopdiff_list[0], 2),
            }
        )

        charge_list = Charge(lobsterdir / "CHARGE.lobster")

        chargedict = getattr(charge_list, charge_type)

        if self.comp.B2 == "Vac":
            lobdata_dict.update(
                {
                    "charge.B1": chargedict[0],
                    "charge.B2": None,
                    "charge.A": (chargedict[1] + chargedict[2]) / 2,
                    "charge.X": (np.sum(chargedict[3:])) / 6,
                    "charge.total": round(np.sum(chargedict), 2),
                    "charge.balanced": np.abs(np.sum(chargedict)) < 0.2,
                }
            )
        else:

            lobdata_dict.update(
                {
                    "popdiff.B2": round(netpopdiff_list[1], 2),
                    "charge.B1": chargedict[0],
                    "charge.B2": chargedict[1],
                    "charge.A": (chargedict[2] + chargedict[3]) / 2,
                    "charge.X": np.sum(chargedict[4:]) / 6,
                    "charge.total": round(np.sum(chargedict), 2),
                    "charge.balanced": np.abs(np.sum(chargedict)) < 0.2,
                }
            )

        icohps = Icohplist(
            filename=lobsterdir / "ICOHPLIST.lobster", is_spin_polarized=True
        ).icohpcollection
        icobis = Icohplist(
            filename=lobsterdir / "ICOBILIST.lobster",
            are_cobis=True,
            is_spin_polarized=True,
        ).icohpcollection

        b1_bondlist = list(icohps.get_icohp_dict_of_site(0))
        # print(b1_bondlist)
        assert (
            len(b1_bondlist) == 6
        ), f"{self.comp.CompID} has incorrect number of bonds {len(b1_bondlist)}"

        lobdata_dict.update(
            {
                "Icohp.B1.sum": icohps.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"], divisor=1, summed_spin_channels=True
                ),
                "Icohp.B1.avg": icohps.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"], divisor=6, summed_spin_channels=True
                ),
                "Icohp.B1.avg_spinup": icohps.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"],
                    divisor=6,
                    summed_spin_channels=False,
                    spin=Spin.up,
                ),
                "Icohp.B1.avg_spindown": icohps.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"],
                    divisor=6,
                    summed_spin_channels=False,
                    spin=Spin.down,
                ),
                "Icohp.B1.x_avg": icohps.get_summed_icohp_by_label_list(
                    ["1", "2"], divisor=2, summed_spin_channels=True
                ),
                "Icohp.B1.x_spinup": icohps.get_summed_icohp_by_label_list(
                    ["1", "2"], divisor=2, summed_spin_channels=False, spin=Spin.up
                ),
                "Icohp.B1.x_spindown": icohps.get_summed_icohp_by_label_list(
                    ["1", "2"], divisor=2, summed_spin_channels=False, spin=Spin.down
                ),
                "Icohp.B1.y_avg": icohps.get_summed_icohp_by_label_list(
                    ["3", "4"], divisor=2, summed_spin_channels=True
                ),
                "Icohp.B1.y_spinup": icohps.get_summed_icohp_by_label_list(
                    ["3", "4"], divisor=2, summed_spin_channels=False, spin=Spin.up
                ),
                "Icohp.B1.y_spindown": icohps.get_summed_icohp_by_label_list(
                    ["3", "4"], divisor=2, summed_spin_channels=False, spin=Spin.down
                ),
                "Icohp.B1.z_avg": icohps.get_summed_icohp_by_label_list(
                    ["5", "6"], divisor=2, summed_spin_channels=True
                ),
                "Icohp.B1.z_spinup": icohps.get_summed_icohp_by_label_list(
                    ["5", "6"], divisor=2, summed_spin_channels=False, spin=Spin.up
                ),
                "Icohp.B1.z_spindown": icohps.get_summed_icohp_by_label_list(
                    ["5", "6"], divisor=2, summed_spin_channels=False, spin=Spin.down
                ),
                "Icobi.B1.sum": icobis.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"], divisor=1, summed_spin_channels=True
                ),
                "Icobi.B1.avg": icobis.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"], divisor=6, summed_spin_channels=True
                ),
                "Icobi.B1.avg_spinup": icobis.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"],
                    divisor=6,
                    summed_spin_channels=False,
                    spin=Spin.up,
                ),
                "Icobi.B1.avg_spindown": icobis.get_summed_icohp_by_label_list(
                    ["1", "2", "3", "4", "5", "6"],
                    divisor=6,
                    summed_spin_channels=False,
                    spin=Spin.down,
                ),
                "Icobi.B1.x_avg": icobis.get_summed_icohp_by_label_list(
                    ["1", "2"], divisor=2, summed_spin_channels=True
                ),
                "Icobi.B1.x_spinup": icobis.get_summed_icohp_by_label_list(
                    ["1", "2"], divisor=2, summed_spin_channels=False, spin=Spin.up
                ),
                "Icobi.B1.x_spindown": icobis.get_summed_icohp_by_label_list(
                    ["1", "2"], divisor=2, summed_spin_channels=False, spin=Spin.down
                ),
                "Icobi.B1.y_avg": icobis.get_summed_icohp_by_label_list(
                    ["3", "4"], divisor=2, summed_spin_channels=True
                ),
                "Icobi.B1.y_spinup": icobis.get_summed_icohp_by_label_list(
                    ["3", "4"], divisor=2, summed_spin_channels=False, spin=Spin.up
                ),
                "Icobi.B1.y_spindown": icobis.get_summed_icohp_by_label_list(
                    ["3", "4"], divisor=2, summed_spin_channels=False, spin=Spin.down
                ),
                "Icobi.B1.z_avg": icobis.get_summed_icohp_by_label_list(
                    ["5", "6"], divisor=2, summed_spin_channels=True
                ),
                "Icobi.B1.z_spinup": icobis.get_summed_icohp_by_label_list(
                    ["5", "6"], divisor=2, summed_spin_channels=False, spin=Spin.up
                ),
                "Icobi.B1.z_spindown": icobis.get_summed_icohp_by_label_list(
                    ["5", "6"], divisor=2, summed_spin_channels=False, spin=Spin.down
                ),
            }
        )

        icohpsfeature = FeaturizeIcoxxlist(
            path_to_icoxxlist=lobsterdir / "ICOHPLIST.lobster",
            path_to_structure=lobsterdir / "POSCAR.lobster.vasp",
        )
        icobisfeature = FeaturizeIcoxxlist(
            path_to_icoxxlist=lobsterdir / "ICOBILIST.lobster",
            are_cobis=True,
            path_to_structure=lobsterdir / "POSCAR.lobster.vasp",
        )

        lobdata_dict.update(
            {
                "Icohp.B1.axial_asym_index": icohpsfeature.calc_site_asymmetry_index(0),
                "Icohp.B1.directional_asym_index": (
                    np.max(
                        np.abs(
                            [
                                lobdata_dict["Icohp.B1.x_avg"],
                                lobdata_dict["Icohp.B1.y_avg"],
                                lobdata_dict["Icohp.B1.z_avg"],
                            ]
                        )
                    )
                    - np.min(
                        np.abs(
                            [
                                lobdata_dict["Icohp.B1.x_avg"],
                                lobdata_dict["Icohp.B1.y_avg"],
                                lobdata_dict["Icohp.B1.z_avg"],
                            ]
                        )
                    )
                )
                / np.abs(lobdata_dict["Icohp.B1.avg"]),
                "Icobi.B1.axial_asym_index": icobisfeature.calc_site_asymmetry_index(0),
                "Icobi.B1.directional_asym_index": (
                    np.max(
                        np.abs(
                            [
                                lobdata_dict["Icobi.B1.x_avg"],
                                lobdata_dict["Icobi.B1.y_avg"],
                                lobdata_dict["Icobi.B1.z_avg"],
                            ]
                        )
                    )
                    - np.min(
                        np.abs(
                            [
                                lobdata_dict["Icobi.B1.x_avg"],
                                lobdata_dict["Icobi.B1.y_avg"],
                                lobdata_dict["Icobi.B1.z_avg"],
                            ]
                        )
                    )
                )
                / np.abs(lobdata_dict["Icobi.B1.avg"]),
            }
        )

        if self.comp.B2 != "Vac":
            b2_bondlist = list(icohps.get_icohp_dict_of_site(1))
            assert (
                len(b2_bondlist) == 6
            ), f"{self.comp.CompID} has incorrect number of bonds {len(b2_bondlist)}"
            lobdata_dict.update(
                {
                    "Icohp.B2.sum": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist, divisor=1, summed_spin_channels=True
                    ),
                    "Icohp.B2.avg": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist, divisor=6, summed_spin_channels=True
                    ),
                    "Icohp.B2.avg_spinup": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist, divisor=6, summed_spin_channels=False, spin=Spin.up
                    ),
                    "Icohp.B2.avg_spindown": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist,
                        divisor=6,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icohp.B2.x_avg": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[:2], divisor=2, summed_spin_channels=True
                    ),
                    "Icohp.B2.x_spinup": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[:2],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.up,
                    ),
                    "Icohp.B2.x_spindown": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[:2],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icohp.B2.y_avg": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[2:4], divisor=2, summed_spin_channels=True
                    ),
                    "Icohp.B2.y_spinup": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[2:4],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.up,
                    ),
                    "Icohp.B2.y_spindown": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[2:4],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icohp.B2.z_avg": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[4:], divisor=2, summed_spin_channels=True
                    ),
                    "Icohp.B2.z_spinup": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[4:],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.up,
                    ),
                    "Icohp.B2.z_spindown": icohps.get_summed_icohp_by_label_list(
                        b2_bondlist[4:],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icobi.B2.sum": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist, divisor=1, summed_spin_channels=True
                    ),
                    "Icobi.B2.avg": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist, divisor=6, summed_spin_channels=True
                    ),
                    "Icobi.B2.avg_spinup": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist, divisor=6, summed_spin_channels=False, spin=Spin.up
                    ),
                    "Icobi.B2.avg_spindown": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist,
                        divisor=6,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icobi.B2.x_avg": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[:2], divisor=2, summed_spin_channels=True
                    ),
                    "Icobi.B2.x_spinup": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[:2],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.up,
                    ),
                    "Icobi.B2.x_spindown": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[:2],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icobi.B2.y_avg": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[2:4], divisor=2, summed_spin_channels=True
                    ),
                    "Icobi.B2.y_spinup": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[2:4],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.up,
                    ),
                    "Icobi.B2.y_spindown": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[2:4],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                    "Icobi.B2.z_avg": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[4:], divisor=2, summed_spin_channels=True
                    ),
                    "Icobi.B2.z_spinup": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[4:],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.up,
                    ),
                    "Icobi.B2.z_spindown": icobis.get_summed_icohp_by_label_list(
                        b2_bondlist[4:],
                        divisor=2,
                        summed_spin_channels=False,
                        spin=Spin.down,
                    ),
                }
            )

            lobdata_dict.update(
                {
                    "Icohp.B2.axial_asym_index": icohpsfeature.calc_site_asymmetry_index(
                        1
                    ),
                    "Icohp.B2.directional_asym_index": (
                        np.max(
                            np.abs(
                                [
                                    lobdata_dict["Icohp.B2.x_avg"],
                                    lobdata_dict["Icohp.B2.y_avg"],
                                    lobdata_dict["Icohp.B2.z_avg"],
                                ]
                            )
                        )
                        - np.min(
                            np.abs(
                                [
                                    lobdata_dict["Icohp.B2.x_avg"],
                                    lobdata_dict["Icohp.B2.y_avg"],
                                    lobdata_dict["Icohp.B2.z_avg"],
                                ]
                            )
                        )
                    )
                    / np.abs(lobdata_dict["Icohp.B2.avg"]),
                    "Icobi.B2.axial_asym_index": icobisfeature.calc_site_asymmetry_index(
                        1
                    ),
                    "Icobi.B2.directional_asym_index": (
                        np.max(
                            np.abs(
                                [
                                    lobdata_dict["Icobi.B2.x_avg"],
                                    lobdata_dict["Icobi.B2.y_avg"],
                                    lobdata_dict["Icobi.B2.z_avg"],
                                ]
                            )
                        )
                        - np.min(
                            np.abs(
                                [
                                    lobdata_dict["Icobi.B2.x_avg"],
                                    lobdata_dict["Icobi.B2.y_avg"],
                                    lobdata_dict["Icobi.B2.z_avg"],
                                ]
                            )
                        )
                    )
                    / np.abs(lobdata_dict["Icobi.B2.avg"]),
                }
            )

        else:

            lobdata_dict.update(
                {
                    "Icohp.B2.sum": None,
                    "Icohp.B2.avg": None,
                    "Icohp.B2.avg_spinup": None,
                    "Icohp.B2.avg_spindown": None,
                    "Icohp.B2.x_avg": None,
                    "Icohp.B2.x_spinup": None,
                    "Icohp.B2.x_spindown": None,
                    "Icohp.B2.y_avg": None,
                    "Icohp.B2.y_spinup": None,
                    "Icohp.B2.y_spindown": None,
                    "Icohp.B2.z_avg": None,
                    "Icohp.B2.z_spinup": None,
                    "Icohp.B2.z_spindown": None,
                    "Icobi.B2.sum": None,
                    "Icobi.B2.avg": None,
                    "Icobi.B2.avg_spinup": None,
                    "Icobi.B2.avg_spindown": None,
                    "Icobi.B2.x_avg": None,
                    "Icobi.B2.x_spinup": None,
                    "Icobi.B2.x_spindown": None,
                    "Icobi.B2.y_avg": None,
                    "Icobi.B2.y_spinup": None,
                    "Icobi.B2.y_spindown": None,
                    "Icobi.B2.z_avg": None,
                    "Icobi.B2.z_spinup": None,
                    "Icobi.B2.z_spindown": None,
                }
            )

            lobdata_dict.update(
                {
                    "Icohp.B2.axial_asym_index": None,
                    "Icohp.B2.directional_asym_index": None,
                    "Icobi.B2.axial_asym_index": None,
                    "Icobi.B2.directional_asym_index": None,
                }
            )

        sitepots = SitePotential(filename=lobsterdir / "SitePotentials.lobster")

        lobdata_dict.update(
            {
                "MadelungEnergy": getattr(
                    sitepots, f"madelungenergies_{charge_type.lower()}"
                ),
                "SitePotential.B1": getattr(
                    sitepots, f"sitepotentials_{charge_type.lower()}"
                )[0],
            }
        )

        if self.comp.B2 == "Vac":
            lobdata_dict.update({"SitePotential.B2": None})

        else:
            lobdata_dict.update(
                {
                    "SitePotential.B2": getattr(
                        sitepots, f"sitepotentials_{charge_type.lower()}"
                    )[1]
                }
            )

        return lobdata_dict

    def get_parsed_coxx(
        self,
        basis_dir: PathLike | str | None = None,
        coxx_type: Literal["COHP", "COBI"] = "COHP",
    ) -> dict:
        """Parse COHPCAR or COBICAR to retrieve the energy resolved COHP/COBI and ICOHP/ICOBI for the B-X bonds.
        The (I)COXX(E) is given for the average of all 6 bonds, and averaged over the 2 bonds aligned with x/y/z axes.

        Args:
            basis_dir (PathLike | str | None, optional): From which basisset (and corresponding subdirectory) the COXXCAR should be parsed.
                                                        If None is given reverts to self.used_basis. Defaults to None.
            coxx_type (Literal[&quot;COHP&quot;, &quot;COBI&quot;], optional): Whether COHP or COBI needs to be parsed. Defaults to "COHP".

        Returns:
            dict: Dictionary of COxx and ICOxx of (B1/B2)-X bonds reported for (average/x/y/z)
        """
        if not basis_dir:
            basis_dir = self.used_basis

        completecoxx = CompleteCohp.from_file(
            fmt="LOBSTER",
            filename=Path(self.jobpaths["5LOB"])
            / basis_dir
            / f"{coxx_type}CAR.lobster",
            structure_file=Path(self.jobpaths["5LOB"])
            / basis_dir
            / "POSCAR.lobster.vasp",
            are_cobis=(coxx_type == "COBI"),
        )

        parsed_coxx = {}

        parsed_coxx.update(
            {
                "efermi": completecoxx.efermi,
                "specie.B1": self.input_dict["input_species"]["B1"],
                "specie.B2": self.input_dict["input_species"]["B2"],
                "specie.A": self.input_dict["input_species"]["A"],
                "specie.X": self.input_dict["input_species"]["X"],
                "energies": completecoxx.energies,
            }
        )

        coxx_b1avg = completecoxx.get_summed_cohp_by_label_list(
            ["1", "2", "3", "4", "5", "6"], divisor=6, summed_spin_channels=False
        ).get_cohp()
        Icoxx_b1avg = completecoxx.get_summed_cohp_by_label_list(
            ["1", "2", "3", "4", "5", "6"], divisor=6, summed_spin_channels=False
        ).get_icohp()
        coxx_b1x = completecoxx.get_summed_cohp_by_label_list(
            ["1", "2"], divisor=2, summed_spin_channels=False
        ).get_cohp()
        Icoxx_b1x = completecoxx.get_summed_cohp_by_label_list(
            ["1", "2"], divisor=2, summed_spin_channels=False
        ).get_icohp()
        coxx_b1y = completecoxx.get_summed_cohp_by_label_list(
            ["3", "4"], divisor=2, summed_spin_channels=False
        ).get_cohp()
        Icoxx_b1y = completecoxx.get_summed_cohp_by_label_list(
            ["3", "4"], divisor=2, summed_spin_channels=False
        ).get_icohp()
        coxx_b1z = completecoxx.get_summed_cohp_by_label_list(
            ["5", "6"], divisor=2, summed_spin_channels=False
        ).get_cohp()
        Icoxx_b1z = completecoxx.get_summed_cohp_by_label_list(
            ["5", "6"], divisor=2, summed_spin_channels=False
        ).get_icohp()

        if self.comp.B2 != "Vac":
            coxx_b2avg = completecoxx.get_summed_cohp_by_label_list(
                ["7", "8", "9", "10", "11", "12"], divisor=6, summed_spin_channels=False
            ).get_cohp()
            Icoxx_b2avg = completecoxx.get_summed_cohp_by_label_list(
                ["7", "8", "9", "10", "11", "12"], divisor=6, summed_spin_channels=False
            ).get_icohp()
            coxx_b2x = completecoxx.get_summed_cohp_by_label_list(
                ["7", "8"], divisor=2, summed_spin_channels=False
            ).get_cohp()
            Icoxx_b2x = completecoxx.get_summed_cohp_by_label_list(
                ["7", "8"], divisor=2, summed_spin_channels=False
            ).get_icohp()
            coxx_b2y = completecoxx.get_summed_cohp_by_label_list(
                ["9", "10"], divisor=2, summed_spin_channels=False
            ).get_cohp()
            Icoxx_b2y = completecoxx.get_summed_cohp_by_label_list(
                ["9", "10"], divisor=2, summed_spin_channels=False
            ).get_icohp()
            coxx_b2z = completecoxx.get_summed_cohp_by_label_list(
                ["11", "12"], divisor=2, summed_spin_channels=False
            ).get_cohp()
            Icoxx_b2z = completecoxx.get_summed_cohp_by_label_list(
                ["11", "12"], divisor=2, summed_spin_channels=False
            ).get_icohp()

            parsed_coxx.update(
                {
                    "coxx_b1avg": coxx_b1avg,
                    "coxx_b2avg": coxx_b2avg,
                    "coxx_b1x": coxx_b1x,
                    "coxx_b1y": coxx_b1y,
                    "coxx_b1z": coxx_b1z,
                    "coxx_b2x": coxx_b2x,
                    "coxx_b2y": coxx_b2y,
                    "coxx_b2z": coxx_b2z,
                    "Icoxx_b1avg": Icoxx_b1avg,
                    "Icoxx_b2avg": Icoxx_b2avg,
                    "Icoxx_b1x": Icoxx_b1x,
                    "Icoxx_b1y": Icoxx_b1y,
                    "Icoxx_b1z": Icoxx_b1z,
                    "Icoxx_b2x": Icoxx_b2x,
                    "Icoxx_b2y": Icoxx_b2y,
                    "Icoxx_b2z": Icoxx_b2z,
                }
            )
        else:
            parsed_coxx.update(
                {
                    "coxx_b1avg": coxx_b1avg,
                    "coxx_b2avg": None,
                    "coxx_b1x": coxx_b1x,
                    "coxx_b1y": coxx_b1y,
                    "coxx_b1z": coxx_b1z,
                    "coxx_b2x": None,
                    "coxx_b2y": None,
                    "coxx_b2z": None,
                    "Icoxx_b1avg": Icoxx_b1avg,
                    "Icoxx_b2avg": None,
                    "Icoxx_b1x": Icoxx_b1x,
                    "Icoxx_b1y": Icoxx_b1y,
                    "Icoxx_b1z": Icoxx_b1z,
                    "Icoxx_b2x": None,
                    "Icoxx_b2y": None,
                    "Icoxx_b2z": None,
                }
            )

        return parsed_coxx


class GroupedAnalysis:
    """This class takes care of collecting the data for the entire HDP database.
    Initializing this class will create SingleHDPanalysis objects for each composition and call its methods and
    collect its data into DataFrames for each.
    Methods that are named 'process_comp_....' are used to collect data from a single composition, which gets used in a
    multiprocess pool by methods called 'get_...' or 'save_...'
    Attributes:
        complist (list[Compound]): list of Compound objects that have completed the required simulation step
        analysislist (list[SingleHDPanalysis]): list of SingleHDPanalysis objects
        basis_selected (bool): indicator if the lobster basisset selection has happened.
    """

    def __init__(
        self, ledger: ProcessLedger, required_step: str = "5LOB", testing: int = 0
    ):
        """Initializes the group of Compounds to be analyzed depending on which step is required to be finished.
        Some functionality may not work if you deviate from requiring the last step ("5LOB") to be finished,
        but you could for example analyze structures based on all comps that completed "1Rel".

        Based on the completion overview from the ProcessLedger, it will initialize SingleHDPanalysis objects
        for each composition that completed the required steps and save it in an attribute.

        Args:
            ledger (ProcessLedger): The ProcessLedger used to track the process for all compositions.
            required_step (str, optional): Which step of the workflow needs to be completed to perform
                                            the intended analysis. Defaults to "5LOB".
            testing (int, optional): If larger than 0 only this number of compositions will be analyzed.
                                    can be used to test on a small subset first. Defaults to 0.
        """
        completionDF = ledger.GetCompletionOverview()
        assert (
            required_step in completionDF.columns.to_list()
        ), f"The given required step for analysis ({required_step}) was not found in LedgerFile"

        if required_step == "5LOB":
            completed_comps = completionDF[
                (
                    (completionDF[required_step] == "-1")
                    | (completionDF[required_step] == "1")
                )
            ].index.to_list()
        else:
            completed_comps = completionDF[
                completionDF[required_step] == "1"
            ].index.to_list()

        self.basis_selected = False
        full_complist = ledger.RestartLedger()
        if testing:

            self.complist = [
                x for x in full_complist if x.CompID in completed_comps[:testing]
            ]
        else:
            self.complist = [x for x in full_complist if x.CompID in completed_comps]
        print(f"{time.strftime('%H:%M:%S')}: starting analysis object construction")
        self.analysis_list = [
            SingleHDPanalysis(x, ledger.LedgerDF[x.CompID]) for x in self.complist
        ]
        return

    def select_lobbasis(
        self,
        saved_qualityoverview: Path | str | None = None,
        save_badprojections: bool = False,
        saving_dir: Path = Path("."),
    ):
        """Collects or reads projection quality for all calculated LOBSTER projection bases.
        If any minimal basisset projection has charge spilling over 3%, it will chose an
        alternative basis (if available) with the lowest charge spilling.

        Args:
            saved_qualityoverview (Path | str | None, optional): If the projection quality
                                    overview was previously generated, it can be read in stead of
                                    re-analyzing all bases. Defaults to None.
            save_badprojections (bool, optional): Saves a separate overview with all bases that
                                        result in charge spilling over 3%. Defaults to False.
            saving_dir (Path, optional): Path to save qualityoverview (and bad projections).
                                        Defaults to Path(".").
        """

        if saved_qualityoverview:
            qual_df = pd.read_csv(saved_qualityoverview, index_col=[0, 1])
        else:
            qual_df = self.get_lobquality_data()
            qual_df.to_csv(
                saving_dir / f'HDP_lobqual_overview_{time.strftime("%y%m%d")}.csv'
            )

        ##Find which comps have charge spilling > 3% in basis0
        b0qual_df = qual_df.xs("basisset0", level=1)
        badqualcomps = b0qual_df[b0qual_df["abs_charge_spilling"] > 3].index
        # save which comps have bad projection quality
        newchoice_df = qual_df.loc[badqualcomps]
        if save_badprojections:
            newchoice_df.to_csv(
                saving_dir / f'HDPs_SpillOver3_{time.strftime("%y%m%d")}.csv'
            )

        for comp in newchoice_df.index.get_level_values(0):
            # print(comp)
            pos_choices = newchoice_df.loc[comp]
            corr_basis = pos_choices[
                pos_choices["abs_charge_spilling"]
                == pos_choices["abs_charge_spilling"].min()
            ].index.to_list()[0]
            for CompAnalysis in self.analysis_list:
                if CompAnalysis.comp.CompID == str(comp):

                    CompAnalysis.update_lobbasis(corr_basis)

        self.basis_selected = True

        return

    def plot_doscar_single(self, comp: SingleHDPanalysis) -> None:
        from pdosplotter_adapted import plot_tdos

        """Plots the total DOS and pDOS for a single composition
        accesses attributes self.dostype and self.output_dir 
        which need to be set by the calling function (plot_doscar_all)
        """

        dos_type = self.dostype
        output_dir = self.output_dir

        dos_dict = comp.get_tdos_per_species(dos_type)
        dos_dict.update({"tdos": comp.get_tdos(dos_type)["densities"]})
        Efermi = comp.analyze_bandedges(dos_type)["combined"]["VBM"]
        filename = str(output_dir / f"{comp.comp.CompID}_{dos_type}dos.png")
        try:
            plot_tdos(dos_dict, comp._get_WFinput_dict(), Efermi, filename)
        except Exception as exc:
            print(f"{comp.comp} generated an exception in plotting:\n{exc}")
            return
        return

    def plot_doscar_all(
        self, dos_type: Literal["vasp", "lso", "lobster"], output_dir: Path = Path(".")
    ):
        """Plots the total DOS and pDOS for all compositions and saves the figures.
        Uses the self.plot_doscar_single() method combined with multiprocess.Pool

        Args:
            dos_type (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;]): which type of DOSCAR to plot from
            output_dir (Path, optional): Path to save the .png's for all plots. Defaults to Path(".").
        """
        self.dostype = dos_type
        self.output_dir = output_dir

        print(
            f"{time.strftime('%H:%M:%S')}: Starting {dos_type}-DOS plotting for {len(self.analysis_list)} Compounds"
        )
        with mp.Pool(processes=8, maxtasksperchild=1) as pool, tqdm(
            total=len(self.analysis_list), desc="Processing compounds"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(self.plot_doscar_single, self.analysis_list)
            ):
                pbar.update()

        return

    def process_comp_quality(self, comp: SingleHDPanalysis):
        """Calls the function to retrieve overview of LOBSTER projection quality for all
        available basissets for a single composition

        Args:
            comp (SingleHDPanalysis): analysis object of specific composition

        Returns:
            dict[compID, dict[projection_quality]]
        """
        return {comp.comp.CompID: comp._get_quality_all_basis()}

    def process_comp_edgecontr(self, comp: SingleHDPanalysis):
        """Calls the function to analyze band edge contributions for a specific composition.
        accesses the attribute self.dostype to determine which DOSCAR type needs to be analyzed
        (is set by calling function)

        Args:
            comp (SingleHDPanalysis): analysis object for specific compound

        Returns:
            dict[compID, dict[bandedge_analysis]]
        """
        return {comp.comp.CompID: comp.analyze_edgecontr(dos_type=self.dostype)}

    def process_comp_completedos(self, comp: SingleHDPanalysis):
        """Retrieves pymatgen CompleteDos object for a specific composition and saves to a json.
        Uses attributes self.dostype and self.output_dir to determine which DOSCAR to parse and
        where to store the json. These attributes are set by the calling function.

        Args:
            comp (SingleHDPanalysis): analysis object of composition to analyze.
        """
        dos_type = self.dostype
        output_dir = self.output_dir
        # print(output_dir)

        compdos = comp.get_CompleteDosObject(dos_type=dos_type).as_dict()

        filename = Path(output_dir) / f"{comp.comp.CompID}_{dos_type}dos.json"
        with open(filename, "w") as f:
            json.dump(jsanitize(compdos), f, indent=4)

        return

    def process_comp_saveparsedcohp(self, comp: SingleHDPanalysis):
        """Calls the function to parse COHPCAR or COBICAR for specific composition and saves to json.
        Uses attributes self.coxx_type and self.output_dir to determine which file to parse
        and where to store the json files. These attributes are set by the calling function.

        Args:
            comp (SingleHDPanalysis): analysis object of specific composition.
        """
        output_dir = self.output_dir
        coxx_type = self.coxx_type
        parescohp = comp.get_parsed_coxx(coxx_type=coxx_type)

        filename = Path(output_dir) / f"{comp.comp.CompID}_parsed{coxx_type}.json"
        with open(filename, "w") as f:
            json.dump(jsanitize(parescohp), f)

        return

    def process_comp_savecohp_smeared(self, comp: SingleHDPanalysis):
        """Calls the function to parse COHPCAR or COBICAR for specific composition,
         Adds additional Gaussian smearing, and saves to json.
        Uses attributes self.coxx_type, self.smearing_width, and self.output_dir to determine which file to parse,
        how much smearing to apply, and where to store the json files. These attributes are set by the calling function.

        Args:
            comp (SingleHDPanalysis): analysis object of specific composition.
        """

        def get_smeared_coxx(coxx, energies, sigma: float) -> dict:
            """Get the the COHPCAR with a Gaussian smearing.

            Args:
                sigma (float): Standard deviation of Gaussian smearing.

            Returns:
                {Spin: NDArray}: Gaussian-smeared DOS by spin.
            """
            diff = [
                energies[idx + 1] - energies[idx] for idx in range(len(energies) - 1)
            ]
            avg_diff = sum(diff) / len(diff)
            return {
                spin: gaussian_filter1d(coxx, sigma / avg_diff)
                for spin, coxx in coxx.items()
            }

        output_dir = self.output_dir
        coxx_type = self.coxx_type
        parescohp = comp.get_parsed_coxx(coxx_type=coxx_type)
        newdict = {}
        newdict.update(
            {
                "efermi": parescohp["efermi"],
                "specie.B1": parescohp["specie.B1"],
                "specie.B2": parescohp["specie.B2"],
                "specie.A": parescohp["specie.A"],
                "specie.X": parescohp["specie.X"],
                "energies": parescohp["energies"],
            }
        )

        for key, val in parescohp.items():
            if key.startswith("coxx") and (val is not None):
                newdict[key] = get_smeared_coxx(
                    val, parescohp["energies"], sigma=self.smearing_width
                )

        filename = Path(output_dir) / f"{comp.comp.CompID}_smeared{coxx_type}.json"
        with open(filename, "w") as f:
            json.dump(jsanitize(newdict), f)

        return

    def process_comp_basic(self, comp: SingleHDPanalysis):
        """Retrieves data on specific composition based on input given to the workflow
        and parameters based on ionic radii.

        Args:
            comp (SingleHDPanalysis): analysis object for specific composition.

        Returns:
            dict[compID, basic_info]: dictionary of input data, tolerance factors, and basis sets.
        """
        inputdic = comp._get_WFinput_dict()
        ion_radius_data = comp.get_ionicradii_based_data()
        basic_info_dict = {
            "compID_num": int(comp.comp.CompID.split("_")[0]),
            "comp_name_simple": str(comp.comp),
            "comp_name_full": f"{comp.comp.A}2{comp.comp.B1}{comp.comp.B2}{comp.comp.X}6",
            "specie.A": inputdic["input_species"]["A"],
            "specie.B1": inputdic["input_species"]["B1"],
            "specie.B2": inputdic["input_species"]["B2"],
            "specie.X": inputdic["input_species"]["X"],
            "element.A": inputdic["input_species"]["A"].element,
            "element.B1": inputdic["input_species"]["B1"].element,
            "element.B2": (
                inputdic["input_species"]["B2"].element
                if comp.comp.B2 != "Vac"
                else None
            ),
            "element.X": inputdic["input_species"]["X"].element,
            "r_ionic.A": ion_radius_data["r_A"],
            "r_ionic.B1": ion_radius_data["r_B1"],
            "r_ionic.B2": ion_radius_data["r_B2"],
            "r_ionic.X": ion_radius_data["r_X"],
            "oct_factor": ion_radius_data["oct_factor"],
            "oct_mismatch": ion_radius_data["oct_mismatch"],
            "gen_t_factor": ion_radius_data["gen_t_factor"],
            "geom_stable": ion_radius_data["geom_stable"],
            "broken_condition": ion_radius_data["broken_condition"],
            "tau_factor": ion_radius_data["tau_factor"],
            "block.B1": inputdic["input_species"]["B1"].element.block,
            "row.B1": inputdic["input_species"]["B1"].row,
            "inputcharge.B1": inputdic["input_species"]["B1"].oxi_state,
            "block.B2": (
                inputdic["input_species"]["B2"].element.block
                if comp.comp.B2 != "Vac"
                else "Vac"
            ),
            "row.B2": (
                inputdic["input_species"]["B2"].row if comp.comp.B2 != "Vac" else None
            ),
            "inputcharge.B2": (
                inputdic["input_species"]["B2"].oxi_state
                if comp.comp.B2 != "Vac"
                else 0
            ),
            "used_lobbasis": inputdic["used_basis_num"],
            "used_lobbasis_func": inputdic["used_basis_funcs"],
        }
        block_list = [basic_info_dict["block.B1"], basic_info_dict["block.B2"]]
        block_list.sort()
        basic_info_dict.update({"block_pairing": f"{block_list[0]}-{block_list[1]}"})

        # print(basic_info_dict)
        return {comp.comp.CompID: basic_info_dict}

    def process_comp_savedos_perspecies(self, comp: SingleHDPanalysis):
        """Parses DOSCAR, resolved on species (A/B1/B2/X) level and saves to json.
        Will save total DOS, total pDOS per species, and orbital resolved pDOS per species.

        Uses attributes self.dostype and self.output_dir to determine which DOSCAR to parse and
        where to store the json. These attributes are set by the calling function.

        Args:
            comp (SingleHDPanalysis): analysis object for specific composition.
        """
        dos_type = self.dostype
        output_dir = self.output_dir
        # print(output_dir)
        output_dos = {}
        output_dos.update({"tdos": comp.get_tdos(dos_type)})
        output_dos.update({"tdos_per_specie": comp.get_tdos_per_species(dos_type)})
        output_dos.update({"spddos_per_specie": comp.get_spddos_per_species(dos_type)})

        filename = Path(output_dir) / f"{comp.comp.CompID}_{dos_type}dos.json"
        with open(filename, "w") as f:
            json.dump(jsanitize(output_dos), f)

    def process_comp_savedos_persite(self, comp: SingleHDPanalysis):
        """Parses DOSCAR to site resolved level and saves to json.
        Will save total DOS, total pDOS per site, and orbital resolved pDOS per site.

        Uses attributes self.dostype and self.output_dir to determine which DOSCAR to parse and
        where to store the json. These attributes are set by the calling function.

        Args:
            comp (SingleHDPanalysis): analysis object for specific composition.
        """
        dos_type = self.dostype
        output_dir = self.output_dir
        # print(output_dir)
        output_dos = {}
        output_dos.update({"tdos": comp.get_tdos(dos_type)})
        output_dos.update({"tdos_per_site": comp.get_tdos_per_site(dos_type)})
        output_dos.update({"spddos_per_site": comp.get_spddos_per_site(dos_type)})

        filename = Path(output_dir) / f"{comp.comp.CompID}_{dos_type}dos_persite.json"
        with open(filename, "w") as f:
            json.dump(jsanitize(output_dos), f)

    def process_comp_save_smeareddos(self, comp: SingleHDPanalysis):
        """Calls fuction that uses Pymatgen CompleteDos object to retrieve DOS, add Gaussian smearing, and save to json.
        Uses attributes self.dostype, self.smearing_width and self.output_dir to determine which DOSCAR to parse,
        how much smearing to ad,d and where to store the json. These attributes are set by the calling function.

        Args:
            comp (SingleHDPanalysis): analysis object for specific composition
        """
        smeared_dos = comp.parse_dos_withsmearing(self.dostype, self.smearing_width)
        filename = (
            Path(self.output_dir)
            / f"{comp.comp.CompID}_{self.dostype}smeareddos_persite.json"
        )
        with open(filename, "w") as f:
            json.dump(jsanitize(smeared_dos), f)

        return

    def process_comp_lobsterdata_basis(self, comp: SingleHDPanalysis):
        """Calls function that returns analysis on lobster outputs for a specific composition.

        Args:
            comp (SingleHDPanalysis): analysis object for specific composition.

        Returns:
            dict[compID, dict[lobster_analysis]]
        """
        return {comp.comp.CompID: comp.analyze_lobster_outputs(self.charge_type)}

    def get_lobquality_data(self):
        """_summary_Get dataframe of lobster projection quality for all calculated basissets for all compositions

        Returns:
            DataFrame[compID, basisset]: quality overview
        """

        print(
            f"{time.strftime('%H:%M:%S')}: Starting lobster quality analysis for {len(self.analysis_list)} compounds"
        )
        tempfile = "lobquality_intermediate.json.tmp"
        recordlist = []

        if os.path.exists(tempfile):
            with open(tempfile, "r") as file:
                recordlist = json.load(file)
                n_recovered = len(recordlist)
        else:
            n_recovered = 0

        with mp.Pool(processes=8, maxtasksperchild=4) as pool, tqdm(
            total=len(self.analysis_list[n_recovered:]), desc="Processing compounds"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_quality, self.analysis_list[n_recovered:]
                )
            ):
                recordlist.append(result)
                # if len(recordlist) >= n_recovered:
                # self.verbose = True
                if len(recordlist) % 50 == 0:
                    with open(tempfile, "w") as f:
                        json.dump(recordlist, f)
                pbar.update()

        records = []
        for compdict in recordlist:
            for outer_k, inner_dict in compdict.items():
                for mid_k, vals in inner_dict.items():
                    row = {"comp": outer_k, "basisset#": mid_k}
                    row.update(vals)
                    records.append(row)
        df = pd.DataFrame(records).set_index(["comp", "basisset#"]).sort_index()
        print(
            f"{time.strftime('%H:%M:%S')}: Finished lobster quality collection for {len(self.analysis_list)} compounds"
        )
        return df

    def get_basic_data(self):
        """Calls process_comp_basic() with MultiProcess for all compositions and gives results in a DataFrame.

        Returns:
            DataFrame[compID, basic_data]
        """
        recordlist = []
        n_recovered = 0

        backup_file = f"basicdata_intermediate.json.tmp"
        if os.path.exists(backup_file):
            with open(backup_file, "r") as f:
                recordlist = json.load(f)
                n_recovered = len(recordlist)

        print(
            f"{time.strftime('%H:%M:%S')}: Collecting BasicInfo for {len(self.analysis_list[n_recovered:])} compounds"
        )
        with mp.Pool(processes=8, maxtasksperchild=4) as pool, tqdm(
            total=len(self.analysis_list[n_recovered:]), desc=f"Parsing BasicInfo"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_basic, self.analysis_list[n_recovered:]
                )
            ):
                recordlist.append(result)
                if len(recordlist) % 50 == 0:
                    with open(backup_file, "w") as f:
                        json.dump(jsanitize(recordlist), f)
                pbar.update()

        records = []
        for compdict in recordlist:
            for key, vals in compdict.items():
                row = {"comp": key}
                row.update(vals)
                records.append(row)

            # for mid_k, vals in inner_dict.items():
            #     row = {'comp': outer_k, 'spin': mid_k}
            #     row.update(vals)
            #     records.append(row)

        df = pd.DataFrame(records).set_index("comp")

        return df

    def get_lobster_data(self, charge_type: Literal["Mulliken", "Loewdin"]):
        """Calls process_comp_lobsterdata_basis with MultiProcess and collects results in a DataFrame.

        Args:
            charge_type (Literal[&quot;Mulliken&quot;, &quot;Loewdin&quot;]): which charge partition scheme to base charge analysis on.

        Returns:
            DataFrame[compID, lobster_data]
        """
        recordlist = []
        n_recovered = 0

        self.charge_type = charge_type

        backup_file = f"basiclobsterdata_intermediate.json.tmp"
        if os.path.exists(backup_file):
            with open(backup_file, "r") as f:
                recordlist = json.load(f)
                n_recovered = len(recordlist)

        print(
            f"{time.strftime('%H:%M:%S')}: Collecting BasicLobsterInfo for {len(self.analysis_list[n_recovered:])} compounds"
        )
        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list[n_recovered:]), desc=f"Parsing Lobster Info"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_lobsterdata_basis,
                    self.analysis_list[n_recovered:],
                )
            ):
                recordlist.append(result)
                if len(recordlist) % 10 == 0:
                    with open(backup_file, "w") as f:
                        json.dump(jsanitize(recordlist), f)
                pbar.update()

        records = []
        for compdict in recordlist:
            for key, vals in compdict.items():
                row = {"comp": key}
                row.update(vals)
                records.append(row)

            # for mid_k, vals in inner_dict.items():
            #     row = {'comp': outer_k, 'spin': mid_k}
            #     row.update(vals)
            #     records.append(row)

        df = pd.DataFrame(records).set_index("comp").sort_index()

        return df

    def get_edgecontr_data(self, dostype: Literal["vasp", "lso", "lobster"] = "vasp"):
        """Calls process_comp_edgecontr() with MultiProcess for all compositions and collects results in a DataFrame.

        Args:
            dostype (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;], optional): which DOSCAR type to use for
                                                                                                analysis. Defaults to "vasp".
        """

        def add_condtype(banddf):
            metallic_index = banddf[banddf["bandgap"] <= 0.1].index
            semicond_index = banddf[
                (banddf["bandgap"] > 0.1) & (banddf["bandgap"] < 4.0)
            ].index
            ins_index = banddf[banddf["bandgap"] >= 4.0].index

            banddf.loc[metallic_index, "cond_type"] = "metallic"
            banddf.loc[semicond_index, "cond_type"] = "semiconductor"
            banddf.loc[ins_index, "cond_type"] = "insulator"
            ##Track down half-metals
            for comp in banddf.index.get_level_values(0):

                if banddf.loc[(comp, "combined")]["cond_type"] == "metallic" and (
                    banddf.loc[(comp, Spin.up)]["cond_type"] != "metallic"
                    or banddf.loc[(comp, Spin.down)]["cond_type"] != "metallic"
                ):
                    banddf.loc[(comp, "combined"), "cond_type"] = "half-metal"

            return banddf

        self.dostype = dostype
        recordlist = []
        n_recovered = 0

        backup_file = f"bandgapdata_{dostype}_intermediate.json.tmp"
        if os.path.exists(backup_file):
            with open(backup_file, "r") as f:
                recordlist = json.load(f)
                n_recovered = len(recordlist)

        print(
            f"{time.strftime('%H:%M:%S')}: Starting edge contribution analysis for {len(self.analysis_list[n_recovered:])} compounds"
        )
        with mp.Pool(processes=8, maxtasksperchild=4) as pool, tqdm(
            total=len(self.analysis_list[n_recovered:]),
            desc=f"Processing bandedges for {dostype}-DOS",
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_edgecontr, self.analysis_list[n_recovered:]
                )
            ):
                recordlist.append(result)
                if len(recordlist) % 50 == 0:
                    with open(backup_file, "w") as f:
                        json.dump(jsanitize(recordlist), f)
                pbar.update()

        records = []
        for compdict in recordlist:
            for outer_k, inner_dict in compdict.items():

                for mid_k, vals in inner_dict.items():
                    row = {"comp": outer_k, "spin": mid_k}
                    row.update(vals)
                    records.append(row)

        df = pd.DataFrame(records).set_index(["comp", "spin"]).sort_index()

        df = add_condtype(df)

        return df

    def save_parsedcohp_all(
        self,
        coxx_type: Literal["COHP", "COBI"] = "COHP",
        output_dir: PathLike | str = ".",
    ):
        """Calls process_comp_saveparsedcohp with MultiProcess. Parses specified COXX-type and save to json files.

        Args:
            coxx_type (Literal[&quot;COHP&quot;, &quot;COBI&quot;], optional): Whether to parse COHPCAR or COBICAR. Defaults to "COHP".
            output_dir (PathLike | str, optional): Where to save the json files. Defaults to ".".
        """
        self.coxx_type = coxx_type
        self.output_dir = output_dir

        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list), desc=f"Saving parsed {coxx_type} to JSON"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_saveparsedcohp, self.analysis_list
                )
            ):
                pbar.update()
        return

    def save_smearedcohp_all(
        self,
        coxx_type: Literal["COHP", "COBI"] = "COHP",
        smearing_width: float = 0.1,
        output_dir: PathLike | str = ".",
    ):
        """Calls process_comp_savecohp_smeared with MultiProcess. Parses specified COXX-type, smeares and save to json files.

        Args:
            coxx_type (Literal[&quot;COHP&quot;, &quot;COBI&quot;], optional): Whether to parse COHPCAR or COBICAR. Defaults to "COHP".
            output_dir (PathLike | str, optional): Where to save the json files. Defaults to ".".
        """
        self.coxx_type = coxx_type
        self.smearing_width = smearing_width
        self.output_dir = output_dir

        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list), desc=f"Saving parsed {coxx_type} to JSON"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_savecohp_smeared, self.analysis_list
                )
            ):
                pbar.update()
        return

    def save_dosjsons_perspecies(
        self,
        dostype: Literal["vasp", "lso", "lobster"] = "vasp",
        output_dir: PathLike | str = ".",
    ):
        """Calls process_comp_savedos_perspecies for all composition using MultiProcess, which saves the
        DOS per (A/B1/B2/X) to json files.

        Args:
            dostype (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;], optional): which DOSCAR type to parse. Defaults to "vasp".
            output_dir (PathLike | str, optional): where to store the json files. Defaults to ".".
        """
        self.dostype = dostype
        self.output_dir = output_dir

        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list), desc=f"Saving DOS to JSON for {dostype}-DOS"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_savedos_perspecies, self.analysis_list
                )
            ):
                # recordlist.append(result)
                # if len(recordlist)%50 == 0:
                #     with open(backup_file,'w') as f:
                #         json.dump(jsanitize(recordlist),f)
                pbar.update()

    def save_dosjsons_persite_all(
        self,
        dostype: Literal["vasp", "lso", "lobster"] = "vasp",
        output_dir: PathLike | str = ".",
    ):
        """Calls process_comp_savedos_persite for all composition using MultiProcess,
        which saves the DOS site resolved to json files.

        Args:
            dostype (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;], optional): which DOSCAR type to parse. Defaults to "vasp".
            output_dir (PathLike | str, optional): where to store the json files. Defaults to ".".
        """
        self.dostype = dostype
        self.output_dir = output_dir

        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list), desc=f"Saving DOS to JSON for {dostype}-DOS"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_savedos_persite, self.analysis_list
                )
            ):
                # recordlist.append(result)
                # if len(recordlist)%50 == 0:
                #     with open(backup_file,'w') as f:
                #         json.dump(jsanitize(recordlist),f)
                pbar.update()

    def save_dosjsons_smeared(
        self,
        dostype: Literal["vasp", "lso", "lobster"] = "vasp",
        smearing_width: float = 0.1,
        output_dir: PathLike | str = ".",
    ):
        """Calls process_comp_save_smeareddos for each composition which parses the DOSCAR and applies additional
        Gaussian smearing and saves the results to json files.

        Args:
            dostype (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;], optional): Which DOSCAR type to parse. Defaults to "vasp".
            smearing_width (float, optional): How much Gaussian smearing to apply. Defaults to 0.1.
            output_dir (PathLike | str, optional): Where to store the json files. Defaults to ".".
        """
        self.dostype = dostype
        self.smearing_width = smearing_width
        self.output_dir = output_dir

        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list), desc=f"Saving DOS to JSON for {dostype}-DOS"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_save_smeareddos, self.analysis_list
                )
            ):
                # recordlist.append(result)
                # if len(recordlist)%50 == 0:
                #     with open(backup_file,'w') as f:
                #         json.dump(jsanitize(recordlist),f)
                pbar.update()

    def save_completedos_all(
        self,
        dostype: Literal["vasp", "lso", "lobster"] = "vasp",
        output_dir: PathLike | str = ".",
    ):
        """Calls the process_comp_completedos for all compositions and stores the CompleteDos objects as json files.

        Args:
            dostype (Literal[&quot;vasp&quot;, &quot;lso&quot;, &quot;lobster&quot;], optional): Which DOSCAR type to parse. Defaults to "vasp".
            output_dir (PathLike | str, optional): Where to store the json files. Defaults to ".".
        """
        self.dostype = dostype
        self.output_dir = output_dir

        with mp.Pool(processes=8, maxtasksperchild=2) as pool, tqdm(
            total=len(self.analysis_list), desc=f"Saving DOS to JSON for {dostype}-DOS"
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(self.process_comp_completedos, self.analysis_list)
            ):
                # recordlist.append(result)
                # if len(recordlist)%50 == 0:
                #     with open(backup_file,'w') as f:
                #         json.dump(jsanitize(recordlist),f)
                pbar.update()

    def process_comp_structuraldata(self, comp: SingleHDPanalysis):
        """Collects data from the relaxed structure for a specific composition.
        Reports lattice parameters, interatomic distances, sizes of the octahedra

        Args:
            comp (SingleHDPanalysis): analysis object for specific composition.

        Returns:
            dict[compID, dictstructural_data]]
        """
        # halide ionic radius in Angstrom to compare octahedral mismatch with ionic radius based one
        r_X = comp.get_ionicradii_based_data()["r_X"] / 100
        assert (
            round(comp.structure.lattice.a, 3)
            == round(comp.structure.lattice.b, 3)
            == round(comp.structure.lattice.c, 3)
        ), f"{comp.comp.CompID} has unequal lattice {comp.structure.lattice.abc};\nOnly cubic structures are supported for structural data extraction"
        structural_dict = {
            "lattice_a_primitive": comp.structure.lattice.a,
            "lattice_a_conventional": comp.structure.to_conventional().lattice.a,
            "distance_B1_X": comp.structure.get_distance(0, 5),
            "distance_B1_A": comp.structure.get_distance(0, 2),
            "size_Oh_B1": comp.structure.get_distance(-2, -1, jimage=(-1, -1, 1)),
            "distance_B2_X": (
                comp.structure.get_distance(1, 5) if comp.comp.B2 != "Vac" else None
            ),
            "distance_B2_A": (
                comp.structure.get_distance(1, 2) if comp.comp.B2 != "Vac" else None
            ),
            "size_Oh_B2": comp.structure.get_distance(-2, -1, jimage=(0, 0, 0)),
        }
        Oh_mismatch = round(
            np.abs(structural_dict["size_Oh_B1"] - structural_dict["size_Oh_B2"])
            / (4 * r_X),
            6,
        )

        structural_dict.update({"calc_oct_mismatch": Oh_mismatch})
        return {comp.comp.CompID: structural_dict}

    def get_structural_data(self):
        """Calls process_comp_structuraldata for all compositions using MultiProcess
        and collects the results in a DataFrame.

        Returns:
            DataFrame[compID,structural_data]
        """
        recordlist = []
        n_recovered = 0
        backup_file = f"structuraldata_intermediate.json.tmp"
        if os.path.exists(backup_file):
            with open(backup_file, "r") as f:
                recordlist = json.load(f)
                n_recovered = len(recordlist)
        print(
            f"{time.strftime('%H:%M:%S')}: Starting structural data analysis for {len(self.analysis_list[n_recovered:])} compounds"
        )
        with mp.Pool(processes=8, maxtasksperchild=4) as pool, tqdm(
            total=len(self.analysis_list[n_recovered:]),
            desc=f"Processing structural data",
        ) as pbar:
            for _, result in enumerate(
                pool.imap_unordered(
                    self.process_comp_structuraldata, self.analysis_list[n_recovered:]
                )
            ):
                recordlist.append(result)
                if len(recordlist) % 50 == 0:
                    with open(backup_file, "w") as f:
                        json.dump(jsanitize(recordlist), f)
                pbar.update()

        records = []
        for compdict in recordlist:
            for key, vals in compdict.items():
                row = {"comp": key}
                row.update(vals)
                records.append(row)
        df = pd.DataFrame(records).set_index("comp").sort_index()
        return df

    @staticmethod
    def combine_dfs(
        banddf: pd.DataFrame,
        lobsterdf: pd.DataFrame,
        structuraldf: pd.DataFrame,
        basicdf: pd.DataFrame,
        add_transitions: bool = True,
        include_x: bool = True,
    ):
        """Combines the DataFrames containing basic, bandedge, lobster output, and structural data into a single DataFrame.
        From the bandedge data only the 'combined' spin channel data is retained.

        Additional functionality isolates data on the main contributing bands from the band edges which was used for plotting purposes.

        Args:
            banddf (pd.DataFrame): DataFrame of bandedgeanalysis
            lobsterdf (pd.DataFrame): DataFrame of lobster output Analysis
            structuraldf (pd.DataFrame): DataFrame of structural data
            basicdf (pd.DataFrame): DataFrame of basic/input data
            add_transitions (bool, optional): Whether to add seperate description of main orbitals involved in transition. Defaults to True.
            include_x (bool, optional): Whether to include the X-site contributions to the transition. Defaults to True.
        """

        def determineTransitionBands(combined_df: pd.DataFrame, include_x: bool = True):
            # figure out which site contributes most to VBM and CBM
            # then get the orbital character of that site
            if include_x:
                vbmcontr = combined_df[
                    [
                        "VBMtotcontr.A",
                        "VBMtotcontr.B1",
                        "VBMtotcontr.B2",
                        "VBMtotcontr.X",
                    ]
                ]
                cbmcontr = combined_df[
                    [
                        "CBMtotcontr.A",
                        "CBMtotcontr.B1",
                        "CBMtotcontr.B2",
                        "CBMtotcontr.X",
                    ]
                ]
            else:
                vbmcontr = combined_df[
                    ["VBMtotcontr.A", "VBMtotcontr.B1", "VBMtotcontr.B2"]
                ]
                cbmcontr = combined_df[
                    ["CBMtotcontr.A", "CBMtotcontr.B1", "CBMtotcontr.B2"]
                ]
            edgedf = pd.DataFrame(index=combined_df.index)
            edgedf["vbmsite"] = [
                x.split(".")[-1] for x in vbmcontr.idxmax(axis=1).astype(str).values
            ]
            edgedf["cbmsite"] = [
                x.split(".")[-1] for x in cbmcontr.idxmax(axis=1).astype(str).values
            ]
            # edgedf['vbmband'] = [combined_df[f'VBMorbital.{x}'] for x in edgedf['vbmsite'].astype(str).values]

            for site, idxs in edgedf.groupby("vbmsite").groups.items():

                edgedf.loc[idxs, "vbmorbital"] = combined_df.loc[idxs][
                    f"VBMorbital.{site}"
                ]
                edgedf.loc[idxs, "vbmorbchar"] = combined_df.loc[idxs][
                    f"VBMorbchar.{site}"
                ]
                edgedf.loc[idxs, "vbmorborder"] = combined_df.loc[idxs][
                    f"VBMorborder.{site}"
                ]
                edgedf.loc[idxs, "vbmorbcontr"] = combined_df.loc[idxs][
                    f"VBMorbcontr.{site}"
                ]

            for site, idxs in edgedf.groupby("cbmsite").groups.items():

                edgedf.loc[idxs, "cbmorbital"] = combined_df.loc[idxs][
                    f"CBMorbital.{site}"
                ]
                edgedf.loc[idxs, "cbmorbchar"] = combined_df.loc[idxs][
                    f"CBMorbchar.{site}"
                ]
                edgedf.loc[idxs, "cbmorborder"] = combined_df.loc[idxs][
                    f"CBMorborder.{site}"
                ]
                edgedf.loc[idxs, "cbmorbcontr"] = combined_df.loc[idxs][
                    f"CBMorbcontr.{site}"
                ]

            # print(edgedf.shape)
            edgedf = edgedf.loc[edgedf[["vbmorbchar", "cbmorbchar"]].dropna().index]
            edgedf["transition_sites"] = [
                x[0] + "-" + y[0] for x, y in zip(edgedf["vbmsite"], edgedf["cbmsite"])
            ]
            # print(edgedf.shape)
            edgedf["transition_bands"] = [
                x + "-" + y for x, y in zip(edgedf["vbmorbchar"], edgedf["cbmorbchar"])
            ]

            return edgedf

        comb_df = pd.concat(
            [basicdf, banddf.xs("combined", level=1), structuraldf, lobsterdf],
            axis=1,
            join="inner",
        )
        if add_transitions:
            trans_df = determineTransitionBands(comb_df, include_x=include_x)
            comb_df = comb_df.join(trans_df)
        # print(comb_df.head())
        return comb_df

    @staticmethod
    def remove_duplicates(
        combineddf: pd.DataFrame,
        banddf: pd.DataFrame,
        lobsterdf: pd.DataFrame,
        structuraldf: pd.DataFrame,
        basicdf: pd.DataFrame,
        to_sim_overview: PathLike | None = None,
        save_dir: PathLike | None = None,
    ):
        """Previous iterations of the composition list featured some duplicate entries (e.g., Cs2HgPdCl6 and Cs2PdHgCl)
        and a few that were assumed stable based on wrong oxidation state assignment (+I,+III in stead of +II, +II).
        These were later removed from the candidate list, but some made it into the data set.
        This function removes duplicate and unwanted entries from all DataFrames

        Args:
            combdf (pd.DataFrame): DataFrame of CombinedInfo
            banddf (pd.DataFrame): DataFrame of band edge analysis
            lobsterdf (pd.DataFrame): DataFrame of lobster output analysis
            structuraldf (pd.DataFrame): DataFrame of structural data
            basicdf (pd.DataFrame): DataFrame of basic/input data
            to_sim_overview (PathLike | None, optional): Path to overview of reworked candidate list
                                                        if None is given this step is fixed and only duplicates are removed.
                                                        Defaults to None.
            save_dir (PathLike | None, optional): Path to where DataFrames without duplicates should be stores. Defaults to None.

        Returns:
            5 * DataFrame: Each input DataFrame with unwanted entries removed.
        """
        combdf = combineddf.copy()
        combdf["element.B2"] = combdf["element.B2"].fillna("Vac")
        combdf["el_set"] = combdf.apply(
            lambda row: {row["element.B1"], row["element.B2"], row["element.X"]}, axis=1
        )
        index_nodupl = combdf["el_set"].drop_duplicates().index
        combdf = combdf.loc[index_nodupl]
        print("Shape without duplicates:", combdf.shape)

        if to_sim_overview is not None:
            to_sim_df = pd.read_csv(to_sim_overview)
            to_sim_df["el_set"] = to_sim_df.apply(
                lambda row: {row["B1"], row["B2"], row["X"]}, axis=1
            )
            combdf["in_sim_list"] = [
                x in to_sim_df["el_set"].to_list() for x in combdf["el_set"]
            ]
            index_nodupl = combdf[combdf["in_sim_list"] == True].index
            combdf = combdf.loc[index_nodupl]

        if save_dir is not None:
            print("Shape after checking to_sim_list", combdf.shape)
            data_output_dir = Path(save_dir)
            basicdf.loc[index_nodupl].to_csv(
                f'{data_output_dir}/HDP_BasicInfo_{time.strftime("%y%m%d")}.csv'
            )
            lobsterdf.loc[index_nodupl].to_csv(
                f'{data_output_dir}/HDP_LobsterInfo_{time.strftime("%y%m%d")}.csv'
            )
            banddf.loc[index_nodupl].to_csv(
                f'{data_output_dir}/HDP_bandedgeInfo_lsodos_{time.strftime("%y%m%d")}.csv'
            )
            structuraldf.loc[index_nodupl].to_csv(
                f'{data_output_dir}/HDP_StructuralInfo_{time.strftime("%y%m%d")}.csv'
            )
            combineddf.loc[index_nodupl].to_csv(
                data_output_dir / f"HDP_CombinedInfo_{time.strftime("%y%m%d")}.csv"
            )

        return (
            combineddf.loc[index_nodupl],
            banddf.loc[index_nodupl],
            lobsterdf.loc[index_nodupl],
            structuraldf.loc[index_nodupl],
            basicdf.loc[index_nodupl],
        )


if __name__ == "__main__":
    data_output_dir = Path("./AnalysisResults")
    g = GroupedAnalysis(p, testing=0)
    #
    print("starting basis selection")
    g.select_lobbasis(
        saved_qualityoverview=data_output_dir / "HDP_lobqual_overview_260509.csv"
    )
    # g.select_lobbasis(save_badprojections=True, saving_dir=data_output_dir)

    # output_dos_base = data_output_dir / 'DOSJSONS'
    # dostypes = ['lso','lobster']
    # output_extras = [ 'LSODOS', 'LOBSTERDOS']
    # for typedos,outputextra in zip(dostypes,output_extras):
    #     outputjson = output_dos_base / outputextra
    #     outputjson.mkdir(exist_ok=True)
    #     g.save_dosjsons_all(typedos,outputjson)

    ### save smeared DOS jsons
    # smeared_dos_output = data_output_dir / "LSO_DOS_SMEARED"
    # smeared_dos_output.mkdir(exist_ok=True)
    # g.save_dosjsons_smeared(dostype='lso',smearing_width=0.05,output_dir=smeared_dos_output)

    # smeared_cohp_output = data_output_dir / "COHP_SMEARED"
    # smeared_cohp_output.mkdir(exist_ok=True)
    # g.save_smearedcohp_all(coxx_type='COHP',smearing_width=0.05, output_dir=smeared_cohp_output)

    # smeared_cobi_output = data_output_dir / "COBI_SMEARED"
    # smeared_cobi_output.mkdir(exist_ok=True)
    # g.save_smearedcohp_all(coxx_type='COBI',smearing_width=0.05, output_dir=smeared_cobi_output)

    # #Recompile dfs
    # dbasic = g.get_basic_data()
    # dbasic.to_csv(f'{data_output_dir}/HDP_BasicInfo_{time.strftime("%y%m%d")}.csv')
    # dlobster = g.get_lobster_data("Loewdin")
    # dlobster.to_csv(f'{data_output_dir}/HDP_LobsterInfo_{time.strftime("%y%m%d")}.csv')
    # dband = g.get_edgecontr_data("lso")
    # dband.to_csv(f'{data_output_dir}/HDP_bandedgeInfo_lsodos_{time.strftime("%y%m%d")}.csv')
    # dstuc = g.get_structural_data()
    # dstuc.to_csv(f'{data_output_dir}/HDP_StructuralInfo_{time.strftime("%y%m%d")}.csv')

    # #read previous data in stead of collecting
    dbasic = pd.read_csv(data_output_dir / "HDP_BasicInfo_260510.csv", index_col=0)
    dlobster = pd.read_csv(data_output_dir / "HDP_LobsterInfo_260510.csv", index_col=0)
    dband = pd.read_csv(
        data_output_dir / "HDP_bandedgeInfo_lsodos_260510.csv", index_col=[0, 1]
    )
    dstuc = pd.read_csv(data_output_dir / "HDP_StructuralInfo_260510.csv", index_col=0)

    # #####combine dataframes into one single df
    dcomb = GroupedAnalysis.combine_dfs(dband, dlobster, dstuc, dbasic)

    # Remove Duplicate entries
    # nodup_output = data_output_dir / "HDP_Data_NoDups"
    # nodup_output.mkdir(exist_ok=True)
    # GroupedAnalysis.remove_duplicates(
    #     combineddf=dcomb,
    #     banddf=dband,
    #     lobsterdf=dlobster,
    #     structuraldf=dstuc,
    #     basicdf=dbasic,
    #     to_sim_overview=Path(
    #         "/home/lwalterb/hdp_project/HDP_WorkFlow_Analysis/ToSimulateComps.csv"
    #     ),
    #     save_dir=nodup_output,
    # )

    # outputcohps = Path(data_output_dir)/"ParsedCOHPs"
    # outputcohps.mkdir(exist_ok=True)
    # g.save_parsedcohp_all('COHP',outputcohps)

    # outputcobis = Path(data_output_dir)/"ParsedCOBIs"
    # outputcobis.mkdir(exist_ok=True)
    # g.save_parsedcohp_all('COBI',outputcobis)
    # # output_dos = Path("/home/lwalterb/hdp_project/NewWF_Analysis/AnalysisResultsPsiK/DOSJSONS/LSODOS")
    # # g.save_completedos_all(dostype='lso',output_dir=output_dos)

    # outputdos = Path(data_output_dir)/"parsedLSODOS"
    # outputdos.mkdir(exist_ok=True)
    # g.save_dosjsons_persite_all('lso',outputdos)
