""" 
HDPCompound.py
Contains the Compound class to create a Halide Double Perovskite Compound with some functionality to writing input files.
"""

import numpy as np
import pandas as pd  # type: ignore
import os
import subprocess
from pathlib import Path

# ============================================================================
# Use python-dotenv to load settings from .env file.
# ============================================================================
from dotenv import load_dotenv
load_dotenv()
# Retrieve to VASP PAW PseudoPotentials
path_to_pseudo = os.environ['PATH_TO_PSEUDO']


chemical_symbols = [
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
]

atomic_numbers = {"A": 1}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z

atomic_spinstates = pd.read_csv(
    "./assets/UnpairedSpins.csv",
    index_col=[0, 1],
    skipinitialspace=True,
)


class Compound:
    """This class contains all compound specific methods. In this case this means that this class contains all basic setups for treating Halide Double Perovskites (HDPs).
    For the initialization we need to know the present ions. This iteration assumes both A and X to be monovalent.

        __init__() for this class requires A(str: element),B1(str: element),nB1(int:charge), B2(str: element),B2(str: element),nB2(int:charge),X(str: element)

    This class contains methods for
        - getting the (electron) spin states of the B-sites
        - writing POTCAR
        - writing KPOINTS grid
        - writing POSCAR
        - Getting basissets for LOBSTER projection

    """

    def __init__(
        self,
        A: str = "Cs",
        B1: str = "Pt",
        nB1: int = 2,
        B2: str = "Am",
        nB2: int = 2,
        X: str = "Cl",
    ):
        self.A = A
        self.B1 = B1
        self.nB1 = nB1
        self.B2 = B2
        self.nB2 = nB2
        self.X = X

    def __repr__(self):
        return self.A + self.B1 + self.B2 + self.X

    def get_spinstate(self, which="B1") -> tuple[float, float]:
        """
        Helper function that returns the atomic spin bohr-magnetron number (Ignoring L and J contributions to the magnetic moment).
        The magnetron numbers are retreived from the file 'UnpairedSpins.csv', and are only defined for the transition metals and
        Lanthanides and Actinides, any other atom will return a magnetic moment of 0
        Inputs:
              which:    Selection of which B-site to return the spinstates for
                        must be string 'B*' with *==1 or 2
        Outputs:
                muHS, muLS  float of mu_s in Bohr magnetrons for the High and Low spinstate
        """
        if which.upper() == "B1":
            Zatom = atomic_numbers[self.B1]
            natom = self.nB1
        elif which.upper() == "B2":
            Zatom = atomic_numbers[self.B2]
            natom = self.nB2
        else:
            raise ValueError(
                "Cannot return spinstates for given B-site, fill in either 'B1' or 'B2'"
            )

        try:
            muHS = atomic_spinstates.loc[Zatom, natom]["mus_HS"]
        except KeyError:
            muHS = 0.0

        try:
            muLS = atomic_spinstates.loc[Zatom, natom]["mus_LS"]
        except KeyError:
            muLS = 0.0

        return float(muHS), float(muLS)

    def write_geninfo(self):
        """writes a short textfile with an overview of the compound information known before any calculations"""
        filename = self.name + "_geninfo.txt"
        b1HS, b1LS = self.get_spinstate("B1")
        # b2HS, b2LS = self.get_spinstate('B2')
        with open(filename, "w") as file:
            file.write("This compound is %s \n" % (self.name))
            file.write(
                "A-site is: %s \t with atomic number: %s \n"
                % (self.A, atomic_numbers[self.A])
            )
            file.write(
                "B1-site is: %s \t with atomic number: %s \n"
                % (self.B1, atomic_numbers[self.B1])
            )
            file.write(
                "\t has charge %s, and spinstates: %s [mu_B] and %s [mu_B] \n"
                % (self.nB1, b1HS, b1LS)
            )
            if self.B2 == "Vac":
                file.write("B2-site is a vacancy \n")
            else:
                b2HS, b2LS = self.get_spinstate("B2")
                file.write(
                    "B2-site is: %s \t with atomic number: %s \n"
                    % (self.B2, atomic_numbers[self.B2])
                )
                file.write(
                    "\t has charge %s, and spinstates: %s [mu_B] and %s [mu_B] \n"
                    % (self.nB2, b2HS, b2LS)
                )
            file.write(
                "X-site is %s \t with atomic number: %s \n"
                % (self.X, atomic_numbers[self.X])
            )
        return

    def write_poscar(self, a=5.5):
        """
        This function writes a standard primitive Fm3m unit cell POSCAR into the current directory.
        It will check if B2-site is a vacancy and adjust the POSCAR accordingly
        Args:
            a (float): the paramater defining the primitive unit cell [[0, a, a], [a, 0, a],[a, a, 0]]
                        the lattice parameter will be sqrt(2*a**2)
        """

        with open("POSCAR", "w") as f:
            f.write(
                "%s \n1.00000000000000 \n%.10f    %.10f    %.10f\n"
                % (repr(self), 0.0, a, a)
            )
            f.write(
                "%.10f    %.10f    %.10f\n%.10f    %.10f    %.10f\n"
                % (a, 0.0, a, a, a, 0.0)
            )

            if self.B2 == "Vac":
                f.write("%s    %s    %s\n1    2    6\n" % (self.B1, self.A, self.X))
                f.write("Direct\n  %.16f  %.16f  %.16f\n" % (0.0, 0.0, 0.0))
            else:
                f.write(
                    "%s    %s    %s    %s\n1    1    2    6\n"
                    % (self.B1, self.B2, self.A, self.X)
                )
                f.write(
                    "Direct\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n"
                    % (0.0, 0.0, 0.0, 0.5, 0.5, 0.5)
                )

            f.write(
                "  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n"
                % (
                    0.25,
                    0.25,
                    0.25,
                    0.75,
                    0.75,
                    0.75,
                    0.75,
                    0.25,
                    0.25,
                    0.25,
                    0.75,
                    0.75,
                )
            )
            f.write(
                "  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n  %.16f  %.16f  %.16f\n"
                % (
                    0.25,
                    0.75,
                    0.25,
                    0.75,
                    0.25,
                    0.75,
                    0.25,
                    0.25,
                    0.75,
                    0.75,
                    0.75,
                    0.25,
                )
            )
        return

    def write_kpoints(self, npoints=4):
        """Creates the file KPOINTS in the current directory.
        This function only supports cubic equidistant gamma-centered kpoint-grids (as this is the relevant grid for this studycase)

        Args:
            npoints (int, optional): number of k-grid points in each direction. Defaults to 4.
        """
        with open("KPOINTS", "w") as file:
            file.write("k-points\n")
            file.write(
                " 0\nGamma \n %i  %i  %i \n 0  0  0" % (npoints, npoints, npoints)
            )
        return

    def write_potcar(self, use_frozen: bool = False, write_file: bool = True, path_to_PAWs: str = path_to_pseudo):
        """This function concatonates the POTCARS for each species defined by the object compound. Which specific POTCAR is used for each
        species is determined by potLUT (potential Look Up Table).
        If any species is given that is not present in the potLUT this function will fail

        Arguments:
            use_frozen(bool):   This specifies whether we want to use potential with extra electrons in the frozen core.
                                This option is specifically used for the relaxation step to improve stability
        """
        if use_frozen:
            LUTpath = "./assets/potLUT"
        else:
            LUTpath = "./assets/potLUT_nofrozen"

        PPpath = path_to_PAWs
        if os.path.isfile("POTCAR") and write_file:
            os.system("rm POTCAR")

        potcar_list = []
        for el in [self.B1, self.B2, self.A, self.X]:
            # print(el)
            if el == "Vac":
                continue
            awk_proc = subprocess.Popen(
                ["awk", "/%s\t/{print $2}" % (el), LUTpath],
                stdout=subprocess.PIPE,
                text=True,
            )
            potname, _ = awk_proc.communicate()
            potcar_list.append(potname.strip())
            if write_file:
                potdir = PPpath + potname.strip() + "/POTCAR"
                os.system("cat %s >> POTCAR" % (potdir))

        return potcar_list

    def get_basis_functions(self, basis="min"):
        """This function will look at BASISPBE.yaml file (adapted from the Lobsterpy package) and return the set of basisfunctions
        that should be used by LOBSTER. There are two sets to chose from, the minimal or maximal basis set.

        Args:
            basis (str, optional): Should be either "min" or "max" and determines which basisset-file to look in. Defaults to "min".

        Returns:
            dict: a Dictonary of 'element':'basis functions' that are to be used in LOBSTER projection
        """

        basisfile = (
            f"./assests/BASIS_PBE_64_{basis}.yaml"
        )
        LUTpath = "./assets/potLUT_nofrozen"

        with open(LUTpath, "rt") as f:
            LUTlines = f.readlines()

        with open(basisfile, "rt") as bases:
            basislines = bases.readlines()

        basis_dict = {}

        for el in [self.B1, self.B2, self.A, self.X]:
            if el == "Vac":
                continue

            potLUTline = [line for line in LUTlines if line.startswith(f"{el}\t")]
            assert (
                len(potLUTline) == 1
            ), f"Looking up the required POTCAR fo r{el} returned multilple options: {potLUTline}"
            potname = potLUTline[0].strip("\n").strip().split("\t")[1]

            basisline = [line for line in basislines if f"{potname}:" in line]
            assert (
                len(basisline) == 1
            ), f"Looking up the required LOBSTER basis for {potname} returned multilple options: {basisline}"
            basisfuncs = basisline[0].strip("\n").strip().split(":")[1]

            basis_dict.update({el: basisfuncs})

        return basis_dict

    def get_all_basis_combos(self):
        """This function will return all possible basisset combinations for the LOBSTER projection. This is based on the basis sets
        available in BASIS_PBE_64_{min/max}.yaml which is ultimately derived from the listed orbitals in the POTCAR files.

        Returns:
            dict: A nested Dictonary of every basisset combination (given iteratively, i.e., basis0:..., basis1:..., etc...) and for each basis 'element':'basis;
        """

        min_basis = self.get_basis_functions(basis="min")
        max_basis = self.get_basis_functions(basis="max")

        all_base_dict = {}
        base_index = 0
        base_name = "basisset" + str(base_index)

        all_base_dict[base_name] = min_basis

        extraB1 = False
        extraB2 = False

        if min_basis[self.B1] != max_basis[self.B1]:
            dummy = min_basis.copy()
            dummy.update({self.B1: max_basis[self.B1]})

            base_index += 1
            base_name = "basisset" + str(base_index)

            all_base_dict[base_name] = dummy
            extraB1 = True

        if self.B2 != "Vac":
            if min_basis[self.B2] != max_basis[self.B2]:
                dummy = min_basis.copy()
                dummy.update({self.B2: max_basis[self.B2]})

                base_index += 1
                base_name = "basisset" + str(base_index)

                all_base_dict[base_name] = dummy

                extraB2 = True

        if extraB1 and extraB2:
            base_index += 1
            base_name = "basisset" + str(base_index)
            all_base_dict[base_name] = max_basis

        return all_base_dict

    def get_maxNBANDS(self):
        """This function will take the max basis functions for a specific HDP, and calculate the number of bands
        required in the DFT calculation for the LOBSTER projection.

        Returns:
            int: The number of bands required for the maximal basis set for this HDP composition
                    This is to be set in the INCAR file as NBANDS

        """
        max_basis = self.get_basis_functions("max")

        nbands = 0

        for el, basis in max_basis.items():
            mult = 0  # multiplier to account for number of occurences in unit cell
            if el.split("_")[0] == self.A:
                mult += 2
            if el.split("_")[0] == self.X:
                mult += 6
            if el.split("_")[0] == self.B1:
                # necessary in case of duplicate entries (e.g. Cs2TlTlCl6 or Cs2CsVBr)
                mult += 1
            if el.split("_")[0] == self.B2:
                mult += 1

            for band in basis.split():

                if band.endswith("s"):
                    nbands += 1 * mult

                elif band.endswith("p"):
                    nbands += 3 * mult

                elif band.endswith("d"):
                    nbands += 5 * mult

                elif band.endswith("f"):
                    nbands += 7 * mult

        return nbands

