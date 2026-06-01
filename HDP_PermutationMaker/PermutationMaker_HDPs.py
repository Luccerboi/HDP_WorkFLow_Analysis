# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 17:09:00 2024

@author: lucwa
"""

from mendeleev.fetch import fetch_ionic_radii
import numpy as np
import pandas as pd
import csv
from os import PathLike
from mendeleev import element


def calculate_tau(nA: int, r_A: float, r_B: float, r_X: float) -> tuple[float]:
    """This function calculates the Bartel tau-factor for a Perovskite composition.
    see DOI:10.1126/sciadv.aav0693 for more info.

    Args:
        nA (int): The oxidation state of the A-cation
        r_A (float): Ionic radius of A-cation
        r_B (float): Ionic radius of B-cation. If dealing with a double perovskite, the average of the two B ionic radii should be given.
        r_X (float): Ionic radius of X-anion. In our case restricted to (-I) because we deal with Halide Double Perovskites.


    Returns:
        float: The calculated tau value for the given ionic radii
    """
    tau = (r_X / r_B) - nA * (nA - (r_A / r_B) / np.log(r_A / r_B))
    return tau


def compare_ens(en1: float, en2: float) -> tuple[int, int]:
    """This function compares two given electronegativities and returns the oxidation states based on their differences.
    The methodology is adapted from Bartel et al., 2019 (DOI: 10.1126/sciadv.aav0693). If the electronegativites differ more
    than 10%, the oxidation state difference is maximized.

    With the restriction of Halide Double Perovskites, there are three options for the oxidation state allocation:
    - (+I, +III)
    - (+II, +II)
    - (+III, +I)

    Args:
        en1 (float): ElectroNegativity for element 1
        en2 (float): ElectroNegativity for element 2

    Returns:
        tuple[int, int]: Oxidation state for element 1, and element 2
    """
    rat = en1 / en2
    if rat < 0.9:
        return 3, 1
    elif rat > 1.1:
        return 1, 3
    else:
        return 2, 2


def determine_simulation_list(
    compositions_input_file: str | PathLike,
    simlist_output_file: str | PathLike,
    tau_cutoff: float = 4.18,
    vacancy_tau_cutoff: float = 5.46,
):
    """Takes in the .csv file with all HDP permutations and applies the tau-factor prediction.
    It also filters out all compositions with Berkellium and elements beyond since VASP is missing PAW pseudopotentials.
    A separate cutoff is used for Vacancy-Ordered HDPs to account for the vacancy.
    The output is saved as csv


    Args:
        compositions_input_file (str | PathLike): filename to the .csv with all HDP permutations
        simlist_output_file (str | PathLike): filename to save all stable and to simulate Cs-HDPs under
        tau_cutoff (float, optional): Cut-off value to use for tau in general case. Defaults to 4.18.
        vacancy_tau_cutoff (float, optional): Specific cut-off value to use for Vacancy-Orderd HDPs. Defaults to 5.46.
    """
    AllCands = pd.read_csv(compositions_input_file)
    NoVacs = AllCands[AllCands["B2"] != "Vac"]
    AllVacs = AllCands[AllCands["B2"] == "Vac"]

    print(
        "All Combs:",
        AllCands.shape,
        "\nWithout Vacs:",
        NoVacs.shape,
        "\nVO-HDPS:",
        AllVacs.shape,
    )
    StableNoVacs = NoVacs[NoVacs["tau"] <= tau_cutoff]
    StableVacs = AllVacs[AllVacs["tau"] <= vacancy_tau_cutoff]
    StableHDPs = pd.concat([StableNoVacs, StableVacs], axis=0, join="inner")
    # print(StableHDPs.tail(15))

    print(
        "\nStable Comps All:",
        StableHDPs.shape,
        "\nWithout Vacs:",
        StableNoVacs.shape,
        "\nVO-HDPs:",
        StableVacs.shape,
    )

    StableCsHDPS = StableHDPs[StableHDPs["A"] == "Cs"]
    ##Filter out the elements without a PAW potential
    ElNoPAWs = ["Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
    StableCsHDPS.loc[:, "NotSim"] = StableCsHDPS.apply(
        lambda row: row["B1"] in ElNoPAWs or row["B2"] in ElNoPAWs, axis=1
    )

    SimCsHDPs = StableCsHDPS[StableCsHDPS["NotSim"] == False]
    SimCsHDPs["comp_name"] = SimCsHDPs.apply(
        lambda row: f"{row['A']}{row['B1']}{row['B2']}{row['X']}", axis=1
    )
    SimCsHDPs.set_index("comp_name").to_csv(simlist_output_file)
    print(
        "\nStable Cs HDPs:",
        StableCsHDPS.shape,
        "\nCs HDPs to simulate:",
        SimCsHDPs.shape,
    )
    return


if __name__ == "__main__":
    """Here we fetch the ionic radii through Mendeleev. Extrapolate the ionic radii to fill potential gaps in coordination number.
    We then loop through all possible charge combinations, making all possible combinations.
    When there are 2 possible oxidation state combinations, we use the ElectroNegativity difference to resolve the issue.
    All charge neutral combinations are written to a .csv file, which is then processed to get a list of predicted stable Cs halide double perovskites.
    """

    # Defining important parameters
    filename = "All_HDP_Permutations.csv"  # filename for writing all combinations
    filename_stableCs = "CsBBX_Stable_Candidates.csv"  # filename for writing the Cs HDPs to simulate (excluding very heavy elements with no PAW PP in VASP)
    en_scale = "mulliken"  # We need to select a ElectroNegativity scale to compare elements on.

    # Set the tau cut-off values for determining the stable compositions
    tau_cutoff = 4.18  # Standard tau cut-off as determined by Bartel et al.
    tau_cutoff_VO = 5.46  # Special cut-off for vacancy ordered HDPs, to match experimentally observed VO-HDPs

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

    # Drop Hydrogen, because it's not of interest
    irs_full = irs_full.drop([1, 1])

    # Seperate periodic table into relevant charges
    chg1 = irs_full.groupby(["charge"]).get_group(1)
    chg2 = irs_full.groupby(["charge"]).get_group(2)
    chg3 = irs_full.groupby(["charge"]).get_group(3)

    # Multiple oxidation states can occur for B-sites with (+I, +II, +III) oxidation states
    # So we treat these as a single group
    bsites = pd.concat([chg1, chg2, chg3]).sort_values(by="atomic_number")

    # +IV oxidation state is only possible with
    chg4 = irs_full.groupby(["charge"]).get_group(4)
    halides = irs_full.groupby(["charge"]).get_group(-1)

    # Safe halide information for easier use
    halide_names = ("F", "Cl", "Br", "I")
    halide_irs = halides["VI"].values
    nX = -1

    # Coordination numbers 6 and 12 are relevant for the A and B sites. We save the indices to use.
    index6 = irs_full.columns.get_loc("VI")
    index12 = irs_full.columns.get_loc("XII")

    # Write the header line for the .csv file
    with open(filename, "w") as file:
        file.write("A2BBX6,A,B1,B2,X,nA,nB1,nB2,nX,rA,rB1,rB2,rBavg,rX,tau\n")

    # Iterator initialization to track progression.
    i_asite = 0
    tot_asites = len(chg1.index.get_level_values("atomic_number").unique()) + len(
        chg2.index.get_level_values("atomic_number").unique()
    )

    # Loop over all potential A-sites with +I oxidation state
    for iter_a1 in chg1.index.get_level_values("atomic_number").unique():
        print(f"treating A-site {i_asite}/{tot_asites}")
        i_asite += 1

        # Store information for this A-site candidate
        ir_a = irs_full.loc[iter_a1, 1]["XII"]
        nA = 1

        # Iterator to track potential B-Site process
        jj = 0
        bsits = len(bsites.index.get_level_values("atomic_number").unique())

        # Loop over all B-sites that have a oxidation state that fits
        for iter_b1 in bsites.index.get_level_values("atomic_number").unique():
            jj += 1
            if jj % 10 == 0:
                print(f"B1-site:{jj}/{bsits}")

            # Loop over all second B-site candidates
            for iter_b2 in (
                bsites.loc[iter_b1:, :].index.get_level_values("atomic_number").unique()
            ):
                # We use these Booleans to track what oxidation state combinations are possible
                perov_bool = False  # Indicates (+II, +II) combination is possible
                elpas_bool = False  # Indicates (+I, +III) combination is possible

                # check available charges
                chgs_b1 = bsites.loc[iter_b1, :].index
                chgs_b2 = bsites.loc[iter_b2, :].index

                # check if perovskite and/or elpasolite is possible
                # if only one combination is possible, the oxidation states will be set
                if 2 in chgs_b1 and 2 in chgs_b2 and iter_b1 != iter_b2:
                    perov_bool = True
                    nB1 = 2
                    nB2 = 2

                if 1 in chgs_b1 and 3 in chgs_b2:
                    elpas_bool = True
                    nB1 = 1
                    nB2 = 3

                elif 3 in chgs_b1 and 1 in chgs_b2:
                    elpas_bool = True
                    nB1 = 3
                    nB2 = 1

                # If both Oxidation State combinations are possible, check ENs
                if perov_bool and elpas_bool:
                    en_b1 = element(iter_b1).electronegativity(en_scale)
                    en_b2 = element(iter_b2).electronegativity(en_scale)

                    nB1, nB2 = compare_ens(en_b1, en_b2)

                # if no charge combination possible continue
                if not perov_bool and not elpas_bool:
                    continue

                # With oxidation states determined, calculate tau based on those irs
                ir_b1 = irs_full.loc[iter_b1, nB1]["VI"]
                ir_b2 = irs_full.loc[iter_b2, nB2]["VI"]

                ir_bavg = (ir_b1 + ir_b2) / 2
                if ir_bavg >= ir_a:
                    # If ionic radius of B is larger than A, the formula for tau fails (and becomes negative)
                    continue

                # We calculate tau for all four halides in one go, to save time.
                taus = calculate_tau(nA, ir_a, ir_bavg, halide_irs)

                nameA = element(iter_a1).symbol
                nameB1 = element(iter_b1).symbol
                nameB2 = element(iter_b2).symbol

                # Write the four tau values (for the four halides) to the file.
                with open(filename, "a") as file:
                    for ii in range(len(taus)):
                        cname = nameA + "2" + nameB1 + nameB2 + halide_names[ii] + "6"
                        file.write(
                            f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n"
                        )

        # Make all Vacancy ordered combinations
        for iter_b1 in chg4.index.get_level_values("atomic_number").unique():
            # calculate vacancy ordered possibilities
            nB1 = 4
            nB2 = 0

            ir_bavg = irs_full.loc[iter_b1, nB1]["VI"]
            if ir_bavg >= ir_a:
                continue

            taus = calculate_tau(nA, ir_a, ir_bavg, halide_irs)

            nameA = element(iter_a1).symbol
            nameB1 = element(iter_b1).symbol
            nameB2 = "Vac"

            with open(filename, "a") as file:
                for ii in range(len(taus)):
                    cname = nameA + "2" + nameB1 + halide_names[ii] + "6"
                    file.write(
                        f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n"
                    )

    # Make all HDPs with A-site having +II oxidation state.
    for iter_a1 in chg2.index.get_level_values("atomic_number").unique():
        # Save A-site ionic radius and start iterator
        ir_a = irs_full.loc[iter_a1, 2]["XII"]
        nA = 2
        print(f"treating A-site {i_asite}/{tot_asites}")
        i_asite += 1

        # Loop over all B-site combinations
        jj = 0
        for iter_b1 in chg1.index.get_level_values("atomic_number").unique():
            nB1 = 1
            jj += 1
            if jj % 10 == 0:
                print(f"B1-site:{jj}/{bsits}")

            for iter_b2 in (
                chg1.loc[iter_b1:, :].index.get_level_values("atomic_number").unique()
            ):
                nB2 = 1
                # Making A(2)B(1)B(1)X(-1) combinations
                if iter_b1 == iter_b2:
                    continue

                ir_b1 = irs_full.loc[iter_b1, nB1]["VI"]
                ir_b2 = irs_full.loc[iter_b2, nB2]["VI"]

                ir_bavg = (ir_b1 + ir_b2) / 2
                if ir_bavg >= ir_a:
                    continue

                taus = calculate_tau(nA, ir_a, ir_bavg, halide_irs)

                nameA = element(iter_a1).symbol
                nameB1 = element(iter_b1).symbol
                nameB2 = element(iter_b2).symbol

                with open(filename, "a") as file:
                    for ii in range(len(taus)):
                        cname = nameA + "2" + nameB1 + nameB2 + halide_names[ii] + "6"
                        file.write(
                            f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n"
                        )

        # Now make Vacancy Ordered combinations for divalent A-sites
        for iter_b1 in chg2.index.get_level_values("atomic_number").unique():
            # calculate vacancy ordered possibilities
            nB1 = 2
            nB2 = 0

            ir_bavg = irs_full.loc[iter_b1, nB1]["VI"]
            if ir_bavg >= ir_a:
                continue

            taus = calculate_tau(nA, ir_a, ir_bavg, halide_irs)

            nameA = element(iter_a1).symbol
            nameB1 = element(iter_b1).symbol
            nameB2 = "Vac"

            with open(filename, "a") as file:
                for ii in range(len(taus)):
                    cname = nameA + "2" + nameB1 + halide_names[ii] + "6"
                    file.write(
                        f"{cname},{nameA},{nameB1},{nameB2},{halide_names[ii]},{nA},{nB1},{nB2},{nX},{ir_a},{ir_b1},{ir_b2},{ir_bavg},{halide_irs[ii]},{taus[ii]}\n"
                    )

    determine_simulation_list(
        compositions_input_file=filename,
        simlist_output_file=filename_stableCs,
        tau_cutoff=tau_cutoff,
        vacancy_tau_cutoff=tau_cutoff_VO,
    )
