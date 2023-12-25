#  Copyright glm_mda_diffusion (C) 2023 Radost Waszkiewicz
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import re  # Parsing sequence descriptions
import numpy as np

import pygrpy
import sarw_spheres

default_aminoacid_masses = {
    "A": 71.08,
    "C": 103.14,
    "D": 115.09,
    "E": 129.12,
    "F": 147.18,
    "G": 57.06,
    "H": 137.15,
    "I": 113.17,
    "K": 128.18,
    "L": 113.17,
    "M": 131.21,
    "N": 114.11,
    "P": 97.12,
    "Q": 128.41,
    "R": 156.2,
    "S": 87.08,
    "T": 101.11,
    "V": 99.14,
    "W": 186.21,
    "Y": 163.18,
    "Z": 0,
    "O": 0,
    "U": 0,
    "J": 0,
    "X": 0,
    "B": 0,
}  # Units: Da


def protein_hydrodynamic_radius(
    sequence,
    steric_radius=1.9025,  # Ang
    hydrodynamic_radius=4.2,  # Ang
    effective_density=0.52,  # Da / Ang^3
    hydration_thickness=3.0,  # Ang
    ensemble_size=30,
    bootstrap_rounds=10,
    aminoacid_masses=default_aminoacid_masses,
    progress_bar=False,
):
    """
    Compute the hydrodynamic radius of a protein based on the provided sequence and parameters.

    Parameters
    ----------
    sequence : str
        Sequence of the protein. Structured domains inside square brackets.
    steric_radius : float, optional
        Steric size of linker beads in Angstroms. Default is 1.9025.
    hydrodynamic_radius : float, optional
        Hydrodynamic size of linker beads in Angstroms. Default is 4.2.
    effective_density : float, optional
        Effective density of large structured domains in Dalton/Angstrom^3. Default is 0.52.
    hydration_thickness : float, optional
        Thickness of hydration shell of domains in Angstroms. Default is 3.0.
    ensemble_size : int, optional
        Number of conformers used in calculations. Default is 30.
    bootstrap_rounds : int or None, optional
        If None no error estimation is done. Number of bootstrap rounds used in Monte Carlo error estimation. Default is 10.
    aminoacid_masses : dict, optional
        Dictionary mapping amino acid symbols to their respective masses in Dalton. Default values provided.
    progress_bar : bool, optional
        If True, display a progress bar for ensemble calculations. Default is False.

    Returns
    -------
    dict
        A dictionary containing:
            - 'protein_rh': Computed GLM-MDA hydrodynamic radius in Angstroms.
            - 'protein_rh_sigma': Monte Carlo uncertainty of the hydrodynamic radius in Angstroms.
                                  Only provided if bootstrap_rounds is not None and > 0.

    """

    bdc = _bead_description_compact(
        annotated_sequence=sequence,
        steric_radius=steric_radius,
        hydrodynamic_radius=hydrodynamic_radius,
        effective_density=effective_density,
        hydration_thickness=hydration_thickness,
        aminoacid_masses=aminoacid_masses,
    )

    bead_steric_radii = _bead_steric_radii(bdc)
    bead_hydrodynamic_radii = _bead_hydrodynamic_radii(bdc)

    conformation_ensemble = []

    if progress_bar:
        import tqdm
        ensemble_iterable = tqdm.tqdm(range(ensemble_size))
    else:
        ensemble_iterable = range(ensemble_size)

    for i in ensemble_iterable:
        conformation_ensemble.append(
            sarw_spheres.generateChain(np.array(bead_steric_radii))
        )

    grand_grand_mu = np.array(
        [
            pygrpy.grpy_tensors.muTT(locations, np.array(bead_hydrodynamic_radii))
            for locations in conformation_ensemble
        ]
    )

    grand_mu = np.mean(grand_grand_mu, axis=0)
    grand_trace = np.trace(grand_mu, axis1=-2, axis2=-1)

    inv_mat = np.linalg.inv(grand_trace)
    total = np.sum(inv_mat)
    effective_rh = total / (2 * np.pi)

    effective_rh_sigma = None

    if bootstrap_rounds != None:
        effective_rh_list = []

        for i in range(bootstrap_rounds):
            ensemble_size = len(grand_grand_mu)
            b_grand_grand_mu = grand_grand_mu[
                np.random.choice(ensemble_size, ensemble_size)
            ]
            b_grand_mu = np.mean(b_grand_grand_mu, axis=0)
            b_grand_trace = np.trace(b_grand_mu, axis1=-2, axis2=-1)
            b_inv_mat = np.linalg.inv(b_grand_trace)
            b_total = np.sum(b_inv_mat)
            b_effective_rh = b_total / (2 * np.pi)

            effective_rh_list.append(b_effective_rh)

        effective_rh_sigma = np.std(effective_rh_list)

    return {"protein_rh": effective_rh, "protein_rh_sigma": effective_rh_sigma}


def _bead_description_compact(
    annotated_sequence,
    hydrodynamic_radius,
    steric_radius,
    effective_density,
    hydration_thickness,
    aminoacid_masses,
):
    blocks = re.split(
        "(\[[A-Z].*?\])", annotated_sequence
    )  # things inside braces are blocks
    bead_description_compact = []
    for block in blocks:
        if len(block) >= 2 and block[0] == "[" and block[-1] == "]":
            block_mass = sum(aminoacid_masses[aa] for aa in block[1:-1])
            block_excluded_volume_radius = (
                block_mass * (3 / (4 * math.pi)) / effective_density
            ) ** (1 / 3)
            block_radius = block_excluded_volume_radius + hydration_thickness
            bead_description_compact.append([block_radius, block_radius, 1])

        elif len(block) > 0 and block[0] != "[" and block[-1] != "]":
            bead_description_compact.append(
                [hydrodynamic_radius, steric_radius, len(block)]
            )
        else:
            raise ValueError("Malformed seqence passed: " + block)

    return bead_description_compact


def _bead_steric_radii(bead_description_compact):
    return sum(
        ([sr] * n for (hr, sr, n) in bead_description_compact), []
    )  # sum is concat


def _bead_hydrodynamic_radii(bead_description_compact):
    return sum(
        ([hr] * n for (hr, sr, n) in bead_description_compact), []
    )  # sum is concat


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Minimum dissipation approximation:\n"
            + "A fast algorithm for the prediction of diffusive properties\n"
            + "of intrinsically disordered proteins.\n"
        )
    )

    parser.add_argument(
        "--sequence",
        type=str,
        required=True,
        help="Sequence of the protein. Structured domains inside square brackets",
    )

    parser.add_argument(
        "--steric_radius",
        type=float,
        default=1.9025,
        help="Steric size of linker beads. Units: Angstrom",
    )
    parser.add_argument(
        "--hydrodynamic_radius",
        type=float,
        default=4.2,
        help="Hydrodynamic size of linker beads. Units: Angstrom",
    )

    parser.add_argument(
        "--effective_density",
        type=float,
        default=0.52,
        help="Effective density of large structured domains. Units: Dalton / Angstrom^3",
    )
    parser.add_argument(
        "--hydration_thickness",
        type=float,
        default=3.0,
        help="Thickness of hydration shell. Units: Angstrom",
    )

    parser.add_argument(
        "--ensemble_size",
        type=int,
        default=30,
        help="Number of conformers used in calculations.",
    )
    parser.add_argument(
        "--bootstrap_rounds",
        type=int,
        default=None,
        help="Number of bootstrap rounds used in MC error estimation.",
    )

    args = parser.parse_args()

    result = protein_hydrodynamic_radius(
        sequence=args.sequence,
        steric_radius=args.steric_radius,
        hydrodynamic_radius=args.hydrodynamic_radius,
        effective_density=args.effective_density,
        hydration_thickness=args.hydration_thickness,
        ensemble_size=args.ensemble_size,
        bootstrap_rounds=args.bootstrap_rounds,
        progress_bar=True,
    )

    print("Computed GLM-MDA hydrodynamic radius [Ang]:")
    print(result["protein_rh"])
    if args.bootstrap_rounds != None and args.bootstrap_rounds > 0:
        print("MC uncertanity [Ang]:")
        print(result["protein_rh_sigma"])


if __name__ == "__main__":
    main()
