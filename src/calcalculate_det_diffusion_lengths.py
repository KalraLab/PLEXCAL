from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist


RESIDUE_TYPES = ("TRP", "TYR", "PHE")
TRANSFER_PARAMETERS = {
    ("TRP", "TRP"): 0.00058,
    ("TYR", "TYR"): 0.0012,
    ("PHE", "PHE"): 0.00098,
    ("TRP", "TYR"): 0.000487,
    ("TYR", "TRP"): 0.0011,
    ("PHE", "TYR"): 0.0012,
    ("PHE", "TRP"): 0.0039,
    ("TRP", "PHE"): 0.00031,
    ("TYR", "PHE"): 0.00098,
}
LIFETIMES_MICROSECONDS = {"TRP": 7.3, "TYR": 5.6, "PHE": 3.1}
EXCITATION_COUNT = 10_000
MAX_TIME_STEPS = 1_000


def _read_pdb_atoms(
    pdb_filename: str | Path,
) -> list[tuple[str, str, int, str, str, float, float, float]]:
    """Read atoms from PDB or mmCIF while preserving residue identity.

    Both PDB and mmCIF inputs retain chain ID and insertion-code identity so
    residues with the same sequence number in different chains are not merged.
    """
    path = Path(pdb_filename)
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return _read_mmcif_atoms(path)

    atoms = []
    saw_model = False
    with path.open("r", encoding="ascii", errors="replace") as pdb_file:
        for line in pdb_file:
            record = line[0:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record == "END":
                break
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip().upper()
                chain_id = line[21:22].strip()
                residue_number = int(line[22:26].strip())
                insertion_code = line[26:27].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except (ValueError, IndexError) as error:
                raise ValueError(f"Invalid PDB atom record: {line.rstrip()}") from error
            # Preserve PDB chain and insertion code, matching the mmCIF reader.
            atoms.append((residue_name, chain_id, residue_number, insertion_code, atom_name, x, y, z))
    return atoms


def _cif_value(value: str) -> str:
    """Normalize CIF null markers to an empty string."""
    return "" if value in {".", "?"} else value


def _read_mmcif_atoms(
    cif_filename: str | Path,
) -> list[tuple[str, str, int, str, str, float, float, float]]:
    """Read the first model from an mmCIF atom_site loop using standard fields."""
    import shlex

    lines = Path(cif_filename).read_text(encoding="utf-8", errors="replace").splitlines()
    atoms: list[tuple[str, str, int, str, str, float, float, float]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip().lower() != "loop_":
            i += 1
            continue
        i += 1
        headers: list[str] = []
        while i < len(lines) and lines[i].lstrip().startswith("_"):
            headers.append(lines[i].strip().split()[0].lower())
            i += 1
        if not headers or not any(h.startswith("_atom_site.") for h in headers):
            continue

        header_index = {name: idx for idx, name in enumerate(headers)}
        def column(*names: str) -> int | None:
            for name in names:
                if name.lower() in header_index:
                    return header_index[name.lower()]
            return None

        idx_group = column("_atom_site.group_pdb")
        idx_atom = column("_atom_site.auth_atom_id", "_atom_site.label_atom_id")
        idx_comp = column("_atom_site.auth_comp_id", "_atom_site.label_comp_id")
        idx_chain = column("_atom_site.auth_asym_id", "_atom_site.label_asym_id")
        idx_seq = column("_atom_site.auth_seq_id", "_atom_site.label_seq_id")
        idx_ins = column("_atom_site.pdbx_pdb_ins_code")
        idx_x = column("_atom_site.cartn_x")
        idx_y = column("_atom_site.cartn_y")
        idx_z = column("_atom_site.cartn_z")
        idx_model = column("_atom_site.pdbx_pdb_model_num")
        required = {"atom": idx_atom, "residue": idx_comp, "sequence": idx_seq,
                    "x": idx_x, "y": idx_y, "z": idx_z}
        missing = [name for name, idx in required.items() if idx is None]
        if missing:
            raise ValueError(f"mmCIF atom_site loop missing required columns: {', '.join(missing)}")

        token_buffer: list[str] = []
        first_model: str | None = None
        while i < len(lines):
            stripped = lines[i].strip()
            lower = stripped.lower()
            if not stripped or stripped.startswith("#"):
                i += 1
                continue
            if lower == "loop_" or lower.startswith("data_") or lower.startswith("save_") or stripped.startswith("_"):
                break
            try:
                row_tokens = shlex.split(stripped, comments=True, posix=True)
            except ValueError as error:
                raise ValueError(f"Could not parse mmCIF row: {stripped}") from error
            token_buffer.extend(row_tokens)
            while len(token_buffer) >= len(headers):
                row = token_buffer[:len(headers)]
                del token_buffer[:len(headers)]
                if idx_group is not None and row[idx_group].upper() not in {"ATOM", "HETATM"}:
                    continue
                model = _cif_value(row[idx_model]) if idx_model is not None else "1"
                if first_model is None:
                    first_model = model
                if model != first_model:
                    continue
                try:
                    residue_name = _cif_value(row[idx_comp]).upper()
                    atom_name = _cif_value(row[idx_atom])
                    chain_id = _cif_value(row[idx_chain]) if idx_chain is not None else ""
                    residue_number = int(_cif_value(row[idx_seq]))
                    insertion_code = _cif_value(row[idx_ins]) if idx_ins is not None else ""
                    x = float(row[idx_x]); y = float(row[idx_y]); z = float(row[idx_z])
                except (ValueError, IndexError) as error:
                    raise ValueError(f"Invalid mmCIF atom_site row: {' '.join(row)}") from error
                atoms.append((residue_name, chain_id, residue_number, insertion_code, atom_name, x, y, z))
            i += 1
        # A partial row indicates malformed atom_site data rather than a valid structure.
        if token_buffer:
            raise ValueError("Incomplete row in mmCIF atom_site loop")
        return atoms
    raise ValueError("No mmCIF atom_site loop found")


def _group_residues(
    atoms: list[tuple[str, str, int, str, str, float, float, float]],
) -> tuple[list[str], list[tuple[list[str], np.ndarray]]]:
    residue_types: list[str] = []
    groups: list[tuple[list[str], np.ndarray]] = []

    for residue_type in RESIDUE_TYPES:
        matching_atoms = [atom for atom in atoms if residue_type in atom[0].upper()]
        # Use chain ID, residue number, and insertion code for both PDB and
        # mmCIF so equivalent residue numbers in different chains stay distinct.
        residue_keys = sorted({(atom[1], atom[2], atom[3]) for atom in matching_atoms},
                              key=lambda key: (key[0], key[1], key[2]))
        for chain_id, residue_number, insertion_code in residue_keys:
            residue_atoms = [
                atom for atom in matching_atoms
                if (atom[1], atom[2], atom[3]) == (chain_id, residue_number, insertion_code)
            ]
            atom_names = [atom[4] for atom in residue_atoms]
            coordinates = np.asarray(
                [[atom[5], atom[6], atom[7]] for atom in residue_atoms], dtype=float
            ).reshape((-1, 3))
            residue_types.append(residue_type)
            groups.append((atom_names, coordinates))
    return residue_types, groups

def _mean_omit_nan(values: np.ndarray) -> np.ndarray:
    means = np.empty(values.shape[1], dtype=float)
    for column in range(values.shape[1]):
        valid_values = values[:, column]
        valid_values = valid_values[~np.isnan(valid_values)]
        means[column] = np.mean(valid_values) if len(valid_values) else np.nan
    return means


def _build_transfer_parameter_matrix(residue_types: list[str]) -> np.ndarray:
    count = len(residue_types)
    matrix = np.zeros((count, count), dtype=float)
    for row, donor_type in enumerate(residue_types):
        for column, acceptor_type in enumerate(residue_types):
            matrix[row, column] = TRANSFER_PARAMETERS[
                (donor_type, acceptor_type)
            ]
    return matrix


def _calculate_transfer_rates(
    distances_matrix_m: np.ndarray,
    transfer_parameter_matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    det_coupling_extrapolated_ev = (3.22e-2 / (1.4 * 1.4)) * np.exp(
        -(1e10 * 2.5579) * (distances_matrix_m - 4.15e-10)
    )
    det_coupling_extrapolated_cm_minus_1 = (
        det_coupling_extrapolated_ev * (1.602e-19) * (5.03 * 10**22)
    )
    det_coupling_extrapolated_cm_minus_1 = (
        det_coupling_extrapolated_cm_minus_1
        - np.diag(np.diag(det_coupling_extrapolated_cm_minus_1))
    )

    hbar_cm_minus_1 = 5.29e-12
    speed_of_light_cm_per_second = 3e10
    det_rates_per_second = (
        1.0 / (speed_of_light_cm_per_second * hbar_cm_minus_1**2)
    ) * (
        transfer_parameter_matrix * det_coupling_extrapolated_cm_minus_1**2
    )

    cutoff_distance_m = math.inf
    distances_matrix_cutoff_m = distances_matrix_m.copy()
    distances_matrix_cutoff_m[
        distances_matrix_cutoff_m > cutoff_distance_m
    ] = np.nan
    det_rates_per_second_cutoff = det_rates_per_second.copy()
    det_rates_per_second_cutoff[np.isnan(distances_matrix_cutoff_m)] = np.nan

    det_rates_per_microsecond = (10**-6) * det_rates_per_second_cutoff
    det_rates_per_microsecond[
        np.isinf(det_rates_per_microsecond)
        | np.isnan(det_rates_per_microsecond)
    ] = 0.0
    return det_rates_per_microsecond, det_coupling_extrapolated_cm_minus_1


def _simulate_diffusion_length(
    distances_matrix_angstrom: np.ndarray,
    det_rates_per_microsecond: np.ndarray,
    residue_types: list[str],
    *,
    random_seed: int | None = None,
) -> float:
    lifetime_chromophores_microseconds = np.asarray(
        [LIFETIMES_MICROSECONDS[residue_type] for residue_type in residue_types],
        dtype=float,
    )
    radiative_rates_per_microsecond = 1.0 / lifetime_chromophores_microseconds
    rates_per_microsecond = np.hstack(
        (
            det_rates_per_microsecond,
            radiative_rates_per_microsecond.reshape((-1, 1)),
        )
    )

    rng = np.random.default_rng(random_seed)
    distances_travelled_as_function_of_time = np.zeros(
        (EXCITATION_COUNT, MAX_TIME_STEPS + 1), dtype=float
    )
    where_it_moved_all_values = np.zeros(
        (EXCITATION_COUNT, MAX_TIME_STEPS + 1), dtype=int
    )
    starting_molecules = np.zeros(EXCITATION_COUNT, dtype=int)

    for excitation_index in range(EXCITATION_COUNT):
        starting_molecule = int(rng.integers(len(residue_types)))
        starting_molecules[excitation_index] = starting_molecule
        where_it_moved_all_values[excitation_index, 0] = starting_molecule + 1
        where_it_moved_location = starting_molecule

        for time_index in range(MAX_TIME_STEPS):
            weights = rates_per_microsecond[where_it_moved_location]
            total_weight = float(np.sum(weights))
            if total_weight <= 0 or not math.isfinite(total_weight):
                raise ValueError(
                    "Transition weights must have a positive finite sum."
                )
            threshold = rng.random() * total_weight
            acceptor_index = int(
                np.searchsorted(np.cumsum(weights), threshold, side="left")
            )

            if acceptor_index == len(weights) - 1:
                break

            where_it_moved_location = acceptor_index
            where_it_moved_all_values[
                excitation_index, time_index + 1
            ] = acceptor_index + 1
            distances_travelled_as_function_of_time[
                excitation_index, time_index + 1
            ] = distances_matrix_angstrom[
                starting_molecule, where_it_moved_location
            ]

    where_it_moved_all_values_last_value = np.max(
        np.where(
            where_it_moved_all_values != 0,
            np.arange(1, where_it_moved_all_values.shape[1] + 1),
            0,
        ),
        axis=1,
    )
    ending_molecules = np.empty(EXCITATION_COUNT, dtype=int)
    for excitation_index, last_column in enumerate(
        where_it_moved_all_values_last_value
    ):
        ending_molecules[excitation_index] = (
            where_it_moved_all_values[excitation_index, last_column - 1] - 1
        )

    distances_travelled_by_each_excitation = np.zeros(EXCITATION_COUNT, dtype=float)
    for excitation_index in range(EXCITATION_COUNT):
        distances_travelled_by_each_excitation[excitation_index] = (
            distances_matrix_angstrom[
                starting_molecules[excitation_index],
                ending_molecules[excitation_index],
            ]
        )

    return float(np.mean(distances_travelled_by_each_excitation))


def process_single_pdb(
    pdb_filename: str | Path,
    # microtubule_diameter_angstrom: float,
    *,
    random_seed: int | None = None,
) -> dict[str, float | int | None]:
    """Return the DET-model displacement and coherent-pair statistics."""
    result: dict[str, float | int | None] = {
        "diffusion_length_angstrom": float("nan"),
        "total_pairs": float("nan"),
        "coherent_pairs": float("nan"),
        "coherent_pair_fraction": None,
    }
    try:
        # microtubule_diameter = float(microtubule_diameter_angstrom)
        # microtubule_radius = microtubule_diameter / 2.0
        # _ = microtubule_radius

        atoms = _read_pdb_atoms(pdb_filename)
        residue_types, groups = _group_residues(atoms)
        number_of_molecules = len(residue_types)
        if number_of_molecules < 2:
            result.update(
                diffusion_length_angstrom=0.0,
                total_pairs=number_of_molecules * (number_of_molecules - 1) // 2,
                coherent_pairs=0,
            )
            return result

        positions_xyz = np.asarray(
            [_mean_omit_nan(coordinates) for _, coordinates in groups],
            dtype=float,
        )
        distances_matrix_angstrom = cdist(
            positions_xyz, positions_xyz, metric="euclidean"
        )
        distances_matrix_m = distances_matrix_angstrom * 10**-10

        transfer_parameter_matrix = _build_transfer_parameter_matrix(residue_types)
        (
            det_rates_per_microsecond,
            det_coupling_extrapolated_cm_minus_1,
        ) = _calculate_transfer_rates(
            distances_matrix_m,
            transfer_parameter_matrix,
        )
        total_pairs = number_of_molecules * (number_of_molecules - 1) // 2
        lower_triangle = np.tril(
            np.ones(det_coupling_extrapolated_cm_minus_1.shape, dtype=bool), k=-1
        )
        coherent_pairs = int(
            np.count_nonzero(
                np.abs(det_coupling_extrapolated_cm_minus_1[lower_triangle]) > 53.052
            )
        )
        result.update(
            total_pairs=total_pairs,
            coherent_pairs=coherent_pairs,
            coherent_pair_fraction=coherent_pairs / total_pairs,
        )
        total_hops = int(np.sum(~np.isnan(det_rates_per_microsecond)))
        _ = total_hops

        result["diffusion_length_angstrom"] = _simulate_diffusion_length(
            distances_matrix_angstrom,
            det_rates_per_microsecond,
            residue_types,
            random_seed=random_seed,
        )
        return result
    except Exception as error:
        print(f"Error processing {pdb_filename}: {error}", file=sys.stderr)
        return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate the DET diffusion length for one PDB or mmCIF structure."
    )
    parser.add_argument("pdb_filename", help="Input PDB (.pdb) or mmCIF (.cif/.mmcif) file")
    # parser.add_argument(
    #     "microtubule_diameter_angstrom",
    #     type=float,
    #     help="Microtubule diameter in Angstroms (unused by the MATLAB algorithm)",
    # )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional Python RNG seed for repeatable Monte Carlo runs",
    )
    args = parser.parse_args()
    result = process_single_pdb(
        args.pdb_filename,
        # args.microtubule_diameter_angstrom,
        random_seed=args.seed,
    )
    print(result)


if __name__ == "__main__":
    main()
