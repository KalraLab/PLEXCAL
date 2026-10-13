# PLEXCAL

# Plexcal: Excitation Energy-Transfer Diffusion Length Calculator

Plexcal is a structure-aware scientific application for estimating excitation energy-transfer diffusion lengths using existing Förster resonance energy transfer (FRET) and Dexter energy transfer (DET) models. It supports uploaded PDB files, RCSB PDB IDs, and AlphaFold structure IDs, with interactive molecular visualization through the PDBe Mol* viewer.

## Installation and Setup

Clone the repository, install the dependencies, and start the server:

```bash
git clone <repository-url>
cd <repository-folder>
python -m pip install -r requirements.txt
python server.py
```

Open the URL displayed in the terminal, usually:

`http://127.0.0.1:8000`

## Features

- **Interactive structure visualization:** Explore molecular structures using the PDBe Mol* viewer, including rotation, sequence, component, representation, selection, and measurement controls.
- **Structure inputs:** Upload a PDB file or load a structure using an RCSB PDB ID or AlphaFold structure ID.
- **FRET and DET calculations:** Estimate excitation energy-transfer diffusion lengths using the existing calculation models.
- **Chromophore visualization:** Optionally highlight tryptophan (TRP), tyrosine (TYR), and phenylalanine (PHE) residues.
- **Results and comparison:** View diffusion lengths, total pairs, coherent pairs, coherent-pair fractions, and comparison charts.
- **Structure metadata:** Retrieve available structural, organism, experimental, and citation information.
- **Batch processing:** Analyze multiple PDB structures and export results to CSV.

## Batch Processing

To process structures in `all_proteins_pdb/`, run:

```bash
python process_all_pdbs.py
```

The script uses up to four worker processes by default and saves results to `all_proteins_pdb/diffusion_lengths.csv`. To specify the number of workers:

```bash
python process_all_pdbs.py --workers 2
```

## Additional Documentation

- `PROCESS_SINGLE_PDB.md` — FRET workflow, equations, physical model, and limitations.
- `DET_CODE_CLEAN.md` — DET workflow and equations.

The uploaded-file FRET implementation is a Python translation of `process_single_pdb_v1.m`. It uses 10,000 excitation trajectories with up to 1,000 transitions per trajectory. Results may differ from MATLAB because the implementations use different random-number streams.

The standalone DET implementation is available in `det_code_clean.py`:

```bash
python det_code_clean.py "1TUB.pdb" 0
```

**Note:** Structure uploads and downloads are limited to 25 MiB. Uploaded structures must be selected again after refreshing the browser, while completed results remain available in the current browser tab.