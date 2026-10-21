
# Plexcal
## Excitation Energy-Transfer Diffusion Length Calculator

Plexcal is a structure-aware scientific application for estimating Photoexcitation energy-transfer diffusion lengths using implemented Förster resonance energy transfer (FRET) and Dexter energy transfer (DET) models.

The application supports uploaded Protein Data Bank (PDB) files, structures identified by RCSB PDB IDs, and AlphaFold structure IDs. Molecular structures can be explored interactively using the PDBe Mol* viewer.

## Features

- **Structure visualization:** Explore molecular structures with rotation, sequence, component, representation, selection, and measurement controls in the PDBe Mol* viewer.
- **Structure inputs:** Upload a PDB file or load a structure using an RCSB PDB ID or AlphaFold structure ID.
- **FRET and DET estimation:** Calculate excitation energy-transfer diffusion lengths using the implemented models.
- **Chromophore visualization:** Optionally highlight tryptophan (TRP), tyrosine (TYR), and phenylalanine (PHE) residues.
- **Results and comparison:** Inspect diffusion lengths, total pairs, coherent pairs, coherent-pair fractions, and comparison charts.
- **Structure metadata:** Retrieve available structural, organism, experimental, and citation information.
- **Batch processing:** Analyze multiple PDB structures and export results to CSV.

## Installation and Setup

### Requirements

- Python version: [Insert the tested Python version]
- Git
- A compatible operating system and the dependencies listed in `requirements.txt`

### Install

Clone the repository and install the dependencies:

```bash
git clone https://github.com/KalraLab/PLEXCAL.git
cd PLEXCAL

python -m venv .venv
```

Activate the environment.

On Linux or macOS:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the dependencies and start the application:

```bash
python -m pip install -r requirements.txt
python server.py
```

Open the local URL displayed in the terminal, usually:

`http://127.0.0.1:8000`

## Batch Processing

To process structures stored in `all_proteins_pdb/`, run:

```bash
python -m src.process_all_pdbs
```

The script uses up to four worker processes by default and saves results to:

`all_proteins_pdb/fret_det_diffusion_lengths.csv`

To specify the number of workers:

```bash
python -m src.process_all_pdbs --workers 2
```

For the standalone DET implementation, run:

```bash
python -m src.calcalculate_det_diffusion_lengths "1TUB.pdb"
```

Check the corresponding script's help output or documentation for any additional arguments.

## Methods and Parameter Provenance

Methods, equations, parameter sources, implementation details, and known limitations are documented in:

[`docs/Plexcal_Methods_Theory_Parameter_Provenance.md`](docs\Plexcal_Theory_Methods_Parameter_Provenance.md)

The uploaded-file FRET implementation is a Python translation of `process_single_pdb_v1.m`. It uses 10,000 excitation trajectories with up to 1,000 transitions per trajectory. Results may differ from MATLAB because the implementations use different random-number streams.

Diffusion lengths are model-derived estimates that depend on the implemented assumptions, parameters, and input structure. They should not be interpreted as experimentally measured transport lengths or as evidence of experimentally validated transfer rates without independent validation.

## Limitations

- Predictions depend on the input structure and the assumptions and parameters used by each calculation model.
- Results from the Python FRET implementation may differ from the corresponding MATLAB implementation.
- Structure metadata availability depends on the relevant external services and structure records.
- Structure uploads and downloads are limited to 25 MiB.
- Uploaded structures must be selected again after refreshing the browser. Completed results remain available in the current browser tab.

Consult the methods documentation before interpreting results or comparing structures.

## Citation

If you use Plexcal in research, please cite the specific software release used in your analysis.

**Software:** Plexcal: Structure-Aware Photoexcitation Energy-Transfer Diffusion lengths claculator in Protein system
**Version:** [Release version]  
**DOI:** [Version-specific DOI]  
**Repository:** https://github.com/KalraLab/PLEXCAL

See `CITATION.cff` for machine-readable citation metadata.

## License

Plexcal is distributed under the terms of the license provided in the [`LICENSE`](LICENSE) file. Consult that file for the permissions and restrictions that apply to use, modification, and redistribution.

## Version and Release History

See [GitHub Releases](https://github.com/KalraLab/PLEXCAL/releases) for version history, release notes, and archived source snapshots.
