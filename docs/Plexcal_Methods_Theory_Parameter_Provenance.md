# Plexcal: Theory, Methods and Parameter Provenance

**Dhanumjaya Reddy Bhavanam¹, Kushal Singh Rawat¹, Aarat P. Kalra¹˒²˒³***

¹ Centre for Biomedical Engineering, Indian Institute of Technology Delhi, New Delhi 110016, India

² Department of Biomedical Engineering, All India Institute of Medical Sciences, New Delhi, India

³ Amar Nath and Shashi Khosla School of Information Technology, Indian Institute of Technology Delhi, New Delhi 110016, India

*Corresponding author:* Aarat P. Kalra

**Software version:** 0.1.0  
**Release date:** October 2026  
**Document:** Software Methods and Technical Documentation

## Abstract

Plexcal is a software tool that estimates how far electronic excitation can migrate through the aromatic amino acids of a protein structure. From a structure in PDB format it locates tryptophan, tyrosine and phenylalanine residues, calculates pairwise electronic couplings and transfer rates for singlet (Förster-type) and triplet (Dexter-type) transfer, and runs 10,000 incoherent kinetic Monte Carlo trajectories in which an excitation hops between residues until it decays. The reported diffusion length is the mean straight-line distance between the residue where an excitation starts and the residue where it ends. The model follows the kinetic-hopping framework developed by Kalra and co-workers for microtubules and other cytoskeletal polymers. This document sets out the theory, the equations and constants as implemented, their relation to the published literature, and the limits on how the outputs should be interpreted. The outputs are model estimates, not measurements.

## 1. Introduction

Tryptophan, tyrosine and phenylalanine absorb ultraviolet light, and an excitation created on one of these side chains can transfer to a neighbouring one. Transfer of this kind between aromatic residues has been recognised in proteins since the 1950s \[1,2\]. When aromatic residues are arranged in a regular polymeric lattice, as in microtubules, repeated transfer events can carry an excitation farther than the extent of a single protein subunit \[4\]. Whether a particular structure supports such migration depends on the separations, orientations and photophysical properties of its aromatic residues, and the first two can be read from atomic coordinates.

Plexcal applies this approach to any protein structure supplied as coordinates. A structure can be uploaded as a PDB file or retrieved from the RCSB Protein Data Bank or the AlphaFold Protein Structure Database, and it is displayed with the Mol\* viewer. Two transfer mechanisms are implemented: singlet transfer in the point-dipole approximation, labelled FRET, and triplet transfer through exchange coupling, labelled DET. For each mechanism Plexcal reports a modelled diffusion length, the number of residue pairs whose coupling exceeds a fixed threshold, and the fraction of all pairs that this number represents. The tool is intended for comparing many structures under one consistent model. It is not a quantum-chemical calculation for any individual structure.

This document describes release 0.1.0. Equations and constants are given as they are implemented. Where a corresponding published value differs from the implemented one, the difference is stated in Section 5.

## 2. Background

### 2.1 Energy transfer between aromatic residues

The aromatic amino acids absorb near 280 nm (tryptophan), 275 nm (tyrosine) and 258 nm (phenylalanine), and the emission of one residue type overlaps the absorption of the next in the order Phe, Tyr, Trp \[3\]. Resonance transfer of excitation among aromatic residues was proposed on this basis by Karreman, Steele and Szent-Györgyi and supported by the depolarisation measurements of Weber \[1,2\]. The theory of resonance transfer in the dipole–dipole limit was given by Förster \[10\] and is reviewed by Scholes \[12\]. Triplet-state transfer, which requires orbital overlap between donor and acceptor, was described by Dexter \[11\].

### 2.2 Energy migration in cytoskeletal polymers

Kalra and co-workers measured tryptophan fluorescence lifetimes of tubulin and microtubules in the presence of an external acceptor, 7-amino-4-methylcoumarin-3-acetic acid (AMCA), and applied Stern–Volmer analysis to obtain a two-dimensional photoexcitation diffusion length of 6.6 ± 0.1 nm in microtubules \[4\]. They also simulated incoherent energy migration by kinetic Monte Carlo on an all-atom model of a 31-dimer microtubule that included tyrosine and tryptophan residues. Förster-type rates gave diffusion lengths of 1.11 nm for tyrosine–tyrosine transfer alone and 1.54 nm when tyrosine–tryptophan transfer was included, well below the experimental value. The authors suggested that the point-dipole approximation may be inaccurate at the shortest separations and that exchange or charge-transfer contributions could raise the coupling; charge-transfer integrals for five close tyrosine dimers reached 419 cm⁻¹ \[4\]. A later Perspective reviews the photophysics of microtubules and actin filaments and the open questions in protein biophotonics \[5\].

Kakati, Karmakar and Kalra extended the approach to triplet states. They computed Dexter transfer rates in microtubules, actin filaments and vimentin filaments and found that interaromatic couplings are mostly weaker than 10⁻³ eV, which limits triplet diffusion lengths to 6.1, 0.5 and 1.0 Å, respectively. Repeated FRET gave singlet diffusion lengths of 12.4 Å in actin filaments and about 8.6 Å in microtubules and vimentin, so singlet migration dominates in these polymers \[6\]. Varshney and Kalra, in a preprint, used kinetic Monte Carlo on idealised one-, two- and three-dimensional lattices of aromatic residues and reported a maximum singlet diffusion length of about 14 nm, reached in tyrosine face-centred cubic lattices at optimal spacing and orientation in a hydrophobic environment \[7\]. Earlier theoretical treatments of excitonic transport in tubulin aromatic networks are given in \[8,9\].

### 2.3 Position of Plexcal

These studies used purpose-built structures and parameter sets for particular polymers. Plexcal adopts the same model class, incoherent hopping with rates obtained from a golden-rule coupling expression and a competing decay channel \[4,6\], and applies it to any structure with a fixed parameter set. Values from Plexcal are therefore comparable between structures processed with the same release. They are not expected to reproduce the published numbers for a given polymer exactly, because the published calculations used different structure preparation (for example, molecular-dynamics ensembles and quantum-chemically computed transition dipoles in \[4\]) and different parameters.

## 3. Theory

### 3.1 Incoherent transfer rates

When the coupling between two chromophores is small compared with the energetic broadening of their transitions, the excitation stays localised on one site and transfer is a rate process. The rate from donor i to acceptor j follows the Fermi golden rule \[4,12,13\]:

$$
k_{ij}=\frac{2\pi}{\hbar}\,\lvert V_{ij}\rvert^{2}\,J_{ij}\qquad(1)
$$

where Vᵢⱼ is the electronic coupling and Jᵢⱼ is the spectral overlap integral of the area-normalised donor emission and acceptor absorption spectra on an energy scale \[4\]. Writing the coupling in wavenumbers, Ṽ = V/(hc), and the overlap on a wavenumber scale, J̃ (in cm), Eq. 1 becomes

$$
k_{ij}=4\pi^{2}c\,\tilde V_{ij}^{2}\,\tilde J_{ij}\qquad(2)
$$

with c in cm s⁻¹, which gives kᵢⱼ in s⁻¹. This is the form used in the implementation (Section 4).

### 3.2 Singlet coupling: point-dipole approximation

For singlet transfer, the coupling is taken as the Coulomb interaction between transition dipoles, screened by the medium \[4,14\]:

$$
V_{ij}=\frac{1}{4\pi\varepsilon_{0}n^{2}}\,\frac{\kappa_{ij}\,\mu_{i}\mu_{j}}{R_{ij}^{3}}\qquad(3)
$$

$$
\kappa_{ij}=\hat{\mathbf u}_{i}\cdot\hat{\mathbf u}_{j}-3\,(\hat{\mathbf u}_{i}\cdot\hat{\mathbf R}_{ij})(\hat{\mathbf u}_{j}\cdot\hat{\mathbf R}_{ij})\qquad(4)
$$

Here μᵢ and μⱼ are the transition-dipole magnitudes, û are unit vectors along the transition dipoles, Rᵢⱼ is the separation, R̂ᵢⱼ is the unit vector joining the sites, ε₀ is the vacuum permittivity and n is the refractive index of the surrounding medium, taken as 1.4 in \[4\]. Because k scales as V², the rate varies as R⁻⁶ for fixed orientation, which is the familiar Förster distance dependence \[10,12\]. In terms of a Förster radius R₀ this is k = (1/τ\_D)(R₀/R)⁶, where τ\_D is the donor lifetime.

### 3.3 Triplet coupling: exchange interaction

Triplet transfer is mediated by exchange and requires overlap of donor and acceptor orbitals. Its coupling therefore falls off much faster with distance than the dipole–dipole coupling and is commonly represented as an exponential \[11,15\]:

$$
V_{ij}=V_{0}\exp\left[-\beta\,(R_{ij}-R_{0})\right]\qquad(5)
$$

The rate is again given by Eq. 1 or 2, with an overlap factor appropriate to the triplet transition in place of Jᵢⱼ \[6\].

### 3.4 Decay and event probabilities

Each site also has a decay channel with rate k\_dec,i = 1/τᵢ, where τᵢ is the excited-state lifetime assigned to the residue type. For an excitation at site i, the probability that the next event is a transfer to site j is

$$
p_{ij}=\frac{k_{ij}}{\sum_{m\neq i}k_{im}+k_{\mathrm{dec},i}}\qquad(6)
$$

and the probability that it decays is k\_dec,i divided by the same sum. Decay ends the trajectory.

### 3.5 Kinetic Monte Carlo and the diffusion length

In a kinetic Monte Carlo scheme the next event is chosen at random with probabilities proportional to the rates, and the event is executed \[4,16,17\]. Because the competing processes are memoryless, the sequence of sites visited, and therefore the final site, has the distribution given by Eq. 6 whether or not waiting times are sampled. Plexcal does not sample waiting times and does not report time-resolved quantities.

For N trajectories, with r₀,i the position of the starting site and r\_f,i the position of the site at which trajectory i ends, the diffusion length reported by Plexcal is

$$
L_{D}=\frac{1}{N}\sum_{i=1}^{N}\left\lVert \mathbf r_{f,i}-\mathbf r_{0,i}\right\rVert\qquad(7)
$$

This is a mean endpoint displacement. It is not the accumulated path length, and a trajectory that wanders away and returns contributes little. It is also a different quantity from the diffusion length L = √(4Dτ) that \[4\] obtains from a diffusion coefficient D and a lifetime τ by Stern–Volmer analysis, although both express how far an excitation travels before it is lost \[18\].

### 3.6 Coupling-threshold descriptor

For each unordered pair of sites Plexcal records whether the magnitude of the coupling exceeds 53.052 cm⁻¹. With nₛ sites there are nₛ(nₛ − 1)/2 pairs. The reported pair count is the number of pairs above the threshold, and the pair fraction is that number divided by the total. The quantity describes how many pairs are strongly coupled. The kinetic Monte Carlo model is incoherent and does not calculate coherence, delocalised excitonic states or superradiance, so a pair above the threshold is not evidence that the pair is quantum coherent. The software and the web interface use the labels "coherent pairs" and "coherent pair fraction" for these two outputs; they should be read as pairs above the coupling threshold.

## 4. Implementation

### 4.1 Input structures

Plexcal reads fixed-column PDB records. Structures can be uploaded, or fetched by identifier from the RCSB PDB \[19,20\] or the AlphaFold Protein Structure Database \[21,22\]. Only the first model of a multi-model file is read, and reading stops at the first ENDMDL or END record. ATOM and HETATM records are both used. A PDB file contains coordinates only; it does not carry excited-state energies, transition densities or spectra, which is why these enter as fixed parameters.

### 4.2 Sites and distances

Residues named TRP, TYR and PHE are selected. Residues are identified by name and residue number. The position of each site is the arithmetic mean of the coordinates of all atoms of that residue present in the file, and the separation between sites is the Euclidean distance

$$
R_{ij}=\lVert\mathbf r_{j}-\mathbf r_{i}\rVert\qquad(8)
$$

A structure with fewer than two aromatic residues returns a diffusion length of zero.

### 4.3 FRET module

The orientation vector of each residue is the normal to its aromatic ring. The ring atoms used are CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3 and CH2 for tryptophan, and CG, CD1, CD2, CE1, CE2 and CZ for tyrosine and phenylalanine. The normal is the singular vector of the mean-centred ring coordinates associated with the smallest singular value, and its sign is chosen to point from the ring centroid towards the residue site position. The sign convention does not affect the rate, which depends on κ² through V². Setting û = n̂ (the ring normal) in Eq. 4 gives the orientation factor κᵢⱼ.

The vacuum dipole–dipole coupling is κᵢⱼ μᵢ μⱼ / (4π ε₀ Rᵢⱼ³), with Rᵢⱼ in metres and μ in C m, and it is converted to wavenumbers with the factor 5.03 × 10²² cm⁻¹ J⁻¹ (approximately 1/hc). The screening factor is applied in the rate, which is

$$
k^{\mathrm{FRET}}_{ij}=\frac{1}{c\,\tilde\hbar^{2}}\,\frac{1}{n^{4}}\,J_{ij}\,\bigl(\tilde V^{\mathrm{vac}}_{ij}\bigr)^{2}\qquad(9)
$$

with c = 3.0 × 10¹⁰ cm s⁻¹, ħ̃ = 5.29 × 10⁻¹² cm⁻¹ s and n = 1.4. The prefactor 1/(cħ̃²) equals 4π²c when ħ̃ = 1/(2πc), and n⁻⁴ Ṽ\_vac² equals (Ṽ\_vac/n²)², so Eq. 9 is Eq. 2 with the coupling of Eq. 3. The rounded value of ħ̃ makes the numerical prefactor 0.6% larger than 4π²c. Rates are converted to ns⁻¹.

### 4.4 DET module

The coupling is evaluated from Eq. 5 with the separation in metres:

$$
V_{ij}\,(\mathrm{eV})=\frac{3.22\times10^{-2}}{n^{2}}\exp\left[-\beta\,(R_{ij}-R_{0})\right],\quad \beta=2.5579\ \text{\AA}^{-1},\ R_{0}=4.15\ \text{\AA},\ n=1.4\qquad(10)
$$

The result is converted to cm⁻¹ with 1.602 × 10⁻¹⁹ J eV⁻¹ and 5.03 × 10²² cm⁻¹ J⁻¹, and the diagonal is set to zero. For this parameter set the coupling is about 1.6 × 10⁻² eV at R₀ and falls below 10⁻³ eV at a separation of about 5.2 Å, consistent in magnitude with the weak couplings described in \[6\]. The rate is

$$
k^{\mathrm{DET}}_{ij}=\frac{1}{c\,\tilde\hbar^{2}}\,P_{ij}\,\tilde V_{ij}^{2}\qquad(11)
$$

where Pᵢⱼ is a donor–acceptor overlap parameter that takes the place of Jᵢⱼ in Eq. 2. Rates are converted to μs⁻¹.

### 4.5 Kinetic Monte Carlo procedure

For each structure the following steps are carried out.

1. Build the site list, the distance matrix and the rate matrix, and append a decay column with rates 1/τᵢ.
2. For each of N = 10,000 trajectories, choose a starting site with equal probability among all sites.
3. At the current site, draw a uniform random number scaled by the total of the transfer and decay rates, and select the event whose cumulative rate interval contains it.
4. If the event is decay, stop. Otherwise move to the selected site and repeat from step 3, up to a maximum of 1,000 events per trajectory.
5. Record the distance between the starting site and the final site, and average over trajectories (Eq. 7).

No distance cutoff is applied to the rates. Random numbers are drawn with NumPy's default generator, and the seed can be set for repeatable runs. A trajectory that reaches the event limit without decaying ends at its current site.

### 4.6 Outputs

| Output | Meaning |
| --- | --- |
| FRET diffusion length (Å) | L\_D from Eq. 7 with the singlet model |
| DET diffusion length (Å) | L\_D from Eq. 7 with the triplet model |
| Total pairs | nₛ(nₛ − 1)/2 |
| Pairs above threshold | Pairs with coupling magnitude above 53.052 cm⁻¹ |
| Pair fraction | Pairs above threshold divided by total pairs |

## 5. Parameters and provenance

All constants below are fixed in release 0.1.0. The literature column states what the published sources define or report; it does not imply that the implemented value was taken from that source.

### 5.1 Common constants

| Quantity | Value | Relation to the literature |
| --- | --- | --- |
| Refractive index n | 1.4 | Same value as assumed in \[4\] |
| Trajectories N | 10,000 | \[4\] reports 3,000 trajectories for its simulations |
| Maximum events per trajectory | 1,000 | Implementation setting |
| c | 3.0 × 10¹⁰ cm s⁻¹ | Physical constant, rounded |
| ħ̃ | 5.29 × 10⁻¹² cm⁻¹ s | Physical constant, rounded (Section 4.3) |
| Wavenumber conversion | 5.03 × 10²² cm⁻¹ J⁻¹ | Approximately 1/hc |
| ε₀ | 8.854 × 10⁻¹² F m⁻¹ | Physical constant |
| Debye | 3.33564 × 10⁻³⁰ C m | Unit conversion |
| Pair threshold | 53.052 cm⁻¹ | Fixed constant of the release (Section 3.6) |

### 5.2 FRET parameters

Transition-dipole magnitudes: tryptophan 2.074 D, tyrosine 1.18 D, phenylalanine 0.28 D. Reference \[4\] used transition dipoles computed at the EOM-CCSD/6-311+G(d) level, 2.77 D for tryptophan and 1.356 D for tyrosine, so the values are not the same.

Singlet lifetimes: tryptophan 3.6 ns, tyrosine 7.5 ns, phenylalanine 3.39 ns. Reference \[4\] used a single dissipation rate of 3.06 × 10⁸ s⁻¹ (a lifetime of about 3.3 ns) for the residues it simulated.

Spectral overlap integrals Jᵢⱼ (rows are donors, columns are acceptors):

| Donor \\ Acceptor | TRP | TYR | PHE |
| --- | --- | --- | --- |
| TRP | 5.5953 × 10⁻⁶ | 2.3247 × 10⁻⁵ | 1.3034 × 10⁻⁷ |
| TYR | 9.2550 × 10⁻⁴ | 4.7933 × 10⁻⁴ | 1.8891 × 10⁻⁵ |
| PHE | 3.9 × 10⁻³ | 4.3 × 10⁻³ | 6.1123 × 10⁻⁴ |

The entries are used with Eq. 9 so that the rate is in s⁻¹, which corresponds to the overlap on a wavenumber scale in cm as defined in Eq. 2. Overlap integrals depend on the spectra, the normalisation and the spectral window used to compute them. Table S11 of \[4\] lists J for tryptophan and tyrosine in units of 10⁻⁶ cm: Trp→Trp 7.53, Tyr→Tyr 16.15, Trp→Tyr 5.63 and Tyr→Trp 30.52. The Plexcal values are not identical to these and should not be assumed to reproduce them.

### 5.3 DET parameters

Triplet lifetimes: tryptophan 7.3 μs, tyrosine 5.6 μs, phenylalanine 3.1 μs. Coupling constants are given in Eq. 10. The overlap parameters Pᵢⱼ (rows are donors, columns are acceptors) are:

| Donor \\ Acceptor | TRP | TYR | PHE |
| --- | --- | --- | --- |
| TRP | 0.00058 | 0.000487 | 0.00031 |
| TYR | 0.0011 | 0.0012 | 0.00098 |
| PHE | 0.0039 | 0.0012 | 0.00098 |

The Supporting Information of \[6\] tabulates Dexter and FRET spectral overlap integrals, molecular-orbital overlap integrals and Dexter couplings for aromatic pairs in tubulin and actin. Quantities of the same kind are used here, but the correspondence between the tabulated entries and the values above has not been checked entry by entry, and the units of Pᵢⱼ are those required for Eq. 11 to return a rate in s⁻¹ with Ṽ in cm⁻¹.

## 6. Interpretation and limitations

1. **Model estimates.** The diffusion lengths are outputs of a simplified model with fixed parameters. They are not measured distances and do not establish that excitation migrates in a given protein.
2. **Incoherent transport only.** Delocalisation, dephasing and superradiance are not treated. The pair threshold in Section 3.6 is a fixed numerical criterion, and a derivation of its value is not part of this release.
3. **Single static structure.** Each calculation uses one set of coordinates. Thermal fluctuations are not sampled, whereas the microtubule simulations in \[4\] used structures drawn from molecular-dynamics ensembles.
4. **Site definition.** The site position is the mean of all atoms of the residue, not the centre of the aromatic ring. The two points differ by a distance that is not negligible compared with the nearest-neighbour separations between aromatic residues.
5. **Chain handling.** Residues are grouped by name and number only. In a structure with several chains, residues of the same type and number in different chains are merged into one site, and alternate locations and insertion codes are not distinguished. Results for multi-chain assemblies should be interpreted with this in mind.
6. **Orientation.** The ring normal is used as the orientation vector in Eq. 4. The π→π\* transition dipoles of the aromatic side chains lie approximately in the ring plane, so the normal is a geometric descriptor of the ring and not the direction of a transition dipole. The factor κ obtained here is therefore a proxy and not the orientation factor of the actual transitions, which \[4\] obtained from quantum-chemical transition dipoles.
7. **Point-dipole approximation.** Eq. 3 is not accurate when the separation is comparable to the size of the chromophores. The Förster-only simulations in \[4\] gave diffusion lengths of 1.11 and 1.54 nm against 6.6 nm measured, and the authors attributed part of the difference to short-range interactions. The FRET module of Plexcal does not include such contributions.
8. **Medium.** A single refractive index and a single set of lifetimes and overlap integrals are used for every residue and structure. Local dielectric heterogeneity and energetic disorder are not represented.
9. **Screening in the pair count.** In the FRET module the screening factor n⁻² is applied through Eq. 9 and not to the coupling used for the threshold count, whereas in the DET module it is part of the coupling in Eq. 10. The pair counts of the two modules are therefore not defined on the same footing.
10. **Predicted structures.** For AlphaFold models the side-chain coordinates are predictions. Plexcal does not weight residues by model confidence.

## 7. Reproducibility and availability

The calculation is stochastic. With N = 10,000 the standard error of L\_D is s/√N, where s is the standard deviation of the displacement distribution. Plexcal reports the mean only, so the spread between runs can be estimated by repeating a calculation with different seeds. To make a result reproducible, record the structure identifier or file checksum, the source and retrieval date, the Plexcal version and the random seed.

The reference implementations are two Python modules, `calculate_fret_diffusion_lengths.py` and `calculate_det_diffusion_lengths.py`. They require Python 3 with NumPy \[24\] and SciPy \[25\] and are run from the command line, for example `python calculate_fret_diffusion_lengths.py structure.pdb --seed 42`. Each returns the diffusion length in ångströms, the total number of pairs, the number of pairs above the threshold and their fraction.

This document is archived together with the software release it describes and should be cited using the DOI assigned on deposition. Calculations that depend on the numerical values in Section 5 should state the release number.

## References

1. Karreman, G.; Steele, R. H.; Szent-Györgyi, A. On Resonance Transfer of Excitation Energy between Aromatic Amino Acids in Proteins. *Proc. Natl. Acad. Sci. U.S.A.* **1958**, *44*, 140–143. https://doi.org/10.1073/pnas.44.2.140
2. Weber, G. Fluorescence-Polarization Spectrum and Electronic-Energy Transfer in Proteins. *Biochem. J.* **1960**, *75*, 345–352. https://doi.org/10.1042/bj0750345
3. Lakowicz, J. R. *Principles of Fluorescence Spectroscopy*, 3rd ed.; Springer: New York, 2006. https://doi.org/10.1007/978-0-387-46312-4
4. Kalra, A. P.; Benny, A.; Travis, S. M.; Zizzi, E. A.; Morales-Sanchez, A.; Oblinsky, D. G.; Craddock, T. J. A.; Hameroff, S. R.; MacIver, M. B.; Tuszyński, J. A.; Petry, S.; Penrose, R.; Scholes, G. D. Electronic Energy Migration in Microtubules. *ACS Cent. Sci.* **2023**, *9*, 352–361. https://doi.org/10.1021/acscentsci.2c01114
5. Kalra, A. P.; Biswas, S.; Mulrain, I.; Wang, M.; Tuszynski, J. A.; Scholes, G. D. All Lit Up: Exploring the Photophysical Properties of Protein Polymers. *J. Phys. Chem. Lett.* **2023**, *14*, 5891–5900. https://doi.org/10.1021/acs.jpclett.3c01032
6. Kakati, A.; Karmakar, T.; Kalra, A. P. Triplet Energy Migration in Cytoskeletal Polymers. *J. Phys. Chem. B* **2025**, *129*, 128–138. https://doi.org/10.1021/acs.jpcb.4c06748
7. Varshney, A.; Kalra, A. P. Incoherent Electronic Energy Migration in Aromatic Amino Acid Lattices Examined by Kinetic Monte Carlo Simulations. *ChemRxiv* **2025**, preprint, version 1. https://chemrxiv.org/engage/chemrxiv/article-details/692c157fa10c9f5ca1733d83
8. Craddock, T. J. A.; Friesen, D.; Mane, J.; Hameroff, S.; Tuszynski, J. A. The Feasibility of Coherent Energy Transfer in Microtubules. *J. R. Soc. Interface* **2014**, *11*, 20140677. https://doi.org/10.1098/rsif.2014.0677
9. Kurian, P.; Obisesan, T. O.; Craddock, T. J. A. Oxidative Species-Induced Excitonic Transport in Tubulin Aromatic Networks: Potential Implications for Neurodegenerative Disease. *J. Photochem. Photobiol. B* **2017**, *175*, 109–124.
10. Förster, T. Zwischenmolekulare Energiewanderung und Fluoreszenz. *Ann. Phys.* **1948**, *437*, 55–75. https://doi.org/10.1002/andp.19484370105
11. Dexter, D. L. A Theory of Sensitized Luminescence in Solids. *J. Chem. Phys.* **1953**, *21*, 836–850. https://doi.org/10.1063/1.1699044
12. Scholes, G. D. Long-Range Resonance Energy Transfer in Molecular Systems. *Annu. Rev. Phys. Chem.* **2003**, *54*, 57–87. https://doi.org/10.1146/annurev.physchem.54.011002.103746
13. Scholes, G. D.; Fleming, G. R. On the Mechanism of Light Harvesting in Photosynthetic Purple Bacteria: B800 to B850 Energy Transfer. *J. Phys. Chem. B* **2000**, *104*, 1854–1868. https://doi.org/10.1021/jp993435l
14. Krueger, B. P.; Scholes, G. D.; Fleming, G. R. Calculation of Couplings and Energy-Transfer Pathways between the Pigments of LH2 by the ab Initio Transition Density Cube Method. *J. Phys. Chem. B* **1998**, *102*, 5378–5386. https://doi.org/10.1021/jp9811171
15. Skourtis, S. S.; Liu, C.; Antoniou, P.; Virshup, A. M.; Beratan, D. N. Dexter Energy Transfer Pathways. *Proc. Natl. Acad. Sci. U.S.A.* **2016**, *113*, 8115–8120. https://doi.org/10.1073/pnas.1517189113
16. Bortz, A. B.; Kalos, M. H.; Lebowitz, J. L. A New Algorithm for Monte Carlo Simulation of Ising Spin Systems. *J. Comput. Phys.* **1975**, *17*, 10–18. https://doi.org/10.1016/0021-9991(75)90060-1
17. Gillespie, D. T. Exact Stochastic Simulation of Coupled Chemical Reactions. *J. Phys. Chem.* **1977**, *81*, 2340–2361. https://doi.org/10.1021/j100540a008
18. Mikhnenko, O. V.; Blom, P. W. M.; Nguyen, T.-Q. Exciton Diffusion in Organic Semiconductors. *Energy Environ. Sci.* **2015**, *8*, 1867–1888. https://doi.org/10.1039/C5EE00925A
19. Berman, H. M.; Westbrook, J.; Feng, Z.; Gilliland, G.; Bhat, T. N.; Weissig, H.; Shindyalov, I. N.; Bourne, P. E. The Protein Data Bank. *Nucleic Acids Res.* **2000**, *28*, 235–242. https://doi.org/10.1093/nar/28.1.235
20. Burley, S. K.; et al. Updated Resources for Exploring Experimentally-Determined PDB Structures and Computed Structure Models at the RCSB Protein Data Bank. *Nucleic Acids Res.* **2025**, *53*, D564–D574. https://doi.org/10.1093/nar/gkae1091
21. Jumper, J.; et al. Highly Accurate Protein Structure Prediction with AlphaFold. *Nature* **2021**, *596*, 583–589. https://doi.org/10.1038/s41586-021-03819-2
22. Varadi, M.; et al. AlphaFold Protein Structure Database: Massively Expanding the Structural Coverage of Protein-Sequence Space with High-Accuracy Models. *Nucleic Acids Res.* **2022**, *50*, D439–D444. https://doi.org/10.1093/nar/gkab1061
23. Sehnal, D.; et al. Mol\* Viewer: Modern Web App for 3D Visualization and Analysis of Large Biomolecular Structures. *Nucleic Acids Res.* **2021**, *49*, W431–W437. https://doi.org/10.1093/nar/gkab314
24. Harris, C. R.; et al. Array Programming with NumPy. *Nature* **2020**, *585*, 357–362. https://doi.org/10.1038/s41586-020-2649-2
25. Virtanen, P.; et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. *Nat. Methods* **2020**, *17*, 261–272. https://doi.org/10.1038/s41592-019-0686-2
