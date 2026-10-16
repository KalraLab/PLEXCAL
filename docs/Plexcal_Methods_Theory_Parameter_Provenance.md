# Plexcal: Methods, Theory, Parameter Provenance, and Reproducibility

**Document type:** Software methods and parameter-provenance record  
**Software:** Plexcal — structure-aware estimation of aromatic-residue excitation-transfer networks and modelled diffusion lengths  
**Version:** 0.1.0 *(replace with the release being archived)*  
**Status:** Methods documentation for a software release; not a claim of experimental validation  
**Primary methods context:** the protein-polymer energy-migration work from A. P. Kalra and collaborators, particularly Kalra et al. (2023), *Electronic Energy Migration in Microtubules*; Kalra et al. (2023), *All Lit Up: Exploring the Photophysical Properties of Protein Polymers*; and Kakati, Karmakar, and Kalra (2025), *Triplet Energy Migration in Cytoskeletal Polymers*. Plexcal adapts this research context into a software implementation; code-specific departures are documented explicitly below.  
**Code reviewed for this document:** `calculate_fret_diffusion_lengths.py` (FRET-labelled implementation) and `calculate_det_diffusion_lengths.py` (DET-labelled implementation), supplied in the project workspace.  

> **Important scope statement.** The scientific foundation is attributed primarily to the published work of Kalra and collaborators and to the Kakati et al. study on triplet energy migration, not to the project team’s disease-protein manuscript. Plexcal is a software adaptation of this established research direction, with implementation-specific changes recorded here. This document distinguishes (i) theory and methods supported by the Kalra/Kakati literature, (ii) operations and constants actually present in the reviewed Python scripts, and (iii) assumptions or provenance details still requiring verification. The supplied PDFs do not include the complete Supplementary Information tables needed to map every hard-coded parameter to an exact table entry and unit. Values without direct traceability are marked **provenance not yet independently verified** rather than assigned to a source by inference.

---

## 1. Purpose and intended use

Plexcal is a computational tool for constructing residue-level aromatic-chromophore networks from protein structures and estimating excitation-transfer behaviour under simplified FRET-like (singlet-state) and Dexter-like (triplet-state) kinetic models. The present implementation identifies tryptophan (TRP/Trp), tyrosine (TYR/Tyr), and phenylalanine (PHE/Phe) residues, calculates residue-centroid separations, assigns model-specific pairwise transfer weights, and simulates stochastic excitation trajectories terminated by a residue-specific decay event.

The primary output currently produced by the reviewed scripts is a mean **initial-to-final Euclidean displacement** over simulated trajectories. This quantity is denoted `diffusion_length_angstrom` in the code. It is not the accumulated path length, and it should not be described as a measured diffusion coefficient or as an experimentally established transport distance. In the Kalra-group energy-migration framework, the corresponding endpoint-displacement statistic is used as the photoexcitation diffusion length, \(L_D\) (Kalra et al., 2023; Kakati et al., 2025).

Plexcal is intended for exploratory structural comparisons, model-based hypothesis generation, and reproducible calculation of the outputs defined by a particular software release. Predictions depend on the chosen structural representation, parameter set, rate equations, decay model, and implementation details. They should not be interpreted as direct measurements of energy transfer in a protein.

## 2. Scientific context

Electronic excitation can migrate between chromophores by different mechanisms. Förster resonance energy transfer (FRET) describes long-range dipole–dipole-mediated transfer between donor and acceptor transitions, subject to spectral overlap and transition-dipole orientation [2–4]. Dexter energy transfer (DET) is an exchange-mediated process whose electronic coupling depends strongly on short-range donor–acceptor orbital overlap and therefore decays rapidly with separation [5,6]. Both can be represented in a kinetic hopping model when transfer events are treated as stochastic transitions between localized sites.

Kalra and collaborators developed this research direction for protein polymers, including microtubules and aromatic-residue networks; Kakati et al. extended the discussion to triplet-state energy migration in cytoskeletal polymers. These studies provide the principal scientific context for Plexcal’s singlet/FRET-like and triplet/DET-like model families. It identifies the residue chromophore sites from TRP, TYR, and PHE, builds distance matrices from their geometric centroids, estimates pairwise coupling/rate matrices, and simulates \(10^4\) photoexcitation trajectories. It defines the reported diffusion length as

$$
L_D
=
\frac{1}{N}
\sum_{i=1}^{N}
\left\|
\mathbf{r}_{f,i}
-
\mathbf{r}_{0,i}
\right\|,
$$

where \(N\) is the number of trajectories, and \(\mathbf r_{0,i}\) and \(\mathbf r_{f,i}\) are the initial and final chromophore positions. This endpoint-displacement definition follows the kinetic-hopping framework used in the cited Kalra-group/Kakati research; the precise code implementation is assessed separately in this document.

This is an **incoherent kinetic Monte Carlo model**. It does not propagate a quantum wavefunction or density matrix and does not by itself calculate quantum coherence, excitonic eigenstates, superradiance, or a coherence length. Any pair-count threshold in the code must therefore be treated as a separate heuristic descriptor, not proof that a pair is quantum coherent.

## 3. Input structures and residue representation

### 3.1 Structural source

The manuscript workflow uses experimentally resolved structures from the RCSB Protein Data Bank (PDB) and describes selection of disease-associated fibril structures, extraction of repeat units, and polymer reconstruction using structure-specific helical geometry [1,7]. Plexcal can analyse a supplied PDB file or a structure obtained by an identifier through its application workflow. Each result should record the exact input structure, source, retrieval date where applicable, model/chain selection, and any preprocessing.

PDB coordinates are structural inputs, not photophysical measurements. A structure alone does not provide all quantities required for rigorous residue-specific transfer calculations, such as environment-dependent excited-state energies, transition densities, spectral overlap integrals, or a direct measurement of the electronic coupling.

### 3.2 Aromatic sites and geometric centres

The reviewed scripts select residues labelled TRP, TYR, and PHE. For each selected residue, the scripts calculate the arithmetic mean of the coordinates of the atoms grouped into that residue and use that mean as the site position. Pairwise site distances are then calculated with the Euclidean metric:

$$
R_{ij}
=
\left\|
\mathbf{r}_j-\mathbf{r}_i
\right\|.
$$

The published methods describe chromophore positions using aromatic-ring centroids (Kalra et al., 2023; Kakati et al., 2025). The code should not be claimed to implement precisely that operation for all calculations without qualification: the site-centroid calculation in the reviewed scripts averages the coordinates of all atoms assigned to the residue, whereas the FRET orientation calculation separately selects a set of aromatic ring atoms. If a ring-system centroid is intended for the distance matrix, the implementation should be changed and versioned accordingly.

### 3.3 Structural parsing limitations in the reviewed scripts

The current parsers group atoms by residue name and residue number but do not retain chain identifier or insertion code in the grouping key. In a multi-chain structure, residues with the same number and residue type can therefore be merged incorrectly. This is especially consequential for fibrils, oligomers, and assemblies with repeated chains. The scripts also read fixed-column PDB records and do not constitute a complete mmCIF parser. Before reporting results for multi-chain assemblies, residue identity should be keyed at minimum by model, chain, residue number, insertion code, and residue name, and alternate locations should be handled explicitly. These are implementation observations, not literature-derived assumptions.

## 4. FRET model: theoretical basis and implementation

### 4.1 Dipole–dipole coupling

In the point-dipole approximation, the electronic coupling between donor \(D\) and acceptor \(A\) has the general distance and orientation dependence

$$
V_{DA}\propto
\frac{\kappa\,\mu_D\mu_A}{n^2R_{DA}^{3}},
$$

where \(R_{DA}\) is donor–acceptor separation, \(\mu_D\) and \(\mu_A\) are transition-dipole magnitudes, \(n\) is the refractive index of the medium, and \(\kappa\) is an orientation factor. The corresponding transfer rate in the dipole–dipole FRET regime scales approximately as \(R^{-6}\), provided the other photophysical terms remain fixed [2–4]. The Kalra-group literature describes dipole–dipole coupling, orientation dependence, and interaction-specific photophysical parameters; the exact mapping of the current script’s constants to the original supplementary tables still needs to be checked against the relevant Supporting Information.

The reviewed FRET script estimates an orientation-related quantity from ring-plane normals. It obtains a normal vector for the selected aromatic ring atoms using singular-value decomposition and computes

$$
\kappa_{\mathrm{proxy}}
=
\mathbf{n}_i\cdot\mathbf{n}_j
-
3\left(\mathbf{n}_i\cdot\hat{\mathbf{R}}_{ij}\right)
\left(\mathbf{n}_j\cdot\hat{\mathbf{R}}_{ij}\right).
$$

It then calculates a dipole–dipole-like coupling from this quantity, fixed residue-type dipole magnitudes, the vacuum permittivity, and the centroid separation converted from Å to metres. The calculated energy in joules is converted to wavenumbers using the numerical factor \(5.03\times10^{22}\), which corresponds approximately to converting joules to cm\(^{-1}\) through \(1/(hc)\) when SI units are used.

**Interpretation caveat:** the normal to an aromatic ring is not generally identical to the electronic transition-dipole direction. The ring-normal construction is therefore a geometric proxy, not a quantum-chemical transition-dipole calculation. The manuscript's use of \(\kappa^2\) and the code's use of a signed ring-normal expression should not be assumed equivalent without an explicit derivation and validation. Because the coupling is subsequently squared in the rate expression, sign may disappear in that step, but the proxy geometry remains an approximation.

### 4.2 FRET rate expression actually implemented

The reviewed script calculates a rate-like matrix using the following expression (with the code's variable names):

$$
k_{ij}^{\mathrm{FRET,code}}
=
\frac{1}{c_{\mathrm{cm/s}}\,\hbar_{\mathrm{cm}}^2}
\frac{1}{n^4}
J_{ij}^{\mathrm{code}}\,
V_{ij,\mathrm{cm}^{-1}}^2,
$$

where the script uses `c = 3.0e10`, `hbar_cm = 5.29e-12`, `n = 1.4`, and the hard-coded pair-specific `j_values` dictionary. It converts the resulting values from s\(^{-1}\) to ns\(^{-1}\) by multiplying by \(10^{-9}\). This is a description of the code, not an endorsement that the equation is dimensionally or physically identical to a conventional Förster-rate expression. The exact units and definitions of the `j_values` must be documented from the original source before the rate can be considered fully traceable.

The conventional Förster relationship is often expressed as

$$
k_{\mathrm{FRET}}
=
\frac{1}{\tau_D}
\left(\frac{R_0}{R_{DA}}\right)^6,
$$

where \(\tau_D\) is the donor excited-state lifetime and \(R_0\) is the Förster radius determined by the donor quantum yield, spectral overlap, orientation factor, and refractive index [2–4]. The reviewed script does **not** directly implement this equation with an explicitly calculated \(R_0\). Do not describe it as a direct implementation of the conventional \(R_0\)-based formula unless the code is revised or an algebraic equivalence is demonstrated.

### 4.3 FRET numerical constants in the current code

The FRET script currently uses the following fixed values:

- **Dipole magnitudes:** TRP = 2.074 D, TYR = 1.18 D, PHE = 0.28 D.
- **Refractive index:** \(n = 1.4\).
- **Pair-specific `j_values`:**
  - TRP–TRP = \(5.5953\times10^{-6}\)
  - TYR–TYR = \(4.7933\times10^{-4}\)
  - PHE–PHE = \(6.1123\times10^{-4}\)
  - TRP–TYR = \(2.3247\times10^{-5}\)
  - TYR–TRP = \(9.2550\times10^{-4}\)
  - PHE–TYR = \(0.0043\)
  - PHE–TRP = \(0.0039\)
  - TRP–PHE = \(1.3034\times10^{-7}\)
  - TYR–PHE = \(1.8891\times10^{-5}\)
- **Decay lifetimes used by the FRET-labelled script:** TRP = 3.6 ns, TYR = 7.5 ns, PHE = 3.39 ns. The script converts these to rates by taking their reciprocals.
- **Trajectory count:** 10,000.
- **Maximum events per trajectory:** 1,000.

**Parameter provenance status:** these exact values are visible in the supplied code. The supplied manuscript says the corresponding FRET parameters are reported in Supplementary Table S3, but that supplementary table was not included in the two supplied PDFs. It is therefore not yet verified that every listed code value, unit, directionality, and lifetime matches Table S3 or the cited source. The asymmetric pair-specific values are especially important to trace: do not silently symmetrise them or attribute them to a particular paper until the original table/source has been checked.

## 5. Dexter energy-transfer model: theoretical basis and implementation

### 5.1 Mechanistic basis

Dexter transfer is an exchange-mediated process requiring sufficient donor–acceptor orbital overlap. Its electronic coupling decreases strongly with separation, commonly represented with an exponential distance dependence [5,6]. This differs from the longer-range dipole–dipole dependence used for FRET. The manuscript states that it uses an empirical Dexter coupling–separation relationship and refers readers to Supplementary Table S4 and the earlier triplet-energy-migration work [1,8].

### 5.2 Distance-dependent coupling implemented

The reviewed DET script computes a coupling proxy using

$$
V_{ij}^{\mathrm{DET,code}}\,(\mathrm{eV})
=
\frac{3.22\times10^{-2}}{1.4\times1.4}
\exp\left[
-(10^{10}\times2.5579)
\left(R_{ij,\mathrm{m}}-4.15\times10^{-10}\right)
\right].
$$

The distance is supplied in metres. The result is then multiplied by
\(1.602\times10^{-19}\) and \(5.03\times10^{22}\) to obtain the code's
wavenumber-like coupling quantity in \(\mathrm{cm}^{-1}\).
The diagonal is set to zero after the exponential calculation.

This expression is reproduced from the current script to make the implementation auditable. Its prefactor, exponent coefficient, reference distance, refractive-index-like divisor, and units need to be matched against the actual Supplementary Table S4 and the source cited for the empirical Dexter relation before they are described as literature-validated constants. The supplied main-text PDFs do not contain enough information to establish the full derivation of each numerical coefficient.

### 5.3 DET rate expression and parameters implemented

The script calculates a rate-like matrix as
$$
k_{ij}^{\mathrm{DET,code}} =
\frac{1}{c_{\mathrm{cm/s}}\,\hbar_{\mathrm{cm}}^2}
\,P_{ij}\,V_{ij,\mathrm{cm}^{-1}}^2,
$$

using \(c = 3.0\times10^{10}\), \(\hbar_{\mathrm{cm}} = 5.29\times10^{-12}\), and a fixed donor–acceptor transfer-parameter matrix \(P_{ij}\). The output is converted from \(\mathrm{s}^{-1}\) to \(\mathrm{\mu s}^{-1}\) by multiplying by \(10^{-6}\). The matrix is constructed from these values in `det_code_clean(1).py`:

- TRP–TRP = 0.00058
- TYR–TYR = 0.0012
- PHE–PHE = 0.00098
- TRP–TYR = 0.000487
- TYR–TRP = 0.0011
- PHE–TYR = 0.0012
- PHE–TRP = 0.0039
- TRP–PHE = 0.00031
- TYR–PHE = 0.00098

The DET-labelled script uses lifetimes TRP = 7.3 µs, TYR = 5.6 µs, and PHE = 3.1 µs. The reciprocal lifetimes form the decay-rate column. It simulates 10,000 trajectories with at most 1,000 event steps per trajectory.

**Parameter provenance status:** the numerical values above are directly verified in the supplied script, but their units, spectral interpretation, and exact correspondence to Supplementary Table S4 are **not verified** because the SI table is absent. They must not be presented as independently validated literature constants until the supplementary material and cited primary source have been checked. The term “transfer parameter” is intentionally used here rather than assigning a physical meaning not established by the supplied code or main text.

## 6. Kinetic Monte Carlo algorithm

For each structure, the scripts perform the following operations:

1. Parse atom coordinates and identify TRP, TYR, and PHE sites.
2. Calculate site positions and all pairwise distances.
3. Build the model-specific transfer-rate matrix.
4. Add a decay channel to each site using the reciprocal of the assigned lifetime.
5. Choose a starting site uniformly from the available sites for each trajectory.
6. At the current site, sum all outgoing transfer rates and the decay rate. Draw a uniform random number over the total rate and select an event according to the cumulative rate weights.
7. If the selected event is decay, terminate the trajectory. Otherwise move the excitation to the selected acceptor and repeat, up to the maximum event count.
8. Calculate the Euclidean distance between the initial site and the final site for each trajectory and average these endpoint distances.

This is the standard event-selection interpretation of a continuous-time kinetic Monte Carlo / kinetic hopping scheme: if all outgoing rates are non-negative, the probability that a particular event is selected at a site is its rate divided by the sum of all rates from that site. The reviewed scripts use the rates as event weights to select the next event. The implementation does not explicitly sample and report physical waiting times between events; the endpoint distribution is therefore conditioned by the relative event rates and terminating decay channel rather than a separately reported trajectory-time distribution.

The referenced protein-polymer energy-migration workflow uses independent stochastic excitation trajectories; the current scripts set the trajectory count to \(10^4\). It is important to preserve the precise output definition: \(L_D\) is the average endpoint displacement, not the sum of distances traversed over all hops. A trajectory that wanders and returns near its starting point can have a small endpoint displacement despite a long path.

## 7. “Coherent pair” count: what the code does and does not establish

Both scripts count unique lower-triangle residue pairs whose absolute coupling-like quantity exceeds 53.052 cm\(^{-1}\). The resulting fields are `total_pairs`, `coherent_pairs`, and `coherent_pair_fraction`.

This is a code-defined threshold statistic. Neither supplied main-text manuscript establishes 53.052 cm\(^{-1}\) as a general criterion for quantum coherence, nor does the kinetic Monte Carlo trajectory algorithm calculate phase coherence, dephasing, exciton delocalisation, or superradiance. The output should therefore be named and described conservatively—for example, **“pairs above the configured coupling threshold”**—until the threshold's physical derivation and applicability are documented. If retained as a legacy field, the software documentation should explicitly state that `coherent_pairs` is a heuristic label and not a quantum-coherence measurement.

## 8. Parameter provenance and evidence hierarchy

Scientific software documentation should let a user trace each numerical value to an identifiable source and to the exact code location that uses it. For Plexcal, provenance should be handled as follows:

- **Directly verified from implementation:** the numerical literals and operations described in Sections 4–7, as found in the two reviewed Python files.
- **Supported by the Kalra/Kakati research context:** aromatic-residue energy-transfer networks, FRET-like singlet transfer, Dexter-like triplet transfer, rate-based kinetic hopping, competition with intrinsic decay, and endpoint-displacement analysis. The exact implementation choices must be checked against the cited papers and their Supporting Information.
- **Cited theoretical background:** FRET distance/orientation/spectral-overlap theory [2–4]; Dexter exchange transfer [5,6]; structural source and coordinate provenance through the PDB [7].
- **Not yet independently traceable:** exact numerical correspondence between the hard-coded `j_values`, `TRANSFER_PARAMETERS`, lifetime choices, and Dexter exponential constants and the original source tables. Kakati et al. (2025) is a key source for triplet-energy-migration methodology, but the current files do not include all relevant Supporting Information needed to verify each numerical mapping. Do not assert exact numerical agreement until this audit is complete.

Before tagging a release intended for DOI citation, obtain and archive the Supporting Information for the relevant Kalra-group and Kakati papers. Compare every constant and its unit against the original table and equation, then update this section with table number, row/column or parameter label, original units, conversion steps, and the code variable that receives the value. Where a constant is derived by the software rather than copied, document the derivation and dimensional analysis.

## 9. References and methodological lineage

The central scientific lineage for Plexcal is the work on electronic energy migration in protein polymers by A. P. Kalra and collaborators, together with Kakati et al.'s treatment of triplet energy migration. These are the principal research-context citations for the model family. The classic FRET and Dexter references below support the underlying physical theory; they are not claimed as the source of every implementation-specific numerical literal.

1. Kalra, A. P.; Biswas, S.; Mulrain, I.; Wang, M.; Tuszynski, J. A.; Scholes, G. D. **All Lit Up: Exploring the Photophysical Properties of Protein Polymers.** *The Journal of Physical Chemistry Letters* **2023**, 14, 5891–5900. https://doi.org/10.1021/acs.jpclett.3c01032.  
   *Core context for protein-polymer photophysics and electronic energy migration.*

2. Kalra, A. P. et al. **Electronic Energy Migration in Microtubules.** *ACS Central Science* **2023**, 9, 352–361. https://doi.org/10.1021/acscentsci.2c01114.  
   *Primary context for aromatic-residue-mediated energy migration in protein polymers and rate-based transport modelling.*

3. Kakati, A.; Karmakar, T.; Kalra, A. P. **Triplet Energy Migration in Cytoskeletal Polymers.** *The Journal of Physical Chemistry B* **2025**, 129, 128–138. https://doi.org/10.1021/acs.jpcb.4c06748.  
   *Principal source to consult for triplet-state migration and the DET-related model; exact parameter-to-code correspondence still requires Supporting Information verification.*

4. Förster, T. **Zwischenmolekulare Energiewanderung und Fluoreszenz.** *Annalen der Physik* **1948**, 2, 55–75. https://doi.org/10.1002/andp.19484370105.  
   *Foundational resonance energy-transfer theory.*

5. Scholes, G. D. **Long-Range Resonance Energy Transfer in Molecular Systems.** *Annual Review of Physical Chemistry* **2003**, 54, 57–87. https://doi.org/10.1146/annurev.physchem.54.011002.103746.  
   *Review of resonance energy-transfer theory.*

6. Dexter, D. L. **A Theory of Sensitized Luminescence in Solids.** *The Journal of Chemical Physics* **1953**, 21, 836–850. https://doi.org/10.1063/1.1699044.  
   *Foundational exchange-mediated energy-transfer theory.*

7. Skourtis, S. S.; Liu, C.; Antoniou, P.; Virshup, A. M.; Beratan, D. N. **Dexter energy transfer pathways.** *Proceedings of the National Academy of Sciences* **2016**, 113, 8115–8120. https://doi.org/10.1073/pnas.1517189113.  
   *Additional mechanistic context for Dexter transfer pathways.*

8. Burley, S. K. et al. **Updated resources for exploring experimentally-determined PDB structures and Computed Structure Models at the RCSB Protein Data Bank.** *Nucleic Acids Research* **2025**, 53, D564–D574. RCSB PDB: https://www.rcsb.org/.  
   *Structural data source; verify final journal metadata and DOI before release.*

9. RCSB Protein Data Bank. **PDB archive and structure records.** https://www.rcsb.org/.