# Zenodo / preprint-server submission fields

## Files to upload
1. `kds_mode_stability.pdf`: the manuscript (main file).
2. `kds_mode_stability.tex`: the LaTeX source.
3. `kds_code_v1.zip`: the code for the computer-assisted steps (Appendix A). The paper's "Data availability" line
   says the code is available on request. If you upload the zip, change that line to "The code is included in
   this record".

## Basic information
- **Resource type:** Publication → Preprint
- **Title:** Mode stability for gravitational perturbations of subextremal Kerr–de Sitter black holes
- **Creator:** Padhi, Jaideep Sai
  - Affiliation: Purdue University
  - ORCID: add yours if you have one
- **Publication date:** 2026-11-08 (or the upload date)
- **Version:** v1
- **Language:** English
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)

## Description (abstract, plain text)
We prove mode stability for the Teukolsky equations of spin s = ±2 on Kerr–de Sitter black hole spacetimes in the
full subextremal range: there are no nontrivial outgoing mode solutions at frequencies σ with Im σ ≥ 0. We then
show that every generalized mode solution of the ungauged linearized Einstein equation D_g Ric(h) = Λh with
Im σ ≥ 0 is pure gauge if σ ≠ 0, and the sum of a linearized Kerr–de Sitter metric and a pure gauge term if σ = 0.
This is the mode stability hypothesis under which Hintz, Petersen and Vasy proved the nonlinear stability of
Kerr–de Sitter black holes in the full subextremal range.

For real frequencies, the Teukolsky part follows Hintz's two-current method for scalar waves, applied to a dual
radial equation obtained from an Euler integral transform. The upper half plane is reached by a continuity argument
in which high-frequency control is borrowed from the gauge-fixed Einstein operator through metric reconstruction.
The passage from the Teukolsky scalars to the metric uses a linearized Goldberg–Sachs argument and the two Killing
fields built from the Killing spinor. Several steps are computer-assisted, using exact arithmetic.

AI use: Claude Opus 5.5 (Anthropic) was used to develop and check proofs, write and run the certificate code,
search the literature, and draft the text. See the "Statement on the use of AI" in the manuscript. The argument has
not yet been checked by an independent human expert.

## Keywords
Kerr–de Sitter; mode stability; Teukolsky equation; linearized Einstein equations; quasinormal modes; black hole
stability; cosmological constant; Hertz potential; computer-assisted proof

## Subjects (MSC 2020)
- 83C57 Black holes
- 83C05 Einstein's equations (general structure, canonical formalism, Cauchy problems)
- 35B40 Asymptotic behavior of solutions to PDEs
- 35P25 Scattering theory for PDEs

## Related works ("Related identifiers", relation: References)
- arXiv:2508.06620 — Hintz, Petersen, Vasy, Conditional non-linear stability of Kerr–de Sitter spacetimes in the
  full subextremal range
- arXiv:2608.25868 — Hintz, Mode stability for the Klein–Gordon equation on subextremal Kerr–de Sitter spacetimes
- doi:10.4310/ACTA.2018.v220.n1.a1 — Hintz, Vasy, Acta Math. 220 (2018)
- arXiv:2112.01355 — Petersen, Vasy, Wave equations in the Kerr–de Sitter spacetime: the full subextremal range

## Notes field (optional)
Preprint v1. In-house adversarial referee reports and the full proof ledger are in the accompanying code archive
(`kds/REFEREE2.md`, `kds/REFEREE3.md`, `kds/LEMMAS.md`).

## If the other server ("Hexagon") asks for different fields
Use the same title, author, abstract, keywords and license. Category: Mathematical Physics / General Relativity
(gr-qc, math.AP, math-ph). For a comments field: "17 pages; computer-assisted proof with exact certificates; code
included."
