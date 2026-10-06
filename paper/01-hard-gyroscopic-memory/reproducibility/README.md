# Paper 01 — Reproducibility Manifest

Status: **PUBLIC RELEASE CANDIDATE / R4_REPRODUCIBILITY_FREEZE**  
Date: 2026-10-05

This manifest separates analytic proof, deterministic floating-point validation, internal audit, and external independent validation. Numerical convergence is not used as a substitute for proof.

## 1. Analytic statements in the manuscript

The manuscript itself proves the following statements.

- With the frozen Bel–Robinson convention,
  `B(k,k,k,k)=1/2 tr(T^2)` and `Q_gamma=1/2 int tr(A^2) dt`.
- The nonlinear Jacobi endpoint has
  `omega=2<a,K0 b>+O(Q^2)`; the antisymmetric endpoint channel contains only even control orders.
- In the unrestricted Jacobi class,
  `kappa_Omega(omega)=|omega|/||K0||+O(omega^2)`.
- Stationary-in/stationary-out TT radiation requires
  `int a = int b = 0`; hence the relevant closed space is
  `U0={f in L2(0,1): int f=0}`.
- With `P0` the orthogonal projection onto `U0`, define
  `K_stat=P0 K0 P0`.  It is compact, skew-adjoint and nonzero.
- Compactness gives an attained top singular pair in `U0`; the exact-endpoint recovery argument then yields
  `kappa_stat(omega)=|omega|/||K_stat||+O(omega^2)`.
- Reconstruction
  `p(t)=2 int_0^t (t-s)a(s) ds`, `q(t)=2 int_0^t (t-s)b(s) ds`
  gives early-zero and stationary-in/out TT strain for every `a,b in U0`.
- In the frozen early shear-free Bondi/radiation-zone convention,
  `omega_gamma^(2)=-2 DeltaPsi_H^(2)`.
- The exact nonlinear stationary Jacobi value law remains
  `kappa_stat(omega)=|omega|/||K_stat||+O(omega^2)`.
- Separately, write `Q_gamma(eps)=eps^2 Q_gamma^(2)+O(eps^3)` and
  `DeltaPsi_H(eps)=eps^2 DeltaPsi_H^(2)+O(eps^3)`.  The superscript `(2)`
  denotes the perturbative coefficient, not a quantity already containing `eps^2`.
- For this leading coefficient problem in the restricted TT/Bondi radiative-data class,
  `kappa_H^lead(psi)=2|psi|/||K_stat||` exactly, where `psi=DeltaPsi_H^(2)`
  and the minimized cost is `Q_gamma^(2)`. This is not a claim about a full finite-amplitude Einstein-solution value function.
- Smooth compactly supported mean-zero tidal controls are dense in `U0`, so the same sharp hard-memory constant is retained for smooth stationary bursts as an infimum.

The decimal coefficient below is validation only; the exact theorem coefficient is `2/||K_stat||`.

## 2. Frozen upstream validation assets

### P6 nonlinear Jacobi / polarization-area validation

Path: `p6_validate_p2_p5.py`

- source commit: `157af997c3a65173bc931e32c8bfbe14051a118d`
- source blob: `c6e64d80edd1b217649112c985c459e9714885da`
- dependencies: Python 3, NumPy, SciPy

Representative frozen outputs:

- `||K0||` at 500 nodes: `0.0910384382384`
- `1/||K0||`: `10.98437121...`
- exact Jacobi `Q/|omega|` at eps `0.001`: `10.98437118`

### E2 Bel–Robinson normalization diagnostic

Path: `e2_bel_robinson_factor_diagnostic.py`

- source commit: `f8c1c1c54694c2cf18afe0d5cbf98a18fbe4e648`
- source blob: `f2bd1c427323fb62223fb6d5e56dea4b9d7aa543`
- dependency: Python 3 standard library

This checks exact rational identities for the frozen `1/4` Bel–Robinson convention, normalization conversion and positive affine-rescaling invariance.

## 3. Paper 01 stationary-sharpness audit

Path: `paper01_stationarity_audit.py`

Dependencies: Python 3, NumPy, SciPy.

The script checks four things:

1. reproduces the unrestricted `K0` norm;
2. diagnoses that the old full top singular subspace is not automatically compatible with the two mean-zero stationarity constraints;
3. computes the Nyström norm of `P0 K0 P0`;
4. reconstructs a continuous mean-zero top pair, evolves the exact Jacobi ODE, checks `Q/|omega| -> 1/||K_stat||`, and verifies the final strain derivatives vanish.

Restricted spectral convergence:

```text
n   full_K0_norm      compressed_mean_zero_norm   2/compressed
50  0.0910384416041   0.0281140528760            71.13880054294
100 0.0910384384527   0.0281140670547            71.13876466581
200 0.0910384382516   0.0281140679578            71.13876238062
300 0.0910384382407   0.0281140680065            71.13876225721
500 0.0910384382384   0.0281140680170            71.13876223067
```

Continuous reconstructed pair / exact Jacobi replay:

```text
s_discrete=0.0281140679892
eta=0.0281140680641
mean_a=+1.388e-17
mean_b=-4.163e-17
1/eta=35.56938105581

stationary exact-Jacobi ratios Q/|omega|
0.03  35.56942171300
0.01  35.56938562391
0.003 35.56938151880
0.001 35.56938115793

stationary_endpoint_derivatives
p1=+2.776e-17
q1=-8.327e-17
```

Interpretation: these numbers validate the implementation and conventions. Restricted sharpness itself is analytic, from compact-operator norm attainment plus exact endpoint recovery.

## 4. Hard-gyro bridge validation retained from P7-2

For `W=int(p qdot-q pdot)du`, the early-zero tests obey `omega^(2)=W/2`, `DeltaPsi_H^(2)=-W/4`.

Representative cases:

- smooth circular: `W=2.35619447`, `omega=1.17809724`, `DeltaPsi=-0.58904862`;
- opposite chirality: signs reverse, ratio remains `-2`;
- ellipse: `W=1.22522112`, `omega=0.61261056`, `DeltaPsi=-0.30630528`;
- fixed axis and chirality cancellation: both channels vanish;
- endpoint memory: `W=-0.62590250`, `omega=-0.31295125`, `DeltaPsi=0.15647563`.

A constant-baseline adversarial test leaves the chord-closed Jacobi area unchanged but changes the hard functional. This is why the manuscript explicitly fixes the early shear-free frame.

## 4.1 R2 local sanity replay and scope checks

The current version of `paper01_stationarity_audit.py` was independently replayed as a non-canonical local sanity check during R2, and the same retained script then passed the governed self-hosted final replay in run `37287597777`. It reproduced the manuscript values to the shown precision:

- `||K0|| = 0.0910384382384` at 500 nodes;
- `||K_stat|| = 0.0281140680170`;
- `2/||K_stat|| = 71.13876223067`;
- continuous-pair `eta = 0.0281140680641`, `1/eta = 35.56938105581`;
- exact-Jacobi ratios `35.56942171300, 35.56938562391, 35.56938151880, 35.56938115793`;
- final-shear magnitude per unit recovery amplitude `0.568784861`.

The endpoint derivatives are at floating-point roundoff and their last digits/signs are environment-sensitive; only the vanishing check is scientifically meaningful. Formal publication evidence is the governed self-hosted run `37287597777`; the local calculation remains only a sanity cross-check.

R2 also records that the theorem uses a fixed local radiation-zone null window and an adapted screen/Bondi polarization identification. Changing the physical window or auxiliary null probe is a different cost problem, not an affine reparameterization of the same one.

## 5. Literature and novelty status

Primary sources used for known structures are listed in the manuscript. R3 reran the collision search against the final post-R2 theorem fingerprint and found no same-condition theorem; historical novelty remains `OPEN_SAME_CONDITION`, not a world-first certification. Adjacent area-constrained bending/elastica literature is now cited explicitly and excluded from the novelty claim.

## 6. R4 clean-environment freeze

Pinned Python dependencies:

- Python observed in the clean self-hosted environment: `3.14.6`;
- `numpy==2.5.3`;
- `scipy==1.18.1`;
- requirements file: `paper/01-hard-gyroscopic-memory/requirements-repro.txt`.

The R4 workflow creates a fresh virtual environment from that requirements file and replays, in order:

1. `p6_validate_p2_p5.py`;
2. `e2_bel_robinson_factor_diagnostic.py`;
3. `paper01_stationarity_audit.py`;
4. manuscript scientific-value markers and stale-number rejection.

The first clean R4 attempt, run `37309548211`, numerically reproduced all three scientific audits but was formally marked failed by an audit-parser bug: the six-column stationarity line was parsed as though it had seven columns. The scientific outputs were already correct. The parser was corrected and the Python dependencies were pinned before the canonical retry.

Canonical R4 retry: self-hosted run `37309724179`, tested SHA `28b16303a080d4301f99fc27ed968574fec3f6ba`.

The clean reproducibility step on that run is **PASS** and reproduces:

- P6 `||K0||` at n=500: `0.0910384382384`;
- P6 `1/||K0||`: `10.98437120904`;
- E2 exact-rational Bel--Robinson diagnostic: `PASS: 4 cases`;
- Paper 01 `||K_stat||`: `0.0281140680170`;
- Paper 01 `2/||K_stat||`: `71.13876223067`;
- continuous-pair `eta=0.0281140680641`, `1/eta=35.56938105581`;
- exact-Jacobi ratios `35.56942171300, 35.56938562392, 35.56938151878, 35.56938115796`;
- final-shear magnitude per unit recovery amplitude `0.568784861`;
- manuscript markers for the signed factor-2 bridge, `K_stat=P0 K0 P0`, and the current numerical coefficient;
- absence of the stale `0.028114068018...` manuscript value.

Endpoint-derivative roundoff signs/last digits are environment-dependent and are not frozen as scientific values.

### Public-export reproducibility subset

The publication-facing reproducibility subset is fixed as:

- `paper/01-hard-gyroscopic-memory/manuscript-ja.tex` and `manuscript-ja/*.tex`;
- `paper01_stationarity_audit.py`;
- `paper/01-hard-gyroscopic-memory/requirements-repro.txt`;
- `paper/01-hard-gyroscopic-memory/REPRODUCIBILITY.md`;
- `p6_validate_p2_p5.py`;
- `e2_bel_robinson_factor_diagnostic.py`.

Internal research roadmaps, review ledgers, branch-management files and unrelated B/G-lane artifacts are not required to reproduce Paper 01 and are not part of this public-export subset.

### External independent replication

No external third-party independent replication of the new restricted coefficient has been obtained as of 2026-10-05. This is explicitly recorded as **not obtained**, not silently treated as satisfied. It is not a mathematical blocker because the proof is analytic and the retained numerical work is validation, but it remains useful post-release evidence if later obtained.

## 7. Prior formal replay provenance

- R2 final formal replay: self-hosted run `37287597777` tested SHA `30d5ff926f00809219a7306cbe91afbfc653af61`; stationarity audit PASS, three-pass LuaLaTeX build PASS, final warning count 0, tracked-tree audit PASS, PDF publish PASS. Published PDF SHA-256: `33367fabe6741dd6b982b44fb7e284ef9e182fd039ac8a5a471279434a1ac12e`.
- R3 final manuscript rebuild: self-hosted run `37296949233` (run #19) PASS; LuaLaTeX build, PDF artifact upload, tracked-tree audit and PDF publish all passed.
- Superseded R2 run `37286957196` had already passed the stationarity audit but failed during isolated TinyTeX preparation with a runner-local TeX Live database I/O error. No manuscript PDF was built or published from that failed run; scientific inputs were unchanged and the later canonical run passed.
