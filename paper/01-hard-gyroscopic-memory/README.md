# Paper 01 — Hard Gyroscopic Gravitational-Wave Memory

**English title:**  
Minimum Fixed-Window Null Curvature Action Associated with Hard Gyroscopic Gravitational-Wave Memory: A Sharp Bel--Robinson Bound

**Target journal:** Physical Review D  
**Article type:** Research Article

## Manuscripts

Both language versions are retained.

- `manuscript-en.tex` / `manuscript-en/` / `manuscript-en.pdf` — English PRD submission manuscript.
- `manuscript-ja.tex` / `manuscript-ja/` / `manuscript-ja.pdf` — Japanese companion manuscript with scientific content aligned to the English version; typesetting differs.

## Main result

For the stationary-in / stationary-out leading-order TT/Bondi radiative-data problem in the fixed convention and fixed auxiliary null window used in the manuscript,

```text
K_stat = P0 K0 P0
||K_stat|| = 0.0281140680170...
2 / ||K_stat|| = 71.13876223067...

omega_gamma^(2) = -2 DeltaPsi_H^(2)
kappa_H^lead(psi) = 2 |psi| / ||K_stat||.
```

The exact nonlinear stationary Jacobi endpoint problem is distinct and obeys

```text
kappa_stat(omega) = |omega| / ||K_stat|| + O(omega^2).
```

## Scope

The theorem is restricted to four dimensions, weak amplitude, leading radiation-zone order, a fixed local radiation-zone patch, a fixed radiation-crossing null window and adapted screen, the standard initially shear-free Bondi frame, and stationary-in / stationary-out radiative data. Final shear is not constrained to vanish.

The cost is a fixed-window auxiliary-null-probe Bel--Robinson curvature action. It is not Bondi energy, Isaacson energy, radiated energy, detector energy, or a direction/window-independent intrinsic exposure.

The paper does not claim a finite-amplitude value law on the full Einstein solution space, a full gyroscopic-memory bound, a spin-memory bound, or a theorem valid in all BMS frames.

## Reproducibility

Public reproduction assets are grouped under `reproducibility/`:

- `README.md` — reproduction manifest and frozen validation values.
- `requirements.txt` — pinned Python dependencies.
- `paper01_stationarity_audit.py` — stationary-spectrum and exact-recovery validation.
- `p6_validate_p2_p5.py` — nonlinear Jacobi / polarization-area validation.
- `e2_bel_robinson_factor_diagnostic.py` — Bel--Robinson normalization diagnostic.

The numerical decimal values validate the analytic theorem; the exact coefficients are defined by operator norms.

## Internal provenance

The research repository additionally retains `internal/` with content-lock, red-team/review, traceability, roadmap, and submission-working records. **`internal/` is not part of the public-release payload.**