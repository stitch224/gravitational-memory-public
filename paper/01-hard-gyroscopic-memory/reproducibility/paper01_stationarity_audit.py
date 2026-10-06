#!/usr/bin/env python3
"""Paper 01 deterministic stationary-sharpness validation.

Numerical validation only.  The analytic sharpness proof uses the compact
skew-adjoint operator K_stat=P0 K0 P0 and exact endpoint recovery.
Dependencies: numpy, scipy.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import svd
from scipy.interpolate import CubicSpline
from scipy.integrate import quad, solve_ivp

SIGMA3 = np.array([[1.0, 0.0], [0.0, -1.0]])
SIGMA1 = np.array([[0.0, 1.0], [1.0, 0.0]])
I2 = np.eye(2)
NMINUS = np.block([[I2, -I2], [np.zeros((2, 2)), I2]])


def quadrature_operator(n):
    x, w = leggauss(n)
    t = (x + 1.0) / 2.0
    w = w / 2.0
    sw = np.sqrt(w)
    d = t[:, None] - t[None, :]
    kernel = 0.5 * d * np.abs(d)
    M = sw[:, None] * kernel * sw[None, :]
    v = sw.copy()                 # ||1||_{L2(0,1)}=1 under quadrature
    v /= np.linalg.norm(v)
    P = np.eye(n) - np.outer(v, v)
    return t, w, M, P


def top_subspace_constraint_det(M, w):
    """Diagnostic: old full top singular plane vs two stationarity conditions.

    A real skew matrix has a two-dimensional top singular subspace.  The
    two constraints are mean(beta)=0 and mean(K0 beta)=0.  Their determinant
    on that plane is nonzero iff no old top vector satisfies both.
    Sign depends on the arbitrary SVD orientation; magnitude is meaningful.
    """
    U, s, _ = svd(M)
    # For a numerically skew matrix, the first two left singular vectors span
    # the real top singular plane.
    E = U[:, :2]
    sw = np.sqrt(w)
    mean_row = sw                  # integral f = <sqrt(w), sqrt(w) f>
    rows = np.vstack([mean_row @ E, mean_row @ (M @ E)])
    # Normalize row 2 by top singular value so determinant is dimensionless.
    rows[1] /= s[0]
    return float(np.linalg.det(rows))


def restricted_norm(n):
    t, w, M, P = quadrature_operator(n)
    return svd(P @ M @ P, compute_uv=False)[0]


def continuous_restricted_pair(n=240):
    t, w, M, P = quadrature_operator(n)
    Mp = P @ M @ P
    U, S, Vh = svd(Mp)
    sw = np.sqrt(w)
    a_vals = U[:, 0] / sw
    b_vals = Vh[0, :] / sw
    a0 = CubicSpline(t, a_vals, bc_type="natural")
    b0 = CubicSpline(t, b_vals, bc_type="natural")

    # Enforce continuous (not only quadrature) mean zero, then normalize.
    ma = quad(lambda z: float(a0(z)), 0.0, 1.0, limit=400)[0]
    mb = quad(lambda z: float(b0(z)), 0.0, 1.0, limit=400)[0]
    na = np.sqrt(quad(lambda z: float((a0(z)-ma)**2), 0.0, 1.0, limit=400)[0])
    nb = np.sqrt(quad(lambda z: float((b0(z)-mb)**2), 0.0, 1.0, limit=400)[0])

    def a(z):
        return (float(a0(z)) - ma) / na

    def b(z):
        return (float(b0(z)) - mb) / nb

    # eta=<a,K0 b>; compute K0 b by nested quadrature.
    def kb(t0):
        return quad(lambda s: 0.5*(t0-s)*abs(t0-s)*b(s), 0.0, 1.0,
                    epsabs=2e-11, epsrel=2e-11, limit=300)[0]

    eta = quad(lambda z: a(z)*kb(z), 0.0, 1.0,
               epsabs=5e-10, epsrel=5e-10, limit=250)[0]
    if eta < 0:
        def a_pos(z, a=a): return -a(z)
        a = a_pos
        eta = -eta
    return S[0], a, b, eta


def exact_omega(rho, a, b, rtol=3e-11, atol=1e-13):
    def rhs(t, y):
        Phi = y.reshape(4, 4)
        A = rho * (a(t)*SIGMA3 + b(t)*SIGMA1)
        G = np.block([[np.zeros((2, 2)), I2], [A, np.zeros((2, 2))]])
        return (G @ Phi).ravel()

    sol = solve_ivp(rhs, (0.0, 1.0), np.eye(4).ravel(), method="DOP853",
                    rtol=rtol, atol=atol, max_step=0.005)
    Phi = sol.y[:, -1].reshape(4, 4)
    S = NMINUS @ Phi
    C = S[:2, :2] + S[2:, 2:]
    O = 0.5*(C-C.T)
    return float(O[1, 0])


def endpoint_derivative(f):
    return 2.0*quad(lambda z: f(z), 0.0, 1.0, limit=400)[0]


def main():
    print("n full_K0_norm top_subspace_constraint_det top_subspace_constant_projection_fraction compressed_mean_zero_norm 2/compressed")
    for n in (50, 100, 200, 300, 500):
        t, w, M, P = quadrature_operator(n)
        full = svd(M, compute_uv=False)[0]
        det = top_subspace_constraint_det(M, w)
        U, _, _ = svd(M)
        const_proj = float(np.sum((np.sqrt(w) @ U[:, :2])**2))
        comp = svd(P @ M @ P, compute_uv=False)[0]
        print(f"{n:3d} {full:.13f} {det:+.12f} {const_proj:.12f} {comp:.13f} {2/comp:.11f}")

    sdisc, a, b, eta = continuous_restricted_pair(240)
    ma = quad(a, 0.0, 1.0, limit=400)[0]
    mb = quad(b, 0.0, 1.0, limit=400)[0]
    print("continuous_spline_pair",
          f"s_discrete={sdisc:.13f}", f"eta={eta:.13f}",
          f"mean_a={ma:+.3e}", f"mean_b={mb:+.3e}", f"1/eta={1/eta:.11f}")
    print("stationary exact-Jacobi ratios Q/|omega| -> 1/s_stat; rho is recovery amplitude")
    for rho in (0.03, 0.01, 0.003, 0.001):
        om = exact_omega(rho, a, b)
        Q = 2.0*rho**2
        print(rho, f"{Q/abs(om):.11f}")
    p1 = 2.0*quad(lambda z: (1.0-z)*a(z), 0.0, 1.0, limit=400)[0]
    q1 = 2.0*quad(lambda z: (1.0-z)*b(z), 0.0, 1.0, limit=400)[0]
    print("stationary_endpoint_derivatives",
          f"pprime1={endpoint_derivative(a):+.3e}",
          f"qprime1={endpoint_derivative(b):+.3e}")
    print("stationary_final_shear_per_unit_rho",
          f"p1={p1:+.9f}", f"q1={q1:+.9f}",
          f"magnitude={np.hypot(p1,q1):.9f}")


if __name__ == "__main__":
    main()
