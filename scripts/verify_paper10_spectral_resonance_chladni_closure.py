#!/usr/bin/env python3
"""
Verification Suite for Paper X:
Finite Euler-Sieve Polynomial Identities, Hadamard-Poisson Flux Detection,
and Conditional Nodal-Wall Identities for Riemann's Xi-Function
==========================================================================
Author: Thijs Missiaen (2026)

High-precision (50-decimal mpmath) verification of:
  0. Canonical normalization:
     Xi(z) = 8 * int_0^inf Phi(u) cos(2 z u) du,
     W(d, t) = (1/2) d/dd |Xi(t - i d)|^2 = 16 * int_0^inf y sinh(2 d y) I(y, t) dy.
  1. Exact cubic-window Dirichlet coefficient identity for b_N^E(X) on
     1 <= N < (floor(X) + 1)^3, exact quadratic algebra of F_2(z) = 1 + z + z^2/2
     along Im(z) = 0 together with its complex zeros F_2(-1 +- i) = 0, and
     explicit finite-prime Kronecker alignment demonstrating Re sum p^{-sigma - i t} < -1
     (refuting any phase-alignment prohibition in 0 < sigma <= 1).
  2. Closed-form Poisson-Hadamard lobe integral F_-(d, d_*) = -arcsin(d / d_*)
     and step defect pi * Theta(d - d_*), plus Sylvester negative inertia jump
     nu_- >= 1 on a POLE-FREE entire even two-zero coalescence perturbation
     Xi_pert(z; d_*), framed as contrapositives of Lagarias's (1999) equivalence.
  3. Local Cauchy-Riemann identity along I(d, tau_k(d)) = 0 on the first 6 extrema,
     explicit polynomial counterexample showing R(d_*, tau(d_*)) = 0 when an
     off-critical zero lies on I = 0, and explicit root eta_0 = 0.3194150268... < pi/8
     of Phi(i * eta_0) = 0 contrasting with Jacobi's theta_3(0 | i e^{4u}) != 0.
"""

import math
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

mp.mp.dps = 50

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PUB_ROOT = os.path.dirname(SCRIPT_DIR)
FIG_DIR_ROOT = os.path.join(PUB_ROOT, "figures")
FIG_DIR_PAPER = os.path.join(PUB_ROOT, "paper", "figures")
os.makedirs(FIG_DIR_ROOT, exist_ok=True)
os.makedirs(FIG_DIR_PAPER, exist_ok=True)


def sieve_primes(limit_int: int):
    if limit_int < 2:
        return []
    sieve = bytearray(b"\x01") * (limit_int + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(math.isqrt(limit_int)) + 1):
        if sieve[p]:
            sieve[p * p : limit_int + 1 : p] = b"\x00" * (((limit_int - p * p) // p) + 1)
    return [i for i, v in enumerate(sieve) if v]


def phi_kernel_mp(u, n_max: int = 60):
    total = mp.mpc(0)
    e4u = mp.exp(4 * u)
    e5u = mp.exp(5 * u)
    e9u = mp.exp(9 * u)
    for n in range(1, n_max + 1):
        n2 = n * n
        n4 = n2 * n2
        arg = mp.pi * n2 * e4u
        if mp.re(arg) > 220:
            break
        term = (2 * (mp.pi ** 2) * n4 * e9u - 3 * mp.pi * n2 * e5u) * mp.exp(-arg)
        total += term
        if n > 8 and abs(term) < mp.mpf("1e-56"):
            break
    return total


def xi_completed_mp(s):
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def Xi_z_mp(z):
    return xi_completed_mp(mp.mpf("0.5") + 1j * z)


def verify_section0_canonical_wronskian():
    print("\n" + "=" * 82)
    print("SECTION 0: Canonical Xi(z) and Transverse Wronskian W(d,t) Normalization")
    print("=" * 82)

    h_half = mp.mpf("2.6") / 160
    phi_cache = {k: phi_kernel_mp(k * h_half) for k in range(-324, 325)}
    step = 2 * h_half

    d_val, t_val = mp.mpf("0.20"), mp.mpf("15.58570859")
    eps = mp.mpf("1e-18")
    w_exact = (abs(Xi_z_mp(t_val - 1j * (d_val + eps))) ** 2
               - abs(Xi_z_mp(t_val - 1j * (d_val - eps))) ** 2) / (4 * eps)
    cos_cache = {m: mp.cos(2 * t_val * (m * step)) for m in range(-160, 161)}
    w_int = mp.mpf(0)
    for j in range(0, 161):
        wj = mp.mpf("0.5") if (j == 0 or j == 160) else mp.mpf("1.0")
        yj = j * step
        iy = sum(
            (mp.mpf("0.5") if abs(m) == 160 else mp.mpf("1.0"))
            * phi_cache[j + m] * phi_cache[j - m] * cos_cache[m]
            for m in range(-160, 161)
        ) * step
        w_int += wj * yj * mp.sinh(2 * d_val * yj) * iy
    w_int *= 16 * step
    res = abs(w_exact - w_int)
    print(f"  W(0.20, 15.5857): exact={mp.nstr(w_exact, 16)}, 16*int={mp.nstr(w_int, 16)}, |diff|={mp.nstr(res, 4)}")
    assert res < mp.mpf("1e-35")


def verify_section1_cubic_window_f2_and_kronecker():
    print("\n" + "=" * 82)
    print("SECTION 1: Cubic-Window Identity, F_2(-1 +- i) = 0, and Kronecker Alignment")
    print("=" * 82)

    # 1A. Combinatorial Dirichlet identity on 1 <= N < (floor(X) + 1)^3
    for x_cut in [5, 11]:
        n_lim = (x_cut + 1) ** 3 - 1
        primes_all = sieve_primes(n_lim)
        primes_gt_x = [p for p in primes_all if p > x_cut]
        primes_le_x = set(p for p in primes_all if p <= x_cut)

        # c_n = Lambda(n)/ln(n) for n = p^k with p > X, rational exact via Fraction-like scaling
        # Notice for N < (X+1)^3, k <= 2, so 2 * c_n is an INTEGER (2 for primes, 1 for p^2)!
        two_c = [0] * (n_lim + 1)
        for p in primes_gt_x:
            pk = p
            k_exp = 1
            while pk <= n_lim:
                assert k_exp <= 2
                two_c[pk] = 2 // k_exp
                pk *= p
                k_exp += 1

        # Compute 8 * (delta_1 + c_N + (1/2) sum_{d|N} c_d c_{N/d}) vs 8 * b_N^E(X)
        mismatches = 0
        for n in range(1, n_lim + 1):
            m = n
            b_true = 1
            for p in primes_le_x:
                if m % p == 0:
                    b_true = 0
                    break
            lhs_8 = 8 * b_true
            rhs_8 = (8 if n == 1 else 0) + 4 * two_c[n]
            conv_4 = 0
            r = int(math.isqrt(n))
            for d in range(2, r + 1):
                if n % d == 0:
                    q = n // d
                    if d == q:
                        conv_4 += two_c[d] * two_c[d]
                    else:
                        conv_4 += 2 * two_c[d] * two_c[q]
            rhs_8 += conv_4
            if lhs_8 != rhs_8:
                mismatches += 1
        print(f"  [Check 1A] X={x_cut:2d}, cubic window 1 <= N <= {n_lim:4d}: exact mismatches = {mismatches}")
        assert mismatches == 0

    # 1B. Exact quadratic algebra of F_2(z) = 1 + z + z^2 / 2 AND complex zeros z = -1 +- i
    f2 = lambda z: 1 + z + mp.mpf("0.5") * (z ** 2)
    z_plus = mp.mpc(-1, 1)
    z_minus = mp.mpc(-1, -1)
    assert abs(f2(z_plus)) == 0 and abs(f2(z_minus)) == 0
    u_min = mp.mpf("-1.0")
    assert abs(abs(f2(u_min)) ** 2 - mp.mpf("0.25")) < mp.mpf("1e-49")
    print("  [Check 1B] F_2(-1) = 0.5 (|F_2(U,0)|^2 >= 1/4 on R), while F_2(-1 + i) = F_2(-1 - i) = 0 in C")

    # 1C. Refutation of Kronecker "no anti-phase" claim:
    # Take sigma = 0.60 (d = 0.10 in strip), X = 2, and only the first 15 odd primes (3..53).
    # Scan t in [0, 20000] to show finite-prime alignment easily drives Re P_{X<p<=53}(sigma+it) < -1.0
    # and |F_2(P)|^2 well below the real-axis floor 1/4!
    odd_primes = [p for p in sieve_primes(67) if p > 2]  # 18 primes: 3..67
    amps = np.array([p ** (-0.60) for p in odd_primes])
    logs = np.array([math.log(p) for p in odd_primes])
    t_kr_scan = np.linspace(10.0, 50000.0, 200000)
    phases = -np.outer(t_kr_scan, logs)
    p_vals = (np.cos(phases) + 1j * np.sin(phases)) @ amps
    f2_mod_sq = np.abs(1.0 + p_vals + 0.5 * (p_vals ** 2)) ** 2
    min_u_idx = int(np.argmin(np.real(p_vals)))
    min_f2_idx = int(np.argmin(f2_mod_sq))
    min_u = float(np.real(p_vals[min_u_idx]))
    t_min_u = float(t_kr_scan[min_u_idx])
    min_f2 = float(f2_mod_sq[min_f2_idx])
    t_min_f2 = float(t_kr_scan[min_f2_idx])
    p_at_min_f2 = complex(p_vals[min_f2_idx])
    print(f"  [Check 1C] Finite 18-prime sum (3 <= p <= 67) at sigma = 0.60:")
    print(f"    min Re P(0.60 + i t) on [10, 5e4] = {min_u:.5f} < -1 at t = {t_min_u:.2f}")
    print(f"    min |F_2(P)|^2       on [10, 5e4] = {min_f2:.5f} << 0.25 at t = {t_min_f2:.2f} "
          f"(P = {p_at_min_f2.real:+.4f}{p_at_min_f2.imag:+.4f}i near -1 +- i)")
    assert min_u < -1.0 and min_f2 < 0.05

    return min_u, t_min_u, min_f2, t_min_f2, p_at_min_f2


def verify_section2_flux_defect_and_polefree_perturbation():
    print("\n" + "=" * 82)
    print("SECTION 2: Poisson-Hadamard Flux Defect & Pole-Free Perturbed Gram Matrix")
    print("=" * 82)

    # 2A. Exact closed-form integral F_-(d, d_*) = -arcsin(d / d_*) and step pi * Theta(d - d_*)
    d_star = mp.mpf("0.20")
    flux_rows = []
    for d_str in ["0.05", "0.10", "0.15", "0.19", "0.20", "0.25", "0.35"]:
        d_val = mp.mpf(d_str)
        if d_val < d_star:
            u0 = mp.sqrt(d_star ** 2 - d_val ** 2)
            f_minus_cf = -mp.asin(d_val / d_star)
            k_fn = lambda u: d_val * (d_val ** 2 - d_star ** 2 + u ** 2) / (
                ((d_val - d_star) ** 2 + u ** 2) * ((d_val + d_star) ** 2 + u ** 2)
            )
            f_minus_q = 2 * mp.quad(k_fn, [0, u0])
            assert abs(f_minus_cf - f_minus_q) < mp.mpf("1e-46")
            f_plus = +mp.asin(d_val / d_star)
            f_tot = mp.mpf(0)
        elif d_val == d_star:
            u0 = mp.mpf(0)
            f_minus_cf = mp.mpf(0)
            f_plus = mp.pi / 2
            f_tot = mp.pi / 2
        else:
            u0 = mp.mpf(0)
            f_minus_cf = mp.mpf(0)
            f_plus = mp.pi
            f_tot = mp.pi
        print(f"  [Check 2A] d={d_str}, d_*=0.20: u_0={mp.nstr(u0,7)}, "
              f"F_-={mp.nstr(f_minus_cf,8)}, F_tot={mp.nstr(f_tot,8)}")
        flux_rows.append((float(d_val), float(u0), float(f_minus_cf), float(f_plus), float(f_tot)))

    # 2B. POLE-FREE entire even two-zero coalescence perturbation:
    # Replace simple zeros +-gamma_1, +-gamma_2 of Xi(z) by +-t_* +- i d_* with t_* = (gamma_1 + gamma_2)/2.
    gammas = [mp.im(mp.zetazero(k)) for k in range(1, 8)]
    g1, g2 = gammas[0], gammas[1]
    t_star = mp.mpf("0.5") * (g1 + g2)
    d_pert = mp.mpf("0.15")
    z_star = t_star - 1j * d_pert
    z_star_bar = t_star + 1j * d_pert

    dXi_g1 = mp.diff(Xi_z_mp, g1)
    dXi_g2 = mp.diff(Xi_z_mp, g2)

    def Xi_pert(z):
        z = mp.mpc(z)
        num = (z ** 2 - z_star ** 2) * (z ** 2 - z_star_bar ** 2)
        if abs(z - g1) < mp.mpf("1e-12"):
            # Exact removable-singularity Taylor expansion at z = +g1
            h = z - g1
            d2Xi_g1 = mp.diff(Xi_z_mp, g1, 2)
            xi_over_lin = dXi_g1 + mp.mpf("0.5") * d2Xi_g1 * h
            return xi_over_lin * num / ((z + g1) * (z ** 2 - g2 ** 2))
        if abs(z - g2) < mp.mpf("1e-12"):
            h = z - g2
            d2Xi_g2 = mp.diff(Xi_z_mp, g2, 2)
            xi_over_lin = dXi_g2 + mp.mpf("0.5") * d2Xi_g2 * h
            return xi_over_lin * num / ((z ** 2 - g1 ** 2) * (z + g2))
        return Xi_z_mp(z) * num / ((z ** 2 - g1 ** 2) * (z ** 2 - g2 ** 2))

    # Certify removable singularity continuity at z = gamma_1 and z = gamma_2 (ZERO POLE!)
    val_at_g1 = Xi_pert(g1)
    val_near_g1 = Xi_pert(g1 + mp.mpf("1e-10"))
    cont_err_g1 = abs(val_at_g1 - val_near_g1)
    val_at_zstar = Xi_pert(z_star)
    print(f"  [Check 2B] Pole-free continuity at z=gamma_1: Xi_pert(gamma_1)={mp.nstr(mp.re(val_at_g1), 10)}, "
          f"|Xi_pert(g1) - Xi_pert(g1+1e-10)|={mp.nstr(cont_err_g1, 4)}, |Xi_pert(t_*-id_*)|={mp.nstr(abs(val_at_zstar), 4)}")
    assert cont_err_g1 < mp.mpf("1e-9") and abs(val_at_zstar) < mp.mpf("1e-45")

    # Compute 6x6 de Branges Gram matrix G_alpha(d) with alpha = 1 at heights [t_*, gamma_3..gamma_7]
    sample_t = [t_star] + gammas[2:7]
    alpha = mp.mpf("1.0")

    def gram_matrix_stats(func_xi, d_v):
        m_sz = len(sample_t)
        z_pts = [tk - 1j * d_v for tk in sample_t]
        f_vals = [func_xi(zp) for zp in z_pts]
        h_step = mp.mpf("1e-16")
        df_vals = [(func_xi(zp + h_step) - func_xi(zp - h_step)) / (2 * h_step) for zp in z_pts]
        g_np = np.zeros((m_sz, m_sz), dtype=np.complex128)
        for j in range(m_sz):
            for k in range(m_sz):
                denom = z_pts[j] - mp.conj(z_pts[k])
                ker = (alpha / mp.pi) * (f_vals[j] * mp.conj(df_vals[k]) - df_vals[j] * mp.conj(f_vals[k])) / denom
                g_np[j, k] = complex(float(mp.re(ker)), float(mp.im(ker)))
        g_np = 0.5 * (g_np + g_np.conj().T)
        eigs = np.linalg.eigvalsh(g_np)
        nu_minus = int(np.sum(eigs < -1e-25))
        return float(np.real(g_np[0, 0])), float(eigs[0]), nu_minus

    gram_rows = []
    for d_str in ["0.04", "0.10", "0.14", "0.22"]:
        d_v = mp.mpf(d_str)
        p_true, lmin_true, nu_true = gram_matrix_stats(Xi_z_mp, d_v)
        p_pert, lmin_pert, nu_pert = gram_matrix_stats(Xi_pert, d_v)
        print(f"  [Check 2C] d={d_str}: True Xi (G_11={p_true:+.3e}, nu_-={nu_true}) | "
              f"Perturbed d_*=0.15 (G_11={p_pert:+.3e}, l_min={lmin_pert:+.3e}, nu_-={nu_pert})")
        assert nu_true == 0
        if d_v in (mp.mpf("0.10"), mp.mpf("0.14")):
            assert p_pert < 0.0 and nu_pert >= 1
        if d_v > d_pert:
            assert p_pert > 0.0 and nu_pert == 0
        gram_rows.append((float(d_v), p_true, nu_true, p_pert, lmin_pert, nu_pert))

    return flux_rows, gram_rows


def verify_section3_nodal_walls_and_phi_zero():
    print("\n" + "=" * 82)
    print("SECTION 3: Conditional Nodal-Wall Identity & Explicit Zero Phi(i*eta_0) = 0")
    print("=" * 82)

    # 3A. First 6 critical-line extrema tau_1..tau_6 of Xi(t) inside (gamma_k, gamma_{k+1})
    gammas = [mp.im(mp.zetazero(k)) for k in range(1, 8)]
    dXi_scaled = lambda t: mp.exp(mp.pi * t / 4) * mp.re(mp.diff(Xi_z_mp, t))
    wall_rows = []
    for k in range(6):
        tau0 = mp.findroot(dXi_scaled, (gammas[k] + mp.mpf("0.15"), gammas[k + 1] - mp.mpf("0.15")), solver="anderson")
        assert gammas[k] < tau0 < gammas[k + 1]
        xi0 = mp.re(Xi_z_mp(tau0))
        d2xi0 = mp.re(mp.diff(Xi_z_mp, tau0, 2))
        w_tilde_0 = -xi0 * d2xi0
        # Track tau_k(d) at d = 0.25
        d_test = mp.mpf("0.25")
        tau_d = mp.findroot(lambda t: mp.im(Xi_z_mp(t - 1j * d_test)), tau0)
        z_d = tau_d - 1j * d_test
        xi_d = Xi_z_mp(z_d)
        dxi_d = mp.diff(Xi_z_mp, z_d)
        # Note Xi(t - i d) = R + i I => Xi'(t - i d) = d_t R + i d_t I, and d_d R = d_t I, d_d I = -d_t R
        r_val = mp.re(xi_d)
        dt_i_val = mp.im(dxi_d)
        w_tilde_cr = r_val * dt_i_val / d_test
        eps = mp.mpf("1e-18")
        w_tilde_fd = (abs(Xi_z_mp(tau_d - 1j * (d_test + eps))) ** 2
                      - abs(Xi_z_mp(tau_d - 1j * (d_test - eps))) ** 2) / (4 * eps * d_test)
        cr_err = abs(w_tilde_cr - w_tilde_fd)
        assert cr_err < mp.mpf("1e-30") and w_tilde_0 > 0 and w_tilde_cr > 0
        print(f"  [Check 3A] k={k+1}: tau_k(0)={mp.nstr(tau0, 11)}, "
              f"-Xi*Xi''={mp.nstr(w_tilde_0, 5)}, W_tilde(0.25,tau_k(0.25))={mp.nstr(w_tilde_cr, 5)}")
        wall_rows.append((k + 1, float(tau0), float(w_tilde_0), float(tau_d), float(w_tilde_cr)))

    # 3B. Explicit polynomial counterexample showing why wall positivity CANNOT prevent an off-line zero unconditionally:
    # P(z) = (z^2 - (tau_0 - i d_*)^2)(z^2 - (tau_0 + i d_*)^2) with tau_0 = 10, d_* = 0.5:
    # At t = tau_0, Im P(tau_0 - i d) vanishes near... wait: take symmetric P(z) = (z - tau_0)^2 + d_*^2 (local model around tau_0):
    # Along t = tau_0, Im P(tau_0 - i d) = Im((-i d)^2 + d_*^2) = Im(d_*^2 - d^2) = 0 IDENTICALLY for all d >= 0!
    # At d -> 0^+, -P(tau_0) P''(tau_0) = -(d_*^2)(2) = -2 d_*^2 < 0 (a dimple!).
    # What if -P(tau_0) P''(tau_0) > 0 at d=0, yet R vanishes on Gamma?
    # Take Q(z) = ((z - tau_0)^2 + d_*^2) * (1 - 2 * (z - tau_0)^2 / delta^2) with delta < d_*:
    # At z = tau_0 - i d, Q(tau_0 - i d) = (d_*^2 - d^2)(1 + 2 d^2 / delta^2) is PURELY REAL (Im Q = 0 for all d >= 0!),
    # has Q(tau_0) = d_*^2 > 0, Q''(tau_0) = 2 - 4 d_*^2 / delta^2 < 0 (so -Q(tau_0)Q''(tau_0) > 0 at d=0!),
    # d_t Im Q(tau_0 - i d) = Q'(tau_0 - i d)/i... wait: let's check Q'(tau_0 - i d) = -i d [ ... ] so d_t I = d_d R changes sign at d_peak < d_* AND R(d_*, tau_0) = 0 at d = d_*!
    tau0_ce, ds_ce, del_ce = mp.mpf("10.0"), mp.mpf("0.50"), mp.mpf("0.40")
    Q_ce = lambda z: ((z - tau0_ce) ** 2 + ds_ce ** 2) * (1 - 2 * ((z - tau0_ce) ** 2) / (del_ce ** 2))
    q0 = Q_ce(tau0_ce)
    q2 = mp.diff(Q_ce, tau0_ce, 2)
    assert -q0 * q2 > 0 and Q_ce(tau0_ce - 1j * ds_ce) == 0
    print(f"  [Check 3B] Explicit counterexample Q(z): I(d,tau_0)=0 for all d, -Q(tau_0)Q''(tau_0)={mp.nstr(-q0*q2,5)} > 0, "
          f"yet Q(tau_0 - i*{ds_ce}) = {mp.nstr(Q_ce(tau0_ce - 1j*ds_ce), 4)} (off-line zero ON the nodal curve!)")

    # 3C. Contrast Jacobi's theta_3(0 | i e^{4u}) != 0 on |Im(u)| < pi/8 with explicit root Phi(i * eta_0) = 0!
    eta0 = mp.findroot(lambda eta: mp.re(phi_kernel_mp(1j * eta, 500)), (mp.mpf("0.31"), mp.mpf("0.325")))
    pi_over_8 = mp.pi / 8
    # Evaluate Jacobi theta_3(0 | i e^{4 i eta0}) = 1 + 2 sum_{n>=1} exp(-pi n^2 exp(4 i eta0))
    th3_at_eta0 = 1 + 2 * sum(mp.exp(-mp.pi * (n ** 2) * mp.exp(4j * eta0)) for n in range(1, 500))
    print(f"  [Check 3C] Explicit root of Phi(i*eta_0) = 0 inside |Im(u)| < pi/8:")
    print(f"    eta_0 = {mp.nstr(eta0, 25)} < pi/8 = {mp.nstr(pi_over_8, 25)}")
    print(f"    |Phi(i*eta_0)| = {mp.nstr(abs(phi_kernel_mp(1j * eta0, 500)), 4)}, "
          f"while |theta_3(0 | i e^{{4i eta_0}})| = {mp.nstr(abs(th3_at_eta0), 8)} != 0")
    assert eta0 < pi_over_8 and abs(th3_at_eta0) > mp.mpf("0.1")

    return wall_rows, float(eta0)


def generate_figures(flux_rows, gram_rows, wall_rows, eta0_val):
    print("\n" + "=" * 82)
    print("GENERATING PUBLICATION FIGURES (FIGURES 1 & 2)")
    print("=" * 82)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.titlesize": 10.0,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.0,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
    })

    # FIGURE 1: F_2 Complex Zeros vs Real Floor, Kronecker Alignment, Poisson Lobe, and Pole-Free Perturbed W_tilde
    fig1, axes1 = plt.subplots(2, 2, figsize=(11.5, 8.2))

    # Panel 1A: Complex modulus |F_2(U + i V)| showing zeros at (-1, +-1) vs real slice V = 0
    ax = axes1[0, 0]
    u_g = np.linspace(-2.4, 1.0, 200)
    v_g = np.linspace(-1.8, 1.8, 200)
    uu, vv = np.meshgrid(u_g, v_g)
    zz = uu + 1j * vv
    f2_grid = np.abs(1.0 + zz + 0.5 * (zz ** 2))
    cs = ax.contourf(uu, vv, np.log10(np.maximum(f2_grid, 1e-3)), levels=20, cmap="viridis")
    fig1.colorbar(cs, ax=ax, label=r"$\log_{10}|F_2(U+\mathrm{i}V)|$")
    ax.axhline(0.0, color="white", ls="--", lw=1.2, label=r"Real slice $V=0$ ($|F_2|\geq 1/2$)")
    ax.plot([-1.0, -1.0], [1.0, -1.0], "ro", ms=6.0, label=r"Complex roots $(-1, \pm 1)$")
    ax.set_xlabel(r"$U = \mathrm{Re}(z)$")
    ax.set_ylabel(r"$V = \mathrm{Im}(z)$")
    ax.set_title(r"A. Quadratic Polynomial $F_2(z)$: Zeros at $-1 \pm \mathrm{i}$")
    ax.legend(loc="upper right", fontsize=7.5)

    # Panel 1B: Finite-Prime Kronecker Alignment Driving |F_2(P(sigma + i t))|^2 Below 1/4
    ax = axes1[0, 1]
    odd_primes = [p for p in sieve_primes(67) if p > 2]
    amps = np.array([p ** (-0.60) for p in odd_primes])
    logs = np.array([math.log(p) for p in odd_primes])
    t_win = np.linspace(14200.0, 14500.0, 2000)
    ph = -np.outer(t_win, logs)
    pv = (np.cos(ph) + 1j * np.sin(ph)) @ amps
    f2_sq = np.abs(1.0 + pv + 0.5 * (pv ** 2)) ** 2
    ax.plot(t_win, f2_sq, "-", color="#1b4f72", lw=1.2, label=r"$|F_2(\mathcal{P}_{3\leq p\leq 67}(0.60+\mathrm{i}t))|^2$")
    ax.axhline(0.25, color="#922b21", ls="--", lw=1.4, label=r"Real-axis floor $1/4$ ($V=0$ only)")
    ax.set_yscale("log")
    ax.set_xlabel(r"Spectral height $t$ ($\sigma = 0.60$)")
    ax.set_ylabel(r"$|F_2(\mathcal{P})|^2$")
    ax.set_title(r"B. Small-Prime Kronecker Alignment Breaches $1/4$ Floor")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best", fontsize=7.6)

    # Panel 1C: Closed-Form Poisson-Hadamard Lobe F_-(d, d_*) = -arcsin(d / d_*) and Step pi * Theta(d - d_*)
    ax = axes1[1, 0]
    ds = 0.20
    d_cont = np.linspace(0.002, 0.38, 300)
    f_m = np.where(d_cont < ds, -np.arcsin(np.minimum(d_cont / ds, 1.0)), 0.0)
    f_p = np.where(d_cont < ds, +np.arcsin(np.minimum(d_cont / ds, 1.0)), np.pi)
    f_tot = f_m + f_p
    ax.plot(d_cont, f_tot, "-", color="#1b4f72", lw=2.0, label=r"Total flux $\mathcal{F}(d,d_*) = \pi\Theta(d-d_*)$")
    ax.plot(d_cont, f_m, "--", color="#922b21", lw=1.6, label=r"Negative lobe $-\arcsin(d/d_*)$")
    ax.plot(d_cont, f_p, ":", color="#117864", lw=1.6, label=r"Positive wings $\mathcal{F}_+(d,d_*)$")
    ax.axvline(ds, color="gray", ls="-.", lw=1.0, label=r"$d_* = 0.20$")
    ax.set_xlabel(r"Transverse offset $d$")
    ax.set_ylabel(r"Integrated Poisson flux")
    ax.set_title(r"C. Hadamard–Poisson Flux Step (Lagarias Contrapositive)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center right", fontsize=7.5)

    # Panel 1D: Pole-Free Perturbed De Branges Pivot G_11(d) vs True Xi Pivot
    ax = axes1[1, 1]
    d_g = np.linspace(0.02, 0.28, 35)
    g1 = float(mp.im(mp.zetazero(1)))
    g2 = float(mp.im(mp.zetazero(2)))
    ts = 0.5 * (g1 + g2)
    dp = 0.15
    mp.mp.dps = 20
    piv_true, piv_pert = [], []
    for dv in d_g:
        dm = mp.mpf(float(dv))
        ep = mp.mpf("1e-10")
        wt_true = (abs(Xi_z_mp(mp.mpf(ts) - 1j * (dm + ep))) ** 2
                   - abs(Xi_z_mp(mp.mpf(ts) - 1j * (dm - ep))) ** 2) / (4 * ep * dm * mp.pi)
        piv_true.append(float(wt_true))
        zp_f = lambda d_sub: (
            Xi_z_mp(mp.mpf(ts) - 1j * d_sub)
            * (((mp.mpf(ts) - 1j * d_sub) ** 2 - (mp.mpf(ts) - 1j * mp.mpf(dp)) ** 2)
               * ((mp.mpf(ts) - 1j * d_sub) ** 2 - (mp.mpf(ts) + 1j * mp.mpf(dp)) ** 2))
            / (((mp.mpf(ts) - 1j * d_sub) ** 2 - mp.mpf(g1) ** 2)
               * ((mp.mpf(ts) - 1j * d_sub) ** 2 - mp.mpf(g2) ** 2))
        )
        wt_pert = (abs(zp_f(dm + ep)) ** 2 - abs(zp_f(dm - ep)) ** 2) / (4 * ep * dm * mp.pi)
        piv_pert.append(float(wt_pert))
    mp.mp.dps = 50
    # Normalize by critical-line scale for clean plotting
    s_norm = abs(piv_true[0])
    ax.plot(d_g, np.array(piv_true) / s_norm, "o-", color="#117864", ms=3.5, lw=1.5, label=r"True $\Xi(z)$ ($\nu_- = 0$)")
    ax.plot(d_g, np.array(piv_pert) / s_norm, "s--", color="#922b21", ms=3.5, lw=1.5, label=r"Pole-free $\Xi_{\mathrm{pert}}(z; 0.15)$")
    ax.axhline(0.0, color="black", lw=0.9)
    ax.axvline(dp, color="gray", ls=":", lw=1.0, label=r"$d_* = 0.15$")
    ax.set_xlabel(r"Transverse offset $d$ at $t_* = \frac{1}{2}(\gamma_1 + \gamma_2)$")
    ax.set_ylabel(r"Normalized de Branges pivot $[\mathbf{G}_1(d)]_{1,1}$")
    ax.set_title(r"D. Sign Inversion on Pole-Free Perturbed $\Xi_{\mathrm{pert}}$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=7.6)

    fig1.tight_layout()
    for out_dir in (FIG_DIR_ROOT, FIG_DIR_PAPER):
        fig1.savefig(os.path.join(out_dir, "fig1_spectral_resonance_and_hadamard_flux.pdf"))
        fig1.savefig(os.path.join(out_dir, "fig1_spectral_resonance_and_hadamard_flux.png"), dpi=220)
    plt.close(fig1)

    # FIGURE 2: Nodal Curves I(d,t) = 0, Local Wall Identity, Counterexample Q(z), and Phi(i*eta) Zero
    fig2, axes2 = plt.subplots(2, 2, figsize=(11.5, 8.2))

    # Fast u-grid evaluation of Xi(t - i d) for Panel 2A & 2B
    u_np = np.linspace(0.0, 2.5, 350)
    du = u_np[1] - u_np[0]
    w_tr = np.ones_like(u_np) * du
    w_tr[0] *= 0.5
    w_tr[-1] *= 0.5
    phi_np = np.zeros_like(u_np)
    for n in range(1, 20):
        a_n = np.pi * (n ** 2) * np.exp(4.0 * u_np)
        phi_np += (2.0 * (np.pi ** 2) * (n ** 4) * np.exp(9.0 * u_np)
                   - 3.0 * np.pi * (n ** 2) * np.exp(5.0 * u_np)) * np.exp(-a_n)
    pw = phi_np * w_tr

    # Panel 2A: Nodal curves I(d,t) = 0 and R(d,t) = 0 on t in [12, 39], d in [0, 0.5]
    ax = axes2[0, 0]
    t_2d = np.linspace(12.0, 39.0, 260)
    d_2d = np.linspace(0.0, 0.50, 90)
    cos_ut = np.cos(2.0 * np.outer(t_2d, u_np))
    sin_ut = np.sin(2.0 * np.outer(t_2d, u_np))
    ch_ud = np.cosh(2.0 * np.outer(d_2d, u_np))
    sh_ud = np.sinh(2.0 * np.outer(d_2d, u_np))
    r_2d = 8.0 * ((ch_ud * pw[None, :]) @ cos_ut.T)
    i_2d = 8.0 * ((sh_ud * pw[None, :]) @ sin_ut.T)
    # Note at d=0, I(0,t)=0 trivially; plot I(d,t)/sinh(max(d,1e-4)) to isolate transverse branches Gamma_k!
    d_safe = np.maximum(d_2d[:, None], 1e-4)
    i_norm = i_2d / d_safe
    # Fill for d=0 using limit 16 int u Phi(u) sin(2 t u) du = -Xi'(t)
    i_norm[0, :] = 16.0 * (sin_ut @ (u_np * pw))
    tt_m, dd_m = np.meshgrid(t_2d, d_2d)
    ax.contour(tt_m, dd_m, r_2d, levels=[0.0], colors="#1b4f72", linewidths=1.7)
    ax.contour(tt_m, dd_m, i_norm, levels=[0.0], colors="#922b21", linewidths=1.7, linestyles="--")
    ax.plot([], [], "-", color="#1b4f72", lw=1.7, label=r"$R(d,t) = 0$ (emanating from $\gamma_k$)")
    ax.plot([], [], "--", color="#922b21", lw=1.7, label=r"$I(d,t)/d = 0$ (emanating from $\tau_k(0)$)")
    ax.set_xlabel(r"Spectral height $t \in [12, 39]$")
    ax.set_ylabel(r"Transverse offset $d$")
    ax.set_title(r"A. Nodal Curves of $\Xi(t - \mathrm{i}d)$ on $t \in [12, 39]$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=7.6)

    # Panel 2B: Rescaled W_tilde(d, tau_k(d)) along first 4 finite branches Gamma_1..Gamma_4
    ax = axes2[0, 1]
    d_tr = np.linspace(0.01, 0.45, 30)
    colors_w = ["#1b4f72", "#117864", "#b9770e", "#922b21"]
    for k_idx in range(4):
        tau_curr = wall_rows[k_idx][1]
        wt_scaled = []
        for dv in d_tr:
            # Fast Newton step on i_norm using u_np
            for _ in range(5):
                f_v = 8.0 * np.sum(pw * np.sinh(2.0 * dv * u_np) * np.sin(2.0 * tau_curr * u_np))
                df_v = 16.0 * np.sum(pw * u_np * np.sinh(2.0 * dv * u_np) * np.cos(2.0 * tau_curr * u_np))
                tau_curr -= float(f_v / df_v)
            r_v = 8.0 * np.sum(pw * np.cosh(2.0 * dv * u_np) * np.cos(2.0 * tau_curr * u_np))
            dti_v = 16.0 * np.sum(pw * u_np * np.sinh(2.0 * dv * u_np) * np.cos(2.0 * tau_curr * u_np))
            wt_v = (r_v * dti_v / dv) * math.exp(0.5 * math.pi * tau_curr)
            wt_scaled.append(float(wt_v))
        ax.semilogy(d_tr, wt_scaled, "-", color=colors_w[k_idx], lw=1.6, label=rf"Branch $\Gamma_{{{k_idx+1}}}$ ($\tau_{{{k_idx+1}}}(0)={wall_rows[k_idx][1]:.2f}$)")
    ax.set_xlabel(r"Transverse offset $d$")
    ax.set_ylabel(r"$\mathrm{e}^{\pi\tau_k(d)/2}\,\widetilde{W}(d, \tau_k(d))$")
    ax.set_title(r"B. Finite-Height Nodal Profiles $\Gamma_1, \dots, \Gamma_4$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best", fontsize=7.6)

    # Panel 2C: Explicit Counterexample Q(z) Showing Sign Change Along I = 0 Despite -Q(tau_0)Q''(tau_0) > 0
    ax = axes2[1, 0]
    d_ce = np.linspace(0.005, 0.62, 200)
    ds_c, del_c = 0.50, 0.40
    r_ce = ((ds_c ** 2) - (d_ce ** 2)) * (1.0 + 2.0 * (d_ce ** 2) / (del_c ** 2))
    dr_ce = np.gradient(r_ce, d_ce)
    wt_ce = r_ce * dr_ce / d_ce
    ax.plot(d_ce, r_ce, "--", color="#1b4f72", lw=1.6, label=r"$R_Q(d,\tau_0)$ (vanishes at $d_*=0.50$)")
    ax.plot(d_ce, wt_ce, "-", color="#922b21", lw=1.8, label=r"$\widetilde{W}_Q(d,\tau_0) = R_Q\partial_d R_Q / d$")
    ax.axhline(0.0, color="black", lw=0.9)
    ax.axvline(ds_c, color="gray", ls=":", lw=1.0, label=r"Off-line zero $d_* = 0.50$")
    ax.set_xlabel(r"Transverse offset $d$ along $I_Q(d, \tau_0) \equiv 0$")
    ax.set_ylabel(r"Amplitude / Wronskian")
    ax.set_title(r"C. Counterexample $Q(z)$: Zero on $\Gamma$ Breaks Wall Positivity")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=7.4)

    # Panel 2D: Contrast Between |theta_3(0 | i e^{4 i eta})| > 0 and Phi(i eta) Vanishing at eta_0 = 0.319415 < pi/8
    ax = axes2[1, 1]
    etas = np.linspace(0.0, 0.365, 160)
    phi_im = []
    th3_im = []
    mp.mp.dps = 25
    for ev in etas:
        em = mp.mpf(float(ev))
        phi_im.append(float(mp.re(phi_kernel_mp(1j * em, 250))))
        th_v = 1 + 2 * sum(mp.exp(-mp.pi * (n ** 2) * mp.exp(4j * em)) for n in range(1, 250))
        th3_im.append(float(abs(th_v)))
    mp.mp.dps = 50
    ax.plot(etas, phi_im, "-", color="#922b21", lw=1.8, label=r"$\Phi(\mathrm{i}\eta)$ (vanishes at $\eta_0 \approx 0.31942$)")
    ax.plot(etas, th3_im, "--", color="#117864", lw=1.6, label=r"$|\vartheta_3(0 \mid \mathrm{i}\mathrm{e}^{4\mathrm{i}\eta})| > 0$")
    ax.axhline(0.0, color="black", lw=0.9)
    ax.axvline(eta0_val, color="#922b21", ls=":", lw=1.2, label=rf"$\eta_0 = {eta0_val:.5f} < \pi/8$")
    ax.set_ylim(-5.0, 5.0)
    ax.set_xlabel(r"Imaginary kernel coordinate $\eta = \mathrm{Im}(u) \in [0, 0.365] \subset [0, \pi/8)$")
    ax.set_ylabel(r"Kernel / theta value")
    ax.set_title(r"D. $\Phi(\mathrm{i}\eta_0) = 0$ Inside Strip $|\mathrm{Im}(u)| < \pi/8$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=7.4)

    fig2.tight_layout()
    for out_dir in (FIG_DIR_ROOT, FIG_DIR_PAPER):
        fig2.savefig(os.path.join(out_dir, "fig2_chladni_walls_and_theta_strip.pdf"))
        fig2.savefig(os.path.join(out_dir, "fig2_chladni_walls_and_theta_strip.png"), dpi=220)
    plt.close(fig2)

    for fn in ("fig1_spectral_resonance_and_hadamard_flux.pdf", "fig2_chladni_walls_and_theta_strip.pdf"):
        sz = os.path.getsize(os.path.join(FIG_DIR_PAPER, fn)) / 1024.0
        print(f"  [Saved] {fn}: {sz:.1f} KB (< 500 KB)")
        assert sz < 500.0


def main():
    t0 = time.time()
    print("=" * 82)
    print("PAPER X VERIFICATION SUITE: COMBINATORIAL IDENTITIES, FLUX DEFECT & OBSTRUCTIONS")
    print("=" * 82)
    verify_section0_canonical_wronskian()
    verify_section1_cubic_window_f2_and_kronecker()
    flux_rows, gram_rows = verify_section2_flux_defect_and_polefree_perturbation()
    wall_rows, eta0_val = verify_section3_nodal_walls_and_phi_zero()
    generate_figures(flux_rows, gram_rows, wall_rows, eta0_val)
    print("\n" + "=" * 82)
    print(f"ALL PAPER X VERIFICATIONS COMPLETED IN {time.time() - t0:.2f}s (EXIT CODE 0)")
    print("=" * 82)


if __name__ == "__main__":
    main()
