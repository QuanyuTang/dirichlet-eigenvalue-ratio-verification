#!/usr/bin/env python3
"""
Exact verification of the arithmetic certificates in

    Counterexamples to a higher-index Dirichlet eigenvalue-ratio conjecture

This script verifies, using exact symbolic/rational arithmetic:

1. The Rayleigh sums in Table "Rayleigh sums used in the proof";
2. The exact integrals for the trial functions f_1 and f_3;
3. The resulting Rayleigh quotient formulas used in the proof;
4. The four elementary constant inequalities, through the same rigorous
   rational certificates used in the manuscript;
5. Several downstream rational comparisons appearing in the Bessel-zero,
   Hardy, and annular-ratio estimates.

No floating-point arithmetic is used.

Requirement:
    sympy

Run:
    python verify_exact_arithmetic.py
"""

import sympy as sp


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def check(name, condition):
    """Raise an error unless condition is exactly true; otherwise print PASS."""
    if condition is not True and condition != sp.true:
        raise AssertionError(f"[FAIL] {name}")
    print(f"[PASS] {name}")


def check_equal(name, lhs, rhs):
    """Check exact symbolic equality."""
    diff = sp.simplify(lhs - rhs)
    if diff != 0:
        raise AssertionError(
            f"[FAIL] {name}\n"
            f"  lhs = {lhs}\n"
            f"  rhs = {rhs}\n"
            f"  lhs-rhs = {diff}"
        )
    print(f"[PASS] {name}")


Q = sp.Rational


# ---------------------------------------------------------------------------
# Part 1. Rayleigh sums
# ---------------------------------------------------------------------------

def rayleigh_sum(nu, p):
    r"""
    Compute sigma_p(nu) recursively from

        sigma_1(nu) = 1 / (4 (nu+1)),

        sigma_p(nu)
          = 1/(nu+p) * sum_{r=1}^{p-1} sigma_r(nu) sigma_{p-r}(nu),
          p >= 2.

    Everything is exact.
    """
    sigma = {1: Q(1, 4 * (nu + 1))}
    for k in range(2, p + 1):
        sigma[k] = sp.simplify(
            sum(sigma[r] * sigma[k-r] for r in range(1, k))
            / Q(nu + k, 1)
        )
    return sigma[p]


def verify_rayleigh_sums():
    print("\n=== Part 1: Rayleigh sums ===")

    expected = {
        (0, 2): Q(1, 32),
        (1, 3): Q(1, 3072),
        (1, 5): Q(13, 8847360),
        (1, 6): Q(11, 110100480),
        (3, 5): Q(1, 110100480),
    }

    for (nu, p), target in expected.items():
        actual = rayleigh_sum(nu, p)
        check_equal(
            rf"sigma_{p}({nu}) = {target}",
            actual,
            target
        )

    # Downstream exact arithmetic used in the Bessel-zero certificates.
    check(
        "3072 > 14^3",
        sp.Integer(3072) > sp.Integer(14)**3
    )

    j3_certificate = (
        sp.Integer(110100480) * sp.Integer(25)**5
        - sp.Integer(1014)**5
    )
    check_equal(
        "110100480*25^5 - 1014^5 = 3212367382176",
        j3_certificate,
        sp.Integer(3212367382176)
    )
    check(
        "110100480*25^5 > 1014^5",
        j3_certificate > 0
    )

    ratio = sp.simplify(rayleigh_sum(1, 5) / rayleigh_sum(1, 6))
    check_equal(
        "sigma_5(1)/sigma_6(1) = 1456/99",
        ratio,
        Q(1456, 99)
    )

    gamma_certificate = (
        sp.Integer(13)**2 * sp.Integer(99)**2 * sp.Integer(32)
        - sp.Integer(5)**2 * sp.Integer(1456)**2
    )
    check_equal(
        "13^2*99^2*32 - 5^2*1456^2 = 5408",
        gamma_certificate,
        sp.Integer(5408)
    )
    check(
        "1456/(99*sqrt(32)) < 13/5 (certified after squaring)",
        gamma_certificate > 0
    )


# ---------------------------------------------------------------------------
# Part 2. Exact trial-function integrals
# ---------------------------------------------------------------------------

def exact_integrals(f, r, a):
    I0 = sp.integrate(r * f**2, (r, a, 1))
    I1 = sp.integrate(r * sp.diff(f, r)**2, (r, a, 1))
    Im1 = sp.integrate(f**2 / r, (r, a, 1))
    return tuple(sp.simplify(x) for x in (I0, I1, Im1))


def verify_trial_integrals():
    print("\n=== Part 2: Trial-function integrals ===")

    r = sp.symbols("r", positive=True)
    a = Q(1, 10)

    f1 = (r - a) * (1 - r) * (1 - r / 2)
    f3 = (r - a) * (1 - r) * (1 - 7*r + 4*r**2)

    I0_f1, I1_f1, Im1_f1 = exact_integrals(f1, r, a)
    I0_f3, I1_f3, Im1_f3 = exact_integrals(f3, r, a)

    # Exact values printed in the appendix.
    check_equal(
        "I_0(f_1)",
        I0_f1,
        Q(23993577, 4480000000)
    )
    check_equal(
        "I_1(f_1)",
        I1_f1,
        Q(2401569, 40000000)
    )
    check_equal(
        "I_{-1}(f_1)",
        Im1_f1,
        sp.log(10) / 100 + Q(27873, 80000000)
    )

    check_equal(
        "I_0(f_3)",
        I0_f3,
        Q(22063678533, 700000000000)
    )
    check_equal(
        "I_1(f_3)",
        I1_f3,
        Q(187046577, 350000000)
    )
    check_equal(
        "I_{-1}(f_3)",
        Im1_f3,
        sp.log(10) / 100 + Q(42655833, 700000000)
    )

    # Verify the two quotient formulas used earlier in the paper:
    #
    # Q_n[f] = (I_1(f) + n^2 I_{-1}(f)) / I_0(f).

    Q1 = sp.simplify((I1_f1 + Im1_f1) / I0_f1)
    Q3 = sp.simplify((I1_f3 + 9*Im1_f3) / I0_f3)

    target_Q1 = (
        56 * (800000 * sp.log(10) + 4831011)
        / sp.Integer(23993577)
    )
    target_Q3 = (
        7000 * (10**6 * sp.log(10) + 12031677)
        / sp.Integer(2451519837)
    )

    check_equal(
        "Q_1[f_1] = 56(800000 log(10)+4831011)/23993577",
        Q1,
        target_Q1
    )
    check_equal(
        "Q_3[f_3] = 7000(10^6 log(10)+12031677)/2451519837",
        Q3,
        target_Q3
    )

    # Downstream rational substitutions using log(10) < 2303/1000.
    log10_upper = Q(2303, 1000)

    Q1_upper = sp.simplify(Q1.subs(sp.log(10), log10_upper))
    Q3_upper = sp.simplify(Q3.subs(sp.log(10), log10_upper))

    check_equal(
        "Q_1 upper certificate",
        Q1_upper,
        Q(373711016, 23993577)
    )
    check(
        "373711016/23993577 < 78/5",
        Q1_upper < Q(78, 5)
    )

    check_equal(
        "Q_3 upper certificate",
        Q3_upper,
        Q(100342739000, 2451519837)
    )
    check(
        "100342739000/2451519837 < 41",
        Q3_upper < 41
    )


# ---------------------------------------------------------------------------
# Part 3. Elementary constant inequalities
# ---------------------------------------------------------------------------

def verify_elementary_constants():
    print("\n=== Part 3: Elementary constant inequalities ===")

    # -----------------------------------------------------------------------
    # (a) log(10) < 2303/1000
    #
    # Since exp(x) = sum_{q>=0} x^q/q! and all terms are positive,
    #
    #   exp(2303/1000)
    #      > sum_{q=0}^9 (2303/1000)^q/q!.
    #
    # Thus it is enough to verify that this finite rational sum is > 10.
    # -----------------------------------------------------------------------

    x = Q(2303, 1000)
    exp_partial = sum(x**q / sp.factorial(q) for q in range(10))
    exp_margin = sp.factor(exp_partial - 10)

    check(
        "sum_{q=0}^9 (2303/1000)^q/q! > 10",
        exp_partial > 10
    )
    print(f"       exact margin = {exp_margin}")
    print("       Hence exp(2303/1000) > 10, so log(10) < 2303/1000.")

    # -----------------------------------------------------------------------
    # (b) pi > 157/50
    #
    # Machin:
    #   pi/4 = 4 atan(1/5) - atan(1/239).
    #
    # Alternating-series estimates:
    #   atan(1/5) > 1/5 - 1/(3*5^3),
    #   atan(1/239) < 1/239.
    #
    # Therefore
    #   pi/4 > 4(1/5 - 1/(3*5^3)) - 1/239.
    # -----------------------------------------------------------------------

    pi4_lower = (
        4 * (Q(1, 5) - Q(1, 3 * 5**3))
        - Q(1, 239)
    )

    check_equal(
        "Machin rational lower bound = 70369/89625",
        pi4_lower,
        Q(70369, 89625)
    )
    check(
        "70369/89625 > 157/200",
        pi4_lower > Q(157, 200)
    )
    print(
        f"       exact margin = "
        f"{sp.factor(pi4_lower - Q(157, 200))}"
    )
    print("       Hence pi/4 > 157/200, so pi > 157/50.")

    # -----------------------------------------------------------------------
    # (c) log(391/100) < 273/200
    #
    # For x>0, t=(x-1)/(x+1):
    #
    #   log x = 2 sum_{q>=0} t^(2q+1)/(2q+1).
    #
    # For x=391/100, t=291/491 > 0.
    # Bound the tail starting at t^7 by
    #
    #   sum_{q>=3} t^(2q+1)/(2q+1)
    #       <= t^7/[7(1-t^2)].
    # -----------------------------------------------------------------------

    xpos = Q(391, 100)
    tpos = sp.simplify((xpos - 1) / (xpos + 1))
    check_equal("t for x=391/100 is 291/491", tpos, Q(291, 491))

    log391_100_upper = 2 * (
        tpos
        + tpos**3 / 3
        + tpos**5 / 5
        + tpos**7 / (7 * (1 - tpos**2))
    )

    check(
        "rigorous upper bound for log(391/100) is < 273/200",
        log391_100_upper < Q(273, 200)
    )
    print(
        f"       exact margin = "
        f"{sp.factor(Q(273, 200) - log391_100_upper)}"
    )

    # -----------------------------------------------------------------------
    # (d) log(391/1000) < -187/200
    #
    # Here t=-609/1391 < 0. Every term
    #
    #   t^(2q+1)/(2q+1)
    #
    # is negative. Hence the infinite series is strictly smaller than
    # its truncation after t^5:
    #
    #   log(391/1000)
    #       < 2(t + t^3/3 + t^5/5).
    #
    # It therefore suffices to verify that this truncation is < -187/200.
    # -----------------------------------------------------------------------

    xneg = Q(391, 1000)
    tneg = sp.simplify((xneg - 1) / (xneg + 1))
    check_equal("t for x=391/1000 is -609/1391", tneg, Q(-609, 1391))

    log391_1000_trunc = 2 * (
        tneg
        + tneg**3 / 3
        + tneg**5 / 5
    )

    check(
        "2(t+t^3/3+t^5/5) < -187/200",
        log391_1000_trunc < Q(-187, 200)
    )
    print(
        f"       exact margin = "
        f"{sp.factor(Q(-187, 200) - log391_1000_trunc)}"
    )
    print(
        "       Since all omitted terms are negative, "
        "log(391/1000) is even smaller."
    )

    # -----------------------------------------------------------------------
    # Downstream Hardy constants in the manuscript.
    # -----------------------------------------------------------------------

    check("100/391 < 32/125", Q(100, 391) < Q(32, 125))
    check("1000/391 < 1279/500", Q(1000, 391) < Q(1279, 500))

    CL_upper = Q(1, 4) * (
        Q(273, 200) + Q(32, 125) - 1
    )
    CR_upper = Q(1, 4) * (
        Q(1279, 500) - Q(187, 200) - 1
    )

    check_equal("C_L rational upper bound = 621/4000", CL_upper, Q(621, 4000))
    check("621/4000 < 39/250", CL_upper < Q(39, 250))

    check_equal("C_R rational upper bound = 623/4000", CR_upper, Q(623, 4000))
    check("623/4000 < 39/250", CR_upper < Q(39, 250))

    # Downstream lower bound for mu_{0,2} using pi > 157/50:
    mu02_lower = Q(211, 250) * (Q(20, 9) * Q(157, 50))**2

    check_equal(
        "(211/250)*(20*(157/50)/9)^2 = 10401878/253125",
        mu02_lower,
        Q(10401878, 253125)
    )
    check(
        "10401878/253125 > 41",
        mu02_lower > 41
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Exact arithmetic verification")
    print("SymPy version:", sp.__version__)
    print("No floating-point arithmetic is used.")

    verify_rayleigh_sums()
    verify_trial_integrals()
    verify_elementary_constants()

    print("\n==============================================")
    print("ALL CHECKS PASSED.")
    print("==============================================")


if __name__ == "__main__":
    main()
