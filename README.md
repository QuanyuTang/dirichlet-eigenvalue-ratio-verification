# Preprint and exact-arithmetic verification code

This repository contains the manuscript and exact-arithmetic verification code accompanying the paper

**Yixin He, Quanyu Tang, and Haiqi Zhang, “Counterexamples to a higher-index Dirichlet eigenvalue-ratio conjecture.”**

## Manuscript

The current public manuscript is available as [`paper.pdf`](paper.pdf), with LaTeX source in [`paper.tex`](paper.tex).

For a fixed snapshot of the first public version, see **Release `preprint-v1` (August 11, 2026)**.

## Preprint status

The manuscript was submitted to arXiv on August 10, 2026 and is currently on hold pending moderation. The arXiv identifier and link will be added here once the submission becomes publicly available.

##

The script `verify_exact_arithmetic.py` verifies the exact computations used in the proof, including:

* the Rayleigh sums appearing in the paper;
* the exact integrals for the trial functions $f_1$ and $f_3$;
* the corresponding Rayleigh quotient formulas;
* the elementary rational bounds for $\log 10$, $\pi$, and the logarithms entering the Hardy estimate;
* the downstream rational comparisons used in the Bessel-zero, Hardy, and annular-ratio estimates.

All computations are performed using symbolic or exact rational arithmetic; no floating-point arithmetic is used.

## Requirements

* Python 3
* SymPy 1.14.0

Install the required package with

```bash
pip install -r requirements.txt
```

## Usage

Run

```bash
python verify_exact_arithmetic.py
```

A complete sample output is provided in `output.txt`. A successful run ends with

```text
ALL CHECKS PASSED.
```
