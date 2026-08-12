# Counterexamples to a higher-index Dirichlet eigenvalue-ratio conjecture

**Yixin He, Quanyu Tang, and Haiqi Zhang**

Preprint, first public version: August 11, 2026.

[Paper (PDF)](paper.pdf) · [LaTeX source](paper.tex) · [Exact-arithmetic verification](verify_exact_arithmetic.py) · [Sample output](output.txt)

> **Main result.** We disprove Ashbaugh's higher-index Payne–Pólya–Weinberger (PPW) conjecture for Dirichlet Laplacian eigenvalue ratios. Counterexamples exist for every integer $m\ge 3$, already in dimension $N=2$.

## The conjecture

Let

$$
0<\lambda_1(\Omega)<\lambda_2(\Omega)\le \lambda_3(\Omega)\le\cdots
$$

be the Dirichlet Laplacian eigenvalues of a bounded domain $\Omega\subset\mathbb R^N$, counted with multiplicity.

The higher-index conjecture, recorded as an open problem by Ashbaugh, asks whether for every integer $m\ge 1$,

$$
\frac{\lambda_{2m}(\Omega)}{\lambda_m(\Omega)}
\le
\frac{j_{N/2,1}^2}{j_{N/2-1,1}^2},
$$

where $j_{\nu,1}$ denotes the first positive zero of the Bessel function $J_\nu$.

For $N=2$, this becomes

$$
\frac{\lambda_{2m}(\Omega)}{\lambda_m(\Omega)}
\le
\frac{j_{1,1}^2}{j_{0,1}^2}.
$$

The cases $m=1,2$ were known. The first open case was therefore

$$
\frac{\lambda_6(\Omega)}{\lambda_3(\Omega)}
\le
\frac{j_{1,1}^2}{j_{0,1}^2}.
$$

## Main result

We prove that the conjecture is false for every $m\ge 3$.

For every integer $m\ge 3$, there exists a bounded planar domain $\Omega_m\subset\mathbb R^2$ with $C^\infty$ boundary such that

$$
\frac{\lambda_{2m}(\Omega_m)}{\lambda_m(\Omega_m)}>\frac{13}{5}>\frac{j_{1,1}^2}{j_{0,1}^2}.
$$

In particular, the first open case $m=3$ already has an explicit counterexample: the circular annulus

$$
x\in\mathbb R^2:\frac{1}{10}<|x|<1.
$$

For this annulus,

$$
\frac{\lambda_6(A)}{\lambda_3(A)}>\frac{13}{5}>\frac{j_{1,1}^2}{j_{0,1}^2}.
$$

Thus the conjectured inequality for $\lambda_6/\lambda_3$ already fails for a circular annulus.

The result also disproves, for $N=2$ and $k=2$, the more general proposed inequality

$$
\frac{\lambda_{km}(\Omega)}{\lambda_m(\Omega)}
\le
\sup_D
\frac{\lambda_k(D)}{\lambda_1(D)}.
$$

## Abstract

Let $\lambda_k(\Omega)$ denote the $k$th Dirichlet Laplacian eigenvalue of a bounded planar domain $\Omega$, with eigenvalues counted with multiplicity. We disprove a conjectured higher-index extension of the Payne–Pólya–Weinberger inequality that was recorded as an open problem by Ashbaugh. For every integer $m\ge 3$, we construct a bounded planar domain $\Omega_m$ with $C^\infty$ boundary such that

$$
\frac{\lambda_{2m}(\Omega_m)}{\lambda_m(\Omega_m)}>\frac{13}{5}>\frac{j_{1,1}^2}{j_{0,1}^2},
$$

where $j_{\nu,1}$ is the first positive zero of $J_\nu$, the Bessel function of the first kind of order $\nu$. For $m=3$, the annulus

$$
x\in\mathbb{R}^2:
\frac{1}{10}<|x|<1
$$

already provides such a counterexample.

## Manuscript

The current public manuscript is available as [`paper.pdf`](paper.pdf), with LaTeX source in [`paper.tex`](paper.tex).

For a fixed snapshot of the first public version, see **Release `preprint-v1` (August 11, 2026)**.

The manuscript was submitted to arXiv on August 10, 2026 and is currently on hold pending moderation. The arXiv identifier and link will be added here once the submission becomes publicly available.

## Exact-arithmetic verification

The proof is explicit, and all numerical comparisons used in the argument are certified using exact symbolic or rational arithmetic.

The script [`verify_exact_arithmetic.py`](verify_exact_arithmetic.py) verifies the exact computations used in the proof, including:

* the Rayleigh sums appearing in the paper;
* the exact integrals for the trial functions $f_1$ and $f_3$;
* the corresponding Rayleigh quotient formulas;
* the elementary rational bounds for $\log 10$, $\pi$, and the logarithms entering the Hardy estimate;
* the downstream rational comparisons used in the Bessel-zero, Hardy, and annular-ratio estimates.

No floating-point arithmetic is used in these verification steps.

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

A complete sample output is provided in [`output.txt`](output.txt).

A successful run ends with

```text
ALL CHECKS PASSED.
```

## Repository contents

* [`paper.pdf`](paper.pdf) — current public manuscript;
* [`paper.tex`](paper.tex) — LaTeX source of the manuscript;
* [`verify_exact_arithmetic.py`](verify_exact_arithmetic.py) — exact-arithmetic verification script;
* [`requirements.txt`](requirements.txt) — Python dependencies;
* [`output.txt`](output.txt) — complete sample output of the verification script.

## Keywords

Dirichlet eigenvalues; Dirichlet Laplacian; eigenvalue ratios; higher-index eigenvalue inequalities; Payne–Pólya–Weinberger inequality; PPW inequality; Ashbaugh conjecture; higher-index PPW conjecture; spectral geometry; counterexamples; annulus; Bessel functions.
