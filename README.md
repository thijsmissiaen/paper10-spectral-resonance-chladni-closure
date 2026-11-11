# Combinatorial Euler-Sieve Identities, Hadamard–Poisson Flux Detection, and Conditional Nodal-Wall Identities for Riemann's Xi-Function

**Author:** Thijs Missiaen (Independent Researcher, Belgium · ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810))  
**DOI (Zenodo / CERN):** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23645893.svg)](https://doi.org/10.5281/zenodo.23645893) [`10.5281/zenodo.23645893`](https://doi.org/10.5281/zenodo.23645893) (`v10.1.0`: [`10.5281/zenodo.24141071`](https://doi.org/10.5281/zenodo.24141071))  
**Manuscript:** [`paper/paper10_spectral_resonance_chladni_closure.pdf`](paper/paper10_spectral_resonance_chladni_closure.pdf) (AMS-LaTeX, 9 pages) · [`paper/paper10_spectral_resonance_chladni_closure.tex`](paper/paper10_spectral_resonance_chladni_closure.tex)  
**Verification Suite:** [`scripts/verify_paper10_spectral_resonance_chladni_closure.py`](scripts/verify_paper10_spectral_resonance_chladni_closure.py) (50-digit `mpmath`)  
**MSC 2020:** Primary `11M06`, `11M26`; Secondary `30D15`, `46C20`

## Overview

Let $\Xi(t - \mathrm{i}d) := \xi(\frac{1}{2} + d + \mathrm{i}t) = R(d, t) + \mathrm{i}I(d, t)$ denote Riemann's completed zeta function, with canonical transverse Wronskian

$$W(d, t) := \frac{1}{2}\frac{\partial}{\partial d}\lvert \Xi(t - \mathrm{i}d)\rvert^2 = 16\int _ 0^{\infty} y\sinh(2dy)\,I(y, t)\,\mathrm{d}y, \qquad \widetilde{W}(d, t) := \frac{W(d, t)}{d}.$$

This note establishes exact combinatorial, calculus, and Cauchy–Riemann identities across three analytical constructions associated with $\Xi$, while proving explicit counterexamples and obstructions that delineate where each construction stops short of unconditional zero-freeness.

## Unconditional Identities vs. Structural Obstructions

1. **Cubic-Window Dirichlet Identity vs. Complex Zeros of $F _ 2$, Strip Divergence, and Kronecker Alignment (Section 2).**
   - **Cubic-window exactness (Proposition 2.1):** Expanding $\zeta(s)\prod _ {p \le X}(1 - p^{-s}) = \exp(\mathcal{P} _ {>X}(s))$ as a formal Dirichlet series (or analytically on $\sigma > 1$) yields the exact quadratic convolution formula for $b _ N^{\mathcal{E}}(X)$ throughout the cubic window $1 \le N < (\lfloor X\rfloor + 1)^3$.
   - **Structural obstructions (Proposition 2.2):** While $\lvert F _ 2(U)\rvert^2 = \frac{1}{4}((1+U)^2 + 1)^2 \ge \frac{1}{4}$ holds along the real slice $V = 0$, the quadratic polynomial $F _ 2(z) = 1 + z + \frac{1}{2}z^2$ vanishes identically in $\mathbb{C}$ at $z = -1 \pm \mathrm{i}$; the prime Dirichlet series $\mathcal{P} _ {>X}(s)$ diverges throughout $0 < \sigma \le 1$; and long simultaneous recurrence times on all primes up to $N(t)$ do not prevent anti-phase alignment, because Kronecker's theorem on a fixed finite subset of small primes drives finite prime sums arbitrarily close to $-1 \pm \mathrm{i}$ (in accord with Bohr–Landau 1914, Montgomery 1977, and Voronin 1975).

2. **Closed-Form Poisson–Hadamard Lobe Integral and Pole-Free Perturbed Gram Matrices (Section 3).**
   - **Closed-form lobe integral (Theorem 3.1):** For an off-critical conjugate quadruple $\pm\gamma _ * \pm \mathrm{i}d _ *$ ($d _ * > 0$), the single-pair Poisson–Hadamard kernel $K(d, d _ *, u)$ has negative central lobe on $\lvert u\rvert < \sqrt{d _ *^2 - d^2}$ ($0 < d < d _ *$) with exact integral

$$\mathcal{F} _ -(d, d _ *) = \int _ {-\sqrt{d _ *^2 - d^2}}^{\sqrt{d _ *^2 - d^2}} K(d, d _ *, u)\,\mathrm{d}u = -\arcsin\left(\frac{d}{d _ *}\right), \qquad \mathcal{F}(d, d _ *) = \pi\,\Theta(d - d _ *).$$

   - **Pole-free Sylvester inertia detection (Proposition 3.2):** Replacing the simple pair $(\pm\gamma _ 1, \pm\gamma _ 2)$ by $\pm t _ * \pm \mathrm{i}d _ *$ yields a pole-free entire even perturbation $\Xi _ {\mathrm{pert}}(z; d _ *)$ whose de Branges Gram matrix exhibits $\nu _ -(\mathbf{G} _ {\alpha}(d)) \ge 1$ below $d _ *$. Both signatures are contrapositives of Lagarias's (1999) equivalence: they detect assumed off-critical zeros rather than preventing them.

3. **Conditional Nodal-Wall Identity and the Pure-Imaginary Zero of $\Phi$ (Section 4).**
   - **Nodal-wall conditions and counterexample (Proposition 4.1):** Along a smooth branch $I(d, \tau(d)) = 0$, Cauchy–Riemann gives $\widetilde{W}(d, \tau(d)) = R(d, \tau(d))\partial _ t I(d, \tau(d))/d$. Positivity requires $-F(\tau _ 0)F''(\tau _ 0) > 0$, $\partial _ t I \neq 0$, and $R \neq 0$; an explicit polynomial counterexample $Q(z)$ satisfies $-Q(\tau _ 0)Q''(\tau _ 0) > 0$ at $d = 0$ yet vanishes on its own $I = 0$ curve at $d = d _ *$.
   - **Pure-imaginary zero inside the Jacobi strip (Proposition 4.2):** Although Jacobi's triple product ensures $\vartheta _ 3(0, \mathrm{i}\mathrm{e}^{4u}) \neq 0$ on $\lvert \mathrm{Im}(u)\rvert < \frac{\pi}{8}$, applying $(\partial _ u^2 - 1)$ does not preserve zero-freeness: $\Phi(\mathrm{i}\eta _ 0) = 0$ at $\eta _ 0 = 0.3194150268\dots < \frac{\pi}{8}$.

## Repository Structure and Reproducibility

- [`paper/paper10_spectral_resonance_chladni_closure.tex`](paper/paper10_spectral_resonance_chladni_closure.tex) and [`paper/paper10_spectral_resonance_chladni_closure.pdf`](paper/paper10_spectral_resonance_chladni_closure.pdf) – AMS-LaTeX source and compiled 9-page PDF.
- [`scripts/verify_paper10_spectral_resonance_chladni_closure.py`](scripts/verify_paper10_spectral_resonance_chladni_closure.py) – 50-decimal `mpmath` and interval-arithmetic verification suite.
- [`figures/`](figures/) – Vector PDF and PNG figures.
- [`.agents/skills/spectral-resonance-chladni-closure/SKILL.md`](.agents/skills/spectral-resonance-chladni-closure/SKILL.md) – Agent mathematical reference skill.

Run the verification suite via:

```bash
python scripts/verify_paper10_spectral_resonance_chladni_closure.py
```

## License and Citation

Released under the [MIT License](LICENSE). Metadata is provided in [`CITATION.cff`](CITATION.cff) and [`.zenodo.json`](.zenodo.json).
