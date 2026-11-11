---
name: spectral-resonance-chladni-closure
description: Mathematical reference and red-team audit guardrails for Paper X on cubic-window Euler-sieve Dirichlet identities, closed-form Poisson-Hadamard lobe integrals (-arcsin(d/d_ast)), pole-free perturbed de Branges Gram matrices, conditional Cauchy-Riemann nodal-wall identities, and the explicit imaginary zero Phi(i eta_0) = 0 inside |Im(u)| < pi/8.
---

# Paper X Skill: Combinatorial Euler-Sieve Identities, Hadamard–Poisson Flux Detection, and Conditional Nodal-Wall Identities

**Author:** Thijs Missiaen (ORCID: `0009-0007-3678-0810`)

## When to Activate

Activate this skill whenever inspecting, extending, or citing Paper X (`paper10-spectral-resonance-chladni-closure`), truncated Euler-product polynomials $F _ 2(z)$, Poisson–Hadamard flux integrals, de Branges reproducing-kernel Gram matrices, or nodal curves $I(d, \tau(d)) = 0$.

## Canonical Normalization and Proven Identities

### 1. Canonical Transverse Wronskian

$$
\Xi(t - \mathrm{i}d) = 8\int _ 0^{\infty} \Phi(u)\cos\bigl(2(t - \mathrm{i}d)u\bigr) \mathrm{d}u, \qquad W(d, t) = \frac{1}{2}\frac{\partial}{\partial d}\lvert \Xi(t - \mathrm{i}d)\rvert^2 = 16\int _ 0^{\infty} y\sinh(2dy) I(y, t) \mathrm{d}y.
$$

### 2. Cubic-Window Dirichlet Coefficient Formula (Proposition 2.1)

For $1 \le N < (\lfloor X\rfloor + 1)^3$, the Dirichlet coefficients of $\zeta(s)\prod _ {p \le X}(1 - p^{-s})$ match those of $1 + \mathcal{P} _ {>X}(s) + \frac{1}{2}\mathcal{P} _ {>X}(s)^2$ identically.

### 3. Closed-Form Poisson–Hadamard Negative Lobe (Theorem 3.1)

For $0 < d < d _ {\ast}$, the symmetric off-critical pair kernel $K(d, d _ {\ast}, u)$ satisfies

$$
\mathcal{F} _ -(d, d _ {\ast}) = \int _ {-\sqrt{d _ {\ast}^2 - d^2}}^{\sqrt{d _ {\ast}^2 - d^2}} K(d, d _ {\ast}, u) \mathrm{d}u = -\arcsin\left(\frac{d}{d _ {\ast}}\right), \qquad \mathcal{F}(d, d _ {\ast}) = \pi \Theta(d - d _ {\ast}).
$$

### 4. Pole-Free Entire Two-Zero Coalescence Perturbation (Equation 3.5)

Moving the simple pair $(\pm\gamma _ 1, \pm\gamma _ 2)$ to $\pm t _ {\ast} \pm \mathrm{i}d _ {\ast}$ ($t _ {\ast} = \frac{1}{2}(\gamma _ 1 + \gamma _ 2){}$) via

$$
\Xi _ {\mathrm{pert}}(z; d _ {\ast}) := \Xi(z) \frac{\bigl(z^2 - (t _ {\ast} - \mathrm{i}d _ {\ast})^2\bigr)\bigl(z^2 - (t _ {\ast} + \mathrm{i}d _ {\ast})^2\bigr)}{\bigl(z^2 - \gamma _ 1^2\bigr)\bigl(z^2 - \gamma _ 2^2\bigr)}
$$

has removable singularities at $\pm\gamma _ 1, \pm\gamma _ 2$ (zero poles), preserves evenness and order $1$, and cleanly demonstrates $\nu _ -(\mathbf{G} _ {\alpha}(d)) \ge 1$ below $d _ {\ast}$.

## Mandatory Guardrails Against Falsified Overclaims

- **Never claim Bohr–Kronecker recurrence prevents anti-phase prime alignment:** In $\mathbb{C}$, $F _ 2(z) = 1 + z + \frac{1}{2}z^2$ vanishes at $z = -1 \pm \mathrm{i}$, $\mathcal{P} _ {>X}(s)$ diverges in $0 < \sigma \le 1$, and Kronecker density on a fixed finite set of small primes already drives finite prime sums arbitrarily close to $-1 \pm \mathrm{i}$ (Proposition 2.2; Bohr–Landau 1914, Voronin 1975, Montgomery 1977).
- **Never present flux steps or Sylvester jumps as zero prevention:** Both are contrapositives of Lagarias's (1999) equivalence (detection of assumed off-critical zeros, holding equally for Davenport–Heilbronn).
- **Never divide $\Xi(z)$ by $(z - \gamma _ 1)^2$:** Since $\gamma _ 1$ is a simple zero, dividing by $(z - \gamma _ 1)^2$ introduces a simple pole at $\gamma _ 1$; always use the removable two-zero factor $\Xi _ {\mathrm{pert}}(z; d _ {\ast})$.
- **Never treat nodal-wall positivity or $\vartheta _ 3 \neq 0$ as unconditional strip proofs:** Positivity along $I(d, \tau(d)) = 0$ requires $-F(\tau _ 0)F''(\tau _ 0) > 0$, $\partial _ t I \neq 0$, and $R \neq 0$ (refuted unconditionally by the polynomial counterexample $Q(z)$ in Proposition 4.1(iii)), and $\Phi(\mathrm{i}\eta _ 0) = 0$ vanishes inside $\lvert \mathrm{Im}(u)\rvert < \frac{\pi}{8}$ at $\eta _ 0 = 0.3194150268\dots$ (Proposition 4.2).

## Verification Command

```bash
python scripts/verify_paper10_spectral_resonance_chladni_closure.py
```
