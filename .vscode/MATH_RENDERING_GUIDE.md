# Math Rendering & Chat Conventions

Purpose: Ensure all future AI/chat interactions in this workspace use properly rendered LaTeX.

## Conventions
- Inline math: wrap expressions with single dollar signs: `$a^2 + b^2 = c^2$`.
- Display math: use double dollar blocks:
  $$
  \int_0^1 x^2 \, dx = \frac{1}{3}
  $$
- Avoid raw LaTeX without delimiters.
- Prefer concise variable names; define them before use.
- When giving a sequence of equations, align with `aligned` inside display math:
  $$
  \begin{aligned}
  y_t &= f(x_t, \theta) \\
  \ell(\theta) &= \sum_t \log p(y_t \mid x_t, \theta)
  \end{aligned}
  $$

## Checklist Before Sending Math
1. Inline variables wrapped in `$...$`.
2. Multi-line derivations in a `$$...$$` block.
3. No stray underscores or backslashes without math delimiters.
4. For probabilities, use `\mathbb{P}` or `p(\cdot)` consistently.

## Common Snippets
- Parameterization: `$x \sim \mathcal{N}(\mu, \sigma^2)$`.
- Logistic: `$\text{logit}^{-1}(x) = 1 / (1 + e^{-x})$`.
- GP prior: `$$ f \sim \mathcal{GP}(0, k(\cdot, \cdot)) $$`.

## Rationale
Consistent formatting avoids rework and ensures downstream renderers (Quarto, Markdown preview, chat UI) display equations correctly.

---
_Last updated: 2025-09-19_
