---
name: tikz-diagram
description: Create publication-quality TikZ diagrams for the Quarto website. Use when the user asks to create a diagram, figure, state diagram, flow chart, or visual for documentation pages.
user-invocable: true
allowed-tools: [Read, Write, Bash, Glob]
---

# TikZ Diagram Creation

Create publication-quality vector diagrams using TikZ/LaTeX, compiled to SVG for use in the Quarto website.

## Workflow

1. Write a `.tex` file in `quarto/website/images/`
2. Compile to PDF with `pdflatex`
3. Convert to SVG with `pdf2svg`
4. Reference the `.svg` in `.qmd` files
5. Clean up build artifacts (`.aux`, `.log`, `.pdf`)

## Prerequisites Check (Run Before First Compile)

Before compiling, verify tools are installed. Run these checks and install anything missing:

```bash
# 1. Check for pdflatex (via TinyTeX)
if [ ! -f "$HOME/.TinyTeX/bin/x86_64-linux/pdflatex" ]; then
  echo "Installing TinyTeX..."
  Rscript -e 'tinytex::install_tinytex()'
fi

# 2. Check for sfmath package (required for sans-serif math)
if ! $HOME/.TinyTeX/bin/x86_64-linux/kpsewhich sfmath.sty > /dev/null 2>&1; then
  echo "Installing sfmath..."
  $HOME/.TinyTeX/bin/x86_64-linux/tlmgr install sfmath
fi

# 3. Check for pdf2svg
if ! command -v pdf2svg > /dev/null 2>&1; then
  echo "Installing pdf2svg..."
  sudo apt-get update -qq && sudo apt-get install -y pdf2svg
fi
```

If TinyTeX is missing entirely and `tinytex` R package is not available:
```bash
# Manual TinyTeX install
wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh
~/.TinyTeX/bin/x86_64-linux/tlmgr install sfmath standalone
```

If pdflatex fails with a missing `.sty` error, install the package:
```bash
~/.TinyTeX/bin/x86_64-linux/tlmgr install <package-name>
```

## Tool Paths

```bash
PDFLATEX="$HOME/.TinyTeX/bin/x86_64-linux/pdflatex"
PDF2SVG="pdf2svg"
OUTDIR="quarto/website/images"
```

## Compile and Convert Commands

Run from the project root (`/mnt/code`):

```bash
PDFLATEX="$HOME/.TinyTeX/bin/x86_64-linux/pdflatex"
NAME="my-diagram"

# Compile .tex to .pdf
$PDFLATEX -output-directory=quarto/website/images quarto/website/images/${NAME}.tex

# Convert .pdf to .svg
pdf2svg quarto/website/images/${NAME}.pdf quarto/website/images/${NAME}.svg

# Clean up build artifacts
rm -f quarto/website/images/${NAME}.{aux,log,pdf}
```

## .tex File Template

Use this standard preamble matching the project's existing diagrams:

```latex
% Compile with: pdflatex DIAGRAM_NAME.tex
% Then convert to SVG: pdf2svg DIAGRAM_NAME.pdf DIAGRAM_NAME.svg
\documentclass[tikz,border=10pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{shapes,arrows.meta,positioning,calc,fit}
\usepackage{amsmath}
\usepackage{helvet}
\renewcommand{\familydefault}{\sfdefault}
\usepackage{sfmath}

\begin{document}

\begin{tikzpicture}[
  % Define styles here
]

% Diagram content here

\end{tikzpicture}
\end{document}
```

Key style conventions:
- **Font**: Helvetica via `helvet` + `sfmath` (sans-serif including math)
- **Document class**: `standalone` with `tikz` option and appropriate `border`
- **Arrow tips**: Use `arrows.meta` library with `Stealth` tips
- **Line widths**: 1.2pt for primary elements, 0.8pt for secondary, 0.5pt for tertiary
- **Rounded corners**: Use `rounded corners=6pt` to `8pt` for boxes

## Referencing in Quarto

In `.qmd` files, reference the SVG with a relative path from the page location:

```markdown
![Caption text describing the diagram.](../images/diagram-name.svg){width=70%}
```

Adjust `width` as needed (60-80% is typical).

## Existing Diagrams

Reference these for style consistency:

| File | Content |
|------|---------|
| `model-diagram.tex` | Data flow architecture (inputs/model/outputs) |
| `multistate-diagram.tex` | Illness-death state diagram (3 states, 3 transitions) |
| `latent-burden-diagram.tex` | Multi-biomarker observation model (latent state to data) |

Read an existing `.tex` file before creating a new diagram to match the visual style.
