# Creating a One-Pager Report Slide

This guide explains how to create a compact, single-slide report using Quarto revealjs format, similar to the PIONEER one-pager (`quarto/reports/pioneer-one-pager-slide.qmd`).

## Overview

A one-pager slide is:
- A **single revealjs slide** (not a full presentation) that presents key findings, plots, and conclusions
- Designed for **printing or sharing as a standalone HTML file** with `embed-resources: true`
- Uses a **multi-column layout** (typically 3 columns: text + plots + plots)
- Applies **corporate styling** from an existing SCSS theme
- Includes **custom CSS** for fine-tuned layout control

## Step 1: Set Up the Quarto File

### Basic YAML Frontmatter

```yaml
---
title: ""  # Empty title - you'll add custom title in body
format:
  revealjs:
    theme: [default, path/to/corporate-theme.scss]
    width: 1600
    height: 1130
    margin: 0.03
    center: false
    controls: false
    progress: false
    slide-number: false
    embed-resources: true  # Critical for standalone HTML
    auto-stretch: false
execute:
  warning: false
  echo: false
  cache: true  # Speed up renders, but see caching notes below
---
```

**Key settings:**
- `width: 1600, height: 1130` - Standard widescreen dimensions
- `embed-resources: true` - Makes HTML file self-contained (no external dependencies)
- `cache: true` - Speeds up rendering BUT requires manual cache clearing when R source files change

### Setup Chunk

```r
```{r setup}
#| include: false
#| cache: false  # Never cache the setup chunk

# Set working directory to project root
knitr::opts_knit$set(root.dir = here::here())
setwd(here::here())
source(here::here("path", "to", "_setup.R"))

# Create a compact theme for report plots
report_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", color = corporate_navy),
    plot.subtitle = element_text(size = 8, color = corporate_grey),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    plot.margin = margin(2, 4, 2, 2, "pt"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 6.5, hjust = 0, color = corporate_grey)
  )

# Load your data using targets or other methods
data_list <- list(
  cutoff1 = tar_read(data_cutoff1),
  cutoff2 = tar_read(data_cutoff2),
  cutoff3 = tar_read(data_cutoff3)
)
```
```

## Step 2: Create Custom CSS

Add inline CSS after the YAML but before content starts. This controls layout and styling.

```markdown
## {.smaller}

<style>
/* Hide default slide footer */
.reveal .slide-footer { display: none; }

/* Custom title styling */
.report-title {
  color: var(--corporate-navy) !important;
  font-size: 1.3em !important;
  margin: 0 !important;
  text-align: center;
  font-weight: 700;
}
.report-subtitle {
  color: var(--corporate-navy);
  font-size: 0.6em;
  text-align: center;
  opacity: 0.85;
  margin: 2px 0 0 0;
}
.report-authors {
  font-size: 0.42em;
  color: var(--corporate-grey);
  text-align: center;
  margin: 2px 0 0 0;
}

/* Gold horizontal rule */
.gold-rule {
  height: 3px;
  background: var(--corporate-gold);
  border: none;
  margin: 4px 0 6px 0;
}

/* Column gap override - default 2em is too wide */
.reveal .columns { gap: 0.8em; }

/* Text column - vertically centered */
.text-col {
  display: flex;
  flex-direction: column;
  justify-content: center;
}

/* Consistent heading and text sizing across columns */
.text-col h3, .plot-col h3 {
  color: var(--corporate-navy) !important;
  font-size: 1.1em !important;
  margin: 6px 0 2px 0 !important;
}
.text-col h3:first-child { margin-top: 0 !important; }

.text-col p, .plot-col p {
  margin: 2px 0;
  font-size: 0.82em;
  line-height: 1.4;
}

.text-col ul { margin: 2px 0; }
.text-col li {
  margin-bottom: 2px;
  font-size: 0.82em;
  line-height: 1.4;
}

/* Plot cell spacing */
.plot-col .cell { margin: 0 0 4px 0; }
.plot-col img { width: 100% !important; border-radius: 2px; }

/* Custom footer bar */
.report-footer-bar {
  position: absolute;
  bottom: 5px;
  left: 3%;
  right: 3%;
  display: flex;
  justify-content: space-between;
  font-size: 0.36em;
  color: var(--corporate-navy);
  border-top: 2px solid var(--corporate-gold);
  padding-top: 3px;
}
</style>
```

**CSS Variables:** Your SCSS theme should export color variables in `:root` (e.g., `--corporate-navy`, `--corporate-gold`). See example SCSS at end of this guide.

## Step 3: Add Title and Content Structure

```markdown
::: {.report-title}
PROJECT: Study Name Analysis
:::

::: {.report-subtitle}
Subtitle Describing the Analysis
:::

::: {.report-authors}
Author Names | Month Year
:::

::: {.gold-rule}
:::

:::: {.columns}

::: {.column width="25%" .text-col}

### Key Findings

- **Finding 1**: Description with emphasis on key points.
- **Finding 2**: Another finding referencing Panel A or Panel B.
- **Finding 3**: More findings.
- **Finding 4**: Final key point.

:::

::: {.column width="36%" .plot-col}

```{r plot-a}
#| fig-width: 7.5
#| fig-height: 5.5

your_plotting_function() +
  labs(title = "A: First Plot Title") +
  report_theme
```

```{r plot-b}
#| fig-width: 7.5
#| fig-height: 2.0

another_plotting_function() +
  labs(title = "B: Second Plot Title") +
  report_theme
```

### Section Title

Text content describing methodology or context. This section appears below the plots in the middle column.

:::

::: {.column width="36%" .plot-col}

```{r plot-c}
#| fig-width: 7.5
#| fig-height: 4.0

third_plotting_function() +
  labs(title = "C: Third Plot Title") +
  report_theme +
  theme(
    # Add plot-specific theme overrides here
    legend.position = "none",
    panel.grid = element_blank()
  )
```

```{r plot-d}
#| fig-width: 7.5
#| fig-height: 3.5

fourth_plotting_function() +
  labs(title = "D: Fourth Plot Title") +
  report_theme
```

### Conclusions

Final summary text highlighting the main takeaways and implications of the analysis.

:::

::::

::: {.report-footer-bar}
Project Name | Analysis Type | Confidentiality Level

Team Name, Department
:::
```

## Step 4: Layout Guidelines

### Column Width Strategy

**Total width must be < 100%** to account for column gaps (`.reveal .columns { gap: 0.8em }`).

Recommended distributions:
- **3 columns:** `25% + 36% + 36% = 97%`
- **2 columns:** `45% + 50% = 95%`
- **Unequal:** `30% + 32% + 32% = 94%`

If content overflows off the right edge, reduce widths or gap further.

### Plot Height Balancing

**Goal:** Align text section headings at the same vertical height across columns.

Example:
```
Column 2: Plot A (5.5) + Plot B (2.0) + Text = 7.5 + text height
Column 3: Plot C (4.0) + Plot D (3.5) + Text = 7.5 + text height
```

If text sections don't align:
1. Adjust plot heights to balance total height per column
2. Check that text sections have similar content length
3. Use flexbox centering in `.text-col` for vertical positioning

### Vertical Centering

The `.text-col` uses flexbox centering:
```css
.text-col {
  display: flex;
  flex-direction: column;
  justify-content: center;
}
```

This vertically centers the Key Findings section. For other columns, text naturally flows from top.

## Step 5: Common Issues and Solutions

### Issue: Plot theme overrides don't work

**Symptom:** You add `theme(panel.grid = element_blank())` in the plot function but grids still appear.

**Cause:** `report_theme` is based on `theme_minimal()`, which includes default grids. When you add `+ report_theme` AFTER your plot function, it re-applies `theme_minimal()` defaults.

**Solution:** Add theme overrides AFTER `report_theme`:
```r
plot_function() +
  labs(title = "Plot Title") +
  report_theme +
  theme(
    panel.grid = element_blank(),  # Add AFTER report_theme
    legend.position = "none"
  )
```

### Issue: Changes not appearing after render

**Symptom:** You modify R source files (e.g., `r/plot_functions.R`) but plots don't update.

**Cause:** `cache: true` in YAML means chunks aren't re-executed unless the chunk itself changes.

**Solution:** Clear the cache directory:
```bash
rm -rf quarto/reports/your-slide_cache
quarto render quarto/reports/your-slide.qmd
```

**Best practice:** Clear cache after ANY R source file modification.

### Issue: X-axis doesn't include zero

**Symptom:** You want to show a zero reference line but the axis is cropped.

**Solution:** Manually expand axis limits:
```r
your_plot +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  scale_x_continuous(
    limits = function(x) c(min(x[1], 0), max(x[2], 0))  # Ensure zero included
  )
```

Or in a plot function with an `add_vline` flag:
```r
if (add_vline) {
  x_min <- min(x_min, 0)
  x_max <- max(x_max, 0)
}
```

### Issue: Content overflows off right edge

**Symptom:** Right column plots or text are cut off or overflow off the slide.

**Solution:**
1. Reduce column widths (e.g., from `26% + 37% + 37% = 100%` to `25% + 36% + 36% = 97%`)
2. Reduce column gap: `.reveal .columns { gap: 0.8em; }` (default is `2em`)
3. Check plot widths - use `fig-width: 7.5` not larger values

### Issue: Confusion matrix shows unwanted grids

**Symptom:** Heatmap/confusion matrix plots show background grids.

**Solution:** Add `panel.grid = element_blank()` in the theme override AND in the base plot function:
```r
plot_confusion_matrix <- function(...) {
  # ... plot code ...
  theme(
    panel.grid = element_blank(),  # Explicitly remove grids
    legend.position = "right"
  )
}
```

Then in the report:
```r
plot_confusion_matrix(...) +
  report_theme +
  theme(panel.grid = element_blank())  # Belt and suspenders
```

## Step 6: Writing Style Guidelines

### Apply Elements of Style

Use `/writing-clearly-and-concisely` skill to review all prose sections (Key Findings, conclusions, methodology text).

**Key principles:**
- **Omit needless words:** "closely track" not "closely follow and track"
- **Active voice:** "Historical trials anchor parameters" not "Parameters are anchored by historical trials"
- **Positive form:** "produces better scores" not "does not produce worse scores"
- **Keep related words together:** "Model predictions closely track observed curves" not "Model predictions track closely the observed curves"
- **Concrete language:** "70-80% sensitivity" not "high sensitivity"

### Panel References

Reference plots as "Panel A", "Panel B", etc. in Key Findings bullets. This helps readers navigate the visual layout.

Example:
```markdown
- **Stable forecasts** (Panel A): Model predictions closely track observed curves.
- **Historical data strengthens predictions** (Panel B): Historical trials produce better scores.
```

## Step 7: Rendering Workflow

### Initial Render

Always render from the project root (where `renv` is configured):

```bash
# From project root
quarto render quarto/reports/your-slide.qmd
```

**Never render from inside the subdirectory** if your `_quarto.yml` uses `execute-dir: project`.

### Iterative Workflow

1. Make changes to `.qmd` file or R source files
2. Clear cache if R sources changed: `rm -rf quarto/reports/your-slide_cache`
3. Render: `quarto render quarto/reports/your-slide.qmd`
4. Open HTML in browser to check layout
5. Adjust column widths, plot heights, or CSS as needed
6. Repeat

### Preview Mode

For faster iteration with live reload:
```bash
quarto preview quarto/reports/your-slide.qmd
```

This opens a browser and auto-renders on save. **Note:** Cache issues still apply - clear cache if needed.

## Step 8: Example SCSS Theme Setup

Your corporate SCSS theme should export CSS variables for use in inline `<style>` blocks.

Example `corporate-theme.scss`:
```scss
/*-- scss:defaults --*/

// Corporate colors
$corporate-navy: #003865;
$corporate-gold: #F0AB00;
$corporate-turquoise: #68D2DF;
$corporate-grey: #666666;

// Revealjs settings
$presentation-heading-color: $corporate-navy;
$presentation-font-size-root: 32px;
$body-bg: #ffffff;
$body-color: $corporate-grey;

/*-- scss:rules --*/

// Export colors as CSS custom properties
:root {
  --corporate-navy: #{$corporate-navy};
  --corporate-gold: #{$corporate-gold};
  --corporate-turquoise: #{$corporate-turquoise};
  --corporate-grey: #{$corporate-grey};
}

// Columns with reasonable gap
.reveal .columns {
  display: flex;
  gap: 2em;  // Can be overridden in inline CSS
}

.reveal .column {
  flex: 0 0 auto;  // Respect width attributes
}
```

## Complete Example Structure

```
quarto/reports/
├── your-slide.qmd              # Main one-pager file
├── your-slide_cache/           # Auto-generated, should be in .gitignore
└── your-slide.html             # Rendered output
```

Typical file size: `.qmd` ~300 lines, `.html` ~500KB-2MB (with embedded plots).

## Troubleshooting Checklist

Before asking for help, verify:

- [ ] Column widths sum to < 100% (recommend 97%)
- [ ] Column gap reduced if needed: `.reveal .columns { gap: 0.8em; }`
- [ ] Cache cleared after R source changes: `rm -rf *_cache`
- [ ] Rendering from project root, not subdirectory
- [ ] Theme overrides placed AFTER `report_theme`
- [ ] Panel references (A, B, C, D) match plot order
- [ ] Prose reviewed with Elements of Style principles
- [ ] `embed-resources: true` for standalone HTML
- [ ] Plot heights balanced across columns for text alignment

## Additional Resources

- Quarto revealjs docs: https://quarto.org/docs/presentations/revealjs/
- Flexbox guide: https://css-tricks.com/snippets/css/a-guide-to-flexbox/
- Elements of Style: Use `/writing-clearly-and-concisely` skill
- Example: `quarto/reports/pioneer-one-pager-slide.qmd` in the SCLC repo
