# Setup plotting for VS Code remote environment
# Run this at the start of your R session

# Set options for better plot display
options(
  device = function(...) {
    png(..., width = 800, height = 600, res = 120, type = "cairo")
  },
  bitmapType = "cairo"
)

# For ggplot2, ensure it uses the right device
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12))
}

# Print confirmation
cat("✓ Plotting setup complete for VS Code remote environment\n")
cat("✓ Using PNG device with Cairo backend\n")
cat("✓ Default plot size: 800x600 at 120 DPI\n")
