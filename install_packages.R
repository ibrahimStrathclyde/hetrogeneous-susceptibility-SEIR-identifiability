# install_packages.R
# Script to install all required packages for this repository

# Function to install packages if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
    cat(paste("Package", pkg, "installed.\n"))
  } else {
    cat(paste("Package", pkg, "already installed.\n"))
  }
}

# List of required packages
required_packages <- c(
  "tidyverse",
  "deSolve",
  "GGally",
  "gridExtra",
  "grid"
)

# Install missing packages
for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat("\nAll required packages are now installed.\n")
