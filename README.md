# Checking the Cox Proportional Hazards Model with Interval-Censored Data

This repository contains the code and data required to reproduce the figures and results from our paper on checking the Cox regression model with interval-censored data.

## 📁 Repository Structure

### `simulation-code/`
Contains the code for generating the simulation results for Figures 1–3.

- **Figure 1**: Run the R scripts in `prototype/`.
- **Figure 2**: Run the R scripts in the following folders:
  - `ff_n200_quad/`, `ff_n400_quad/`, `ff_n800_quad/`
  - `ff_n200_twisted/`, `ff_n400_twisted/`, `ff_n800_twisted/`
  - `omnibus_n200_quad/`, `omnibus_n400_quad/`, `omnibus_n800_quad/`
  - `omnibus_n200_twisted/`, `omnibus_n400_twisted/`, `omnibus_n800_twisted/`
- **Figure 3**: Run the R scripts in:
  - `prop_n200_mono/`, `prop_n400_mono/`, `prop_n800_mono/`
  - `prop_n200_quad/`, `prop_n400_quad/`, `prop_n800_quad/`

Each R script must be run multiple times to obtain enough simulation replicates.

---

### `application-code/`
Contains the code for generating the results for Figures 4–5 and numerical summaries in Section 4 of the paper.

- **Figures 4**: R scripts in `prop_diabt_glucose/` and `ff_diabt_bmi_logtGlucose/`
- **Figures 5**: R scripts in `ff_hyper_SysBP/` and `ff_hyper_DiaBP/`
- Additional R scripts produce outputs referenced in Section 4.

---

### `figures-code/`
Contains all R scripts used to generate each figure and its sub-figures.

Run the scripts in the respective subfolders:
- `figures-code/Figure1/`
- `figures-code/Figure2/`
- `figures-code/Figure3/`
- `figures-code/Figure4/`
- `figures-code/Figure5/`
- `figures-code/FigureS1/`
- `figures-code/FigureS2/`

---

### `figures/`
Contains all the figures included in the paper.

---

### `application-synthetic-data/`
Contains the synthetic datasets used in the application section, along with descriptions.

---

### `output/`
Contains summarized results required to reproduce all the figures.

---

### `wrapper-scripts/`
Contains R scripts that iterate through all sub-folders in the correct order, eliminating the need for manual navigation and ensuring that the regenerated outputs match those in `figures`:

- `wrapper-scripts/run_simulations.R`
- `wrapper-scripts/render_figures.R`

---

## 🔧 Setup and Installation

### System Requirements

- **OS**: Linux, macOS, or Windows
- **R Version**: ≥ 4.3.2
- **Compiler**: Ensure R is configured with a C/C++ compiler (e.g., `g++`, `clang`)

### Required R Packages

Install the required R packages using the script below:

```r
# Use 'remotes' to install specific versions
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_version("tidyverse", version = "2.0.0")
remotes::install_version("Rcpp", version = "1.0.13")
remotes::install_version("RcppArmadillo", version = "0.12.8.4.0")
remotes::install_version("ggpubr", version = "0.6.0")
remotes::install_version("latex2exp", version = "0.9.6")
remotes::install_version("dplyr", version = "1.1.4")
```

Or use the provided R script `install_packages.R` if available.

---

## 🚀 Reproducing Results

1. **Clone the repository**

```bash
git clone https://github.com/yangjc-xu/check-cox-int.git
cd check-cox-int
```

2. **Open R or RStudio**

3. **Run the relevant scripts as described above**

   For example, to reproduce Figure 1:

   ```r
   source("simulation-code/prototype/run_simulation.R")
   ```

4. **Repeat simulations as necessary for stable results**

