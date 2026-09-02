# ci-ma-re

This repository is for the paper **"Random Effects Meta-Analysis When the Studies Use Multiple Regression Models With Different Sets of Covariates"** by Jaewoong Joo, Bikram Karmakar, and Hani Doss.

## Contents

- `ci-ma-re.R`: R script for the simulation studies. It runs 5,000 simulation replicates for each simulation scenario and saves `ci-ma-re-simulation-results.rds` in the `results` directory.
- `results/ci-ma-re-simulation-results.rds`: simulation results used to produce the numerical summaries and figures reported in the paper.
- `ci-ma-re-figure1.R`: R script that generates the three panels of Figure 1 from the simulation results in the `results` directory.
- `ci-ma-re-figure2.R`: R script that generates the three panels of Figure 2 from the simulation results in the `results` directory.
- `ci-ma-re-figure3-supp-figures-s1-s4.R`: R script that generates Figure 3 in the main paper and Figures S1--S4 in the supplementary material from the simulation results in the `results` directory.
- `ci-ma-re-real-data.R`: R script for the real-data analysis. Running the script saves the result tables used to construct the real-data figures.
- `ci-ma-re-figures5-6-supp-figure-s5.R`: R script that generates Figures 5 and 6 in the main paper and Figure S5 in the supplementary material. Run `ci-ma-re-real-data.R` first to create the required result tables.

Each script installs any missing R packages, sets its working directory to the directory containing the script, and writes its output beside the code or in the indicated output directory.

## Reproducing the Simulation Figures

The final simulation output is included in `results/ci-ma-re-simulation-results.rds`, so the simulation does not need to be rerun to reproduce the figures. Run `ci-ma-re-figure1.R`, `ci-ma-re-figure2.R`, and `ci-ma-re-figure3-supp-figures-s1-s4.R` to generate Figures 1--3 in the main paper and Figures S1--S4 in the supplementary material.

## Reproducing the Real-Data Figures

Run ci-ma-re-real-data.R and then ci-ma-re-figures5-6-supp-figure-s5.R to generate Figures 5 and 6 in the main paper and Figure S5 in the supplementary material. Because the UNOS data are subject to restrictions on redistribution, they are not included in this repository.
