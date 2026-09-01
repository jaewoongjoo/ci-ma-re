# ci-ma-re

This repository is for the paper **"Random Effects Meta-Analysis When the Studies Use Multiple Regression Models With Different Sets of Covariates"** by Jaewoong Joo, Bikram Karmakar, and Hani Doss.

## Contents

- `ci-ma-re.R`: R script for the simulation studies. It runs 5,000 simulation replicates and saves `ci-ma-re-simulation-results.rds` in the `results` directory.
- `ci-ma-re-figure1.R`: R script that generates the three panels of Figure 1 from the simulation results in the `results` directory.
- `ci-ma-re-figure2.R`: R script that generates the three panels of Figure 2 from the simulation results in the `results` directory.
- `ci-ma-re-figure3-supp-figures-s1-s4.R`: R script that generates Figure 3 in the main paper and Figures S1--S4 in the supplementary material from the simulation results in the `results` directory.
- `ci-ma-re-real-data.R`: R script for thereal-data analysis. Running the script saves the result tables used to construct the real-data figures.
- `ci-ma-re-figures5-6-supp-figure-s5.R`: R script that generates Figures 5 and 6 in the main paper and Figure S5 in the supplementary material. Run `ci-ma-re-real-data.R` first to create the required result tables.

Each script installs any missing R packages, sets its working directory to the directory containing the script, and writes its output beside the code or in the indicated output directory.

Running `ci-ma-re.R` creates the simulation files required by the Figure 1--3 and Figure S1--S4 scripts. Running `ci-ma-re-real-data.R` creates the tables required by the Figure 5, Figure 6, and Figure S5.
