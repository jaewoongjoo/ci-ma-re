packages <- c("stats", "base", "rstudioapi")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
}

invisible(lapply(packages, library, character.only = TRUE))
loc.current <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
    }
    else if (Sys.getenv("RSTUDIO") == "1") {
        return(dirname(rstudioapi::getSourceEditorContext()$path))
    }
}

code.dir <- loc.current()
setwd(code.dir)
results.dir <- file.path(code.dir, "results")
figures.dir <- file.path(code.dir, "figures")
dir.create(figures.dir, showWarnings = FALSE, recursive = TRUE)
simulation.file <- file.path(results.dir, "ci-ma-re-simulation-results.rds")
simulation.results <- readRDS(simulation.file)
rho.tags <- c("corr0.3", "corr0.6", "corr0.8", "corr0.9")
rho.vals <- c("0.3", "0.6", "0.8", "0.9")
p <- 8
k <- 20
n <- p + 1
summ <- setNames(vector("list", length(rho.tags)), rho.tags)
for (tag in rho.tags) {
    scenario.results <- simulation.results[[tag]]
    nsim.ok <- scenario.results$nsim_ok
    coverage.beta <- scenario.results$coverage_beta
    rownames(coverage.beta) <- paste0("beta[", 0:p, "]")
    colnames(coverage.beta) <- paste0("study_", 1:k)
    summ[[tag]] <- list(nsim_ok = nsim.ok, coverage_beta = coverage.beta)
}

idx.noint <- 2:n
pch.values <- c(0, 1, 2, 5)
pt.cex <- 1.4
rho.legend <- as.expression(lapply(rho.vals, function(r) {
    bquote(rho[X] == .(as.numeric(r)))
}))

plot.coverage <- function(study.index, file.out) {
    coverage.by.rho <- sapply(seq_along(rho.tags), function(i) {
        as.numeric(summ[[rho.tags[i]]]$coverage_beta[idx.noint, study.index])
    })
    beta.lab <- as.expression(lapply(1:p, function(j) {
        bquote(beta[.(j)]^{
            (.(study.index))
        })
    }))
    pdf(file.out, width = 5, height = 5)
    par(mar = c(4, 3, 3, 1))
    plot(NULL, xlim = c(1, p), ylim = c(0.87, 1.04), xlab = "", ylab = "", xaxt = "n")
    axis(1, at = 1:p, labels = beta.lab, cex.axis = 1)
    abline(h = 0.95, lty = 2, col = "red")
    for (i in seq_along(rho.tags)) {
        points(1:p, coverage.by.rho[, i], pch = pch.values[i], cex = pt.cex)
    }
    legend("topleft", legend = rho.legend, pch = pch.values, pt.cex = pt.cex, bty = "o")
    dev.off()
}

for (study.index in 1:k) {
    plot.coverage(study.index, file.path(figures.dir, sprintf("fig3_beta_cov_study%d.pdf", 
        study.index)))
}
