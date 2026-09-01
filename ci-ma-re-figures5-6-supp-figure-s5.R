rm(list = ls())
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
load(file.path(code.dir, "ci-ma-re-real-data-results.Rdata"))
model_order_mu <- c("Naive", "DL", "Meta")

plot_data_mu <- mu_comparison[match(model_order_mu, mu_comparison$Model), , drop = FALSE]
pdf("mu-hyper.pdf", width = 10, height = 3.5)
par(mfrow = c(1, 1))
plot(plot_data_mu$Estimate, 1:nrow(plot_data_mu), xlim = range(c(plot_data_mu$Lower_CI, plot_data_mu$Upper_CI), 
    na.rm = TRUE), ylim = c(0.5, nrow(plot_data_mu) + 0.5), xaxt = "n", yaxt = "n", xlab = expression(mu[ht]), 
    ylab = "", pch = 16, cex = 1.5, col = "black")

arrows(plot_data_mu$Lower_CI, 1:nrow(plot_data_mu), plot_data_mu$Upper_CI, 1:nrow(plot_data_mu), 
    angle = 90, code = 3, length = 0.1, col = "black")

axis(2, at = 1:nrow(plot_data_mu), labels = plot_data_mu$Model, cex.axis = 1)
axis(1, cex.axis = 1)
abline(h = 1:nrow(plot_data_mu), col = "lightgray", lty = "dotted")
dev.off()
plot_beta_region <- function(plot_data_k, region_k) {
    model_order <- c("Naive", "Full", "Meta")
    plot_data_k <- plot_data_k[match(model_order, plot_data_k$Model), , drop = FALSE]
    finite_rows <- is.finite(plot_data_k$Estimate) & is.finite(plot_data_k$Lower_CI) & is.finite(plot_data_k$Upper_CI)
    x_rng <- range(c(plot_data_k$Lower_CI[finite_rows], plot_data_k$Upper_CI[finite_rows]), 
        na.rm = TRUE)
    if (!is.finite(diff(x_rng)) || diff(x_rng) == 0) {
        x_rng <- x_rng + c(-0.1, 0.1)
    }
    region_idx <- as.integer(region_k)
    plot(plot_data_k$Estimate[finite_rows], which(finite_rows), xlim = x_rng, ylim = c(0.5, 
        nrow(plot_data_k) + 0.5), xaxt = "n", yaxt = "n", xlab = bquote(beta[ht]^.(paste0("(", 
        region_idx, ")"))), ylab = "", main = "", pch = 16, cex = 1.3, col = "black")
    arrows(plot_data_k$Lower_CI[finite_rows], which(finite_rows), plot_data_k$Upper_CI[finite_rows], 
        which(finite_rows), angle = 90, code = 3, length = 0.08, col = "black")
    axis(2, at = 1:nrow(plot_data_k), labels = plot_data_k$Model, cex.axis = 0.9)
    axis(1, cex.axis = 0.85)
    abline(h = 1:nrow(plot_data_k), col = "lightgray", lty = "dotted")
}

beta_pdf_dir <- file.path(code.dir, "beta-hyper")
if (dir.exists(beta_pdf_dir)) {
    unlink(beta_pdf_dir, recursive = TRUE)
}

dir.create(beta_pdf_dir)
region_levels <- sort(unique(as.character(beta_all_comparison$region)))
for (region_k in region_levels) {
    plot_data_k <- beta_all_comparison[as.character(beta_all_comparison$region) == region_k, 
        , drop = FALSE]
    file_region <- gsub("[^A-Za-z0-9_]+", "_", region_k)
    out_file <- file.path(beta_pdf_dir, paste0("beta-hyper-study", file_region, ".pdf"))
    pdf(file = out_file, width = 4, height = 4)
    par(mar = c(4, 2.5, 1, 1))
    plot_beta_region(plot_data_k, region_k)
    dev.off()
}
