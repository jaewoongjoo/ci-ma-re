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
mu.true <- c(-3, 3, 1, -1, -3, 3, 1, -1, -3)
p <- 8
n <- p + 1
summ <- setNames(vector("list", length(rho.tags)), rho.tags)
for (tag in rho.tags) {
    scenario.results <- simulation.results[[tag]]
    mu.hat <- scenario.results$mu_hat_mat
    omega.opt <- scenario.results$omega_opt
    replicate.status <- scenario.results$replicate_status
    ok <- replicate.status == "ok" & is.finite(omega.opt)
    nsim.ok <- sum(ok)
    mu.hat.ok <- mu.hat[ok, , drop = FALSE]
    summ[[tag]] <- list(nsim_ok = nsim.ok, mean_mu = colMeans(mu.hat.ok), var_mu = apply(mu.hat.ok, 
        2, var), mse_mu = colMeans((mu.hat.ok - matrix(rep(mu.true, nsim.ok), nrow = nsim.ok, 
        ncol = n, byrow = TRUE))^2), mu_hat_mat = mu.hat.ok)
}

idx.noint <- 2:n
pch.values <- c(0, 1, 2, 5)
pch.values.filled <- c(15, 16, 17, 18)
pt.cex <- 1.4
mu.lab <- as.expression(lapply(1:p, function(i) bquote(mu[.(i)])))
rho.legend <- as.expression(lapply(rho.vals, function(r) {
    bquote(rho[X] == .(as.numeric(r)))
}))

pdf(file.path(figures.dir, "fig1_mu_mean.pdf"), width = 5, height = 5)
par(mar = c(4, 3, 3, 1))
plot(NULL, xlim = c(1, p), ylim = c(-6, 6), xlab = "", ylab = "", xaxt = "n")
axis(1, at = 1:p, labels = mu.lab, cex.axis = 1.25)
points(1:p, mu.true[idx.noint], pch = 3, col = "red", cex = 1.6)
for (i in seq_along(rho.tags)) {
    points(1:p, summ[[rho.tags[i]]]$mean_mu[idx.noint], pch = pch.values[i], cex = pt.cex)
}

legend("bottomleft", legend = c("True", rho.legend), pch = c(3, pch.values), col = c("red", 
    rep("black", length(rho.tags))), pt.cex = c(1.6, rep(pt.cex, length(rho.tags))))

dev.off()
pdf(file.path(figures.dir, "fig1_mu_sqrtmse.pdf"), width = 5, height = 5)
par(mar = c(4, 3, 3, 1))
plot(NULL, xlim = c(1, p), ylim = c(0, 35), xlab = "", ylab = "", xaxt = "n")
axis(1, at = 1:p, labels = mu.lab, cex.axis = 1.25)
for (i in seq_along(rho.tags)) {
    points(1:p, sqrt(summ[[rho.tags[i]]]$mse_mu[idx.noint]), pch = pch.values[i], cex = pt.cex)
    points(1:p, sqrt(summ[[rho.tags[i]]]$var_mu[idx.noint]), pch = pch.values.filled[i], 
        col = "red", cex = 0.9)
}

legend("topright", legend = rho.legend, pch = pch.values, pt.cex = pt.cex)
dev.off()
mu.component <- 2
pdf(file.path(figures.dir, "fig1_mu_density.pdf"), width = 5, height = 5)
par(mar = c(4, 3, 3, 1))
plot(NULL, xlim = c(-40, 40), ylim = c(0, 0.35), xlab = "", ylab = "")
for (i in seq_along(rho.tags)) {
    values <- summ[[rho.tags[i]]]$mu_hat_mat[, mu.component]
    values <- values[is.finite(values)]
    density.fit <- density(values)
    lines(density.fit$x, density.fit$y, lwd = 2, lty = i)
}

abline(v = mu.true[mu.component], lty = 2, col = "red")
legend("topright", legend = rho.legend, lty = seq_along(rho.tags), bty = "o")
dev.off()
