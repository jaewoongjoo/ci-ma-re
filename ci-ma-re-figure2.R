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
n <- p + 1
c1 <- 0.3
Sigma.true <- 5 * matrix(c(1, c1, c1, c1, c1, c1, c1, c1, c1, c1, 1, c1, c1, c1, c1, c1, 
    c1, c1, c1, c1, 1, c1, c1, c1, c1, c1, c1, c1, c1, c1, 1, c1, c1, c1, c1, c1, c1, c1, 
    c1, c1, 1, c1, c1, c1, c1, c1, c1, c1, c1, c1, 1, c1, c1, c1, c1, c1, c1, c1, c1, c1, 
    1, c1, c1, c1, c1, c1, c1, c1, c1, c1, 1, c1, c1, c1, c1, c1, c1, c1, c1, c1, 1), nrow = n, 
    ncol = n)

summ <- setNames(vector("list", length(rho.tags)), rho.tags)
for (tag in rho.tags) {
    scenario.results <- simulation.results[[tag]]
    Sigma.hat <- scenario.results$Sigma_hat_arr
    omega.opt <- scenario.results$omega_opt
    replicate.status <- scenario.results$replicate_status
    ok <- replicate.status == "ok" & is.finite(omega.opt)
    nsim.ok <- sum(ok)
    Sigma.hat.ok <- Sigma.hat[, , ok, drop = FALSE]
    mean.Sigma <- apply(Sigma.hat.ok, c(1, 2), mean)
    summ[[tag]] <- list(nsim_ok = nsim.ok, mean_sigma = mean.Sigma, mse_sigma = apply((Sigma.hat.ok - 
        array(rep(Sigma.true, nsim.ok), dim = dim(Sigma.hat.ok)))^2, c(1, 2), mean), var_sigma = apply((Sigma.hat.ok - 
        array(rep(mean.Sigma, nsim.ok), dim = dim(Sigma.hat.ok)))^2, c(1, 2), mean), Sigma_hat_arr = Sigma.hat.ok)
}

pch.values <- c(0, 1, 2, 5)
pch.values.filled <- c(15, 16, 17, 18)
pt.cex <- 1.4
sig.diag.lab <- as.expression(lapply(1:p, function(i) bquote(sigma[.(i)]^2)))
sig.off.lab <- as.expression(list(bquote(sigma[23]), bquote(sigma[35])))
sig.all.lab <- c(sig.diag.lab, sig.off.lab)
rho.legend <- as.expression(lapply(rho.vals, function(r) {
    bquote(rho[X] == .(as.numeric(r)))
}))

diag.index <- 2:n
off.index <- list(c(3, 4), c(4, 6))
pick.Sigma <- function(M) {
    c(vapply(diag.index, function(i) M[i, i], numeric(1)), vapply(off.index, function(ij) M[ij[1], 
        ij[2]], numeric(1)))
}

true.Sigma.selected <- pick.Sigma(Sigma.true)
x.position <- seq_along(true.Sigma.selected)
Sigma.means <- do.call(rbind, lapply(rho.tags, function(tag) {
    pick.Sigma(summ[[tag]]$mean_sigma)
}))

pdf(file.path(figures.dir, "fig2_Sigma_mean.pdf"), width = 5, height = 5)
par(mar = c(4, 3, 3, 1), mgp = c(0, 1.5, 0))
plot(NULL, xlim = range(x.position), ylim = c(0, 10), xlab = "", ylab = "", xaxt = "n")
axis(1, at = x.position, labels = sig.all.lab, cex.axis = 1.25)
points(x.position, true.Sigma.selected, pch = 3, col = "red", cex = 1.6)
for (i in seq_along(rho.tags)) {
    points(x.position, Sigma.means[i, ], pch = pch.values[i], cex = pt.cex)
}

legend("topright", legend = c("True", rho.legend), pch = c(3, pch.values), col = c("red", 
    rep("black", length(rho.tags))), pt.cex = c(1.6, rep(pt.cex, length(rho.tags))))

dev.off()
Sigma.sqrt.MSE <- do.call(rbind, lapply(rho.tags, function(tag) {
    sqrt(pmax(pick.Sigma(summ[[tag]]$mse_sigma), 0))
}))

Sigma.SD <- do.call(rbind, lapply(rho.tags, function(tag) {
    sqrt(pmax(pick.Sigma(summ[[tag]]$var_sigma), 0))
}))

pdf(file.path(figures.dir, "fig2_Sigma_mse.pdf"), width = 5, height = 5)
par(mar = c(4, 3, 3, 1), mgp = c(0, 1.5, 0))
plot(NULL, xlim = range(x.position), ylim = c(0, 8), xlab = "", ylab = "", xaxt = "n")
axis(1, at = x.position, labels = sig.all.lab, cex.axis = 1.25)
for (i in seq_along(rho.tags)) {
    points(x.position, Sigma.sqrt.MSE[i, ], pch = pch.values[i], cex = pt.cex)
    points(x.position, Sigma.SD[i, ], pch = pch.values.filled[i], col = "red", cex = 0.9)
}

legend("topright", legend = rho.legend, pch = pch.values, pt.cex = pt.cex)
dev.off()
Sigma.component <- 2
pdf(file.path(figures.dir, "fig2_logSigma_density_mm.pdf"), width = 5, height = 5)
par(mar = c(4, 3, 3, 1))
plot(NULL, xlim = c(-1, 4), ylim = c(0, 0.9), xlab = expression(log ~ sigma[1]^2), ylab = "")
for (i in seq_along(rho.tags)) {
    values <- summ[[rho.tags[i]]]$Sigma_hat_arr[Sigma.component, Sigma.component, ]
    values <- log(values[is.finite(values) & values > 0])
    density.fit <- density(values)
    lines(density.fit$x, density.fit$y, lwd = 2, lty = i)
}

abline(v = log(Sigma.true[Sigma.component, Sigma.component]), lty = 2, col = "red")
legend("topright", legend = rho.legend, lty = seq_along(rho.tags), bty = "o")
dev.off()
