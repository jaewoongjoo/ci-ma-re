packages <- c("MASS", "stats", "base", "mvtnorm", "foreach", "doParallel", "Matrix", "matrixcalc", 
    "expm", "rstudioapi")

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
output.dir <- file.path(code.dir, "results")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
num_cores <- max(1, parallel::detectCores() - 1)
k <- 20
n_external <- 50000
p <- 8
kap <- 10
d_alpha <- 0.5
c_armijo <- 0.5
step_init <- 0.001
step_min <- 1e-15
gamma_init_scale <- 1
gamma_init_offdiag <- 0
gamma_init_structure <- "dense"
alpha_level <- 0.05
zcrit <- qnorm(1 - alpha_level/2)
mu_true <- c(-3, 3, 1, -1, -3, 3, 1, -1, -3)
c1 <- 0.3
c2 <- c1
c3 <- c1
c4 <- c1
c5 <- c1
c6 <- c1
c7 <- c1
c8 <- c1
Sigma_true <- 5 * matrix(c(1, c1, c2, c3, c4, c5, c6, c7, c8, c1, 1, c1, c2, c3, c4, c5, 
    c6, c7, c2, c1, 1, c1, c2, c3, c4, c5, c6, c3, c2, c1, 1, c1, c2, c3, c4, c5, c4, c3, 
    c2, c1, 1, c1, c2, c3, c4, c5, c4, c3, c2, c1, 1, c1, c2, c3, c6, c5, c4, c3, c2, c1, 
    1, c1, c2, c7, c6, c5, c4, c3, c2, c1, 1, c1, c8, c7, c6, c5, c4, c3, c2, c1, 1), nrow = p + 
    1, ncol = p + 1)

mu_x <- c(-1, 0, 1, -1, 1, -3, 1, 3)
num_study <- c(rep(400, k/4), rep(500, k/4), rep(600, k/4), rep(700, k/4))
patterns <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(1, 2, 7, 8), c(3, 4, 5, 6))
num_patterns <- length(patterns)
A_list <- vector("list", k)
for (kk in seq_len(k)) {
    A_list[[kk]] <- patterns[[((kk - 1)%%num_patterns) + 1]][patterns[[((kk - 1)%%num_patterns) + 
        1]] <= p]
}

nsim <- 5000
base_seed <- 20260716
omega_grid <- c(0.01, seq(0.1, 3, 0.1))
rho_values <- c(0.3, 0.6, 0.8, 0.9)
rho_tags <- c("corr0.3", "corr0.6", "corr0.8", "corr0.9")
commutation_matrix <- function(n) {
    K <- matrix(0, n * n, n * n)
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            K[(j - 1) * n + i, (i - 1) * n + j] <- 1
        }
    }
    K
}

Sym <- function(A) (A + t(A))/2
safe_log_determinant <- function(A, jitter = 0) {
    A <- Sym(A)
    n <- nrow(A)
    tryCatch({
        R <- chol(A + diag(jitter, n))
        2 * sum(log(diag(R)))
    }, error = function(e) {
        as.numeric(determinant(A + diag(jitter, n), logarithm = TRUE)$modulus)
    })
}

safe_inverse <- function(A, jitter = 0) {
    A <- Sym(A)
    n <- nrow(A)
    tryCatch({
        R <- chol(A + diag(jitter, n))
        chol2inv(R)
    }, error = function(e) {
        solve(A + diag(jitter, n))
    })
}

vech <- function(A) A[lower.tri(A, diag = TRUE)]
vech_to_mat <- function(v) {
    n <- (-1 + sqrt(1 + 8 * length(v)))/2
    mat <- matrix(0, nrow = n, ncol = n)
    mat[lower.tri(mat, diag = TRUE)] <- v
    mat
}

make_gamma_init <- function(n) {
    g <- matrix(gamma_init_offdiag, nrow = n, ncol = n)
    diag(g) <- gamma_init_scale
    g
}

diag_vech_indices <- function(n) {
    j <- seq_len(n)
    as.integer(1 + (j - 1) * (n + 1) - (j - 1) * j/2)
}

log_T_func <- function(mu_vec, vech_gam, M, c_vec, W, theta) {
    gam <- vech_to_mat(vech_gam)
    G <- gam %*% t(gam)
    log_det_G <- safe_log_determinant(G)
    val <- numeric(dim(M)[3])
    for (kk in seq_len(dim(M)[3])) {
        A <- M[, , kk] + G
        v <- c_vec[, kk] + G %*% mu_vec
        val[kk] <- log_det_G - safe_log_determinant(A) - t(mu_vec) %*% G %*% mu_vec + t(v) %*% 
            safe_inverse(A) %*% v - t(theta[[kk]]) %*% W[[kk]] %*% theta[[kk]]
    }
    0.5 * sum(val)
}

lag_func <- function(mu_vec, vech_gam, M, c_vec, W, theta, rho_vec, kap) {
    gam <- vech_to_mat(vech_gam)
    pen <- 0
    for (pp in seq_along(rho_vec)) {
        pen <- pen + rho_vec[pp] * (1/kap) * log(1 + exp(kap * (-gam[pp, pp])))
    }
    -log_T_func(mu_vec, vech_gam, M, c_vec, W, theta) + pen
}

vec_gradient <- function(vech_gam, mu_vec, M, c_vec, rho_vec, kap, p, k) {
    gam <- vech_to_mat(vech_gam)
    n <- p + 1
    G <- gam %*% t(gam)
    Ginv <- safe_inverse(G)
    I <- diag(n)
    I2 <- diag(n^2)
    Kc <- commutation_matrix(n)
    make_P_lt <- function(n) {
        m <- n * (n + 1)/2
        P <- matrix(0, n * n, m)
        idx <- 1
        for (j in 1:n) {
            for (i in j:n) {
                e <- matrix(0, n, n)
                e[i, j] <- 1
                P[, idx] <- as.vector(e)
                idx <- idx + 1
            }
        }
        P
    }
    P_lt <- make_P_lt(n)
    tmp_sum <- numeric(n^2)
    for (kk in 1:k) {
        A <- M[, , kk] + G
        Ainv <- safe_inverse(A)
        v <- c_vec[, kk] + G %*% mu_vec
        Av <- Ainv %*% v
        tmp_sum <- tmp_sum + as.vector(Ginv) - as.vector(Ainv) - as.vector(kronecker(mu_vec, 
            mu_vec)) + as.vector(kronecker(mu_vec, Av)) + as.vector(kronecker(Av, mu_vec)) - 
            as.vector(kronecker(Av, Av))
    }
    S <- 0.5 * ((I2 + Kc) %*% (kronecker(gam, I)) %*% P_lt)
    d_log_T_d_vech_gam <- t(S) %*% tmp_sum
    d_penalty <- numeric(length(vech_gam))
    diag_idx <- diag_vech_indices(n)
    for (pp in 1:n) {
        idx <- diag_idx[pp]
        d_penalty[idx] <- rho_vec[pp] * (plogis(kap * gam[pp, pp]) - 1)
    }
    as.vector(-d_log_T_d_vech_gam + d_penalty)
}

optimize_rho <- function(rho_vec, mu_vec, vech_gam, M, c_vec, W, theta, kap, p, eta = 0.1, 
    lower = 1, upper = 3) {
    gam <- vech_to_mat(vech_gam)
    n <- p + 1
    violation <- (1/kap) * log(1 + exp(kap * (-diag(gam))))
    out <- tryCatch(optim(par = rho_vec, fn = function(rho) {
        -sum(rho * violation) + (1/(2 * eta)) * sum((rho - rho_vec)^2)
    }, method = "L-BFGS-B", lower = rep(lower, n), upper = rep(upper, n)), error = function(e) {
        list(par = rho_vec)
    })
    out$par
}

build_data_object <- function(seed_i, sigma2_x) {
    set.seed(seed_i)
    n <- p + 1
    bet_true <- mvrnorm(n = k, mu = mu_true, Sigma = Sigma_true)
    theta <- vector("list", k)
    S_inv <- vector("list", k)
    x_full <- vector("list", k)
    x_red <- vector("list", k)
    x_full_ext <- vector("list", k)
    x_red_ext <- vector("list", k)
    y <- vector("list", k)
    for (kk in 1:k) {
        x_full[[kk]] <- mvrnorm(n = num_study[kk], mu = mu_x, Sigma = sigma2_x)
        y[[kk]] <- cbind(1, x_full[[kk]]) %*% bet_true[kk, ] + 3 * rt(n = num_study[kk], 
            df = 4.1)
        x_red[[kk]] <- cbind(1, x_full[[kk]][, A_list[[kk]]])
        theta[[kk]] <- safe_inverse(t(x_red[[kk]]) %*% x_red[[kk]]) %*% t(x_red[[kk]]) %*% 
            y[[kk]]
        res <- as.numeric(y[[kk]] - x_red[[kk]] %*% theta[[kk]])
        sigma_hat_k_2 <- sum(res^2)/(num_study[kk] - ncol(x_red[[kk]]))
        S_inv[[kk]] <- (1/sigma_hat_k_2) * (1/num_study[kk]) * t(x_red[[kk]]) %*% x_red[[kk]]
        x_full_ext[[kk]] <- cbind(1, mvrnorm(n = n_external, mu = mu_x, Sigma = sigma2_x))
        x_red_ext[[kk]] <- x_full_ext[[kk]][, c(1, (A_list[[kk]] + 1))]
    }
    list(bet_true = bet_true, theta = theta, S_inv = S_inv, x_full_ext = x_full_ext, x_red_ext = x_red_ext)
}

build_M_c_W_given_omega <- function(obj, omega_value) {
    n <- p + 1
    M <- array(0, c(n, n, k))
    c_vec <- array(0, c(n, k))
    W <- vector("list", k)
    for (kk in 1:k) {
        S_inv_adj <- num_study[kk] * obj$S_inv[[kk]]
        M[, , kk] <- (8 * omega_value)/(n_external^2) * t(obj$x_full_ext[[kk]]) %*% obj$x_red_ext[[kk]] %*% 
            S_inv_adj %*% t(obj$x_red_ext[[kk]]) %*% obj$x_full_ext[[kk]]
        c_vec[, kk] <- (8 * omega_value)/(n_external^2) * t(obj$x_full_ext[[kk]]) %*% obj$x_red_ext[[kk]] %*% 
            S_inv_adj %*% t(obj$x_red_ext[[kk]]) %*% obj$x_red_ext[[kk]] %*% obj$theta[[kk]]
        W[[kk]] <- (8 * omega_value)/(n_external^2) * t(obj$x_red_ext[[kk]]) %*% obj$x_red_ext[[kk]] %*% 
            S_inv_adj %*% t(obj$x_red_ext[[kk]]) %*% obj$x_red_ext[[kk]]
    }
    list(M = M, c_vec = c_vec, W = W)
}

predvar_objective_given_fit <- function(obj, hat_mu, hat_gamma, omega_value) {
    n <- p + 1
    m <- n * (n + 1)/2
    MCW <- build_M_c_W_given_omega(obj, omega_value)
    M <- MCW$M
    c_vec <- MCW$c_vec
    W <- MCW$W
    theta <- obj$theta
    Gam <- hat_gamma
    G <- Gam %*% t(Gam)
    I_n <- diag(n)
    I2 <- diag(n^2)
    Kc <- commutation_matrix(n)
    Ssym <- 0.5 * (I2 + Kc)
    make_P_lt <- function(n) {
        m <- n * (n + 1)/2
        P <- matrix(0, n * n, m)
        idx <- 1
        for (j in 1:n) {
            for (i in j:n) {
                e <- matrix(0, n, n)
                e[i, j] <- 1
                P[, idx] <- as.vector(e)
                idx <- idx + 1
            }
        }
        P
    }
    P_lt <- make_P_lt(n)
    S_Gamma_fix <- (I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt
    rho_hat_local <- optimize_rho(rep(1, n), hat_mu, vech(Gam), M, c_vec, W, theta, kap, 
        p)
    Vk_list <- vector("list", k)
    mk_list <- vector("list", k)
    for (kk in 1:k) {
        Ak <- M[, , kk] + G
        Vk <- safe_inverse(Ak)
        Vk_list[[kk]] <- Vk
        mk_list[[kk]] <- Vk %*% (c_vec[, kk] + G %*% hat_mu)
    }
    J <- matrix(0, n + m, n + m)
    for (kk in 1:k) {
        Vk <- Vk_list[[kk]]
        vk <- c_vec[, kk] + G %*% hat_mu
        ak <- Vk %*% vk
        rk <- hat_mu - ak
        J_mm_k <- G %*% Vk %*% G - G
        J_mG_vec_k <- (t(ak - hat_mu) %x% diag(n)) + (t(hat_mu - ak) %x% (G %*% Vk))
        J_mg_k <- J_mG_vec_k %*% S_Gamma_fix
        M_dr_mu <- (diag(n) - Vk %*% G)
        d_rrT_mu <- (kronecker(I_n, rk) + kronecker(rk, I_n)) %*% M_dr_mu
        J_gm_k <- t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% d_rrT_mu)
        Ginv <- safe_inverse(G)
        tmp_k <- as.vector(Ginv) - as.vector(Vk) - as.vector(hat_mu %*% t(hat_mu)) + as.vector(hat_mu %*% 
            t(ak)) + as.vector(ak %*% t(hat_mu)) - as.vector(ak %*% t(ak))
        z <- Ssym %*% tmp_k
        Zmat <- matrix(z, n, n)
        TermA <- -t(P_lt) %*% (diag(n) %x% Zmat) %*% P_lt
        M_dr_vec <- -(t(rk) %x% Vk)
        d_rrT_vec <- (kronecker(I_n, rk) + kronecker(rk, I_n)) %*% M_dr_vec
        M_tmp_vec <- (-(t(Ginv) %x% Ginv) + (t(Vk) %x% Vk) - d_rrT_vec)
        TermB <- -t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% (M_tmp_vec %*% S_Gamma_fix))
        J_gg_k <- TermA + TermB
        J[1:n, 1:n] <- J[1:n, 1:n] + J_mm_k
        J[1:n, (n + 1):(n + m)] <- J[1:n, (n + 1):(n + m)] + J_mg_k
        J[(n + 1):(n + m), 1:n] <- J[(n + 1):(n + m), 1:n] + J_gm_k
        J[(n + 1):(n + m), (n + 1):(n + m)] <- J[(n + 1):(n + m), (n + 1):(n + m)] + J_gg_k
    }
    H_pen <- matrix(0, m, m)
    diag_idx <- diag_vech_indices(n)
    for (ii in 1:n) {
        idx <- diag_idx[ii]
        pii <- plogis(kap * Gam[ii, ii])
        H_pen[idx, idx] <- rho_hat_local[ii] * kap * pii * (1 - pii)
    }
    J[(n + 1):(n + m), (n + 1):(n + m)] <- J[(n + 1):(n + m), (n + 1):(n + m)] + H_pen
    B_meat <- matrix(0, n + m, n + m)
    for (kk in 1:k) {
        Vk <- Vk_list[[kk]]
        mk <- mk_list[[kk]]
        s_mu_k <- G %*% Vk %*% (c_vec[, kk] + G %*% hat_mu) - G %*% hat_mu
        Ginv <- safe_inverse(G)
        tmp_k <- as.vector(Ginv) - as.vector(Vk) - as.vector(hat_mu %*% t(hat_mu)) + as.vector(hat_mu %*% 
            t(mk)) + as.vector(mk %*% t(hat_mu)) - as.vector(mk %*% t(mk))
        s_vGam_k <- -as.vector(t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% tmp_k))
        s_k <- c(s_mu_k, s_vGam_k)
        B_meat <- B_meat + s_k %*% t(s_k)
    }
    J_inv <- solve(J)
    Cov_theta_hat_Gam <- J_inv %*% B_meat %*% t(J_inv)
    total_pluginvar <- 0
    for (kk in 1:k) {
        Vk <- Vk_list[[kk]]
        mk <- mk_list[[kk]]
        J_mu_k <- Vk %*% G
        b <- as.numeric(hat_mu - mk)
        J_vGam_k <- (t(b) %x% Vk) %*% S_Gamma_fix
        J_theta_k <- cbind(J_mu_k, J_vGam_k)
        PluginVar_k <- J_theta_k %*% Cov_theta_hat_Gam %*% t(J_theta_k)
        total_pluginvar <- total_pluginvar + sum(diag(PluginVar_k))
    }
    return(total_pluginvar)
}

run_one_with_obj <- function(obj, omega_val, current_niter) {
    n <- p + 1
    bet_true <- obj$bet_true
    theta <- obj$theta
    S_inv <- obj$S_inv
    x_full_ext <- obj$x_full_ext
    x_red_ext <- obj$x_red_ext
    M <- array(0, c(n, n, k))
    c_vec <- array(0, c(n, k))
    W <- vector("list", k)
    for (kk in 1:k) {
        S_inv_adj <- num_study[kk] * S_inv[[kk]]
        M[, , kk] <- (8 * omega_val)/(n_external^2) * t(x_full_ext[[kk]]) %*% x_red_ext[[kk]] %*% S_inv_adj %*% 
            t(x_red_ext[[kk]]) %*% x_full_ext[[kk]]
        c_vec[, kk] <- (8 * omega_val)/(n_external^2) * t(x_full_ext[[kk]]) %*% x_red_ext[[kk]] %*% 
            S_inv_adj %*% t(x_red_ext[[kk]]) %*% x_red_ext[[kk]] %*% theta[[kk]]
        W[[kk]] <- (8 * omega_val)/(n_external^2) * t(x_red_ext[[kk]]) %*% x_red_ext[[kk]] %*% S_inv_adj %*% 
            t(x_red_ext[[kk]]) %*% x_red_ext[[kk]]
    }
    M_w_optimal <- M
    c_vec_w_optimal <- c_vec
    tilde_mu <- matrix(0, nrow = current_niter, ncol = n)
    tilde_gamma <- array(0, c(current_niter, n, n))
    rho_mat <- matrix(0, nrow = current_niter, ncol = n)
    rho_mat[1, ] <- rep(1, n)
    tilde_mu[1, ] <- rep(1, n)
    tilde_gamma[1, , ] <- make_gamma_init(n)
    j_final <- current_niter
    converged <- FALSE
    line_search_failed <- FALSE
    for (j in 1:(current_niter - 1)) {
        gam <- tilde_gamma[j, , ]
        G <- gam %*% t(gam)
        temp_mu_1 <- matrix(0, n, n)
        temp_mu_2 <- rep(0, n)
        for (kk in 1:k) {
            A <- M[, , kk] + G
            Ainv <- safe_inverse(A)
            temp_mu_1 <- temp_mu_1 + (G - G %*% Ainv %*% G)
            temp_mu_2 <- temp_mu_2 + (G %*% Ainv %*% c_vec[, kk])
        }
        tilde_mu[j + 1, ] <- as.vector(safe_inverse(temp_mu_1) %*% temp_mu_2)
        z_old <- vech(tilde_gamma[j, , ])
        grad_val <- vec_gradient(z_old, tilde_mu[j + 1, ], M, c_vec, rho_mat[j, ], kap, p, 
            k)
        t_k <- step_init
        func_old <- lag_func(tilde_mu[j + 1, ], z_old, M, c_vec, W, theta, rho_mat[j, ], 
            kap)
        line_search_failed_iter <- FALSE
        repeat {
            z_try <- z_old - t_k * grad_val
            func_new <- lag_func(tilde_mu[j + 1, ], z_try, M, c_vec, W, theta, rho_mat[j, 
                ], kap)
            if (!is.finite(func_new)) {
                t_k <- d_alpha * t_k
                if (t_k < step_min) {
                  t_k <- 0
                  line_search_failed_iter <- TRUE
                  break
                }
                next
            }
            if (func_new <= (func_old - c_armijo * t_k * sum(grad_val^2))) 
                break
            t_k <- d_alpha * t_k
            if (t_k < step_min) {
                t_k <- 0
                line_search_failed_iter <- TRUE
                break
            }
        }
        if (line_search_failed_iter) {
            tilde_gamma[j + 1, , ] <- tilde_gamma[j, , ]
            z_cand <- z_old
        }
        else {
            z_cand <- z_old - t_k * grad_val
            tilde_gamma[j + 1, , ] <- vech_to_mat(z_cand)
        }
        rho_mat[j + 1, ] <- optimize_rho(rho_mat[j, ], tilde_mu[j + 1, ], z_cand, M, c_vec, 
            W, theta, kap, p)
        delta_mu <- tilde_mu[j + 1, -1] - tilde_mu[j, -1]
        delta_gamma <- as.vector(tilde_gamma[j + 1, -1, -1] - tilde_gamma[j, -1, -1])
        max_diff_mu <- max(abs(delta_mu))
        max_diff_gamma <- max(abs(delta_gamma))
        max_diff <- max(max_diff_mu, max_diff_gamma)
        conv_ok <- max_diff < 0.001
        if (conv_ok) {
            j_final <- j + 1
            converged <- TRUE
            break
        }
    }
    hat_mu <- tilde_mu[j_final, ]
    hat_gamma <- tilde_gamma[j_final, , ]
    hat_sigma <- safe_inverse(hat_gamma %*% t(hat_gamma))
    rho_hat <- optimize_rho(rho_mat[j_final, ], hat_mu, vech(hat_gamma), M, c_vec, W, theta, 
        kap, p)
    n <- p + 1
    m <- n * (n + 1)/2
    Gam <- hat_gamma
    G <- Gam %*% t(Gam)
    I_n <- diag(n)
    I2 <- diag(n^2)
    Kc <- commutation_matrix(n)
    Ssym <- 0.5 * (I2 + Kc)
    make_P_lt <- function(n) {
        m <- n * (n + 1)/2
        P <- matrix(0, n * n, m)
        idx <- 1
        for (j in 1:n) for (i in j:n) {
            e <- matrix(0, n, n)
            e[i, j] <- 1
            P[, idx] <- as.vector(e)
            idx <- idx + 1
        }
        P
    }
    P_lt <- make_P_lt(n)
    S_Gamma_fix <- (I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt
    Vk_list <- vector("list", k)
    mk_list <- vector("list", k)
    for (kk in 1:k) {
        Ak <- M_w_optimal[, , kk] + G
        Vk <- safe_inverse(Ak)
        mk <- Vk %*% (c_vec_w_optimal[, kk] + G %*% hat_mu)
        Vk_list[[kk]] <- Vk
        mk_list[[kk]] <- mk
    }
    J <- matrix(0, n + m, n + m)
    for (kk in 1:k) {
        Ak <- M_w_optimal[, , kk] + G
        Vk <- Vk_list[[kk]]
        vk <- c_vec_w_optimal[, kk] + G %*% hat_mu
        ak <- Vk %*% vk
        rk <- hat_mu - ak
        J_mm_k <- G %*% Vk %*% G - G
        J_mG_vec_k <- (t(ak - hat_mu) %x% diag(n)) + (t(hat_mu - ak) %x% (G %*% Vk))
        J_mg_k <- J_mG_vec_k %*% S_Gamma_fix
        M_dr_mu <- (diag(n) - Vk %*% G)
        d_rrT_mu <- (kronecker(I_n, rk) + kronecker(rk, I_n)) %*% M_dr_mu
        J_gm_k <- t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% d_rrT_mu)
        tmp_k <- as.vector(safe_inverse(G)) - as.vector(Vk) - as.vector(hat_mu %*% t(hat_mu)) + 
            as.vector(hat_mu %*% t(ak)) + as.vector(ak %*% t(hat_mu)) - as.vector(ak %*% 
            t(ak))
        z <- Ssym %*% tmp_k
        Zmat <- matrix(z, n, n)
        TermA <- -t(P_lt) %*% (diag(n) %x% Zmat) %*% P_lt
        Ginv <- safe_inverse(G)
        M_dr_vec <- -(t(rk) %x% Vk)
        d_rrT_vec <- (kronecker(I_n, rk) + kronecker(rk, I_n)) %*% M_dr_vec
        M_tmp_vec <- (-(t(Ginv) %x% Ginv) + (t(Vk) %x% Vk) - d_rrT_vec)
        TermB <- -t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% (M_tmp_vec %*% S_Gamma_fix))
        J_gg_k <- TermA + TermB
        J[1:n, 1:n] <- J[1:n, 1:n] + J_mm_k
        J[1:n, (n + 1):(n + m)] <- J[1:n, (n + 1):(n + m)] + J_mg_k
        J[(n + 1):(n + m), 1:n] <- J[(n + 1):(n + m), 1:n] + J_gm_k
        J[(n + 1):(n + m), (n + 1):(n + m)] <- J[(n + 1):(n + m), (n + 1):(n + m)] + J_gg_k
    }
    H_pen <- matrix(0, m, m)
    diag_idx <- diag_vech_indices(n)
    for (ii in 1:n) {
        idx <- diag_idx[ii]
        pii <- plogis(kap * Gam[ii, ii])
        H_pen[idx, idx] <- rho_hat[ii] * kap * pii * (1 - pii)
    }
    J[(n + 1):(n + m), (n + 1):(n + m)] <- J[(n + 1):(n + m), (n + 1):(n + m)] + H_pen
    B <- matrix(0, n + m, n + m)
    for (kk in 1:k) {
        Ak <- M_w_optimal[, , kk] + G
        Vk <- Vk_list[[kk]]
        vk <- c_vec_w_optimal[, kk] + G %*% hat_mu
        ak <- Vk %*% vk
        s_mu_k <- G %*% Vk %*% vk - G %*% hat_mu
        Ginv <- safe_inverse(G)
        tmp_k <- as.vector(Ginv) - as.vector(Vk) - as.vector(hat_mu %*% t(hat_mu)) + as.vector(hat_mu %*% 
            t(mk)) + as.vector(mk %*% t(hat_mu)) - as.vector(mk %*% t(mk))
        s_vGam_k <- -as.vector(t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% tmp_k))
        s_k <- c(s_mu_k, s_vGam_k)
        B <- B + s_k %*% t(s_k)
    }
    J_inv <- solve(J)
    Cov_theta_hat_Gam <- J_inv %*% B %*% t(J_inv)
    hat_beta <- array(0, c(n, k))
    beta_covered_pred <- matrix(0, nrow = n, ncol = k)
    for (kk in 1:k) {
        Vk <- Vk_list[[kk]]
        mk <- mk_list[[kk]]
        hat_beta[, kk] <- mk
        J_mu_k <- Vk %*% G
        b <- as.numeric(hat_mu - mk)
        J_vGam_k <- (t(b) %x% Vk) %*% ((I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt)
        J_theta_k <- cbind(J_mu_k, J_vGam_k)
        PluginVar_k <- J_theta_k %*% Cov_theta_hat_Gam %*% t(J_theta_k)
        PredVar_k <- Vk + PluginVar_k
        se_pred <- sqrt(pmax(diag(PredVar_k), 0))
        pi_low <- mk - zcrit * se_pred
        pi_high <- mk + zcrit * se_pred
        beta_true_vec <- bet_true[kk, ]
        inside <- (beta_true_vec >= pi_low) & (beta_true_vec <= pi_high)
        beta_covered_pred[, kk] <- as.integer(inside)
    }
    list(hat_mu = hat_mu, hat_beta = hat_beta, hat_gamma = hat_gamma, hat_sigma = hat_sigma, 
        beta_covered = beta_covered_pred, bet_true = bet_true, j_final = j_final, converged = converged, 
        rho_mat = rho_mat)
}

run_two_stage <- function(seed_i, grid = omega_grid, current_niter, sigma2_x) {
    obj_stage1 <- build_data_object(seed_i, sigma2_x)
    fit_stage1 <- tryCatch(run_one_with_obj(obj_stage1, omega_val = 1, current_niter), error = function(e) NULL)
    if (is.null(fit_stage1)) 
        return(list(fit_stage2 = NULL, omega_opt = NA))
    best_obj <- Inf
    best_omega <- 1
    for (w in grid) {
        obj_val <- tryCatch(predvar_objective_given_fit(obj_stage1, fit_stage1$hat_mu, fit_stage1$hat_gamma, 
            w), error = function(e) Inf)
        if (is.finite(obj_val) && obj_val < best_obj) {
            best_obj <- obj_val
            best_omega <- w
        }
    }
    fit_stage2 <- tryCatch(run_one_with_obj(obj_stage1, omega_val = best_omega, current_niter), 
        error = function(e) NULL)
    if (is.null(fit_stage2)) 
        fit_stage2 <- fit_stage1
    list(fit_stage2 = fit_stage2, omega_opt = best_omega)
}

run_simulation <- function() {
    n <- p + 1
    all_chunks <- list()
    for (rho_idx in seq_along(rho_values)) {
        rho <- rho_values[rho_idx]
        current_niter <- 1500
        r1 <- rho
        r2 <- r1
        r3 <- r1
        r4 <- r1
        r5 <- r1
        r6 <- r1
        r7 <- r1
        sigma2_x <- 5 * matrix(c(1, r1, r2, r3, r4, r5, r6, r7, r1, 1, r1, r2, r3, r4, r5, 
            r6, r2, r1, 1, r1, r2, r3, r4, r5, r3, r2, r1, 1, r1, r2, r3, r4, r4, r3, r2, 
            r1, 1, r1, r2, r3, r5, r4, r3, r2, r1, 1, r1, r2, r6, r5, r4, r3, r2, r1, 1, 
            r1, r7, r6, r5, r4, r3, r2, r1, 1), nrow = p, ncol = p)
        mu_hat_mat <- matrix(NA_real_, nrow = nsim, ncol = n)
        Sigma_hat_arr <- array(NA_real_, dim = c(n, n, nsim))
        coverage_sum <- matrix(0, nrow = n, ncol = k)
        omega_opt_vec <- rep(NA_real_, nsim)
        replicate_status <- rep(NA_character_, nsim)
        replicate_error <- rep(NA_character_, nsim)
        jfinal_vec <- rep(NA_real_, nsim)
        converged_vec <- rep(NA, nsim)
        results <- foreach::foreach(ii = seq_len(nsim), .packages = c("MASS", "Matrix"), 
            .inorder = TRUE) %dopar% {
            seed_i <- base_seed + ceiling(ii/50) * 1000000 + 1000 + ((ii - 1)%%50 + 1)
            tryCatch({
                res <- run_two_stage(seed_i, grid = omega_grid, current_niter = current_niter, 
                  sigma2_x = sigma2_x)
                fit2 <- res$fit_stage2
                list(ok = TRUE, mu_hat = fit2$hat_mu, sigma_hat = fit2$hat_sigma, beta_covered = fit2$beta_covered, 
                  omega_opt = res$omega_opt, j_final = fit2$j_final, converged = fit2$converged)
            }, error = function(e) {
                list(ok = FALSE, msg = conditionMessage(e))
            })
        }
        for (ii in seq_len(nsim)) {
            r <- results[[ii]]
            if (isTRUE(r$ok)) {
                mu_hat_mat[ii, ] <- r$mu_hat
                Sigma_hat_arr[, , ii] <- r$sigma_hat
                coverage_sum <- coverage_sum + r$beta_covered
                omega_opt_vec[ii] <- r$omega_opt
                jfinal_vec[ii] <- r$j_final
                converged_vec[ii] <- r$converged
                replicate_status[ii] <- "ok"
            }
            else {
                replicate_status[ii] <- "error"
                replicate_error[ii] <- r$msg
            }
        }
        ok <- replicate_status == "ok" & is.finite(omega_opt_vec)
        nsim_ok <- sum(ok)
        mu_hat_ok <- mu_hat_mat[ok, , drop = FALSE]
        Sigma_hat_ok <- Sigma_hat_arr[, , ok, drop = FALSE]
        mean_mu <- colMeans(mu_hat_ok)
        var_mu <- apply(mu_hat_ok, 2, var)
        mse_mu <- colMeans((mu_hat_ok - matrix(rep(mu_true, nsim_ok), nrow = nsim_ok, 
            ncol = n, byrow = TRUE))^2)
        mean_Sigma <- apply(Sigma_hat_ok, c(1, 2), mean)
        mse_sigma <- apply((Sigma_hat_ok - array(rep(Sigma_true, nsim_ok), dim = dim(Sigma_hat_ok)))^2, 
            c(1, 2), mean)
        var_sigma <- apply((Sigma_hat_ok - array(rep(mean_Sigma, nsim_ok), dim = dim(Sigma_hat_ok)))^2, 
            c(1, 2), mean)
        coverage_beta_matrix <- coverage_sum/nsim_ok
        rownames(coverage_beta_matrix) <- paste0("beta[", 0:p, "]")
        colnames(coverage_beta_matrix) <- paste0("study_", 1:k)
        all_chunks[[rho_tags[rho_idx]]] <- list(scenario = rho_tags[rho_idx], rho_x = rho, 
            nsim_ok = nsim_ok, mean_mu = mean_mu, var_mu = var_mu, mse_mu = mse_mu, mean_sigma = mean_Sigma, 
            mse_sigma = mse_sigma, var_sigma = var_sigma, coverage_beta = coverage_beta_matrix, 
            omega_opt = omega_opt_vec, jfinal_vec = jfinal_vec, converged_vec = converged_vec, 
            mu_hat_mat = mu_hat_mat, Sigma_hat_arr = Sigma_hat_arr, replicate_status = replicate_status, 
            replicate_error = replicate_error)
    }
    summ <- all_chunks
    omega_summary_df <- do.call(rbind, lapply(seq_along(rho_tags), function(i) {
        x <- summ[[rho_tags[i]]]$omega_opt
        data.frame(rho_x = rho_values[i], min = min(x, na.rm = TRUE), q25 = quantile(x, 0.25, 
            na.rm = TRUE), median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE), 
            q75 = quantile(x, 0.75, na.rm = TRUE), max = max(x, na.rm = TRUE))
    }))
    output_dir <- output.dir
    out_rds <- file.path(output_dir, "ci-ma-re-simulation-results.rds")
    saveRDS(all_chunks, file = out_rds)
    list(results = all_chunks, omega_summary = omega_summary_df)
}

cl <- NULL
if (num_cores > 1) {
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    objects_to_export <- setdiff(ls(envir = .GlobalEnv), "cl")
    parallel::clusterExport(cl, objects_to_export, envir = .GlobalEnv)
} else {
    foreach::registerDoSEQ()
}

simulation_output <- tryCatch(run_simulation(), finally = if (!is.null(cl)) parallel::stopCluster(cl))
