packages <- c("MASS", "stats", "base", "Matrix", "matrixcalc", "rstudioapi",
              "vroom", "rvest", "tibble", "xml2")

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

set.seed(123)

kap                  <- 10
d_alpha              <- 0.5
c_armijo             <- 1e-4
step_init            <- 5
step_min             <- 1e-14
gamma_init_scale     <- 0.5
gamma_init_offdiag   <- 0.05
gamma_init_structure <- "dense"
GAMMA_INIT_GRID      <- c(0.5)
GAMMA_INIT_OFFDIAG_GRID <- c(0.05)
alpha_level          <- 0.05
zcrit                <- qnorm(1 - alpha_level / 2)

Niter       <- 1500
OMEGA_GRID <- sort(unique(c(
  seq(0.01, 0.09, 0.01),
  seq(0.1, 3.0, 0.1),
  seq(2.5, 2.9, 0.05)
)))

commutation_matrix <- function(n) {
  K <- matrix(0, n * n, n * n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      K[(j - 1) * n + i, (i - 1) * n + j] <- 1
    }
  }
  K
}
Sym <- function(A) (A + t(A)) / 2
safe_log_determinant <- function(A, jitter = 0) {
  n <- nrow(A)
  as.numeric(determinant(A + diag(jitter, n), logarithm = TRUE)$modulus)
}

safe_inverse <- function(A, jitter = 0) {
  n <- nrow(A)
  solve(A + diag(jitter, n))
}
vech <- function(A) A[lower.tri(A, diag = TRUE)]
vech_to_mat <- function(v) {
  n <- (-1 + sqrt(1 + 8 * length(v))) / 2
  mat <- matrix(0, nrow = n, ncol = n)
  mat[lower.tri(mat, diag = TRUE)] <- v
  mat
}
make_gamma_init <- function(n) {
  g <- matrix(0, nrow = n, ncol = n)
  g[lower.tri(g)] <- gamma_init_offdiag
  diag(g) <- gamma_init_scale
  g
}
diag_vech_indices <- function(n) {
  j <- seq_len(n)
  as.integer(1 + (j - 1) * (n + 1) - (j - 1) * j / 2)
}
omega_for_study <- function(omega_value, kk) { omega_value }
objective_study_weight <- function(kk) { 1 }

log_T_func <- function(mu_vec, vech_gam, M, c_vec, W, theta){
  gam <- vech_to_mat(vech_gam); G <- gam %*% t(gam); log_det_G <- safe_log_determinant(G); val <- numeric(dim(M)[3])
  for (kk in seq_len(dim(M)[3])) { A <- M[,,kk] + G; v <- c_vec[,kk] + G %*% mu_vec; val[kk] <- log_det_G - safe_log_determinant(A) - t(mu_vec) %*% G %*% mu_vec + t(v) %*% safe_inverse(A) %*% v - t(theta[[kk]]) %*% W[[kk]] %*% theta[[kk]] }
  0.5 * sum(val)
}

lag_func <- function(mu_vec, vech_gam, M, c_vec, W, theta, rho_vec, kap){
  gam <- vech_to_mat(vech_gam); pen <- 0; for (pp in seq_along(rho_vec)) { pen <- pen + rho_vec[pp] * (1/kap) * log(1 + exp(kap * (-gam[pp, pp]))) }
  -log_T_func(mu_vec, vech_gam, M, c_vec, W, theta) + pen
}

vec_gradient <- function(vech_gam, mu_vec, M, c_vec, rho_vec, kap, p, k) {
  gam <- vech_to_mat(vech_gam); n <- p + 1; G <- gam %*% t(gam); Ginv <- safe_inverse(G)
  I <- diag(n); I2 <- diag(n^2); Kc <- commutation_matrix(n)
  make_P_lt <- function(n){ m <- n*(n+1)/2; P <- matrix(0, n*n, m); idx <- 1; for (j in 1:n) { for (i in j:n) { e <- matrix(0, n, n); e[i, j] <- 1; P[, idx] <- as.vector(e); idx <- idx + 1 } }; P }
  P_lt <- make_P_lt(n)
  tmp_sum <- numeric(n^2)
  for (kk in 1:k) { A <- M[,,kk] + G; Ainv <- safe_inverse(A); v <- c_vec[,kk] + G %*% mu_vec; Av <- Ainv %*% v; tmp_sum <- tmp_sum + as.vector(Ginv) - as.vector(Ainv) - as.vector(kronecker(mu_vec, mu_vec)) + as.vector(kronecker(mu_vec, Av)) + as.vector(kronecker(Av, mu_vec)) - as.vector(kronecker(Av, Av)) }
  S <- 0.5 * ((I2 + Kc) %*% (kronecker(gam, I)) %*% P_lt)
  d_log_T_d_vech_gam <- t(S) %*% tmp_sum
  d_penalty <- numeric(length(vech_gam)); diag_idx <- diag_vech_indices(n)
  for (pp in 1:n) { idx <- diag_idx[pp]; d_penalty[idx] <- rho_vec[pp] * ( plogis(kap * gam[pp, pp]) - 1 ) }
  as.vector(- d_log_T_d_vech_gam + d_penalty)
}

optimize_rho <- function(rho_vec, mu_vec, vech_gam, M, c_vec, W, theta, kap, p,
                         eta = 0.1, lower = 1, upper = 3) {
  gam <- vech_to_mat(vech_gam)
  n <- p + 1
  violation <- (1 / kap) * log(1 + exp(kap * (-diag(gam))))
  out <- tryCatch(
    optim(
      par = rho_vec,
      fn = function(rho) {
        -sum(rho * violation) + (1 / (2 * eta)) * sum((rho - rho_vec)^2)
      },
      method = "L-BFGS-B",
      lower = rep(lower, n),
      upper = rep(upper, n)
    ),
    error = function(e) {
      list(par = rho_vec)
    }
  )
  out$par
}


file_path <- file.path(code.dir, "DECEASED_DONOR_DATA.DAT")
htm_path <- file.path(code.dir, "DECEASED_DONOR_DATA.htm")

dat <- vroom(
  file_path,
  delim = "\t",
  col_names = FALSE,
  na = ".",
  trim_ws = TRUE,
  col_types = vroom::cols(.default = vroom::col_character()),
  show_col_types = FALSE
)
dat <- as.data.frame(dat)

htm_data <- read_html(htm_path) %>% html_table()
layout_info <- htm_data[[1]]
new_colnames <- layout_info$LABEL
colnames(dat) <- new_colnames


to_numeric_na <- function(x) suppressWarnings(as.numeric(x))

recode_yn01 <- function(x) {
  x <- trimws(as.character(x)); out <- rep(NA_real_, length(x))
  out[x == "Y"] <- 1; out[x == "N"] <- 0; out
}

recode_diabetes01 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x)); out <- rep(NA_real_, length(x_num))
  out[x_num == 1] <- 0; out[x_num %in% c(2, 3, 4, 5)] <- 1; out[x_num == 998] <- NA; out
}

recode_region <- function(x) factor(suppressWarnings(as.numeric(x)))

df_raw <- data.frame(
  creat        = to_numeric_na(dat$CREAT_DON),
  hypertension = recode_yn01(dat$HIST_HYPERTENS_DON),
  age_raw      = to_numeric_na(dat$AGE_DON),
  bmi_raw      = to_numeric_na(dat$BMI_DON_CALC),
  diabetes     = recode_diabetes01(dat$HIST_DIABETES_DON),
  bun_raw      = to_numeric_na(dat$BUN_DON),
  cig          = recode_yn01(dat$HIST_CIG_DON),
  region       = recode_region(dat$REGION)
)

model_vars_raw <- c("creat", "hypertension", "age_raw", "bmi_raw",
                    "diabetes", "bun_raw", "cig", "region")

missing_table <- data.frame(
  Variable = model_vars_raw,
  N_missing = sapply(model_vars_raw, function(v) sum(is.na(df_raw[[v]]))),
  Pct_missing = sapply(model_vars_raw, function(v) mean(is.na(df_raw[[v]])) * 100)
)

df <- df_raw[complete.cases(df_raw[, model_vars_raw]), ]
df$region <- droplevels(df$region)

df <- df[, c("creat", "hypertension", "age_raw", "bmi_raw",
             "diabetes", "bun_raw", "cig", "region")]
names(df) <- c("creat", "hypertension", "age", "bmi",
               "diabetes", "bun", "cig", "region")

region_levels <- levels(df$region)
k <- length(region_levels)

raw_region_n <- table(df_raw$region, useNA = "ifany")
cc_region_n  <- table(df$region)

region_complete_table <- data.frame(
  region = region_levels,
  n_raw = as.integer(raw_region_n[region_levels]),
  n_complete = as.integer(cc_region_n[region_levels])
)

study_size_by_region <- c(
  "1" = 400, "2" = 400, "3" = 400, "4" = 500, "5" = 500, "6" = 500,
  "7" = 600, "8" = 600, "9" = 600, "10" = 700, "11" = 700
)

region_complete_table$n_dropped <- region_complete_table$n_raw - region_complete_table$n_complete
region_complete_table$pct_complete <- round(100 * region_complete_table$n_complete / region_complete_table$n_raw, 2)
region_complete_table$n_study <- as.integer(study_size_by_region[as.character(region_complete_table$region)])
region_complete_table$n_external <- region_complete_table$n_complete - region_complete_table$n_study

X_full <- model.matrix(~ hypertension + age + bmi + diabetes + bun + cig, data = df)
y <- df$creat
region <- df$region

x_names <- colnames(X_full)[-1]
p <- ncol(X_full) - 1
n <- p + 1

concept_to_cols <- list(
  grep("^hypertension$", x_names), grep("^age$", x_names), grep("^bmi$", x_names),
  grep("^diabetes$", x_names), grep("^bun$", x_names), grep("^cig$", x_names)
)

expand_pattern <- function(concept_idx) sort(unique(unlist(concept_to_cols[concept_idx])))

X_corr <- X_full[, -1, drop = FALSE]
cor_X <- cor(X_corr)

build_real_objects <- function(seed_i = 123, x_formula, A_list_use, x_names_use, study_size_by_region) {
  set.seed(seed_i)
  k_use <- length(region_levels); p_use <- length(x_names_use); n_use <- p_use + 1
  theta <- vector("list", k_use); S_inv <- vector("list", k_use)
  x_full <- vector("list", k_use); x_red <- vector("list", k_use)
  x_full_ext <- vector("list", k_use); x_red_ext <- vector("list", k_use)
  y_study <- vector("list", k_use); study_data <- vector("list", k_use); external_data <- vector("list", k_use)
  num_study <- integer(k_use); num_ext <- integer(k_use)
  
  for (kk in 1:k_use) {
    dat_k <- df[df$region == region_levels[kk], , drop = FALSE]
    study_size_k <- as.integer(study_size_by_region[as.character(region_levels[kk])])
    idx <- sample(seq_len(nrow(dat_k)), replace = FALSE)
    stu_idx <- idx[1:study_size_k]; ext_idx <- idx[(study_size_k + 1):nrow(dat_k)]
    stu_dat <- dat_k[stu_idx, , drop = FALSE]; ext_dat <- dat_k[ext_idx, , drop = FALSE]
    study_data[[kk]] <- stu_dat; external_data[[kk]] <- ext_dat
    num_study[kk] <- nrow(stu_dat); num_ext[kk] <- nrow(ext_dat)
    
    X_stu_full <- model.matrix(x_formula, data = stu_dat)
    X_ext_full <- model.matrix(x_formula, data = ext_dat)
    x_full[[kk]] <- X_stu_full
    x_red[[kk]]  <- X_stu_full[, c(1, A_list_use[[kk]] + 1), drop = FALSE]
    y_study[[kk]] <- stu_dat$creat
    
    theta[[kk]] <- safe_inverse(t(x_red[[kk]]) %*% x_red[[kk]]) %*% t(x_red[[kk]]) %*% y_study[[kk]]
    res <- as.numeric(y_study[[kk]] - x_red[[kk]] %*% theta[[kk]])
    dk <- ncol(x_red[[kk]])
    sigma_hat_k_2 <- sum(res^2) / (num_study[kk] - dk)

    S_inv[[kk]] <- Sym((1 / num_study[kk]) * (1 / sigma_hat_k_2) * crossprod(x_red[[kk]]))
    
    x_full_ext[[kk]] <- X_ext_full
    x_red_ext[[kk]]  <- X_ext_full[, c(1, A_list_use[[kk]] + 1), drop = FALSE]
  }
  
  list(k_use = k_use, p_use = p_use, n_use = n_use, theta = theta, S_inv = S_inv,
       x_full = x_full, x_red = x_red, x_full_ext = x_full_ext, x_red_ext = x_red_ext,
       y_study = y_study, study_data = study_data, external_data = external_data,
       num_study = num_study, num_ext = num_ext, A_list_use = A_list_use, x_names_use = x_names_use)
}

build_M_c_W_real <- function(obj, omega_value) {
  k_use <- obj$k_use; n_use <- obj$n_use
  M <- array(0, c(n_use, n_use, k_use)); c_vec <- array(0, c(n_use, k_use)); W <- vector("list", k_use)
  
  for (kk in 1:k_use) {
    omega_k <- omega_for_study(omega_value, kk)
    m_k <- obj$num_ext[kk]
    S_inv_adj <- obj$num_study[kk] * obj$S_inv[[kk]]
    B_k <- crossprod(obj$x_full_ext[[kk]], obj$x_red_ext[[kk]])
    C_k <- crossprod(obj$x_red_ext[[kk]])
    
    M[,,kk] <- Sym((8 * omega_k) / (m_k^2) * B_k %*% S_inv_adj %*% t(B_k))
    
    c_vec[,kk] <- as.vector((8 * omega_k) / (m_k^2) *
      B_k %*% S_inv_adj %*% C_k %*% obj$theta[[kk]])
    
    W[[kk]] <- Sym((8 * omega_k) / (m_k^2) * C_k %*% S_inv_adj %*% C_k)
  }
  list(M = M, c_vec = c_vec, W = W)
}

compute_sandwich_predvar <- function(M, c_vec, W, theta, hat_mu, hat_gamma, rho_hat, p_use, k_use) {
  n <- p_use + 1; m <- n * (n + 1) / 2
  Gam <- hat_gamma; G <- Gam %*% t(Gam)
  I_n <- diag(n); I2 <- diag(n^2); Kc <- commutation_matrix(n); Ssym <- 0.5 * (I2 + Kc)
  make_P_lt <- function(n){ m <- n*(n+1)/2; P <- matrix(0, n*n, m); idx <- 1; for (j in 1:n) for (i in j:n) { e <- matrix(0,n,n); e[i,j] <- 1; P[,idx] <- as.vector(e); idx <- idx+1 }; P }
  P_lt <- make_P_lt(n); S_Gamma_fix <- (I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt
  
  Vk_list <- vector("list", k_use); mk_list <- vector("list", k_use)
  for (kk in 1:k_use) { Ak <- M[,,kk] + G; Vk <- safe_inverse(Ak); mk <- Vk %*% (c_vec[,kk] + G %*% hat_mu); Vk_list[[kk]] <- Vk; mk_list[[kk]] <- mk }
  
  J <- matrix(0, n + m, n + m)
  for (kk in 1:k_use) {
    Vk <- Vk_list[[kk]]; vk <- c_vec[,kk] + G %*% hat_mu; ak <- Vk %*% vk; rk <- hat_mu - ak
    J_mm_k <- G %*% Vk %*% G - G
    J_mG_vec_k <- (t(ak - hat_mu) %x% diag(n)) + (t(hat_mu - ak) %x% (G %*% Vk)); J_mg_k <- J_mG_vec_k %*% S_Gamma_fix
    M_dr_mu <- (diag(n) - Vk %*% G); d_rrT_mu <- (kronecker(I_n, rk) + kronecker(rk, I_n)) %*% M_dr_mu; J_gm_k <- t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% d_rrT_mu)
    Ginv <- safe_inverse(G)
    tmp_k <- as.vector(Ginv) - as.vector(Vk) - as.vector(hat_mu %*% t(hat_mu)) + as.vector(hat_mu %*% t(ak)) + as.vector(ak %*% t(hat_mu)) - as.vector(ak %*% t(ak))
    z <- Ssym %*% tmp_k; Zmat <- matrix(z, n, n); TermA <- - t(P_lt) %*% (diag(n) %x% Zmat) %*% P_lt
    M_dr_vec <- -(t(rk) %x% Vk); d_rrT_vec <- (kronecker(I_n, rk) + kronecker(rk, I_n)) %*% M_dr_vec; M_tmp_vec <- (-(t(Ginv) %x% Ginv) + (t(Vk) %x% Vk) - d_rrT_vec)
    TermB <- - t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% (M_tmp_vec %*% S_Gamma_fix)); J_gg_k <- TermA + TermB
    J[1:n, 1:n] <- J[1:n, 1:n] + J_mm_k; J[1:n, (n+1):(n+m)] <- J[1:n, (n+1):(n+m)] + J_mg_k; J[(n+1):(n+m), 1:n] <- J[(n+1):(n+m), 1:n] + J_gm_k; J[(n+1):(n+m), (n+1):(n+m)] <- J[(n+1):(n+m), (n+1):(n+m)] + J_gg_k
  }
  
  H_pen <- matrix(0, m, m); diag_idx <- diag_vech_indices(n)
  for (ii in 1:n) { idx <- diag_idx[ii]; pii <- plogis(kap * Gam[ii, ii]); H_pen[idx, idx] <- rho_hat[ii] * kap * pii * (1 - pii) }
  J[(n+1):(n+m), (n+1):(n+m)] <- J[(n+1):(n+m), (n+1):(n+m)] + H_pen
  
  B <- matrix(0, n + m, n + m)
  for (kk in 1:k_use) {
    Vk <- Vk_list[[kk]]; vk <- c_vec[,kk] + G %*% hat_mu; ak <- Vk %*% vk; mk <- mk_list[[kk]]
    s_mu_k <- G %*% Vk %*% vk - G %*% hat_mu
    Ginv <- safe_inverse(G)
    tmp_k <- as.vector(Ginv) - as.vector(Vk) - as.vector(hat_mu %*% t(hat_mu)) + as.vector(hat_mu %*% t(mk)) + as.vector(mk %*% t(hat_mu)) - as.vector(mk %*% t(mk))
    s_vGam_k <- - as.vector(t(P_lt) %*% t(kronecker(Gam, I_n)) %*% (Ssym %*% tmp_k)); s_k <- c(s_mu_k, s_vGam_k); B <- B + s_k %*% t(s_k)
  }
  
  J_inv <- solve(J); Cov_theta_hat_Gam <- J_inv %*% B %*% t(J_inv)
  
  hat_beta <- array(0, c(n, k_use)); PredVar_arr <- array(0, c(n, n, k_use))
  for (kk in 1:k_use) {
    Vk <- Vk_list[[kk]]; mk <- mk_list[[kk]]; hat_beta[, kk] <- mk
    J_mu_k <- Vk %*% G; b <- as.numeric(hat_mu - mk); J_vGam_k <- (t(b) %x% Vk) %*% ((I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt); J_theta_k <- cbind(J_mu_k, J_vGam_k)
    PredVar_k <- Vk + J_theta_k %*% Cov_theta_hat_Gam %*% t(J_theta_k); PredVar_arr[,,kk] <- PredVar_k
  }
  list(Cov_theta_hat_Gam = Cov_theta_hat_Gam, hat_beta = hat_beta, PredVar_arr = PredVar_arr)
}

fit_real_from_objects <- function(obj, omega_value) {
  k_use <- obj$k_use; p_use <- obj$p_use; n_use <- obj$n_use; theta <- obj$theta
  MCW <- build_M_c_W_real(obj, omega_value); M <- MCW$M; c_vec <- MCW$c_vec; W <- MCW$W
  M_w_optimal <- M; c_vec_w_optimal <- c_vec
  
  tilde_mu <- matrix(0, nrow = Niter, ncol = n_use); tilde_gamma <- array(0, c(Niter, n_use, n_use)); rho_mat <- matrix(0, nrow = Niter, ncol = n_use)
  rho_mat[1,] <- rep(1, n_use)
  tilde_mu[1,] <- rep(1, n_use); tilde_mu[1, 1] <- mean(df$creat)
  tilde_gamma[1,,] <- make_gamma_init(n_use)
  
  diag_history <- data.frame(
    iter = integer(0),
    func_old = numeric(0),
    func_new = numeric(0),
    grad_norm = numeric(0),
    step = numeric(0),
    max_diff = numeric(0),
    line_search_failed = logical(0)
  )
  
  j_final <- Niter; converged <- FALSE
  for (j in 1:(Niter - 1)) {
    gam <- tilde_gamma[j,,]; G <- gam %*% t(gam)
    temp_mu_1 <- matrix(0, n_use, n_use); temp_mu_2 <- rep(0, n_use)
    for (kk in 1:k_use) { A <- M[,,kk] + G; Ainv <- safe_inverse(A); temp_mu_1 <- temp_mu_1 + (G - G %*% Ainv %*% G); temp_mu_2 <- temp_mu_2 + (G %*% Ainv %*% c_vec[,kk]) }
    tilde_mu[j+1, ] <- as.vector(safe_inverse(temp_mu_1) %*% temp_mu_2)
    
    z_old <- vech(tilde_gamma[j,,])
    grad_val <- vec_gradient(z_old, tilde_mu[j+1,], M, c_vec, rho_mat[j,], kap, p_use, k_use)
    grad_norm <- sqrt(sum(grad_val^2))
    t_k <- step_init
    func_old <- lag_func(tilde_mu[j+1,], z_old, M, c_vec, W, theta, rho_mat[j,], kap)
    
    if (!is.finite(grad_norm) || !is.finite(func_old)) {
      warning(sprintf("Non-finite gradient or objective at iteration %d.", j))
      j_final <- j
      converged <- FALSE
      break
    }
    
    line_search_failed_iter <- FALSE
    func_new <- NA_real_
    
    repeat {
      z_try <- z_old - t_k * grad_val
      func_new <- lag_func(tilde_mu[j+1,], z_try, M, c_vec, W, theta, rho_mat[j,], kap)
      
      if (!is.finite(func_new)) {
        t_k <- d_alpha * t_k
        if (t_k < step_min) {
          line_search_failed_iter <- TRUE
          break
        }
        next
      }
      
      if (func_new <= (func_old - c_armijo * t_k * sum(grad_val^2))) break
      
      t_k <- d_alpha * t_k
      if (t_k < step_min) {
        line_search_failed_iter <- TRUE
        break
      }
    }
    
    if (line_search_failed_iter) {
      warning(sprintf(
        "Gamma line search failed at iteration %d: func_old = %.6f, grad_norm = %.3e.",
        j, func_old, grad_norm
      ))
      tilde_gamma[j+1,,] <- tilde_gamma[j,,]
      z_cand <- z_old
    } else {
      z_cand <- z_old - t_k * grad_val
      tilde_gamma[j+1,,] <- vech_to_mat(z_cand)
    }
    
    rho_mat[j+1,] <- optimize_rho(rho_mat[j,], tilde_mu[j+1,], z_cand, M, c_vec, W, theta, kap, p_use)
    
    delta_mu <- tilde_mu[j+1,-1] - tilde_mu[j,-1]
    delta_gamma <- as.vector(tilde_gamma[j+1,-1,-1] - tilde_gamma[j,-1,-1])
    max_diff <- max(max(abs(delta_mu)), max(abs(delta_gamma)))
    
    diag_history <- rbind(
      diag_history,
      data.frame(
        iter = j,
        func_old = func_old,
        func_new = func_new,
        grad_norm = grad_norm,
        step = t_k,
        max_diff = max_diff,
        line_search_failed = line_search_failed_iter
      )
    )
    
    if (line_search_failed_iter) {
      j_final <- j + 1
      converged <- FALSE
      break
    }
    
    conv_ok <- max_diff < 1e-5
    
    if (conv_ok) {
      j_final <- j + 1
      converged <- TRUE
      break
    }
  }
  
  hat_mu <- tilde_mu[j_final, ]; hat_gamma <- tilde_gamma[j_final,,]
  rho_hat <- optimize_rho(rho_mat[j_final,], hat_mu, vech(hat_gamma), M, c_vec, W, theta, kap, p_use)
  
  post <- compute_sandwich_predvar(M = M_w_optimal, c_vec = c_vec_w_optimal, W = W, theta = theta,
                                   hat_mu = hat_mu, hat_gamma = hat_gamma, rho_hat = rho_hat, p_use = p_use, k_use = k_use)
  hat_beta <- post$hat_beta; PredVar_arr <- post$PredVar_arr; Cov_theta_hat_Gam <- post$Cov_theta_hat_Gam
  m <- n_use * (n_use + 1) / 2
  
  rownames(hat_beta) <- c("(Intercept)", obj$x_names_use); colnames(hat_beta) <- region_levels
  names(hat_mu) <- c("(Intercept)", obj$x_names_use)
  dimnames(PredVar_arr) <- list(c("(Intercept)", obj$x_names_use), c("(Intercept)", obj$x_names_use), region_levels)
  rownames(Cov_theta_hat_Gam) <- c(names(hat_mu), paste0("gamma_", seq_len(m))); colnames(Cov_theta_hat_Gam) <- c(names(hat_mu), paste0("gamma_", seq_len(m)))
  
  list(hat_mu = hat_mu, hat_beta = hat_beta, hat_gamma = hat_gamma,
       hat_sigma = safe_inverse(hat_gamma %*% t(hat_gamma)),
       Cov_theta_hat_Gam = Cov_theta_hat_Gam, PredVar_arr = PredVar_arr, theta = theta,
       S_inv = obj$S_inv, A_list = obj$A_list_use, x_names = obj$x_names_use, regions = region_levels,
       num_study = obj$num_study, num_ext = obj$num_ext, study_data = obj$study_data, external_data = obj$external_data,
       optimal_w = rep(omega_value, k_use), omega_value = omega_value, j_final = j_final,
       gamma_init_scale = gamma_init_scale, gamma_init_offdiag = gamma_init_offdiag,
       converged = converged, diag_history = diag_history,
       tilde_mu = tilde_mu[seq_len(j_final), , drop = FALSE],
       tilde_gamma = tilde_gamma[seq_len(j_final), , , drop = FALSE],
       rho_mat = rho_mat[seq_len(j_final), , drop = FALSE])
}

run_real <- function(seed_i = 123, x_formula, A_list_use, x_names_use, omega_value = 1, study_size_by_region) {
  obj <- build_real_objects(seed_i = seed_i, x_formula = x_formula, A_list_use = A_list_use, x_names_use = x_names_use, study_size_by_region = study_size_by_region)
  fit_real_from_objects(obj = obj, omega_value = omega_value)
}

patterns_concept_cycle <- list(
  c(1), c(1), c(2, 3, 4), c(2, 3, 4), c(1, 5, 6), c(1, 5, 6),
  c(1, 2), c(1, 2), c(3, 4, 5, 6), c(3, 4, 5, 6), c(3, 4, 5, 6)
)
patterns_concept_meta <- patterns_concept_cycle[((seq_len(k) - 1) %% length(patterns_concept_cycle)) + 1]

study_design_concept_meta <- setNames(patterns_concept_meta, region_levels)

covered_concepts <- sort(unique(unlist(study_design_concept_meta)))

concept_inclusion_count <- sapply(seq_len(p), function(jj) sum(sapply(study_design_concept_meta, function(a) jj %in% a)))
names(concept_inclusion_count) <- paste0("X", seq_len(p))

A_list_meta <- vector("list", k)
for (kk in seq_len(k)) A_list_meta[[kk]] <- expand_pattern(study_design_concept_meta[[region_levels[kk]]])

A_table_meta <- data.frame(
  region = region_levels,
  n_complete = as.integer(cc_region_n[region_levels]),
  conceptual_indices = sapply(region_levels, function(rr) paste(study_design_concept_meta[[rr]], collapse = ", ")),
  actual_indices = sapply(A_list_meta, function(a) paste(a, collapse = ", ")),
  actual_columns = sapply(A_list_meta, function(a) paste(x_names[a], collapse = ", "))
)

predvar_objective_given_fit_real <- function(obj, hat_mu, hat_gamma, omega_value, target_name = NULL) {
  MCW <- build_M_c_W_real(obj = obj, omega_value = omega_value)
  M <- MCW$M; c_vec <- MCW$c_vec; W <- MCW$W; theta <- obj$theta
  p_use <- obj$p_use; k_use <- obj$k_use; n_use <- obj$n_use
  
  rho_hat_local <- optimize_rho(rep(1, n_use), hat_mu, vech(hat_gamma), M, c_vec, W, theta, kap, p_use)
  post <- compute_sandwich_predvar(M = M, c_vec = c_vec, W = W, theta = theta,
                                   hat_mu = hat_mu, hat_gamma = hat_gamma, rho_hat = rho_hat_local, p_use = p_use, k_use = k_use)
  Cov_theta_hat_Gam <- post$Cov_theta_hat_Gam
  
  Gam <- hat_gamma; G <- Gam %*% t(Gam)
  I_n <- diag(n_use); I2 <- diag(n_use^2); Kc <- commutation_matrix(n_use)
  make_P_lt <- function(n){ m <- n*(n+1)/2; P <- matrix(0, n*n, m); idx <- 1; for (j in 1:n) for (i in j:n) { e <- matrix(0,n,n); e[i,j] <- 1; P[,idx] <- as.vector(e); idx <- idx+1 }; P }
  P_lt <- make_P_lt(n_use)
  
  all_beta_names <- c("(Intercept)", obj$x_names_use)
  if (!is.null(target_name)) {
    idx_target <- match(target_name, all_beta_names)
  }
  
  total_plugin_var <- 0
  for (kk in seq_len(k_use)) {
    Vk <- safe_inverse(M[,,kk] + G); mk <- Vk %*% (c_vec[, kk] + G %*% hat_mu)
    J_mu_k <- Vk %*% G; b <- as.numeric(hat_mu - mk); J_vGam_k <- (t(b) %x% Vk) %*% ((I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt); J_theta_k <- cbind(J_mu_k, J_vGam_k)
    PluginVar_k <- J_theta_k %*% Cov_theta_hat_Gam %*% t(J_theta_k)
    sw <- objective_study_weight(kk)
    if (is.null(target_name)) total_plugin_var <- total_plugin_var + sw * sum(diag(PluginVar_k))
    else total_plugin_var <- total_plugin_var + sw * PluginVar_k[idx_target, idx_target]
  }
  total_plugin_var <- as.numeric(total_plugin_var)
  if (!is.finite(total_plugin_var)) total_plugin_var <- .Machine$double.xmax
  total_plugin_var
}

total_pred_objective_given_fit_real <- function(fit, target_name = NULL) {
  all_beta_names <- rownames(fit$hat_beta)
  if (is.null(target_name)) {
    return(sum(sapply(seq_along(fit$regions), function(kk) sum(diag(fit$PredVar_arr[,,kk])))))
  }
  idx_target <- match(target_name, all_beta_names)
  sum(fit$PredVar_arr[idx_target, idx_target, ])
}

target_pred_se_range_from_fit_real <- function(fit, target_name = "hypertension") {
  idx_target <- which(rownames(fit$hat_beta) == target_name)
  if (length(idx_target) != 1) return(c(NA_real_, NA_real_))
  se <- sqrt(pmax(fit$PredVar_arr[idx_target, idx_target, ], 0))
  c(min(se, na.rm = TRUE), max(se, na.rm = TRUE))
}

run_real_two_stage <- function(seed_i = 123, x_formula, A_list_use, x_names_use,
                               grid = OMEGA_GRID, study_size_by_region, target_name = NULL,
                               gamma_init_grid = GAMMA_INIT_GRID,
                               gamma_init_offdiag_grid = GAMMA_INIT_OFFDIAG_GRID) {
  obj <- build_real_objects(seed_i = seed_i, x_formula = x_formula, A_list_use = A_list_use, x_names_use = x_names_use, study_size_by_region = study_size_by_region)
  
  old_gamma_init_scale <- gamma_init_scale
  old_gamma_init_offdiag <- gamma_init_offdiag
  on.exit(gamma_init_scale <<- old_gamma_init_scale, add = TRUE)
  on.exit(gamma_init_offdiag <<- old_gamma_init_offdiag, add = TRUE)
  
  gamma_init_scale <<- old_gamma_init_scale
  gamma_init_offdiag <<- old_gamma_init_offdiag
  fit_stage1 <- fit_real_from_objects(obj = obj, omega_value = 1)
  
  best_obj <- Inf; best_omega <- 1; best_fit <- NULL
  omega_history <- data.frame()
  for (w in grid) {
    for (ginit in gamma_init_grid) {
      for (ginit_offdiag in gamma_init_offdiag_grid) {
        gamma_init_scale <<- ginit
        gamma_init_offdiag <<- ginit_offdiag
        fit_w <- tryCatch(fit_real_from_objects(obj = obj, omega_value = w), error = function(e) NULL)
        obj_val <- Inf
        total_pred_obj <- Inf
        se_rng <- c(NA_real_, NA_real_)
        converged_w <- FALSE
        j_final_w <- NA_integer_
        gamma_move_fro <- NA_real_
        gamma_move_max <- NA_real_
        
        if (!is.null(fit_w)) {
          converged_w <- isTRUE(fit_w$converged)
          j_final_w <- fit_w$j_final
          gamma_start <- fit_w$tilde_gamma[1,,]
          gamma_move_fro <- sqrt(sum((fit_w$hat_gamma - gamma_start)^2))
          gamma_move_max <- max(abs(fit_w$hat_gamma - gamma_start))
          obj_val <- tryCatch(
            predvar_objective_given_fit_real(obj, fit_w$hat_mu, fit_w$hat_gamma, w, target_name = target_name),
            error = function(e) Inf
          )
          total_pred_obj <- tryCatch(total_pred_objective_given_fit_real(fit_w, target_name = target_name), error = function(e) Inf)
          se_rng <- target_pred_se_range_from_fit_real(fit_w, target_name = "hypertension")
        }
        
        omega_history <- rbind(
          omega_history,
          data.frame(
            omega_init = 1,
            omega = w,
            gamma_init_scale = ginit,
            gamma_init_offdiag = ginit_offdiag,
            gamma_move_fro = gamma_move_fro,
            gamma_move_max = gamma_move_max,
            plugin_objective = obj_val,
            total_pred_objective = total_pred_obj,
            converged = converged_w,
            j_final = j_final_w,
            hypertension_pred_se_min = se_rng[1],
            hypertension_pred_se_max = se_rng[2]
          )
        )
        
        if (converged_w && is.finite(obj_val) && obj_val < best_obj) {
          best_obj <- obj_val
          best_omega <- w
          best_fit <- fit_w
        }
      }
    }
  }
  
  if (is.null(best_fit)) {
    warning("No converged omega-grid fit had a finite plugin objective; using omega = 1.")
    best_omega <- 1
    best_obj <- predvar_objective_given_fit_real(obj, fit_stage1$hat_mu, fit_stage1$hat_gamma, 1, target_name = target_name)
    best_fit <- fit_stage1
  }
  
  list(fit_stage1 = fit_stage1, omega_opt = best_omega, objective = best_obj,
       fit_stage2 = best_fit, omega_history = omega_history, target_name = target_name)
}

naive_formula_lm <- creat ~ hypertension
full_formula_lm  <- creat ~ hypertension + age + bmi + diabetes + bun + cig
x_formula_meta   <- ~ hypertension + age + bmi + diabetes + bun + cig
x_names_meta     <- colnames(model.matrix(x_formula_meta, data = df))[-1]

fit_naive_lm <- lm(naive_formula_lm, data = df)
fit_full_lm  <- lm(full_formula_lm,  data = df)

res_meta <- run_real_two_stage(
  seed_i = 123, x_formula = x_formula_meta, A_list_use = A_list_meta, x_names_use = x_names_meta,
  grid = OMEGA_GRID, study_size_by_region = study_size_by_region, target_name = NULL
)
fit_meta <- res_meta$fit_stage2

omega_summary <- data.frame(Model = "Meta", omega_init = 1, omega_opt = res_meta$omega_opt,
                            gamma_init_scale = fit_meta$gamma_init_scale,
                            gamma_init_offdiag = fit_meta$gamma_init_offdiag,
                            objective = res_meta$objective)

region_n <- table(df$region)
split_summary <- data.frame(
  region = fit_meta$regions,
  n_complete = as.integer(region_n[fit_meta$regions]),
  n_study_used = fit_meta$num_study,
  n_external_used = fit_meta$num_ext
)

zcrit_compare <- qnorm(1 - 0.025 / 3)

get_mu_hypertension_ci_meta <- function(fit, label) {
  idx <- match("hypertension", names(fit$hat_mu))
  est <- fit$hat_mu[idx]; se <- sqrt(pmax(fit$Cov_theta_hat_Gam[idx, idx], 0))
  data.frame(Model = label, Estimate = as.numeric(est),
             Lower_CI = as.numeric(est - zcrit_compare * se), Upper_CI = as.numeric(est + zcrit_compare * se),
             SE = as.numeric(se), N_used = NA_integer_, tau2 = NA_real_)
}

get_hypertension_ci_lm <- function(fit, label, N_used = NA_integer_) {
  cf <- coef(fit); vc <- vcov(fit); idx <- match("hypertension", names(cf))
  est <- cf[idx]; se <- sqrt(vc[idx, idx])
  data.frame(Model = label, Estimate = as.numeric(est),
             Lower_CI = as.numeric(est - zcrit_compare * se), Upper_CI = as.numeric(est + zcrit_compare * se),
             SE = as.numeric(se), N_used = as.integer(N_used))
}

get_region_hypertension_ci_lm <- function(fit, label, region_label, N_used) {
  cf <- coef(fit); vc <- vcov(fit); idx <- which(names(cf) == "hypertension")
  if (length(idx) != 1) {
    return(data.frame(Model = label, region = region_label, Estimate = NA_real_, Lower_CI = NA_real_,
                      Upper_CI = NA_real_, SE = NA_real_, N_used = as.integer(N_used), Var_type = "lm_wald"))
  }
  est <- cf[idx]; se <- sqrt(vc[idx, idx])
  data.frame(Model = label, region = region_label, Estimate = as.numeric(est),
             Lower_CI = as.numeric(est - zcrit_compare * se), Upper_CI = as.numeric(est + zcrit_compare * se),
             SE = as.numeric(se), N_used = as.integer(N_used), Var_type = "lm_wald")
}

meta_pool_hypertension <- function(formula, label) {
  betas <- numeric(0)
  ses <- numeric(0)
  for (rk in region_levels) {
    dk <- df[df$region == rk, , drop = FALSE]
    fk <- lm(formula, data = dk)
    cf <- coef(fk)
    vc <- vcov(fk)
    idx <- which(names(cf) == "hypertension")
    if (length(idx) != 1 || !is.finite(cf[idx]) || !is.finite(vc[idx, idx]) || vc[idx, idx] <= 0) next
    betas <- c(betas, cf[idx])
    ses <- c(ses, sqrt(vc[idx, idx]))
  }
  K <- length(betas)
  v <- ses^2
  wf <- 1 / v
  mu_f <- sum(wf * betas) / sum(wf)
  Q <- sum(wf * (betas - mu_f)^2)
  C <- sum(wf) - sum(wf^2) / sum(wf)
  tau2 <- if (K > 1 && C > 0) max(0, (Q - (K - 1)) / C) else 0
  wr <- 1 / (v + tau2)
  mu_hat <- sum(wr * betas) / sum(wr)
  se_hat <- sqrt(1 / sum(wr))
  data.frame(Model = label, Estimate = as.numeric(mu_hat),
             Lower_CI = as.numeric(mu_hat - zcrit_compare * se_hat),
             Upper_CI = as.numeric(mu_hat + zcrit_compare * se_hat),
             SE = as.numeric(se_hat), N_used = K, tau2 = as.numeric(tau2))
}

mu_comparison <- rbind(
  meta_pool_hypertension(creat ~ hypertension, "Naive"),
  meta_pool_hypertension(creat ~ hypertension + age + bmi + diabetes + bun + cig, "DL"),
  get_mu_hypertension_ci_meta(fit_meta, "Meta")
)

target_name <- "hypertension"
idx_target <- match(target_name, rownames(fit_meta$hat_beta))

obj_dbg <- build_real_objects(seed_i = 123, x_formula = x_formula_meta, A_list_use = A_list_meta, x_names_use = x_names_meta, study_size_by_region = study_size_by_region)
MCW_dbg <- build_M_c_W_real(obj = obj_dbg, omega_value = fit_meta$omega_value)
M_dbg <- MCW_dbg$M; c_dbg <- MCW_dbg$c_vec

Gam <- fit_meta$hat_gamma; G <- Gam %*% t(Gam)
I_n <- diag(n); I2 <- diag(n^2); Kc <- commutation_matrix(n)
make_P_lt <- function(n){ m <- n*(n+1)/2; P <- matrix(0, n*n, m); idx <- 1; for (j in 1:n) for (i in j:n) { e <- matrix(0,n,n); e[i,j] <- 1; P[,idx] <- as.vector(e); idx <- idx+1 }; P }
P_lt <- make_P_lt(n)

beta_all_comparison <- data.frame(); meta_plugin_decomp <- data.frame()
for (kk in seq_along(fit_meta$regions)) {
  region_k <- as.character(fit_meta$regions[kk])
  region_complete_data_k <- df[df$region == region_k, , drop = FALSE]; n_region_k <- nrow(region_complete_data_k)
  
  fit_naive_region_k <- lm(creat ~ hypertension, data = region_complete_data_k)
  fit_full_region_k  <- lm(creat ~ hypertension + age + bmi + diabetes + bun + cig, data = region_complete_data_k)
  
  naive_row <- get_region_hypertension_ci_lm(fit_naive_region_k, "Naive", region_k, n_region_k)
  full_row  <- get_region_hypertension_ci_lm(fit_full_region_k, "Full", region_k, n_region_k)
  
  Vk <- safe_inverse(M_dbg[,,kk] + G); mk <- Vk %*% (c_dbg[, kk] + G %*% fit_meta$hat_mu)
  est_meta_saved <- as.numeric(fit_meta$hat_beta[idx_target, kk]); est_meta_recomputed <- as.numeric(mk[idx_target])
  
  J_mu_k <- Vk %*% G; b <- as.numeric(fit_meta$hat_mu - mk)
  J_vGam_k <- (t(b) %x% Vk) %*% ((I2 + Kc) %*% (kronecker(Gam, I_n)) %*% P_lt); J_theta_k <- cbind(J_mu_k, J_vGam_k)
  
  plugin_var <- (J_theta_k %*% fit_meta$Cov_theta_hat_Gam %*% t(J_theta_k))[idx_target, idx_target]
  pred_var   <- Vk[idx_target, idx_target] + plugin_var
  se_meta    <- sqrt(pmax(pred_var, 0))
  
  meta_row <- data.frame(Model = "Meta", region = region_k, Estimate = est_meta_recomputed,
                         Lower_CI = as.numeric(est_meta_recomputed - zcrit_compare * se_meta),
                         Upper_CI = as.numeric(est_meta_recomputed + zcrit_compare * se_meta),
                         SE = as.numeric(se_meta), N_used = NA_integer_, Var_type = "prediction")
  
  beta_all_comparison <- rbind(beta_all_comparison, naive_row, full_row, meta_row)
  meta_plugin_decomp <- rbind(meta_plugin_decomp,
                              data.frame(region = region_k, beta_hat = est_meta_recomputed,
                                         Vk_var = as.numeric(Vk[idx_target, idx_target]),
                                         plugin_var = as.numeric(plugin_var),
                                         pred_var = as.numeric(pred_var),
                                         pred_se = as.numeric(se_meta)))
}

save(
  missing_table,
  region_complete_table,
  cor_X,
  A_table_meta,
  omega_summary,
  split_summary,
  mu_comparison,
  beta_all_comparison,
  meta_plugin_decomp,
  file = file.path(code.dir, "ci-ma-re-real-data-results.Rdata")
)
