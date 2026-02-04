#' Asymptotic critical values for multivariate normal mixture MEM test
#'
#' Computes asymptotic critical values and p-values for the modified EM test
#' of Kasahara and Shimotsu for multivariate normal mixtures.
#'
#' @export
#' @title mvnmixCrit
#' @name mvnmixCrit
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#'   in the form of (alpha = (alpha_1,...,alpha_m), mu = (mu_1',...,mu_m'),
#'   sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param values Vector of test statistic values at which p-values are computed
#' @param nrep Number of simulation replications (default 10000)
#' @param ninits.crit Number of random starts for cone projection optimization (default 25)
#' @return A list with:
#' \item{crit}{Critical values at 10\%, 5\%, 1\% levels}
#' \item{pvals}{P-values corresponding to \code{values}}
mvnmixCrit <- function(y, parlist, values = NULL, nrep = 10000, ninits.crit = 25) {

  y <- as.matrix(y)
  n <- nrow(y)
  d <- ncol(y)

  alpha <- parlist$alpha
  m     <- length(alpha)

  pvals <- NULL

  # d=1: the two cones cover R^2 so the limit is chi-sq(2) per component
  if (d == 1) {
    if (m == 1) {
      crit <- qchisq(c(0.1, 0.05, 0.01), 2, lower.tail = FALSE)
      if (!is.null(values)) {
        pvals <- pchisq(values, 2, lower.tail = FALSE)
      }
      return(list(crit = crit, pvals = pvals))
    } else {
      return(mvnmixCrit_simulate(y, parlist, values, nrep, ninits.crit = 1))
    }
  }

  # d >= 2: need cone projection
  return(mvnmixCrit_simulate(y, parlist, values, nrep, ninits.crit))
}


# ============================================================
# Core simulation engine
# ============================================================

mvnmixCrit_simulate <- function(y, parlist, values, nrep, ninits.crit = 25) {

  set.seed(123456)

  n <- nrow(y)
  d <- ncol(y)
  dsig <- d * (d + 1) / 2

  alpha <- parlist$alpha
  mu    <- parlist$mu
  sigma <- parlist$sigma
  m     <- length(alpha)

  # Precompute index tuples and permutation tables (once, not per draw)
  tup3 <- mvn_tuples3(d)
  tup4 <- mvn_tuples4(d)
  d_muv <- nrow(tup3)
  d_mu4 <- nrow(tup4)
  n_lam <- d_muv + d_mu4

  perm12_list <- vector("list", d_muv)
  for (s in 1:d_muv)
    perm12_list[[s]] <- perm12(tup3[s, 1], tup3[s, 2], tup3[s, 3])

  perm22_list <- vector("list", d_mu4)
  mc4_vec     <- numeric(d_mu4)
  for (s in 1:d_mu4) {
    perm22_list[[s]] <- perm22(tup4[s, 1], tup4[s, 2], tup4[s, 3], tup4[s, 4])
    mc4_vec[s] <- perm4_count(tup4[s, 1], tup4[s, 2], tup4[s, 3], tup4[s, 4])
  }

  # Rebuild mu and Sigma matrices
  mu_mat <- matrix(mu, nrow = d, ncol = m)
  Sigma_list <- vector("list", m)
  for (j in 1:m) {
    sv <- sigma[((j - 1) * dsig + 1):(j * dsig)]
    Sigma_list[[j]] <- sigmavec2mat(sv, d)
  }

  # Compute scores
  scores <- mvn_scores(y, alpha, mu_mat, Sigma_list, m, d, tup3, tup4)
  S_eta <- scores$S_eta
  S_lam <- scores$S_lam

  # Information matrices
  I_eta <- crossprod(S_eta) / n
  I_lam <- crossprod(S_lam) / n
  I_el  <- crossprod(S_eta, S_lam) / n

  # Regularize if ill-conditioned
  Iall <- rbind(cbind(I_eta, I_el), cbind(t(I_el), I_lam))
  if (rcond(Iall) < .Machine$double.eps) {
    eig <- eigen(Iall, symmetric = TRUE)
    tol2 <- 1e-14 * eig$values[1]
    vals <- pmax(eig$values, tol2)
    Iall <- eig$vectors %*% (vals * t(eig$vectors))
    dim_eta <- ncol(S_eta)
    I_eta <- Iall[1:dim_eta, 1:dim_eta]
    I_el  <- Iall[1:dim_eta, (dim_eta + 1):ncol(Iall)]
    I_lam <- Iall[(dim_eta + 1):nrow(Iall), (dim_eta + 1):ncol(Iall)]
  }

  # Schur complement
  I_lam_eta <- I_lam - t(I_el) %*% solve(I_eta, I_el)

  # Generate draws from N(0, I_lam.eta)
  e <- eigen(I_lam_eta, symmetric = TRUE)
  e$values <- pmax(e$values, 0)
  sqrtI <- e$vectors %*% (sqrt(e$values) * t(e$vectors))
  u <- matrix(rnorm(nrep * m * n_lam), nrep, m * n_lam) %*% sqrtI

  # For each component, compute LR
  EM <- matrix(0, nrow = nrep, ncol = m)

  for (jj in 1:m) {
    idx <- ((jj - 1) * n_lam + 1):(jj * n_lam)
    I_jj <- I_lam_eta[idx, idx]
    I_jj_inv <- solve(I_jj)
    u_jj <- u[, idx]

    if (d == 1) {
      # d=1: chi-sq(2) per component
      EM[, jj] <- rowSums((u_jj %*% I_jj_inv) * u_jj)
    } else {
      # d >= 2: cone projection needed
      Z_jj <- u_jj %*% I_jj_inv   # nrep x n_lam
      R_chol <- chol(I_jj)

      # Precompute initial value grid (deterministic, done once)
      init_grid <- make_init_grid(d, dsig, ninits.crit)

      for (rr in 1:nrep) {
        Z_r <- Z_jj[rr, ]
        EM[rr, jj] <- cone_project_maxLR(Z_r, R_chol, d, dsig, d_muv, d_mu4,
                                          perm12_list, perm22_list, mc4_vec,
                                          tup4, init_grid)
      }
    }
  }

  max_EM <- apply(EM, 1, max)
  max_EM_sort <- sort(max_EM)

  q <- ceiling(nrep * c(0.90, 0.95, 0.99))
  crit <- max_EM_sort[q]

  pvals <- NULL
  if (!is.null(values)) {
    k <- length(values)
    pvals <- rowMeans(t(matrix(rep.int(max_EM_sort, k), ncol = k)) >= values)
  }

  return(list(crit = crit, pvals = pvals))
}


# ============================================================
# Index tuple generators
# ============================================================

mvn_tuples3 <- function(d) {
  out <- matrix(0L, nrow = 0, ncol = 3)
  for (i in 1:d) for (j in i:d) for (k in j:d)
    out <- rbind(out, c(i, j, k))
  out
}

mvn_tuples4 <- function(d) {
  out <- matrix(0L, nrow = 0, ncol = 4)
  for (i in 1:d) for (j in i:d) for (k in j:d) for (l in k:d)
    out <- rbind(out, c(i, j, k, l))
  out
}


# ============================================================
# Permutation helpers (precomputed)
# ============================================================

perm12 <- function(i, j, k) {
  perms <- rbind(c(i, j, k), c(j, i, k), c(k, i, j))
  perms[, 2:3] <- t(apply(perms[, 2:3, drop = FALSE], 1, sort))
  unique(perms)
}

perm22 <- function(i, j, k, l) {
  idx <- c(i, j, k, l)
  all_perms <- .permutations4(idx)
  out <- matrix(0L, nrow = 0, ncol = 4)
  for (r in 1:nrow(all_perms)) {
    p <- all_perms[r, ]
    a <- sort(p[1:2])
    b <- sort(p[3:4])
    if (a[1] > b[1] || (a[1] == b[1] && a[2] > b[2])) {
      tmp <- a; a <- b; b <- tmp
    }
    out <- rbind(out, c(a, b))
  }
  unique(out)
}

perm4_count <- function(i, j, k, l) {
  tab <- table(c(i, j, k, l))
  factorial(4) / prod(factorial(tab))
}

.permutations4 <- function(x) {
  n <- length(x)
  if (n == 1) return(matrix(x, 1, 1))
  out <- matrix(0L, nrow = 0, ncol = n)
  for (i in 1:n) {
    rest <- .permutations4(x[-i])
    out <- rbind(out, cbind(x[i], rest))
  }
  unique(out)
}


# ============================================================
# Fast cone mapping using precomputed permutation tables
# ============================================================

cone1_map_fast <- function(lam_mu, lam_v_mat, d_muv, d_mu4, perm12_list, perm22_list) {
  t_muv <- numeric(d_muv)
  for (s in 1:d_muv) {
    perms <- perm12_list[[s]]
    val <- 0
    for (r in 1:nrow(perms))
      val <- val + lam_mu[perms[r, 1]] * lam_v_mat[perms[r, 2], perms[r, 3]]
    t_muv[s] <- val
  }

  t_mu4 <- numeric(d_mu4)
  for (s in 1:d_mu4) {
    perms <- perm22_list[[s]]
    val <- 0
    for (r in 1:nrow(perms))
      val <- val + lam_v_mat[perms[r, 1], perms[r, 2]] * lam_v_mat[perms[r, 3], perms[r, 4]]
    t_mu4[s] <- val
  }

  c(t_muv, t_mu4)
}

cone2_map_fast <- function(lam_mu, lam_v_mat, d_muv, d_mu4, perm12_list, mc4_vec, tup4) {
  t_muv <- numeric(d_muv)
  for (s in 1:d_muv) {
    perms <- perm12_list[[s]]
    val <- 0
    for (r in 1:nrow(perms))
      val <- val + lam_mu[perms[r, 1]] * lam_v_mat[perms[r, 2], perms[r, 3]]
    t_muv[s] <- val
  }

  t_mu4 <- numeric(d_mu4)
  for (s in 1:d_mu4) {
    t_mu4[s] <- -mc4_vec[s] * lam_mu[tup4[s, 1]] * lam_mu[tup4[s, 2]] *
                               lam_mu[tup4[s, 3]] * lam_mu[tup4[s, 4]]
  }

  c(t_muv, t_mu4)
}


# ============================================================
# Initial value grid for cone projection
# ============================================================

make_init_grid <- function(d, dsig, ninits) {
  n_par <- d + dsig
  # Systematic grid: different scales and sign patterns
  grid <- matrix(0, nrow = ninits, ncol = n_par)
  for (i in 1:ninits) {
    scale <- c(0.1, 0.5, 1.0, 2.0)[(i %% 4) + 1]
    grid[i, ] <- rnorm(n_par, sd = scale)
  }
  grid
}


# ============================================================
# Cone projection: max of LR over both cones
# ============================================================

cone_project_maxLR <- function(Z, R_chol, d, dsig, d_muv, d_mu4,
                                perm12_list, perm22_list, mc4_vec,
                                tup4, init_grid) {
  n_par <- d + dsig

  # Quadratic form matrices
  W <- crossprod(R_chol)  # = I_jj (the information submatrix)

  # Unconstrained LR = Z' W Z
  LR_full <- sum(Z * (W %*% Z))

  # --- Cone 1 objective ---
  obj1 <- function(par) {
    lam_mu <- par[1:d]
    lam_v_mat <- sigmavec2mat(par[(d + 1):n_par], d)
    t_vec <- cone1_map_fast(lam_mu, lam_v_mat, d_muv, d_mu4, perm12_list, perm22_list)
    diff <- t_vec - Z
    sum(diff * (W %*% diff))
  }

  # --- Cone 2 objective ---
  obj2 <- function(par) {
    lam_mu <- par[1:d]
    lam_v_mat <- sigmavec2mat(par[(d + 1):n_par], d)
    t_vec <- cone2_map_fast(lam_mu, lam_v_mat, d_muv, d_mu4, perm12_list, mc4_vec, tup4)
    diff <- t_vec - Z
    sum(diff * (W %*% diff))
  }

  best1 <- obj1(rep(0, n_par))
  best2 <- obj2(rep(0, n_par))

  ninits <- nrow(init_grid)
  for (i in 1:ninits) {
    par0 <- init_grid[i, ]

    res1 <- tryCatch(
      nlminb(par0, obj1, control = list(iter.max = 100, rel.tol = 1e-8)),
      error = function(e) list(objective = Inf)
    )
    if (res1$objective < best1) best1 <- res1$objective

    res2 <- tryCatch(
      nlminb(par0, obj2, control = list(iter.max = 100, rel.tol = 1e-8)),
      error = function(e) list(objective = Inf)
    )
    if (res2$objective < best2) best2 <- res2$objective
  }

  LR1 <- max(LR_full - best1, 0)
  LR2 <- max(LR_full - best2, 0)
  max(LR1, LR2)
}


# ============================================================
# Score computation
# ============================================================

mvn_scores <- function(y, alpha, mu_mat, Sigma_list, m, d, tup3, tup4) {
  n <- nrow(y)
  dsig <- d * (d + 1) / 2
  d_muv <- nrow(tup3)
  d_mu4 <- nrow(tup4)
  n_lam <- d_muv + d_mu4

  P_list <- vector("list", m)
  f_mat <- matrix(0, n, m)
  z_list <- vector("list", m)

  for (j in 1:m) {
    Sigma_j <- Sigma_list[[j]]
    P_j <- solve(Sigma_j)
    P_list[[j]] <- P_j

    resid <- t(t(y) - mu_mat[, j])
    z_j <- resid %*% P_j
    z_list[[j]] <- z_j

    log_det <- determinant(Sigma_j, logarithm = TRUE)$modulus
    mahal <- rowSums(z_j * resid)
    f_mat[, j] <- exp(-0.5 * (d * log(2 * pi) + log_det + mahal))
  }

  f0 <- f_mat %*% alpha
  w_mat <- t(t(f_mat) * alpha) / as.vector(f0)

  # --- Identifiable scores ---
  if (m >= 2) {
    S_alpha <- (f_mat[, 1:(m - 1), drop = FALSE] - f_mat[, m]) / as.vector(f0)
  } else {
    S_alpha <- matrix(0, n, 0)
  }

  S_mu <- matrix(0, n, m * d)
  for (j in 1:m)
    S_mu[, ((j - 1) * d + 1):(j * d)] <- w_mat[, j] * z_list[[j]]

  S_v <- matrix(0, n, m * dsig)
  for (j in 1:m) {
    P_j <- P_list[[j]]
    z_j <- z_list[[j]]
    sv_j <- matrix(0, n, dsig)
    idx <- 0
    for (b in 1:d) {
      for (a in b:d) {
        idx <- idx + 1
        c_ab <- ifelse(a == b, 1, 2)
        sv_j[, idx] <- w_mat[, j] * 0.5 * c_ab * (z_j[, a] * z_j[, b] - P_j[a, b])
      }
    }
    S_v[, ((j - 1) * dsig + 1):(j * dsig)] <- sv_j
  }

  S_eta <- cbind(S_alpha, S_mu, S_v)

  # --- Unidentifiable scores ---
  S_lam <- matrix(0, n, m * n_lam)

  for (j in 1:m) {
    P_j <- P_list[[j]]
    z_j <- z_list[[j]]

    s_muv <- matrix(0, n, d_muv)
    for (s in 1:d_muv) {
      ii <- tup3[s, 1]; jj <- tup3[s, 2]; kk <- tup3[s, 3]
      s_muv[, s] <- w_mat[, j] * (
        z_j[, ii] * z_j[, jj] * z_j[, kk]
        - P_j[ii, jj] * z_j[, kk]
        - P_j[ii, kk] * z_j[, jj]
        - P_j[jj, kk] * z_j[, ii]
      ) / 6
    }

    s_mu4 <- matrix(0, n, d_mu4)
    for (s in 1:d_mu4) {
      ii <- tup4[s, 1]; jj <- tup4[s, 2]; kk <- tup4[s, 3]; ll <- tup4[s, 4]
      s_mu4[, s] <- w_mat[, j] * (
        z_j[, ii] * z_j[, jj] * z_j[, kk] * z_j[, ll]
        - P_j[ii, jj] * z_j[, kk] * z_j[, ll]
        - P_j[ii, kk] * z_j[, jj] * z_j[, ll]
        - P_j[ii, ll] * z_j[, jj] * z_j[, kk]
        - P_j[jj, kk] * z_j[, ii] * z_j[, ll]
        - P_j[jj, ll] * z_j[, ii] * z_j[, kk]
        - P_j[kk, ll] * z_j[, ii] * z_j[, jj]
        + P_j[ii, jj] * P_j[kk, ll]
        + P_j[ii, kk] * P_j[jj, ll]
        + P_j[ii, ll] * P_j[jj, kk]
      ) / 24
    }

    col_start <- (j - 1) * n_lam + 1
    S_lam[, col_start:(col_start + d_muv - 1)] <- s_muv
    S_lam[, (col_start + d_muv):(col_start + n_lam - 1)] <- s_mu4
  }

  list(S_eta = S_eta, S_lam = S_lam)
}
