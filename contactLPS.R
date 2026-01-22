contactLPS <- function(data, 
                       K.age,
                       K.disp,
                       WAIC = FALSE,
                       L = 100,
                       penorder = 2){
  
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  # Data preprocessing
  data$group <- factor(data$group)
  tpart0 <- which(data$tpart == 0)
  keep_idx <- setdiff(seq_len(nrow(data)), tpart0)
  data_sub <- data[keep_idx, ]
  data_sub$time <- factor(data_sub$time)
  times <- levels(data_sub$time)
  
  min_age_part <- min(data$age_part)
  max_age_part <- max(data$age_part)
  min_age_cont <- min(data$age_cont)
  max_age_cont <- max(data$age_cont)
  
  age_p <- seq(min_age_part, max_age_part)
  age_c <- seq(min_age_cont, max_age_cont)
  
  build_symmetric_basis <- function(Bi, Bj, K) {
    Bgrid <- kronecker(Bj, Bi)
    Bmat <- matrix(1:(K^2), nrow = K, ncol = K)
    ind1 <- Bmat[lower.tri(Bmat)]
    ind2 <- t(Bmat)[lower.tri(t(Bmat))]
    Bgrid[,ind1] <- Bgrid[,ind1] + Bgrid[,ind2]
    Bsym <- Bgrid[,-c(Bmat[upper.tri(Bmat)])]
    Bsym
  }
  ###############################################################################
  ### Design matrices for contact rate (mu)
  Bi   <- cubicbs(age_p, lower = min_age_part, upper = max_age_part, K = K.age)$Bmatrix
  Bj   <- cubicbs(age_c, lower = min_age_cont, upper = max_age_cont, K = K.age)$Bmatrix
  Bsym <- build_symmetric_basis(Bi, Bj, K.age)
  Bfull <- do.call(rbind, replicate(length(times), Bsym, simplify = FALSE))
  B <- Matrix::Matrix(Bfull, sparse = TRUE)[keep_idx, , drop = FALSE]
  
  ### Model matrices
  
  # Design matrix for overdispersion (phi)
  Zi    <- cubicbs(age_p, lower = min_age_part, upper = max_age_part, K = K.disp)$Bmatrix
  Zj    <- cubicbs(age_c, lower = min_age_cont, upper = max_age_cont, K = K.disp)$Bmatrix
  Z1    <- kronecker(Zj, Zi)
  Zfull <- do.call(rbind, replicate(length(times), Z1, simplify = FALSE))
  Z     <- Matrix::Matrix(Zfull, sparse = TRUE)[keep_idx, , drop = FALSE]
  
  # Group covariate matrix
  Wfull <- model.matrix(~ group, data = data)
  W     <- Matrix::Matrix(Wfull, sparse = TRUE)[keep_idx, , drop = FALSE]
  
  W_no_int_full <- Wfull[, -1, drop = FALSE]
  
  # Compute interaction design matrix for all interactions
  BW_list <- lapply(1:ncol(W_no_int_full), function(j) Bfull * W_no_int_full[, j])
  
  # Combine all interaction matrices column-wise
  BW_full <- do.call(cbind, BW_list)
  
  BW <- Matrix::Matrix(BW_full, sparse = TRUE)[keep_idx, , drop = FALSE]
  
  X <- cbind(W, B, BW)
  
  y <- data_sub$y
  offset <- data_sub$offset
  pB <- ncol(B)
  pZ <- ncol(Z)
  pW    <- ncol(W)
  
  ###########################################################################
  
  construct_penalty <- function(K, order = 2, sym = TRUE) {
    
    # build the 'order'-th difference operator D of size (n-order) x n
    D <- diag(K)
    for (k in seq_len(order)) {
      D <- diff(D)
    }
    
    if (!sym) {
      # non‐symmetric penalty: P = D'D
      DTD <- crossprod(D)
      P <- kronecker(diag(K), DTD) + kronecker(DTD, diag(K))
      
    } else {
      # symmetric penalty: restrict (I ⊗ D) and (D ⊗ I) to lower‐triangle entries
     
      # create an index matrix for an n×n theta matrix
      Theta <- matrix(seq_len(K^2), nrow = K, ncol = K)
      colD <- Theta[lower.tri(Theta, diag = TRUE)]
      
      # row‐indices for (I ⊗ D):
      # block‐diagonal D’s; each block has (n-order) rows and n columns
      D1mat   <- matrix(seq_len((K - order) * K),
                        nrow = (K - order), ncol = K)
      rowD1 <- D1mat[lower.tri(D1mat, diag = TRUE)]
      
      D2mat_full <- matrix(seq_len(K * (K - order)),
                           nrow = K, ncol = (K - order))
      D2mat      <- D2mat_full[-seq_len(order), , drop = FALSE]
      rowD2   <- D2mat[lower.tri(D2mat, diag = TRUE)]
      
      # restrict each to the lower‐triangle entries and sum their cross‐products
      DTD1 <- crossprod(kronecker(diag(K), D)[rowD1, colD, drop = FALSE] )
      DTD2 <- crossprod(kronecker(D, diag(K))[rowD2, colD, drop = FALSE])
      
      P <- DTD1 + DTD2
    }
    
    return(P)
  }
  
  
  P_B <- construct_penalty(K.age, penorder, sym = TRUE)
  P_Z <- construct_penalty(K.disp, penorder, sym = FALSE)
  
  QvB <- function(v) (exp(v[1]) * P_B) + diag(1e-5, pB)
  QvBW <- function(v) Matrix::bdiag(replicate(pW - 1, QvB(v), simplify = FALSE))

  Qv.mu <- function(v) Matrix::bdiag(diag(1e-5, pW), QvB(v), QvBW(v))
  Qv.phi <- function(v) Matrix::Matrix((exp(v[2]) * P_Z) + diag(1e-5, pZ), sparse = TRUE)
  Qv.xi <- function(v) Matrix::bdiag(Qv.mu(v), Qv.phi(v))
  
  
  
  # Mean function for NB
  mu.nb <- function(xi_mu) as.numeric(exp(X %*% xi_mu) * offset)
  
  # Log-posterior of xi given v
  log_pxi <- function(xi_mu, xi_phi, v) {
    xi     <- c(xi_mu, xi_phi)
    Qv.xi  <- Qv.xi(v)
    mu     <- mu.nb(xi_mu)
    phi    <- as.numeric(exp(Z %*% xi_phi))
    ll     <- sum(lgamma(y + phi) - lgamma(phi) + (y * (as.numeric(X %*% xi_mu) + log(offset))) + 
                    phi * as.numeric(Z %*% xi_phi) - (y + phi) * log(mu + phi))
    prior  <- - 0.5 * sum((xi * Qv.xi) %*% xi)
    ll + prior
  }  
  
  
  # Gradient with respect to theta
  Grad.theta <- function(xi_mu, xi_phi) {
    phi <- exp(as.numeric(Z %*% xi_phi))
    mu <- mu.nb(xi_mu)
    
    result <- colSums((y - ((y + phi) * mu )/ (mu + phi)) * X)
    return(result)
  }
  
  Hess.theta <- function(xi_mu, xi_phi) {
    phi <- exp(as.numeric(Z %*% xi_phi))
    mu <- mu.nb(xi_mu)
    result <- (-1) * (t(X) %*% ((phi * (y + phi) * (mu / (mu + phi) ^ 2)) * X))
    return(result)
  }
  
  Grad.phi <- function(xi_phi, xi_mu) {
    phi <- exp(as.numeric(Z %*% xi_phi))
    mu <- mu.nb(xi_mu)
    result <- colSums((digamma(y + phi) * phi - digamma(phi) * phi + 
                         phi * (1 + as.numeric(Z %*% xi_phi)) - 
                         phi * (log(mu + phi) + ((y + phi) / (mu + phi))))* Z)
    return(result)
  }
  
  Hess.phi <- function(xi_phi, xi_mu) {
    phi <- exp(as.numeric(Z %*% xi_phi))
    mu <- mu.nb(xi_mu)
    term1 <- t(Z) %*% (((digamma(y + phi) + (trigamma(y + phi) * phi)) * phi) * Z)
    term2 <- (-1) * (t(Z) %*% (((digamma(phi) + (trigamma(phi) * phi)) * phi) * Z))
    term3 <- t(Z) %*% ((phi * (2 + as.numeric(Z %*% xi_phi))) * Z)
    term4 <- (-1) * (t(Z) %*% ((log(mu + phi) * phi + phi^2 / (mu + phi) + 
                                  ((y * mu + 2 * mu * phi + phi^2) / (mu + phi)^2) * phi) * Z))
    result <- term1 + term2 + term3 + term4
    
    return(result)
  }
  
  
  Hess.theta.phi <- function(xi_mu, xi_phi) {
    phi <- exp(as.numeric(Z %*% xi_phi))
    mu <- mu.nb(xi_mu)
    term1 <- (-1) * t(X) %*% (((phi * mu) / (mu + phi)) * Z)
    term2 <- t(X) %*% ((((y + phi) * (phi * mu)) / (mu + phi)^2) * Z)
    result <- term1 + term2
    return(result)
  }
  
  
  Grad.logptheta <- function(xi_mu, xi_phi, v){
    
    Qv.mu <- Qv.mu(v)
    result <- as.numeric(c(Grad.theta(xi_mu = xi_mu, xi_phi = xi_phi)) - Qv.mu%*%xi_mu)
    return(as.numeric(result))
  }
  
  Hess.logptheta <- function(xi_mu, xi_phi, v){
    Qv.mu <- Qv.mu(v)
    result <- Hess.theta(xi_mu = xi_mu, xi_phi = xi_phi) - Qv.mu
    return(result)
  }
  
  Grad.logpphi <- function(xi_phi, xi_mu, v){
    Qv.phi <- Qv.phi(v)
    result <- as.numeric(Grad.phi(xi_phi = xi_phi, xi_mu = xi_mu)) - Qv.phi%*%xi_phi
    return(as.numeric(result))
  }
  
  Hess.logpphi <- function(xi_phi, xi_mu, v){
    Qv.phi <- Qv.phi(v)
    result <- Hess.phi(xi_phi = xi_phi, xi_mu = xi_mu) - Qv.phi
    return(result)
  }
  
  # Initial values for log-penalty and log-overdispersion parameter
  v_init <- c(1, 1)
  
  NR_ximu <- function(xi_mu0, xi_phi, v){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dximu <- as.numeric((-1) * solve(Hess.logptheta(xi_mu = xi_mu0, xi_phi = xi_phi, v),
                                       Grad.logptheta(xi_mu = xi_mu0, xi_phi = xi_phi, v)))
      
      ximu.new <- xi_mu0 + dximu
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi_mu = xi_mu0, xi_phi = xi_phi, v)
      while (log_pxi(xi_mu = ximu.new, xi_phi = xi_phi, v) <= logpxi.current) {
        step <- step * .5
        ximu.new <- xi_mu0 + (step * dximu)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((ximu.new - xi_mu0) ^ 2))
      iter <- iter + 1
      xi_mu0 <- ximu.new
      if(dist < epsilon) break
    }
    
    ximu.star <- xi_mu0
    return(ximu.star)
  }
  
  
  NR_xiphi <- function(xi_phi0, xi_mu, v){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxiphi <- as.numeric((-1) * solve(Hess.logpphi(xi_phi = xi_phi0, xi_mu = xi_mu, v),
                                        Grad.logpphi(xi_phi = xi_phi0, xi_mu = xi_mu, v)))
      xiphi.new <- xi_phi0 + dxiphi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi_mu = xi_mu, xi_phi = xi_phi0, v)
      while (log_pxi(xi_mu = xi_mu, xi_phi = xiphi.new, v) <= logpxi.current) {
        step <- step * .5
        xiphi.new <- xi_phi0 + (step * dxiphi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xiphi.new - xi_phi0) ^ 2))
      iter <- iter + 1
      xi_phi0 <- xiphi.new
      if(dist < epsilon) break
    }
    
    xiphi.star <- xi_phi0
    return(xiphi.star)
  }
  
  # Initialize
  v_init      <- c(1, 1)
  xi_mu0  <- rep(0.1, (pB * pW) + pW)
  xi_phi0 <- rep(0.1, pZ)
  ximu_init   <- NR_ximu(xi_mu0, xi_phi0, v_init)
  xiphi_init  <- NR_xiphi(xi_phi0, ximu_init, v_init)
  xi_init     <- c(ximu_init, xiphi_init)
  
  # Hessian at init
  H_mm <- Hess.theta(ximu_init, xiphi_init)
  H_mp <- Hess.theta.phi(ximu_init, xiphi_init)
  H_pp <- Hess.phi(xiphi_init, ximu_init)
  Hess.xi_init <- rbind(cbind(H_mm, H_mp), 
                        cbind(t(H_mp), H_pp))
  
  # Initialize mu and phi
  mu_init <- mu.nb(ximu_init)
  phi_init <- exp(as.numeric(Z %*% xiphi_init))
  
  # Log-posterior for v
  a.delta <- 1e-5; b.delta <- 1e-5; nu <- 3
  log_pv <- function(v) {
    Qv.xi <- Qv.xi(v)
    eigs  <- eigen(-(Hess.xi_init - Qv.xi), only.values = TRUE)$values
    a1    <- -0.5 * sum(xi_init * (Qv.xi %*% xi_init))
    a2    <- sum(lgamma(y + phi_init) - lgamma(phi_init) +
                   y * (X %*% ximu_init + log(offset)) +
                   phi_init * (Z %*% xiphi_init) -
                   (y + phi_init) * log(mu_init + phi_init))
    a3    <- -0.5 * sum(log(eigs[eigs > 0]))
    a4    <- 0.5 * (pB*pW + nu) * v[1] + 0.5 * (pZ + nu) * v[2]
    a5    <- - (0.5 * nu + a.delta) * (log(b.delta + 0.5 * nu * exp(v[1])) +
                                         log(b.delta + 0.5 * nu * exp(v[2])))
    a1 + a2 + a3 + a4 + a5
  }
  
  
  # Optimize v
  v_mode <- stats::optim(par = c(1, 1), fn = log_pv, method = "Nelder-Mead",
                         control = list(fnscale = -1, reltol = 1e-5))$par
  
  
  # Re-fit xi at v_mode
  ximu_mode  <- NR_ximu(ximu_init, xiphi_init, v_mode)
  xiphi_mode <- NR_xiphi(xiphi_init, ximu_mode, v_mode)
  
  # Compute Hessian at mode
  H_mm_m <- Hess.theta(ximu_mode, xiphi_mode)
  H_mp_m <- Hess.theta.phi(ximu_mode, xiphi_mode)
  H_pp_m <- Hess.phi(xiphi_mode, ximu_mode)
  Hess.xi_m <- rbind(cbind(H_mm_m, H_mp_m),
                     cbind(t(H_mp_m),H_pp_m))
  Qv.xi_m <- Qv.xi(v_mode)
  Hess.logpxi_m <- Hess.xi_m - Qv.xi_m
  
  # Covariance of posterior
  Covar <- -solve(Hess.logpxi_m)
  
  ##########################################################################
  # Dispersion surface
  Zi_grid <- cubicbs(min_age_part:max_age_part, lower = min_age_part,
                     upper = max_age_part, K = K.disp)$Bmatrix
  Zj_grid <- cubicbs(min_age_cont:max_age_cont, lower = min_age_cont,
                     upper = max_age_cont, K = K.disp)$Bmatrix
  Z_grid  <- kronecker(Zi_grid, Zj_grid)
  disp    <- exp(Z_grid %*% xiphi_mode)
  disp_mat <- matrix(disp, nrow = length(min_age_part:max_age_part),
                     ncol = length(min_age_cont:max_age_cont), byrow = FALSE)
  
  ##########################################################################
  # Information criteria
  log_lik <- function(xi_mu, xi_phi) {
    mu  <- mu.nb(xi_mu)
    phi <- exp(as.numeric(Z %*% xi_phi))
    sum(lgamma(y + phi) - lgamma(phi) +
          y * (X %*% xi_mu + log(offset)) +
          phi * (Z %*% xi_phi) -
          (y + phi) * log(mu + phi))
  }

  
  ED <- diag(Covar %*% (-Hess.xi_m))  
  sumED <- sum(ED)
  BIC     <- -2*log_lik(ximu_mode, xiphi_mode) + sumED * log(length(y))
  twologL <- 2 * log_lik(ximu_mode, xiphi_mode)
  n       <- length(y)
  
  
  # WAIC
  if (WAIC) {
    L <- L  # number of posterior samples
    xi_mode <- c(ximu_mode, xiphi_mode)
    nxi <- length(ximu_mode) + length(xiphi_mode)
    xi_phi_index <- (nxi - (length(xiphi_mode) - 1)):nxi
    xi_mu_index <- 1:(nxi - length(xiphi_mode))
    
    # Sample from approximate posterior
    xi_samples <- MASS::mvrnorm(n = L, mu = xi_mode, Sigma = Covar)
    
    # Log conditional posterior of xi given v
    log_lik_per_obs <- function(xi_mu, xi_phi) {
      mu <- mu.nb(xi_mu)
      phi <- exp(as.numeric(Z %*% xi_phi))
      
      result <- lgamma(y + phi) - lgamma(phi) +
        (y * (as.numeric(X %*% xi_mu) + log(offset))) +
        phi * as.numeric(Z %*% xi_phi) -
        (y + phi) * log(mu + phi)
      
      return(as.numeric(result))
    }
    loglik_matrix <- matrix(NA, nrow = L, ncol = n)
    for (ll in 1:L) {
      loglik_matrix[ll, ] <- log_lik_per_obs(
        xi_mu = xi_samples[ll, xi_mu_index],
        xi_phi = xi_samples[ll, xi_phi_index]
      )
    }
    
    # WAIC components
    lppd <- sum(log(colMeans(exp(loglik_matrix))))  # log pointwise predictive density
    p_waic <- sum(apply(loglik_matrix, 2, var))      # effective number of parameters
    waic <- -2 * (lppd - p_waic)
  } else {
    lppd <- p_waic <- waic <- NULL
  }
  
  # Return results
  list(
    data_sub      = data_sub,
    ximu_mode     = ximu_mode,
    xiphi_mode    = xiphi_mode,
    lambda        = exp(v_mode),
    disp_surface  = disp_mat,
    Covariance    = Covar,
    BIC           = BIC,
    ED            = sumED,
    twologL       = twologL,
    n             = n,
    WAIC          = waic,
    lppd          = lppd,
    p_waic        = p_waic,
    pB            = pB,
    pZ            = pZ,
    pW            = pW,
    W             = W,
    K.age         = K.age,
    K.disp         = K.disp
  )
}
