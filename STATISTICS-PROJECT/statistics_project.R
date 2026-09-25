# -------------------------------
# Task 1 (Fisher scoring)
# -------------------------------
makexy <- function(df, yname, drop = character()) {
  xnames <- setdiff(names(df), c(yname, drop))
  for (nm in xnames) if (is.character(df[[nm]])) df[[nm]] <- factor(df[[nm]])
  form_text <- paste(yname, "~", paste(xnames, collapse = " + "))
  X <- model.matrix(as.formula(form_text), df)
  y <- df[[yname]]
  list(X = X, y = y)
}

softplus <- function(x) ifelse(x > 30, x, log1p(exp(x)))

fisher_scoring <- function(X, y, model = c("logistic", "poisson"),
                           tol = 1e-6, maxit = 200) {
  
  model <- match.arg(model)
  p <- ncol(X)
  beta <- rep(0, p)
  
  for (iter in 1:maxit) {
    
    eta <- as.vector(X %*% beta)
    
    if (model == "logistic") {
      mu <- 1 / (1 + exp(-eta))
      w <- mu * (1 - mu)
      w <- pmax(w, 1e-10)
      
      score <- t(X) %*% (y - mu)
      info_mat <- (t(X) * w) %*% X
      ll_now <- sum(y * eta - softplus(eta))
    } else {
      mu <- exp(eta)
      w <- mu
      w <- pmax(w, 1e-10)
      
      score <- t(X) %*% (y - mu)
      info_mat <- (t(X) * w) %*% X
      ll_now <- sum(y * eta - mu)
    }
    
    step <- solve(info_mat, score)
    
    step_s <- 1
    b_new <- beta
    new_ll <- ll_now
    
    for (k in 1:12) {
      b_try <- beta + step_s * step
      eta_try <- as.vector(X %*% b_try)
      
      try_ll <- if (model == "logistic") {
        sum(y * eta_try - softplus(eta_try))
      } else {
        sum(y * eta_try - exp(eta_try))
      }
      
      if (is.finite(try_ll) && try_ll >= ll_now) {
        b_new <- b_try
        new_ll <- try_ll
        break
      }
      step_s <- step_s / 2
    }
    
    if (max(abs(b_new - beta)) < tol) {
      beta <- b_new
      break
    }
    beta <- b_new
  }
  
  eta_final <- as.vector(X %*% beta)
  
  if (model == "logistic") {
    mu_final <- 1 / (1 + exp(-eta_final))
    w_final <- mu_final * (1 - mu_final)
    w_final <- pmax(w_final, 1e-10)
  } else {
    mu_final <- exp(eta_final)
    w_final <- mu_final
    w_final <- pmax(w_final, 1e-10)
  }
  
  info_mat_final <- (t(X) * w_final) %*% X
  cov_mat <- solve(info_mat_final)
  se <- sqrt(diag(cov_mat))
  
  list(
    beta = beta,
    se = se,
    iter = iter,
    logLik = new_ll,
    info = info_mat_final
  )
}

# -------------------------------
# CreditCard (logistic)
# -------------------------------
credit_data <- read.csv("CreditCard.csv", stringsAsFactors = FALSE)
credit_data <- na.omit(credit_data)

credit_data$card <- as.integer(credit_data$card == "yes")

if ("owner" %in% names(credit_data)) credit_data$owner <- as.integer(credit_data$owner == "yes")
if ("selfemp" %in% names(credit_data)) credit_data$selfemp <- as.integer(credit_data$selfemp == "yes")

drop_credit <- intersect(c("share", "expenditure"), names(credit_data))
credit_xy <- makexy(credit_data, "card", drop_credit)

credit_fit <- fisher_scoring(credit_xy$X, credit_xy$y, model = "logistic")

credit_table <- cbind(
  Estimate = credit_fit$beta,
  SE = credit_fit$se,
  z = credit_fit$beta / credit_fit$se
)
rownames(credit_table) <- colnames(credit_xy$X)

cat("CreditCard Data – Logistic Model\n")
print(round(credit_table, 4))
cat("LogLik:", round(credit_fit$logLik, 4), " Iter:", credit_fit$iter, "\n\n")

# -------------------------------
# azcabgptca (Poisson)
# -------------------------------
los_data <- read.csv("azcabgptca.csv", stringsAsFactors = FALSE)
los_data <- na.omit(los_data)

los_xy <- makexy(los_data, "los")
los_fit <- fisher_scoring(los_xy$X, los_xy$y, model = "poisson")

los_table <- cbind(
  Estimate = los_fit$beta,
  SE = los_fit$se,
  z = los_fit$beta / los_fit$se
)
rownames(los_table) <- colnames(los_xy$X)

cat("\n==============================\n")
cat("Hospital LOS Data – Poisson Model\n")
print(round(los_table, 4))
cat("LogLik:", round(los_fit$logLik, 4), " Iter:", los_fit$iter, "\n")

# -------------------------------
# Task 2 (Monte Carlo Markov Chain)
# -------------------------------
logpost_logistic <- function(beta, X, y, sigma2) {
  eta <- as.vector(X %*% beta)
  ll <- sum(y * eta - softplus(eta))
  lp <- -sum(beta^2) / (2 * sigma2)
  ll + lp
}

logpost_poisson <- function(beta, X, y, sigma2) {
  eta <- as.vector(X %*% beta)
  mu <- exp(eta)
  ll <- sum(y * eta - mu - lgamma(y + 1))
  lp <- -sum(beta^2) / (2 * sigma2)
  ll + lp
}

rw_metropolis_info_mv <- function(logpost, X, y,
                                  n_iter, burn,
                                  s, sigma2,
                                  info_mat,
                                  beta_init = NULL,
                                  tune = TRUE, target = 0.25) {
  
  p <- ncol(X)
  beta <- if (is.null(beta_init)) rep(0, p) else beta_init
  
  chain <- matrix(NA, n_iter, p)
  colnames(chain) <- colnames(X)
  
  # Make info_mat PD if needed (jitter)
  info_mat <- as.matrix(info_mat)
  eps <- 1e-8
  R <- NULL
  for (k in 1:10) {
    try_R <- try(chol(info_mat + diag(eps, p)), silent = TRUE)
    if (!inherits(try_R, "try-error")) {
      R <- try_R
      break
    }
    eps <- eps * 10
  }
  if (is.null(R)) stop("chol() failed: info matrix not positive definite even after jitter.")
  
  logp_curr <- logpost(beta, X, y, sigma2)
  acc <- 0
  
  for (i in 1:n_iter) {
    # step = (info_mat)^(-1/2) z  via backsolve(R, z)
    z <- rnorm(p)
    step <- backsolve(R, z, upper.tri = TRUE)
    beta_prop <- beta + s * step
    
    logp_prop <- logpost(beta_prop, X, y, sigma2)
    loga <- logp_prop - logp_curr
    
    accepted <- 0
    if (is.finite(loga) && log(runif(1)) < loga) {
      beta <- beta_prop
      logp_curr <- logp_prop
      acc <- acc + 1
      accepted <- 1
    }
    
    if (tune && i <= burn) {
      gamma <- 1 / sqrt(i)
      s <- s * exp(gamma * (accepted - target))
    }
    
    chain[i, ] <- beta
  }
  
  list(
    chain = chain,
    post = chain[(burn + 1):n_iter, , drop = FALSE],
    acc_rate = acc / n_iter,
    final_s = s
  )
}

set.seed(1)

mcmc_credit <- rw_metropolis_info_mv(
  logpost = logpost_logistic,
  X = credit_xy$X,
  y = credit_xy$y,
  n_iter = 20000,
  burn = 6000,
  s = 1.0,
  sigma2 = 100,
  info_mat = credit_fit$info,
  beta_init = credit_fit$beta,
  tune = TRUE
)

mcmc_los <- rw_metropolis_info_mv(
  logpost = logpost_poisson,
  X = los_xy$X,
  y = los_xy$y,
  n_iter = 20000,
  burn = 6000,
  s = 1.0,
  sigma2 = 100,
  info_mat = los_fit$info,
  beta_init = los_fit$beta,
  tune = TRUE
)

apply(mcmc_credit$post, 2, mean)
apply(mcmc_credit$post, 2, sd)
mcmc_credit$acc_rate

apply(mcmc_los$post, 2, mean)
apply(mcmc_los$post, 2, sd)
mcmc_los$acc_rate

j_income <- which(colnames(credit_xy$X) == "income")
j_proc <- which(colnames(los_xy$X) == "procedure")

plot(mcmc_credit$post[, j_income], type = "l",
     main = "Trace plot of income coefficient (logistic regression)",
     xlab = "Post-burn iteration",
     ylab = "Coefficient value")

acf(mcmc_credit$post[, j_income],
    main = "Autocorrelation of income coefficient (logistic regression)")

plot(mcmc_los$post[, j_proc], type = "l",
     main = "Trace plot of procedure coefficient (Poisson regression)",
     xlab = "Post-burn iteration",
     ylab = "Coefficient value")

acf(mcmc_los$post[, j_proc],
    main = "Autocorrelation of procedure coefficient (Poisson regression)")
