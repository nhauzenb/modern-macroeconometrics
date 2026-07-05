#### --------------------------------------------------------------------- ####
#### --------                Main estimation script              --------- ####
#### --------   Machine-learning Growth-at-Risk: BART / GP / BNN --------- ####
#### --------------------------------------------------------------------- ####
rm(list = ls())
set.seed(2026)   
fig.dir <- "figs/"; dir.create(fig.dir, showWarnings = FALSE)

# colours of the slide deck
bblue <- rgb( 32,  66, 133, maxColorValue = 255)
rred  <- rgb(194,  40,  46, maxColorValue = 255)

#### --------------------------------------------------------------------- ####
#### ------------------    Main MCMC functions     ----------------------- ####
#### --------------------------------------------------------------------- ####
# (1) Bayesian additive regression trees via dbarts package
mcmc.bart <- function(y, X, nburn = 2500, nsave = 2500, ntree = 250,
                     keeptrees = FALSE){
  library(dbarts)
  X  <- data.frame(NFCI = X)
  bf <- bart(x.train = X, y.train = y, x.test = X, ntree = ntree,
             ndpost = nsave, nskip = nburn, keeptrees = keeptrees, verbose = FALSE)
  sig <- tail(bf$sigma, nsave)            # align sigma draws with the saved f draws
  list(f = bf$yhat.test, sig2 = sig^2,    # yhat.test: ndpost x N (original y-scale)
       bf = if (keeptrees) bf else NULL)
}

# (2) Gaussian process (GP) regression
mcmc.gp <- function(y, X, ell, nburn = 2500, nsave = 2500){
  N  <- length(y)
  D2 <- as.matrix(dist(X))^2; sig.f <- as.numeric(var(y))/2
  Kn <- sig.f * exp(- 0.5*D2 / ell^2)       # signal-variance scaling
  sig2_draw <- 0.1; ntot <- nburn + nsave; si <- 0
  f_store <- matrix(NA, nsave, N)
  sig2_store <- rep(NA, nsave)
  pb <- txtProgressBar(1, ntot, style = 3)
  for (irep in 1:ntot){
    # posterior mean (fhat) and covariance (Vhat) of the function values,
    # given the current error variance
    Dinv <- solve(Kn + diag(sig2_draw, N))
    Vhat <- Kn - Kn %*% Dinv %*% Kn
    fhat <- Kn %*% Dinv %*% y
    # draw the function values (robust fallback if the Cholesky fails)
    f_draw <- try(fhat + t(chol(Vhat + diag(1e-8, N))) %*% rnorm(N), silent = TRUE)
    if (is(f_draw, "try-error"))
      f_draw <- t(mvtnorm::rmvnorm(1, fhat, Matrix::forceSymmetric(Vhat)))
    # update the error variance from its inverse-gamma full conditional
    eps_draw <- y - f_draw
    sig2_draw   <- 1/rgamma(1, 10 + N/2, 0.1 + sum(eps_draw^2)/2)   # IG(10, 0.1) prior
    if (irep > nburn){ si <- si + 1
      f_store[si, ] <- as.numeric(f_draw); sig2_store[si] <- sig2_draw }
    setTxtProgressBar(pb, irep)
  }
  close(pb); list(f = f_store, sig2 = sig2_store)
}

# (3) Bayesian neural network 

## Auxiliary functions
## Activation functions
act.fc.set.all <- list(
  "tanh"      = list("func" = function(x) tanh(x),
                     "grad" = function(x) 1 - (tanh(x))^2),
  "relu"      = list("func" = function(x) (abs(x) + x)/2,
                     "grad" = function(x) ifelse(x >= 0, 1, 0)),
  "sigmoid"   = list("func" = function(x) 1/(1 + exp(-x)),
                     "grad" = function(x) (1/(1 + exp(-x)))*(1 - (1/(1 + exp(-x))))),
  "leakyrelu" = list("func" = function(x) ifelse(x >= 0, x, 0.01*x),
                     "grad" = function(x) ifelse(x >= 0, 1, 0.01)))

##  Horseshoe prior for BNN parameters
get.hs.bnn <- function(bdraw, lambda.hs, nu.hs, tau.hs, zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1
  }else{
    tau.hs <- 1/rgamma(1, shape = (k+1)/2, rate = 1/zeta.hs + sum(bdraw^2/lambda.hs)/2)
  }
  lambda.hs <- 1/rgamma(k, shape = 1, rate = 1/nu.hs + bdraw^2/(2*tau.hs))
  nu.hs     <- 1/rgamma(k, shape = 1, rate = 1 + 1/lambda.hs)
  zeta.hs   <- 1/rgamma(1, shape = 1, rate = 1 + 1/tau.hs)
  list("psi" = lambda.hs*tau.hs, "lambda" = lambda.hs, "tau" = tau.hs,
       "nu" = nu.hs, "zeta" = zeta.hs)
}

## Hamiltonian Monte Carlo (HMC) step for weights 
hmc_deep <- function(theta, f, grad_f, f_list, epsilon = .1, L = 10){
  p <- rnorm(length(theta), 0, 1)
  theta_tilde <- theta
  p_tilde <- p + epsilon*grad_f(theta = theta_tilde, f_list = f_list)/2
  for (i in 1:L){
    theta_tilde <- theta_tilde + epsilon*p_tilde
    if (i != L) p_tilde <- p_tilde + epsilon*grad_f(theta = theta_tilde, f_list = f_list)
  }
  p_tilde <- p_tilde + epsilon*grad_f(theta = theta_tilde, f_list = f_list)/2
  p_tilde <- -p_tilde
  log.prob_acc  <- f(theta = theta, f_list = f_list)
  log.corr_acc  <- sum(p^2)/2
  log.prob_prop <- f(theta = theta_tilde, f_list = f_list)
  log.corr_prop <- sum(p_tilde^2)/2
  if (isTRUE(log(runif(1)) < log.prob_prop - log.prob_acc + log.corr_acc - log.corr_prop)){
    theta <- theta_tilde
  }
  theta
}

## Log conditional posterior
get.post_k <- function(theta, f_list){
  k_draw <- f_list$k_draw.nr;  k.V <- f_list$k.V
  nr1 <- f_list$nr1; nr2 <- f_list$nr2; QQ <- f_list$QQ; Q <- f_list$Q
  MM <- f_list$MM; acf_draw <- f_list$acf_draw
  y <- f_list$y; X.hat.nr <- f_list$X.hat.nr; X.hat.wonr <- f_list$X.hat.wonr
  b_draw.nr <- f_list$b_draw.nr; b_draw.wonr <- f_list$b_draw.wonr
  sig2_draw <- f_list$sig2_draw; acf_set <- f_list$acf_set

  k_draw[1:MM[nr1],,nr1] <- theta
  fit.wonr <- X.hat.wonr[,,QQ] %*% b_draw.wonr
  for (nn in 1:Q){                       # matrix() guards against dropped dims (K = 1)
    Mlay <- MM[nn]
    X.hat.nr[,nr2,nn+1] <- acf_set[[acf_draw[nn]]][["func"]](matrix(X.hat.nr[,1:Mlay,nn], ncol = Mlay) %*% matrix(k_draw[1:Mlay,,nn], nrow = Mlay))
  }
  fit.nr <- as.matrix(X.hat.nr[,nr2,QQ]) %*% b_draw.nr

  loglik  <- sum(dnorm(y, fit.wonr + fit.nr, sqrt(sig2_draw), log = TRUE))
  logpr   <- sum(dnorm(theta, 0, sqrt(k.V), log = TRUE))
  loglik + logpr
}

## Log conditional posterior gradient
get.post.grad_k <- function(theta, f_list){
  k_draw <- f_list$k_draw.nr;  k.V <- f_list$k.V
  nr1 <- f_list$nr1; nr2 <- f_list$nr2; QQ <- f_list$QQ; Q <- f_list$Q
  MM <- f_list$MM; acf_draw <- f_list$acf_draw
  y <- f_list$y; X.hat.nr <- f_list$X.hat.nr; X.hat.wonr <- f_list$X.hat.wonr
  b_draw.nr <- f_list$b_draw.nr; b_draw.wonr <- f_list$b_draw.wonr
  sig2_draw <- f_list$sig2_draw; acf_set <- f_list$acf_set

  k_draw[1:MM[nr1],,nr1] <- theta
  normalizer <- 1/sqrt(sig2_draw)
  fit.wonr <- (X.hat.wonr[,,QQ]*normalizer) %*% b_draw.wonr

  chain.grad <- matrix(1, length(y), 1)
  for (nn in 1:Q){                     
    Mlay <- MM[nn]
    z.nn <- matrix(X.hat.nr[,1:Mlay,nn], ncol = Mlay) %*% matrix(k_draw[1:Mlay,,nn], nrow = Mlay)
    X.hat.nr[,nr2,nn+1] <- acf_set[[acf_draw[nn]]][["func"]](z.nn)
    if (nn > nr1) chain.grad <- chain.grad * acf_set[[acf_draw[nn]]][["grad"]](z.nn) * k_draw[nr2,1,nn]
  }
  fit.nr <- as.matrix(X.hat.nr[,nr2,QQ]*normalizer) %*% b_draw.nr

  yy <- (y*normalizer - fit.wonr)
  X.in1   <- matrix(X.hat.nr[,1:MM[nr1],nr1], ncol = MM[nr1])
  dhqdkq  <- as.numeric(b_draw.nr)*chain.grad*acf_set[[acf_draw[nr1]]][["grad"]](X.in1 %*% matrix(theta, ncol = 1))*normalizer
  dloglik <- crossprod(X.in1, (yy - fit.nr)*dhqdkq)
  dlogpr  <- -1/k.V*theta
  dloglik + dlogpr
}

# Main BNN sampler
mcmc.bnn <- function(y, X, M = 20, act = "tanh",
                    nburn = 2500, nsave = 2500, nthin = 2,
                    eps = 0.01, L = 20, s_pr = 5, S_pr = 0.5){
  
  library(MASS)
  library(Matrix)
  
  Q  <- 1                                    # Single hidden layer
  
  y.mu <- mean(y); y.sd <- sd(y)
  y <- matrix((y - y.mu)/y.sd, ncol = 1)
  X <- matrix((X - mean(X))/sd(X), ncol = 1)
  X <- cbind(X, "cons" = 1)                  # constant column: first-layer bias
  
  acf_set  <- act.fc.set.all[act]            # activation function (fixed)
  acf_draw <- rep(1, Q)                      # activation indicator
  
  ## Key dimensions
  K  <- ncol(X); N <- nrow(X)
  QQ <- Q + 1
  MM <- c(K, rep(M+1, Q))                    # layer input dims (incl. bias input)
  Kmax <- max(K, M+1)
  
  ## Design matrices, priors, and starting values
  k_draw <- k.V <- array(0, dim = c(Kmax, M, Q))
  X.hat  <- array(0, dim = c(N, Kmax, QQ))
  X.hat[,1:K,1] <- X
  
  # 1st layer: K x M (input x neurons)
  k_draw[1:K,,1] <- t(matrix(runif(K*M, 0, 1), M, K))
  X.hat[,1:M,2]  <- matrix(X.hat[,1:K,1], ncol = K) %*% matrix(k_draw[1:K,,1], nrow = K)/K
  X.hat[,M+1,2]  <- 1                        # bias input for the next layer
  k.V[1:K,,1]    <- matrix(1e-10, K, M)
  
  # Output layer (beta) and its HS prior
  b_draw <- matrix(0, M, 1)
  b.v <- rep(1, M); b.v.inv <- 1/b.v
  lambda.beta.mat <- matrix(0.1, M, 1); nu.beta.mat <- matrix(0.1, M, 1)
  tau.beta <- 0.1; zeta.beta <- 0.1
  
  # Horseshoe prior on the hidden-layer weights (one column per neuron)
  lam.mat <- nu.mat <- tau.mat <- zeta.mat <- list()
  lam.mat[[1]] <- matrix(0.1, K, M); nu.mat[[1]] <- matrix(0.1, K, M)
  tau.mat[[1]] <- matrix(0.1, M, 1); zeta.mat[[1]] <- matrix(0.1, M, 1)
  
  g_draw <-  matrix(0,K,1)                   # linear-part coefficients (gamma)
  g.v <- rep(1, K); g.v.inv <- 1/g.v 

  # Horseshoe prior on the linear-part coefficients
  g.lam.mat <- matrix(0.1, K, 1)
  g.nu.mat  <- matrix(0.1,  K, 1)
  g.tau     <- 0.1
  g.zeta    <- 0.1
  
  sig2_draw <- 0.1
  acc.k <- matrix(0, M, Q)

  ## Storage
  ntot <- nburn + nsave*nthin
  save.set <- seq(nthin, nsave*nthin, nthin) + nburn
  save.ind <- 0
  f_store  <- matrix(NA, nsave, N)
  sig2_store <- rep(NA, nsave)

  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  for (irep in seq_len(ntot)){
    ## Step 1: jointly draw the linear coefficients (gamma) and the output
    ##         weights (beta) from their Gaussian full conditional
    norm.sig <- 1/sqrt(sig2_draw)
    y.lin <- y*norm.sig
    x.lin <- cbind(X, X.hat[,1:M,QQ])*norm.sig
    
    gb.V_po <- try(solve(crossprod(x.lin) + diag(c(g.v.inv, b.v.inv))), silent=F) # Conditional posterior variance-covariance of beta and gamma
    if (is(gb.V_po,"try-error")) gb.V_po <- ginv(crossprod(x.lin) + diag(c(g.v.inv, b.v.inv)))
    gb.m_po <- gb.V_po%*%crossprod(x.lin, y.lin) # Conditional posterior mean of beta and gamma
    gb_draw <- try(gb.m_po + t(chol(gb.V_po))%*%rnorm(K+M), silent=F) # Simulate from multivariate normal distribution
    if (is(gb_draw, "try-error")) gb_draw <- matrix(as.numeric(mvtnorm::rmvnorm(1, gb.m_po, as.matrix(forceSymmetric(gb.V_po)))), K+M,1)
    
    g_draw <- gb_draw[1:K,,drop = F]
    b_draw <- gb_draw[(K+1):(K+M),, drop = F]
    
    # Fit of linear part:
    fit_lin <- X%*%g_draw 
    y.nolin <- y - fit_lin
    
    ## Step 2.1: HS prior variances for gamma
    g_hs <- get.hs.bnn(g_draw,lambda.hs = g.lam.mat, nu.hs = g.nu.mat, tau.hs = g.tau,zeta.hs=g.zeta)
    g.lam.mat <- g_hs$lambda  # Local scales
    g.nu.mat  <- g_hs$nu
    g.tau     <- g_hs$tau     # Global scales
    g.zeta    <- g_hs$zeta
    g.v       <- g_hs$psi     # Diagonal elements of prior variance matrix
    g.v[g.v > 1e-3] <- 1e-3   # numerical ceiling on the prior variance
    g.v.inv <- 1/g.v          # Diagonal elements prior precision matrix
    
    
    ## Step 2.2: HS prior variances for beta
    hs.beta <- get.hs.bnn(b_draw, lambda.hs = lambda.beta.mat, nu.hs = nu.beta.mat,
                          tau.hs = tau.beta, zeta.hs = zeta.beta)
    lambda.beta.mat <- hs.beta$lambda; nu.beta.mat <- hs.beta$nu
    tau.beta <- hs.beta$tau; zeta.beta <- hs.beta$zeta
    b.v <- hs.beta$psi
    b.v.inv <- 1/b.v

    
    ## Step 3: update the hidden-layer weights one neuron at a time (HMC)
    nr1 <- 1
    for (nr2 in seq_len(M)){
    theta <- k_draw[1:MM[nr1],nr2,nr1]
    wonr.slct   <- setdiff(1:M, nr2)
    X.hat.wonr  <- X.hat[,wonr.slct,, drop = FALSE]
    b_draw.nr   <- b_draw[nr2,, drop = FALSE]
    b_draw.wonr <- b_draw[wonr.slct,, drop = FALSE]
    k_draw.nr   <- k_draw[,nr2,, drop = FALSE]
    k_draw.nr[1:MM[nr1],,nr1] <- theta
    k_star <- hmc_deep(theta   = theta,
                           f       = get.post_k,
                           grad_f  = get.post.grad_k,
                           f_list  = list(k_draw.nr = k_draw.nr, y = y.nolin,
                                          X.hat.nr = X.hat, X.hat.wonr = X.hat.wonr,
                                          k.V = k.V[1:MM[nr1],nr2,nr1],
                                          nr1 = nr1, nr2 = nr2, QQ = QQ, Q = Q, MM = MM,
                                          acf_draw = acf_draw,
                                          b_draw.nr = b_draw.nr, b_draw.wonr = b_draw.wonr,
                                          sig2_draw = sig2_draw, acf_set = acf_set),
                           epsilon = eps, L = L)
      accept <- !identical(as.vector(k_star), as.vector(theta))
      if (accept){
        k_draw[1:MM[nr1],nr2,nr1] <- k_star
        acc.k[nr2,nr1] <- acc.k[nr2,nr1] + 1
        X.hat[,nr2,nr1+1] <- acf_set[[acf_draw[nr1]]][["func"]](matrix(X.hat[,1:MM[nr1],nr1], ncol = MM[nr1]) %*% matrix(k_star, ncol = 1))
      }
    }

     ## Step 4: horseshoe shrinkage variances for the hidden-layer weights
    for (j in 1:M){
        k_hs.j <- get.hs.bnn(bdraw     = k_draw[1:MM[nr1],j,nr1],
                             lambda.hs = lam.mat[[nr1]][,j],
                             nu.hs     = nu.mat[[nr1]][,j],
                             tau.hs    = tau.mat[[nr1]][j,1],
                             zeta.hs   = zeta.mat[[nr1]][j,1])
      lam.mat[[nr1]][,j]   <- k_hs.j$lambda
      nu.mat[[nr1]][,j]    <- k_hs.j$nu
      tau.mat[[nr1]][j,1]  <- k_hs.j$tau
      zeta.mat[[nr1]][j,1] <- k_hs.j$zeta
      k.V[1:MM[nr1],j,nr1] <- k_hs.j$psi
    }
    k.V[,,nr1][k.V[,,nr1] > 10]   <- 10
    k.V[,,nr1][k.V[,,nr1] < 1e-5] <- 1e-5

    ## Step 5: recompute the hidden-layer outputs with the updated weights
    for (nn in 1:Q){
      X.hat[,1:M,nn+1]  <- acf_set[[acf_draw[nn]]][["func"]](matrix(X.hat[,1:MM[nn],nn], ncol = MM[nn]) %*% matrix(k_draw[1:MM[nn],,nn], nrow = MM[nn]))
      X.hat[,M+1,nn+1]  <- 1                 # bias input for the next layer
    }
    

    # Fit of the neural network
    fit_nn   <- X.hat[,1:M,QQ]%*%b_draw
    fit_nn <- fit_nn - mean(fit_nn)          # demean: level is identified via the linear part
    fit  <- fit_lin + fit_nn
    
    ## Step 6: draw the error variance from its inverse-gamma full conditional
    s_po <- s_pr + N
    S_po <- S_pr + as.numeric(crossprod(y - fit))
    sig2_draw <- 1/rgamma(1, s_po/2, S_po/2)

    ## Storage (on the original y-scale)
    if (irep %in% save.set){
      save.ind <- save.ind + 1
      f_store[save.ind, ]  <- as.numeric(fit)*y.sd + y.mu
      sig2_store[save.ind]   <- sig2_draw * y.sd^2
    }
    setTxtProgressBar(pb, irep)
  }
  close(pb)
  cat(sprintf("\n[BNN: Q = %d, M = %d, %s] mean HMC acceptance: %.1f%%\n",
              Q, M, act, 100*mean(acc.k/ntot)))
  list(f = f_store, sig2 = sig2_store, acc = acc.k/ntot)
}

#### --------------------------------------------------------------------- ####
#### -------------------  Plot helpers (per specification)  -------------- ####
#### --------------------------------------------------------------------- ####
plot.curve <- function(fit, x, y, h, file){
  ord  <- order(x)
  f.lo <- apply(fit$f, 2, quantile, 0.05)
  f.md <- apply(fit$f, 2, quantile, 0.50)
  f.hi <- apply(fit$f, 2, quantile, 0.95)
  pdf(file, width = 8.5, height = 5.4)
  par(mar = c(3.6, 3.8, 1.0, 0.8), mgp = c(2.4, 0.7, 0))
  plot(x, y, pch = 16, col = adjustcolor("grey55", 0.6), cex = 0.8,
       xlab = "NFCI",
       ylab = if (h == 4) "One-year-ahead GDP growth" else paste0(h, "-quarter-ahead GDP growth"))
  abline(h = 0, v = 0, lty = 3, col = "grey50")
  polygon(c(x[ord], rev(x[ord])), c(f.lo[ord], rev(f.hi[ord])),
          col = adjustcolor(rred, 0.20), border = NA)
  lines(x[ord], f.md[ord], col = rred, lwd = 2.6)
  abline(lm(y ~ x), col = bblue, lwd = 2, lty = 2)          # linear (OLS) reference
  legend("topright", bty = "n", cex = 0.85,
         legend = c("Realisations", "Posterior median", "90% credible band", "OLS"),
         pch = c(16, NA, 15, NA), lty = c(NA, 1, NA, 2),
         col = c("grey55", rred, adjustcolor(rred, 0.4), bblue), lwd = c(NA, 2.6, 6, 2))
  dev.off()
}

plot.fan <- function(fit, dts, y, file){
  fan.tau <- seq(0.05, 0.95, 0.05)
  ypred   <- fit$f + sqrt(fit$sig2) * matrix(rnorm(length(fit$f)), nrow(fit$f), ncol(fit$f))
  qfan    <- t(apply(ypred, 2, quantile, probs = fan.tau))
  P       <- length(fan.tau)
  pal     <- colorRampPalette(c(bblue, "white"))(ceiling(P/2) + 1)
  pdf(file, width = 9, height = 5)
  par(mar = c(3.4, 3.8, 2, 0.8), mgp = c(2.4, 0.7, 0))
  plot(dts, y, type = "n", ylim = range(qfan, y),
       xlab = "Year", ylab = expression(y[t]))
  for (j in 1:floor(P/2))                                # nested symmetric quantile bands
    polygon(c(dts, rev(dts)), c(qfan[, j], rev(qfan[, P + 1 - j])), col = pal[j], border = NA)
  points(dts, y, pch = 16, cex = 0.5, col = rred)        # realised growth
  abline(h = 0, col = "grey55", lty = 2)
  legend("topleft", bty = "n", cex = 0.8,
         legend = c("Realised growth", "90% credible set"),
         pch = c(16, 15), col = c(rred, pal[2]))
  box(col = "grey60")
  dev.off()
}


#### --------------------------------------------------------------------- ####
#### ------------  Data and settings (as in Session 5)  ------------------ ####
#### --------------------------------------------------------------------- ####
str.smp <- 1973; end.smp <- 2025 + 3/4; h <- 4   # h = horizon (quarters); 4 = one-year-ahead GaR
Y.slct  <- "GDPC1_PCH"

# Growth-at-Risk (Adrian, Boyarchenko & Giannone, 2019): annualised h-quarter
# average GDP growth y_t = 4*mean(g_{t-h+1..t}), paired with the NFCI dated t-h
gar <- read.csv("data/GaR-data.csv", stringsAsFactors = FALSE)
gar <- ts(gar[, c("GDPC1_PCH", "NFCI")], start = c(1970,1), frequency = 4)
gar <- window(gar, start = c(1971,1))

g    <- as.numeric(gar[, "GDPC1_PCH"])
nfci <- as.numeric(gar[, "NFCI"])
n    <- length(g)
y4   <- rep(NA_real_, n)
for (t in h:n) y4[t] <- 4*mean(g[(t-h+1):t])       # annualised trailing h-qtr avg
gar[, "GDPC1_PCH"] <- y4                           # GaR target
gar[, "NFCI"]      <- c(rep(NA,h), head(nfci, -h)) # predictor dated t-h
gar <- window(gar, start = str.smp, end = end.smp) # drop leading NA quarters

data <- list(y = gar[, Y.slct], X = gar[, "NFCI", drop = FALSE],
             dates = time(gar), target = Y.slct, h = h)
dts <- as.numeric(data$dates)
y   <- as.numeric(data$y)
X  <- as.numeric(data$X)
N   <- length(y)

#### --------------------------------------------------------------------- ####
#### ---------  Estimate one model per specification  -------------------- ####
#### --------------------------------------------------------------------- ####
spec.list <- list(
  list(model = "bart", ntree = 1),
  list(model = "bart", ntree = 25),
  list(model = "bart", ntree = 250),
  list(model = "gp",   ell = 4),
  list(model = "gp",   ell = 1),
  list(model = "gp",   ell = 0.25),
  list(model = "bnn", M = 20, act = "tanh"),      
  list(model = "bnn", M = 20, act = "relu"),
  list(model = "bnn", M = 20, act = "sigmoid"))

res.spec <- list()

for (sp in spec.list){
  tag <- switch(sp$model,                        
                bart = paste0("bart_S", sp$ntree),
                gp   = paste0("gp_ell", sp$ell),
                bnn  = paste0("bnn_shallow_", sp$act))
  fit <- switch(sp$model,
                bart = mcmc.bart(y, X, ntree = sp$ntree, keeptrees = (sp$ntree == 1)),
                gp   = mcmc.gp(y, X, ell = sp$ell),
                bnn  = mcmc.bnn(y, X, M = sp$M, act = sp$act))
  res.spec[[tag]] <- fit
  pref <- paste0(fig.dir, "S06_", tag, "_GaR_NFCI_")
  plot.curve(fit, X, y, h, paste0(pref, "curve.pdf"))
  plot.fan(fit, dts, y, paste0(pref, "fan.pdf"))
}

