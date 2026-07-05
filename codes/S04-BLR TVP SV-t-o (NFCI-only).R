#### --------------------------------------------------------------------- ####
#### --------                Main estimation script              --------- ####
#### ----  Bayesian TVP regression: homosk. / t / SV / t-SV / o-SV ------- ####
#### ----  Non-centered parameterisation (NCP) + FFBS  ------------------- ####
#### --------------------------------------------------------------------- ####
rm(list = ls())
set.seed(2026)
plot.dir <- "figs/"; dir.create(plot.dir, showWarnings = FALSE)

library(Rcpp)
library(stochvol)        # required for the SV / SVt / SVo error models
sourceCpp("ffbs.cpp")    # forward-filter backward-sampler for the RW states

# colours of the slide deck
bblue <- rgb( 32,  66, 133, maxColorValue = 255)
rred  <- rgb(194,  40,  46, maxColorValue = 255)

#### --------------------------------------------------------------------- ####
#### ------------------    Auxiliary functions     ----------------------- ####
#### --------------------------------------------------------------------- ####

# Horseshoe update (Makalic & Schmidt, 2016)
get.hs <- function(bdraw, lam.hs, nu.hs, tau.hs, zeta.hs){
  k <- length(bdraw)
  lam.hs  <- 1/rgamma(k, shape = 1,       rate = 1/nu.hs + bdraw^2/(2*tau.hs))
  nu.hs   <- 1/rgamma(k, shape = 1,       rate = 1 + 1/lam.hs)
  tau.hs  <- 1/rgamma(1, shape = (k+1)/2, rate = 1/zeta.hs + sum(bdraw^2/lam.hs)/2)
  zeta.hs <- 1/rgamma(1, shape = 1,       rate = 1 + 1/tau.hs)
  list(psi = lam.hs*tau.hs, lam = lam.hs, tau = tau.hs, nu = nu.hs, zeta = zeta.hs)
}

# Log full-conditional of the t degrees of freedom nu
log.post.nu <- function(nu, sum.log.lam, sum.lam, n, nu0){
  if (nu <= 2) return(-Inf)                          # finite-variance support
  n * ((nu/2)*log(nu/2) - lgamma(nu/2)) +
    (nu/2 - 1)*sum.log.lam - (nu/2)*sum.lam - nu/nu0
}

#### --------------------------------------------------------------------- ####
#### ------------------    Main MCMC sampler       ----------------------- ####
#### --------------------------------------------------------------------- ####
mcmc.sampler <- function(y, X, err.var = c("homo", "t", "SV", "SVt", "SVo"),
                         nburn = 3000, nsave = 6000){

  # Add the intercept as the LAST column; K counts ALL regressors
  X <- cbind(X, "cons" = 1)
  N <- nrow(X); K <- ncol(X); KK <- 2*K     # K constant parts + K NCP scales
  shrink.slct <- 1:(K-1)                     # shrink everything but the intercept

  ## Prior mean (shrink standardised predictors and NCP scales towards zero)
  b_pr <- rep(0, K); o_pr <- rep(0, K); obseq_pr <- c(b_pr, o_pr)

  ## Horseshoe scalings for the constant part b and the NCP scales sqrt(omega)
  lam_b <- rep(1, K-1); nu_b <- rep(1, K-1); tau_b <- 0.5; zeta_b <- 0.5
  lam_o <- rep(1, K-1); nu_o <- rep(1, K-1); tau_o <- 0.5; zeta_o <- 0.5

  v_pr     <- rep(1, KK)
  v_pr[K]  <- 1e4      # constant intercept: essentially unshrunk
  v_pr[KK] <- 1e-10    # TVP intercept switched off (constant mean)
  v_pr.inv <- 1/v_pr

  ## Inverse-Gamma prior for the error variance
  s_pr <- 3; S_pr <- 0.3; s_po <- s_pr + N/2

  ## SV prior (err.var in SV / SVt / SVo); "SVt" gives t-distributed SV innovations
  sv_nu     <- if (err.var == "SVt") sv_exponential(0.1) else sv_infinity()
  sv_priors <- specify_priors(mu = sv_normal(0, 1), phi = sv_beta(25, 1.5),
                              sigma2 = sv_gamma(0.5, 100), nu = sv_nu, rho = sv_constant(0))
  svdraw    <- list(mu = 0, phi = 0.99, sigma = 0.01,
                    nu = if (err.var == "SVt") 10 else Inf, rho = 0, beta = NA, latent0 = 0)
  sv_latent <- rep(0, N)

  ## Stock-Watson outlier component o_t in {1,...,20} (err.var == "SVo")
  ot_draw   <- rep(1, N)
  ot_grid   <- seq(1, 20, by = 1); ot_grid_n <- length(ot_grid)
  ot_p      <- 1e-4
  ot_p_grid <- c(1 - ot_p, rep(ot_p/(ot_grid_n - 1), ot_grid_n - 1))
  ot_Ba <- 1; ot_Bb <- 99

  ## Storage
  ntot <- nburn + nsave; save.set <- 1:nsave + nburn; save.ind <- 0
  b_store    <- matrix(NA, nsave, K)
  o_store    <- matrix(NA, nsave, K)
  bt_store   <- array (NA, c(nsave, N, K))
  sig2_store <- matrix(NA, nsave, N)        # per-period variance
  ot_store   <- matrix(NA, nsave, N)        # outlier scalings o_t (1 if none)
  nu_store   <- rep(NA, nsave)              # t degrees of freedom (t / SVt)
  ypred_store <- matrix(NA, nsave, N)       # posterior predictive draws (in-sample fan)

  ## Initialisation
  b.ols     <- qr.solve(crossprod(X), crossprod(X, y))
  sig2_sc   <- as.numeric(crossprod(y - X %*% b.ols)/(N-K))   # scalar sigma^2
  sig2_draw <- rep(sig2_sc, N)
  lam_t     <- rep(1, N)                                      # t mixing weights
  nu_t <- 10; nu0 <- 10; sd_nu <- 2; acc_nu <- 0             # t dof, prior mean, MH step, acc.
  lat_draw  <- matrix(0, N, K); o_draw <- rep(0, K); bt_draw <- matrix(0, N, K)

  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  for (irep in 1:ntot){
    # Total observation std. dev. = outlier scaling x sqrt(variance)
    norm <- 1/(ot_draw * sqrt(sig2_draw))

    # Constant coefficients and non-centered scales (stacked regression)
    yy <- y*norm
    xx <- cbind(X, X*lat_draw)*norm
    V_po       <- solve(crossprod(xx) + diag(v_pr.inv))
    obseq_po   <- V_po %*% (crossprod(xx, yy) + diag(v_pr.inv) %*% obseq_pr)
    obseq_draw <- obseq_po + t(chol(V_po)) %*% rnorm(KK)
    b_draw <- obseq_draw[1:K]
    o_draw <- obseq_draw[(K+1):KK]
    O_draw <- if (K == 1) o_draw else diag(o_draw)

    # TVP latent states via FFBS (unit-variance random walk, scaled observations)
    yy <- (y - X %*% b_draw)*norm
    xx <- (X %*% O_draw)*norm
    lat_draw <- t(ffbs(t(yy), xx, matrix(1, N, 1), t(matrix(1, K, N)),
                       K, 1, N, matrix(0, K, 1), diag(K)))
    for (tt in 1:N) bt_draw[tt, ] <- b_draw + o_draw*lat_draw[tt, ]

    # Horseshoe prior variances: constant part and NCP scales
    hs_b <- get.hs((b_draw - b_pr)[shrink.slct], lam_b, nu_b, tau_b, zeta_b)
    lam_b <- hs_b$lam
    nu_b <- hs_b$nu
    tau_b <- hs_b$tau
    zeta_b <- hs_b$zeta
    v_pr[shrink.slct] <- hs_b$psi
    
    hs_o <- get.hs((o_draw - o_pr)[shrink.slct], lam_o, nu_o, tau_o, zeta_o)
    lam_o <- hs_o$lam
    nu_o <- hs_o$nu
    tau_o <- hs_o$tau
    zeta_o <- hs_o$zeta
    v_pr[K + shrink.slct] <- hs_o$psi
    v_pr.inv <- 1/v_pr

    # Error variance (residuals de-scaled by the outlier component)
    f_draw   <- rowSums(X * bt_draw)            # fitted conditional mean
    eps_draw <- y - f_draw
    ssr      <- as.numeric((eps_draw/ot_draw)^2)
    nu_cur   <- NA
    if (err.var == "homo"){
      sig2_sc   <- 1/rgamma(1, s_po, S_pr + sum(ssr)/2)
      sig2_draw <- rep(sig2_sc, N)
    } else if (err.var == "t"){
      # t-distributed errors: eps_t ~ N(0, sigma^2 / lam_t), lam_t ~ Gamma(nu/2, nu/2)
      lam_t   <- rgamma(N, shape = (nu_t + 1)/2, rate = (nu_t + ssr/sig2_sc)/2)
      sig2_sc <- 1/rgamma(1, s_po, S_pr + sum(lam_t * ssr)/2)
      # nu via random-walk Metropolis-Hastings
      slog <- sum(log(lam_t)); ssum <- sum(lam_t)
      nu_prop <- nu_t + rnorm(1, 0, sd_nu)
      if (log.post.nu(nu_prop, slog, ssum, N, nu0) -
          log.post.nu(nu_t,    slog, ssum, N, nu0) > log(runif(1))){
        nu_t <- nu_prop; acc_nu <- acc_nu + 1
      }
      sig2_draw <- sig2_sc / lam_t
      nu_cur    <- nu_t
    } else {  # SV / SVt / SVo
      svdraw <- svsample_fast_cpp(eps_draw/ot_draw, startpara = svdraw,
                                  startlatent = sv_latent, priorspec = sv_priors)
      svdraw[c("mu","phi","sigma","nu","rho")] <-
        as.list(svdraw$para[, c("mu","phi","sigma","nu","rho")])
      sv_latent <- svdraw$latent
      sig2_draw <- exp(as.numeric(sv_latent))
      if (err.var == "SVt"){
        sig2_draw <- sig2_draw * as.numeric(svdraw$tau)   # t-scaling (heavy tails)
        nu_cur    <- svdraw$nu
      } else if (err.var == "SVo"){
        loglik <- sapply(ot_grid, function(g) dnorm(eps_draw, 0, g*sqrt(sig2_draw), log = TRUE))
        logpr  <- sweep(loglik, 2, log(ot_p_grid), "+")
        logpr  <- logpr - apply(logpr, 1, max)          # log-sum-exp guard against underflow
        probs  <- exp(logpr); probs <- probs/rowSums(probs)
        for (tt in 1:N) ot_draw[tt] <- sample(ot_grid, 1, prob = probs[tt, ])
        ot_sum    <- sum(ot_draw != 1)
        ot_p      <- rbeta(1, ot_Ba + ot_sum, ot_Bb + N - ot_sum)
        ot_p_grid <- c(1 - ot_p, rep(ot_p/(ot_grid_n - 1), ot_grid_n - 1))
      }
    }

    # Storage
    if (irep %in% save.set){
      save.ind <- save.ind + 1
      b_store[save.ind, ]    <- b_draw
      o_store[save.ind, ]    <- o_draw
      bt_store[save.ind, , ] <- bt_draw
      sig2_store[save.ind, ] <- sig2_draw
      ot_store[save.ind, ]   <- ot_draw
      nu_store[save.ind]     <- nu_cur
      # posterior predictive draw: y_t = x_t' b_t + o_t * sqrt(sig2_t) * z
      ypred_store[save.ind, ] <- f_draw + ot_draw * sqrt(sig2_draw) * rnorm(N)
    }
    setTxtProgressBar(pb, irep)
  }
  close(pb)
  list(b = b_store, o = o_store, bt = bt_store, sig2 = sig2_store,
       ot = ot_store, nu = nu_store, ypred = ypred_store, b.ols = b.ols)
}

#### --------------------------------------------------------------------- ####
#### ------------------    Data and settings    -------------------------- ####
#### --------------------------------------------------------------------- ####
err.set <- c("homo", "t", "SV", "SVt", "SVo")   # error models to loop over
str.smp <- 1973; end.smp <- 2025 + 3/4; h <- 4   # h = cumulative-ahead horizon (quarters); 4 = one-year-ahead GaR (main)
nburn   <- 3000; nsave <- 6000
Y.slct  <- "GDPC1_PCH"

# Growth-at-Risk (Adrian, Boyarchenko & Giannone, 2019): annualised h-quarter
# average GDP growth y_t = 4*mean(g_{t-h+1..t}), paired with the NFCI dated t-h
gar <- read.csv("data/GaR-data.csv", stringsAsFactors = FALSE)
gar <- ts(gar[, c("GDPC1_PCH", "NFCI")], start = c(1970,1), frequency = 4)
gar <- window(gar, start = c(1971,1))

# GaR target: annualised trailing h-quarter average growth
g    <- as.numeric(gar[, "GDPC1_PCH"])
nfci <- as.numeric(gar[, "NFCI"])
n    <- length(g)
y4   <- rep(NA_real_, n)
for (t in h:n) y4[t] <- 4*mean(g[(t-h+1):t])   # annualised trailing h-qtr avg

gar[, "GDPC1_PCH"] <- y4                            # one-year-ahead target
gar[, "NFCI"]      <- c(rep(NA,h), head(nfci, -h))         # predictor dated t-h

# window drops the leading h quarters (NA target / NA predictor)
gar <- window(gar, start = str.smp, end = end.smp)

data <- list(y    = gar[, Y.slct],
             X     = gar[, "NFCI", drop = F],
             dates = time(gar),
             target = "GDPC1_PCH", h = h)


dts <- data$dates

# y-axis limits, fixed by horizon (comparable across error models)
if(h == 1){       lims <- c(-0.5, 0.5); lims.v <- c(0, 10)
} else if(h == 4){ lims <- c(-5, 5);    lims.v <- c(0, 10) }

#### --------------------------------------------------------------------- ####
#### ----------------------- Run MCMC sampler ---------------------------- ####
#### --------------------------------------------------------------------- ####
for (err.var in err.set){            # produce all figures for every error model
fit <- mcmc.sampler(data$y, data$X, err.var = err.var, nburn = nburn, nsave = nsave)

pref  <- paste0(plot.dir, "S04_", err.var, "_TVP_HS_", Y.slct, "_h", h, "_")
bt.lo <- apply(fit$bt[, , 1], 2, quantile, 0.05)
bt.md <- apply(fit$bt[, , 1], 2, quantile, 0.50)
bt.hi <- apply(fit$bt[, , 1], 2, quantile, 0.95)
bt.q25 <- apply(fit$bt[, , 1], 2, quantile, 0.25)
bt.q75 <- apply(fit$bt[, , 1], 2, quantile, 0.75)

#### ----  Figure 1: time-varying NFCI coefficient  --------------------- ####
pdf(paste0(pref, "NFCI.pdf"), width = 8, height = 5)
par(mar = c(3.4, 3.6, 1.2, 0.8), mgp = c(2.2, 0.6, 0))
plot(dts, bt.md, col = rred, lwd = 2.6, type = "l", ylim = lims,
     xlab = "Year", ylab = expression(beta[t]))
polygon(c(dts, rev(dts)), c(bt.lo,  rev(bt.hi)),  col = adjustcolor(rred, 0.15), border = NA)
polygon(c(dts, rev(dts)), c(bt.q25, rev(bt.q75)), col = adjustcolor(rred, 0.30), border = NA)
abline(h = fit$b.ols[1],       col = bblue, lwd = 1.6, lty = 2)
abline(h = 0, lwd = 0.8)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("Posterior median", "50% / 90% credible sets", "OLS"),
       col = c(rred, adjustcolor(rred, 0.4), bblue),
       lwd = c(2.6, 6, 1.6), lty = c(1, 1, 2))
dev.off()

#### --------         Figure 2: stochastic volatility paths      --------- ####

vola <- fit$ot * sqrt(fit$sig2)
v.lo <- apply(vola, 2, quantile, 0.05)
v.md <- apply(vola, 2, quantile, 0.50)
v.hi <- apply(vola, 2, quantile, 0.95)
pdf(paste0(pref, "vola.pdf"), width = 8, height = 5)
par(mar = c(3.4, 3.8, 2.4, 0.8), mgp = c(2.2, 0.6, 0))
plot(dts, v.md, col = bblue, lwd = 2.4, type = "l", ylim = lims.v, xlab = "Year", ylab = "")
polygon(c(dts, rev(dts)), c(v.lo, rev(v.hi)), col = adjustcolor(bblue, 0.20), border = NA)
legend("topleft", bty = "n", cex = 0.8, legend = c("Posterior median", "90% credible set"),
       col = c(bblue, adjustcolor(bblue, 0.4)), lwd = c(2.4, 6))
dev.off()

#### --------        Figure 3: In-sample predictive fan            -------- ####
grid.p <- seq(0.05, 0.95, by = 0.05)                # quantile grid for the fan
qfan   <- t(apply(fit$ypred, 2, quantile, probs = grid.p))   # N x P conditional quantiles
P    <- length(grid.p)
pal  <- colorRampPalette(c(bblue, "white"))(ceiling(P/2) + 1)

pdf(paste0(pref, "fan.pdf"), width = 9, height = 5)
par(mar = c(3.4, 3.8, 2, 0.8), mgp = c(2.4, 0.7, 0))
plot(dts, data$y, type = "n", ylim = range(-12,12),
     xlab = "Year", ylab = expression(y[t]),
     main = "")
for (j in 1:floor(P/2))                              # nested symmetric quantile bands
  polygon(c(dts, rev(dts)), c(qfan[, j], rev(qfan[, P + 1 - j])), col = pal[j], border = NA)
points(dts, data$y, pch = 16, cex = 0.5, col = rred)   # realized growth
abline(h = 0, col = "grey55", lty = 2)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("Realised growth", "90% credible set"),
       pch = c(16, 15), col = c(rred, pal[2]))
box(col = "grey60")
dev.off()
}

