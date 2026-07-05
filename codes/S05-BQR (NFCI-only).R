#### --------------------------------------------------------------------- ####
#### --------                Main estimation script              --------- ####
#### --------     Bayesian quantile regression with HS prior     --------- ####
#### --------------------------------------------------------------------- ####
rm(list = ls())
set.seed(2026)

library(GIGrvg)   # rgig(): draws of the AL specific auxiliary variable 
library(sn)       # qst():  skew-t quantile function for fitting a predictive density

# colors of the slide deck
bblue <- rgb( 32,  66, 133, maxColorValue = 255)
rred  <- rgb(194,  40,  46, maxColorValue = 255)


#### --------------------------------------------------------------------- ####
#### ------------------    Auxiliary functions     ----------------------- ####
#### --------------------------------------------------------------------- ####

# Horseshoe update
get.hs <- function(bdraw, lam.hs, nu.hs, tau.hs, zeta.hs){
  k <- length(bdraw)
  lam.hs <- 1/rgamma(k, shape = 1, rate = 1/nu.hs + bdraw^2/(2*tau.hs))
  nu.hs  <- 1/rgamma(k, shape = 1, rate = 1 + 1/lam.hs)
  tau.hs  <- 1/rgamma(1, shape = (k+1)/2, rate = 1/zeta.hs + sum(bdraw^2/lam.hs)/2)
  zeta.hs <- 1/rgamma(1, shape = 1, rate = 1 + 1/tau.hs)
  list(psi = lam.hs*tau.hs, lam = lam.hs, tau = tau.hs, nu = nu.hs, zeta = zeta.hs)
}

# Skew-t quantiles at the levels tau.out
fit.skewt <- function(q.emp, tau.in, tau.out, start = NULL){
  obj <- function(par){                       # par = (xi, log omega, alpha, log(nu-1))
    qhat <- tryCatch(qst(tau.in, xi = par[1], omega = exp(par[2]),
                         alpha = par[3], nu = 1 + exp(min(par[4], 13.8))),
                     error = function(e) rep(NA_real_, length(tau.in)))
    if (!all(is.finite(qhat))) return(1e10)
    sum((q.emp - qhat)^2)
  }
  fit1 <- function(st){
    op <- optim(st, obj, method = "Nelder-Mead", control = list(maxit = 1000))
    if (op$convergence != 0)                  # hit maxit: continue from incumbent
      op <- optim(op$par, obj, method = "Nelder-Mead", control = list(maxit = 1000))
    op
  }
  op <- fit1(c(median(q.emp), log(diff(range(q.emp))/4 + 1e-6), 0, log(3)))
  if (!is.null(start)){                       # warm start from the previous quarter
    op.w <- fit1(start)
    if (op.w$value < op$value) op <- op.w
  }
  list(q = qst(tau.out, xi = op$par[1], omega = exp(op$par[2]),
               alpha = op$par[3], nu = 1 + exp(min(op$par[4], 13.8))),
       par = op$par, sse = op$value, conv = op$convergence)
}


#### --------------------------------------------------------------------- ####
#### ------------------    Main MCMC functions     ----------------------- ####
#### --------------------------------------------------------------------- ####
# tau = quantile level in (0,1)
# tau = 0.5 reduces to median regression
mcmc.sampler <- function(y, X, tau = 0.5, nburn = 2000, nsave = 3000){

  ## Add the intercept as the last column; K counts all regressors
  X <- cbind(X, "cons" = 1)
  N <- nrow(X)
  K <- ncol(X)
  shrink.slct <- 1:(K-1)  # shrink all coefficients except the intercept

  ## Asymmetric-Laplace scale parameters 
  theta <- (1 - 2*tau)/(tau*(1 - tau))
  psi2  <- 2/(tau*(1 - tau))
  v_draw     <- rep(1, N) # AL auxiliary variable

  ## Prior mean for regression coefficients
  b_pr <- rep(0, K)

  ## HS prior setup
  # Local scalings 
  lam_b <- rep(1, K-1)
  nu_b  <- rep(1, K-1)
  # Global scalings 
  tau_b  <- 1
  zeta_b <- 1

  v_pr <- rep(1, K)
  v_pr[K] <- 1e4                    # intercept, no shrinkage
  v_pr.inv <- 1/v_pr

  ## Inverse-Gamma prior for the error variance
  s_pr <- 3; S_pr <- 0.3
  sig2_draw  <- as.numeric(var(y))
  
  ## Storage
  ntot     <- nburn + nsave
  save.set <- 1:nsave + nburn
  save.ind <- 0

  b_store     <- matrix(NA, nsave, K)
  sig2_store  <- matrix(NA, nsave, 1)
  v_store     <- matrix(NA, nsave, N)
  
  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  for (irep in 1:ntot){
    # Per-observation normalization from the AL mixture
    norm <- as.numeric(1/sqrt(psi2*sig2_draw*v_draw))
    xx   <- X*norm                          # normalized X
    yy   <- (y - theta*v_draw)*norm         # normalized (location-shifted) y

    # Step 1: Sample regression coefficients (weighted ridge with shrinkage prior)
    V_po <- solve(crossprod(xx) + diag(v_pr.inv))
    b_po <- V_po %*% (crossprod(xx, yy) + diag(v_pr.inv) %*% b_pr)
    b_draw <- b_po + t(chol(V_po)) %*% rnorm(K)

    # Step 2: Sample HS
    hs_b <- get.hs(bdraw   = (b_draw - b_pr)[shrink.slct],
                     lam.hs  = lam_b,
                     nu.hs   = nu_b,
                     tau.hs  = tau_b,
                     zeta.hs = zeta_b)
    lam_b  <- hs_b$lam
    nu_b   <- hs_b$nu
    tau_b  <- hs_b$tau
    zeta_b <- hs_b$zeta

    v_pr[shrink.slct] <- hs_b$psi
    v_pr.inv <- 1/v_pr

    # Step 3: Sample the AL auxiliary variable from its GIG full conditional
    res_draw <- as.numeric(y - X %*% b_draw)
    gamma2   <- 2/sig2_draw + theta^2/(psi2*sig2_draw)
    for (tt in 1:N){
      delta2 <- res_draw[tt]^2/(psi2*sig2_draw)
      v_draw[tt] <- rgig(n = 1, lambda = 1/2, chi = delta2, psi = gamma2)
    }

    # Step 4: Sample the error variance from its IG full conditional
    eps_draw  <- res_draw - theta*v_draw
    s_po   <- s_pr + 3*N
    S_po   <- S_pr + 2*sum(v_draw) + sum(eps_draw^2/(psi2*v_draw))
    sig2_draw <- 1/rgamma(1, s_po/2, S_po/2)

    # Storage
    if (irep %in% save.set){
      save.ind <- save.ind + 1
      b_store[save.ind, ]     <- b_draw
      sig2_store[save.ind, ]  <- sig2_draw
      v_store[save.ind, ]     <- v_draw
    }
    setTxtProgressBar(pb, irep)
  }
  close(pb)

  list(b = b_store, sig2 = sig2_store, v = v_store, tau = tau)
}


#### --------------------------------------------------------------------- ####
#### ------------------    Data and settings    -------------------------- ####
#### --------------------------------------------------------------------- ####
tau.grid  <- seq(0.05, 0.95, 0.05)   # quantile levels to estimate (the QR grid)
str.smp <- 1973; end.smp <- 2025 + 3/4; h <- 4   # h = horizon (quarters); 4 = one-year-ahead GaR
nburn   <- 3000; nsave <- 6000
Y.slct  <- "GDPC1_PCH"
fig.dir <- "figs/"; dir.create(fig.dir, showWarnings = FALSE)

# Growth at risk (Adrian, Boyarchenko & Giannone, 2019): annualised h-quarter average GDP growth 
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
dts <- data$dates

#### --------------------------------------------------------------------- ####
#### ------------------    Estimate one QR per quantile level    --------- ####
#### --------------------------------------------------------------------- ####
res.q <- vector("list", length(tau.grid))
names(res.q) <- sprintf("tau%02d", round(tau.grid*100))
for (j in seq_along(tau.grid)){
  message("Estimating quantile tau = ", tau.grid[j])
  res.q[[j]] <- mcmc.sampler(data$y, data$X, tau = tau.grid[j],
                             nburn = nburn, nsave = nsave)
}

# Posterior-median fitted conditional quantiles  q_tau(t) = x_t' b_tau
Xfull <- cbind(as.numeric(data$X), 1)                                 # T x 2 (NFCI, intercept)
qfit  <- sapply(res.q, function(r) as.numeric(Xfull %*% apply(r$b, 2, median)))  # T x P
qfit  <- t(apply(qfit, 1, sort))                                      # enforce monotonicity

#### --------------------------------------------------------------------- ####
#### -------  Figure 1: NFCI coefficient across quantiles  --------------- ####
#### --------------------------------------------------------------------- ####
bq       <- sapply(res.q, function(r) quantile(r$b[, 1], c(0.05, 0.25, 0.5, 0.75, 0.95)))
ols.nfci <- coef(lm(data$y ~ as.numeric(data$X)))[2]    # OLS (mean-regression) slope

pdf(paste0(fig.dir, "S05_GaR_QR_coefficients.pdf"), width = 8, height = 5)
par(mar = c(3.4, 3.6, 1.2, 0.8), mgp = c(2.2, 0.6, 0))
plot(tau.grid, bq[3, ], type = "l", col = rred, lwd = 2.6, ylim = range(bq),
     xlab = expression(tau), ylab = expression(beta[tau]))
polygon(c(tau.grid, rev(tau.grid)), c(bq[1, ], rev(bq[5, ])),   # 90% credible set
        col = adjustcolor(rred, 0.15), border = NA)
polygon(c(tau.grid, rev(tau.grid)), c(bq[2, ], rev(bq[4, ])),   # 50% credible set
        col = adjustcolor(rred, 0.30), border = NA)
abline(h = ols.nfci, col = bblue, lwd = 1.6, lty = 2)      # OLS (constant) slope
abline(h = 0, lwd = 0.8)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("Posterior median", "50% / 90% credible sets", "OLS"),
       col = c(rred, adjustcolor(rred, 0.4), bblue), lwd = c(2.6, 6, 1.6), lty = c(1, 1, 2))
dev.off()

#### --------------------------------------------------------------------- ####
#### -----------  Figure 2: Quantile-specific estimates  ----------------- ####
#### --------------------------------------------------------------------- ####
qs   <- c(0.05, 0.25, 0.50, 0.75, 0.95)              # Quantiles to draw
qlab <- c("5th", "25th", "50th (median)", "75th", "95th")
qcol <- colorRampPalette(c(rred, "grey40", bblue))(length(qs))
xt   <- as.numeric(data$X); yt <- as.numeric(data$y)
xg   <- seq(min(xt), max(xt), length.out = 100)
Xg   <- cbind(xg, 1)                                 # design grid (NFCI, intercept)

fit.med <- fit.lo <- fit.hi <- matrix(NA, length(xg), length(qs))
for (j in seq_along(qs)){
  bdr <- res.q[[sprintf("tau%02d", round(qs[j]*100))]]$b   # nsave x 2 (slope, intercept)
  fdr <- bdr %*% t(Xg)                                     # nsave x ng posterior fits
  fit.med[, j] <- apply(fdr, 2, median)
  fit.lo[, j]  <- apply(fdr, 2, quantile, 0.05)
  fit.hi[, j]  <- apply(fdr, 2, quantile, 0.95)
}

pdf(paste0(fig.dir, "S05_GaR_QR_NFCI.pdf"), width = 8.5, height = 5.4)
par(mar = c(3.6, 3.8, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
plot(xt, yt, type = "n", ylim = c(-15, 10),
     xlab = "NFCI", ylab = "One-year ahead GDP growth")
abline(h = 0, v = 0, lty = 3, col = "grey60")
for (j in seq_along(qs))                             # 90% credible bands (drawn first)
  polygon(c(xg, rev(xg)), c(fit.lo[, j], rev(fit.hi[, j])),
          col = adjustcolor(qcol[j], 0.15), border = NA)
points(xt, yt, pch = 16, col = adjustcolor("grey55", 0.5), cex = 0.7)
for (j in seq_along(qs))                             # posterior-median lines
  lines(xg, fit.med[, j], col = qcol[j], lwd = 2.6)
legend("bottomleft", bty = "n", cex = 0.85, seg.len = 1.4,
       title = expression("Quantile " * tau * " (90% credible set)"),
       legend = rev(qlab), col = rev(qcol), lwd = 2.6)
box(col = "grey60")
dev.off()

#### --------------------------------------------------------------------- ####
#### ----  Figure 3: Predictive fan plot via a skew-t approximation  ----- ####
#### --------------------------------------------------------------------- ####
fan.tau <- seq(0.05, 0.95, 0.05)
qfan   <- matrix(NA_real_, nrow(qfit), length(fan.tau))
sse.st <- rep(NA_real_, nrow(qfit)); st.par <- NULL
for (t in seq_len(nrow(qfit))){
  ft <- fit.skewt(qfit[t, ], tau.grid, fan.tau, start = st.par)
  qfan[t, ] <- ft$q; st.par <- ft$par; sse.st[t] <- ft$sse
}
P     <- length(fan.tau)
pal   <- colorRampPalette(c(bblue, "white"))(ceiling(P/2) + 1)

pdf(paste0(fig.dir, "S05_GaR_QR_fan.pdf"), width = 9, height = 5)
par(mar = c(3.4, 3.8, 2, 0.8), mgp = c(2.4, 0.7, 0))
plot(dts, data$y, type = "n", ylim = range(qfan, data$y),
     xlab = "Year", ylab = expression(y[t]), main = "")
for (j in 1:floor(P/2))                                  # nested symmetric quantile bands
  polygon(c(dts, rev(dts)), c(qfan[, j], rev(qfan[, P + 1 - j])), col = pal[j], border = NA)
points(dts, data$y, pch = 16, cex = 0.5, col = rred)   # realised growth
abline(h = 0, col = "grey55", lty = 2)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("Realised growth", "90% credible set"),
       pch = c(16, 15), col = c(rred, pal[2]))
box(col = "grey60")
dev.off()

