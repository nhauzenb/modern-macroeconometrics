#### --------------------------------------------------------------------- ####
#### --------                Main estimation script              --------- ####
#### ----  Bayesian TVP regression with Horseshoe shrinkage prior    ----- ####
#### ----        Non-centered parameterisation (NCP) + FFBS          ----- ####
#### --------------------------------------------------------------------- ####

rm(list = ls())
set.seed(2026)
plot.dir <- "figs/"; dir.create(plot.dir, showWarnings = FALSE)

#### --------------------------------------------------------------------- ####
#### ----    Application inspired by "When is Growth at Risk?"       ----- ####
#### ----  Plagborg-Moller, Reichlin, Ricco, and Hasenzagl (2020)    ----- ####
#### --------------------------------------------------------------------- ####

library(Rcpp)
library(stochvol)       
sourceCpp("ffbs.cpp")    # forward filtering backward sampling for the RW states

# colors of the slide deck
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

#### --------------------------------------------------------------------- ####
#### ------------------    Main MCMC sampler       ----------------------- ####
#### --------------------------------------------------------------------- ####
mcmc.sampler <- function(y, X, err.var = c("homo", "SV"), nburn = 3000, nsave = 6000){
  err.var <- match.arg(err.var)

  # Data: add the intercept as the LAST column; K counts ALL regressors
  X <- cbind(X, "cons" = 1)
  N  <- nrow(X)
  K  <- ncol(X)
  KK <- 2*K                          # K constant parts + K non-centered scales
  shrink.slct <- 1:(K-1)             # shrink all coefficients except the intercept

  ###---------------- Prior mean for regression coefficients ----------------###
  b_pr <- rep(0, K)                  # standardized predictors -> shrink to zero
  o_pr <- rep(0, K)                  # NCP scales shrunk to zero (constant coef.)
  obseq_pr <- c(b_pr, o_pr)

  ###-------------- Prior variances for regression coefficients -------------###
  # Local and global horseshoe scalings for the constant part b ...
  lam_b <- rep(1, K-1); nu_b <- rep(1, K-1); tau_b <- 0.5; zeta_b <- 0.5
  # ... and for the non-centered TVP scales sqrt(omega)
  lam_o <- rep(1, K-1); nu_o <- rep(1, K-1); tau_o <- 0.5; zeta_o <- 0.5

  v_pr     <- rep(1, KK)
  v_pr[K]  <- 1e4                    # constant intercept:  essentially unshrunk
  v_pr[KK] <- 1e-4                   # TVP intercept:       heavily shrunk (constant mean)
  v_pr.inv <- 1/v_pr

  ###------------------ Prior moments for error variance --------------------###
  s_pr <- 3    # Prior DoF
  S_pr <- 0.3  # Prior scaling
  s_po <- s_pr + N/2

  # Stochastic-volatility priors (used only if err.var == "SV")
  sv_priors <- specify_priors(
    mu = sv_normal(0, 1), phi = sv_beta(25, 1.5),
    sigma2 = sv_gamma(0.5, 10), nu = sv_infinity(), rho = sv_constant(0))
  svdraw    <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0,
                    beta = NA, latent0 = 0)
  sv_latent <- rep(0, N)

  ###----------------------------- MCMC setup -------------------------------###
  ntot     <- nburn + nsave
  save.set <- 1:nsave + nburn
  save.ind <- 0

  b_store    <- matrix(NA, nsave, K)            # constant part
  o_store    <- matrix(NA, nsave, K)            # sqrt(omega) (NCP scales)
  bt_store   <- array (NA, c(nsave, N, K))      # time-varying coefficients
  sig2_store <- matrix(NA, nsave, N)

  # Initialization
  b.ols     <- qr.solve(crossprod(X), crossprod(X, y))
  sig2_draw <- rep(as.numeric(crossprod(y - X %*% b.ols)/(N-K)), N)
  lat_draw  <- matrix(0, N, K)
  o_draw    <- rep(0, K)
  bt_draw   <- matrix(0, N, K)

  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  for (irep in 1:ntot){
    # Sample constant coefficients and non-centered scales (stacked regression)
    yy <- y/sqrt(sig2_draw)
    xx <- cbind(X, X*lat_draw)/sqrt(sig2_draw)
    V_po       <- solve(crossprod(xx) + diag(v_pr.inv))
    obseq_po   <- V_po %*% (crossprod(xx, yy) + diag(v_pr.inv) %*% obseq_pr)
    obseq_draw <- obseq_po + t(chol(V_po)) %*% rnorm(KK)
    b_draw <- obseq_draw[1:K]
    o_draw <- obseq_draw[(K+1):KK]
    O_draw <- if (K == 1) o_draw else diag(o_draw)
    for (tt in 1:N) bt_draw[tt, ] <- b_draw + o_draw*lat_draw[tt, ]

    # Sample TVP latent states with FFBS (unit-variance random walk, scaled obs.)
    yy <- (y - X %*% b_draw)/sqrt(sig2_draw)
    xx <- (X %*% O_draw)/sqrt(sig2_draw)
    lat_draw <- t(ffbs(t(yy), xx, matrix(1, N, 1), t(matrix(1, K, N)),
                       K, 1, N, matrix(0, K, 1), diag(K)))
    for (tt in 1:N) bt_draw[tt, ] <- b_draw + o_draw*lat_draw[tt, ]

    # Sample prior variances of the constant part (all coefficients but the intercept)
    hs_b <- get.hs(bdraw   = (b_draw - b_pr)[shrink.slct],
                   lam.hs  = lam_b,
                   nu.hs   = nu_b,
                   tau.hs  = tau_b,
                   zeta.hs = zeta_b)
    lam_b <- hs_b$lam; nu_b <- hs_b$nu; tau_b <- hs_b$tau; zeta_b <- hs_b$zeta
    v_pr[shrink.slct] <- hs_b$psi

    # Sample prior variances of the NCP scales (all coefficients but the intercept)
    hs_o <- get.hs(bdraw   = (o_draw - o_pr)[shrink.slct],
                   lam.hs  = lam_o,
                   nu.hs   = nu_o,
                   tau.hs  = tau_o,
                   zeta.hs = zeta_o)
    lam_o <- hs_o$lam; nu_o <- hs_o$nu; tau_o <- hs_o$tau; zeta_o <- hs_o$zeta
    v_pr[K + shrink.slct] <- hs_o$psi
    v_pr.inv <- 1/v_pr

    # Sample error variance
    f_draw   <- rowSums(X * bt_draw)
    eps_draw <- y - f_draw
    if (err.var == "homo"){
      S_po      <- S_pr + crossprod(eps_draw)/2
      sig2_draw <- rep(1/rgamma(1, s_po, S_po), N)
    }else if (err.var == "SV"){
      svdraw <- svsample_fast_cpp(eps_draw, startpara = svdraw,
                                  startlatent = sv_latent, priorspec = sv_priors)
      svdraw[c("mu","phi","sigma","nu","rho")] <-
        as.list(svdraw$para[, c("mu","phi","sigma","nu","rho")])
      sv_latent <- pmax(t(svdraw$latent), -8)   # floor on the log-volatility (numerical safeguard)
      sig2_draw <- exp(as.numeric(sv_latent))
    }

    # Storage
    if (irep %in% save.set){
      save.ind <- save.ind + 1
      b_store[save.ind, ]    <- b_draw
      o_store[save.ind, ]    <- o_draw
      bt_store[save.ind, , ] <- bt_draw
      sig2_store[save.ind, ] <- sig2_draw
    }
    setTxtProgressBar(pb, irep)
  }
  close(pb)

  list(b = b_store, o = o_store, bt = bt_store, sig2 = sig2_store)
}

#### --------------------------------------------------------------------- ####
#### ------------------    Application    -------------------------------- ####
#### --------------------------------------------------------------------- ####
# GaR-style FRED-QD panel: one-quarter-ahead GDP growth on a focused macro +
# financial predictor block
csv_path <- "data/2026-04-QD.csv"   # stored FRED-QD vintage
str.smp  <- 1973               # Sample start
end.smp  <- 2025 + 3/4         # Sample end
h        <- 1                  # Forecast horizon (quarters)
Y.slct   <- "GDPC1"            # Target variable
p        <- 1                  # Predictor lags: 1 = period t-h information only

X.slct <- c(# --- real activity / nominal (blue) 
            "INDPRO",      # INDPRO     industrial production
            "UNRATE",      # UNRATE     unemployment rate
            "PAYEMS",      # EMPL       nonfarm employment
            "DPIC96",      # DISPINC    real disposable income
            "INVCQRMTSPL", # INVENTO    change in private inventories (top mean predictor)
            "GPDIC1",      # INVESTM    gross private domestic investment
            "IMPGSC1",     # IMPORT     real imports
            "HOUST",       # HOUSESTART housing starts
            "PERMIT",      # HOUSEPERMIT new housing permits
            "CUMFNS",      # CAPUTIL    capacity utilisation
            "UMCSENTx",    # CONSSENT   consumer sentiment
            "PCECTPI",     # PCEPRICE   PCE price index
            # --- financial block (red) 
            "FEDFUNDS",    # FEDFUNDS   federal funds rate
            "GS10TB3Mx",   # TERMSPR    term spread (10y - 3m)
            "BAA10YM",     # BAASPR     BAA corporate - 10y Treasury spread
            "COMPAPFF",    # CPAPERSPR  commercial paper - fed funds spread
            "S&P 500",     # STOCKPRICE S&P 500
            "S&P div yield",# DIVYIELD  S&P 500 dividend yield
            "VIXCLSx")     # VXO        (implied) stock-market volatility
source("fred_setup.R")


fit <- mcmc.sampler(data$y, data$X, err.var = "SV", nburn = 3000, nsave = 6000)

lbl  <- c(colnames(data$X), "cons")
dts  <- data$dates
N    <- length(dts); K <- length(lbl)

## Real vs financial classification (by FRED-QD mnemonic)
fin.mnem <- c("TB3MS","TB6MS","GS1","GS5","GS10","MORTGAGE30US","AAA","BAA",
              "BAA10YM","MORTG10YRx","TB6M3Mx","GS1TB3Mx","GS10TB3Mx",
              "CPF3MTB3Mx","TB3SMFFM","T5YFFM","AAAFFM","COMPAPFF","M2REAL",
              "BUSLOANSx","NONREVSLx","REALLNx","TOTALSLx","DTCOLNVHFNM",
              "DTCTHFNM","CONSPIx","INVEST","VIXCLSx","NIKKEI225","NASDAQCOM",
              "S&P 500","S&P div yield","S&P PE ratio",
              "EXSZUSx","EXJPUSx","EXUSUKx","EXCAUSx")
mnem  <- sub("\\(t(-[0-9]+)?\\)$", "", lbl)      # strip (t-1) suffix
is.fin <- mnem %in% fin.mnem
col.k  <- ifelse(is.fin, rred, bblue); col.k[K] <- "black"   # intercept black

## Posterior summaries of the TVP paths
bt.lo <- apply(fit$bt, c(2,3), quantile, 0.05)
bt.md <- apply(fit$bt, c(2,3), quantile, 0.50)
bt.hi <- apply(fit$bt, c(2,3), quantile, 0.95)
tv    <- apply(bt.md, 2, sd)                     # time variation per coefficient

#### --------------------------------------------------------------------- ####
#### ----  Figure 1: time-varying coefficients (median + 90% band)  ------ ####
#### --------------------------------------------------------------------- ####
pdf(paste0(plot.dir, "S03_TVP_HS_", Y.slct, "_beta.pdf"), width = 11, height = 8)
par(mfrow = c(ceiling(K/3), 3), mar = c(2.4, 3.6, 2, 0.8), mgp = c(2.2, 0.6, 0))
for (kk in 1:K){
  ttl <- if (kk == K) expression(alpha[t]~"(intercept)") else mnem[kk]
  plot(dts, bt.md[, kk], type = "n", ylim = range(c(bt.lo[, kk], bt.hi[, kk])),
       xlab = "", ylab = "", main = ttl, col.main = col.k[kk])
  polygon(c(dts, rev(dts)), c(bt.lo[, kk], rev(bt.hi[, kk])),
          col = adjustcolor(col.k[kk], 0.20), border = NA)
  lines(dts, bt.md[, kk], col = col.k[kk], lwd = 2)
  abline(h = 0, col = "black", lwd = 0.8)
}
dev.off()

#### --------------------------------------------------------------------- ####
#### ----           Figure 2: TVP scales sqrt(omega) (NCP)            ---- ####
#### --------------------------------------------------------------------- ####
pdf(paste0(plot.dir, "S03_TVP_HS_", Y.slct, "_omega.pdf"), width = 11, height = 8)
par(mfrow = c(ceiling(K/3), 3), mar = c(2.4, 3.6, 2, 0.8), mgp = c(2.2, 0.6, 0))
for (kk in 1:K){
  ttl <- bquote(sqrt(omega)[.(kk)] ~ .(if (kk==K) "(intercept)" else mnem[kk]))
  dd  <- density(c(fit$o[, kk], -fit$o[, kk]))
  plot(dd, main = ttl, col.main = col.k[kk], xlab = "", ylab = "", yaxt = "n")
  polygon(dd, col = adjustcolor(col.k[kk], 0.25), border = col.k[kk])
  abline(v = 0, col = "black", lwd = 0.8)
}
dev.off()

#### --------------------------------------------------------------------- ####
#### ----  Figure 3: the single most time-varying coefficient  ----------- ####
#### --------------------------------------------------------------------- ####
star     <- which.max(tv[1:(K-1)])          # most time-varying predictor (excl. intercept)
col.star <- col.k[star]
bt.q25 <- apply(fit$bt[, , star], 2, quantile, 0.25)
bt.q75 <- apply(fit$bt[, , star], 2, quantile, 0.75)
pdf(paste0(plot.dir, "S03_TVP_HS_", Y.slct, "_mostTVP.pdf"), width = 9, height = 4.8)
par(mar = c(3.4, 3.6, 3, 0.8), mgp = c(2.2, 0.6, 0))
plot(dts, bt.md[, star], type = "n", ylim = range(c(bt.lo[, star], bt.hi[, star])),
     xlab = "Year", ylab = expression(beta[t]),
     main = paste0("Time-varying coefficient of GDP growth on ", mnem[star],
                   "\n(most time-varying coefficient)"))
polygon(c(dts, rev(dts)), c(bt.lo[, star], rev(bt.hi[, star])),
        col = adjustcolor(col.star, 0.15), border = NA)
polygon(c(dts, rev(dts)), c(bt.q25, rev(bt.q75)),
        col = adjustcolor(col.star, 0.30), border = NA)
lines(dts, bt.md[, star], col = col.star, lwd = 2.6)
abline(h = median(fit$b[, star]), col = "black", lwd = 1.6, lty = 2)
abline(h = 0, col = "black", lwd = 0.8)
legend("topright", bty = "n", cex = 0.8,
       legend = c("Posterior median", "50% / 90% credible sets"),
       col = c(col.star, adjustcolor(col.star, 0.4)),
       lwd = c(2.6, 6), lty = c(1, 1))
dev.off()

#### --------------------------------------------------------------------- ####
#### ---- Figure 4: Stochastic volatility path exp(h_t/2)            ----- ####
#### --------------------------------------------------------------------- ####
if (ncol(fit$sig2) > 1 && length(unique(fit$sig2[1, ])) > 1){
  sd.t  <- sqrt(fit$sig2)
  s.lo  <- apply(sd.t, 2, quantile, 0.05)
  s.md  <- apply(sd.t, 2, quantile, 0.50)
  s.hi  <- apply(sd.t, 2, quantile, 0.95)
  pdf(paste0(plot.dir, "S03_TVP_HS_", Y.slct, "_SV.pdf"), width = 9, height = 4.4)
  par(mar = c(3.4, 3.8, 2.4, 0.8), mgp = c(2.2, 0.6, 0))
  plot(dts, s.md, type = "n", ylim = range(c(s.lo, s.hi)),
       xlab = "Year", ylab = "",
       main = expression("Stochastic volatility evolution  " *
                         sigma[t] == exp(h[t]/2)))
  polygon(c(dts, rev(dts)), c(s.lo, rev(s.hi)),
          col = adjustcolor(bblue, 0.20), border = NA)
  lines(dts, s.md, col = bblue, lwd = 2.4)
  legend("topleft", bty = "n", cex = 0.8,
         legend = c("Posterior median", "90% credible set"),
         col = c(bblue, adjustcolor(bblue, 0.4)), lwd = c(2.4, 6))
  dev.off()
}

cat("Most time-varying coefficient:", mnem[star], "\n")
cat("Figures written to", plot.dir, "\n")
