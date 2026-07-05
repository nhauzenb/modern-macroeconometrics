#### --------------------------------------------------------------------- ####
#### --------                Main estimation script              --------- ####
#### --------  Bayesian linear regression with shrinkage priors  --------- ####
#### --------------------------------------------------------------------- ####
rm(list = ls())
set.seed(2026)

#### --------------------------------------------------------------------- ####
#### ----    Application inspired by "When is Growth at Risk?"       ----- ####
#### ----  Plagborg-Moller, Reichlin, Ricco, and Hasenzagl (2020)    ----- ####
#### ----  Non-conjugate shrinkage with SSVS and horseshoe priors    ----- ####
#### --------------------------------------------------------------------- ####

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

#### --------------------------------------------------------------------- ####
#### ------------------    Main MCMC functions     ----------------------- ####
#### --------------------------------------------------------------------- ####
mcmc.sampler <- function(y, X, prior = c("HS", "SSVS"), nburn = 5000, nsave = 10000){
  prior <- match.arg(prior)

  # Data: add the intercept as the LAST column; K counts ALL regressors
  X <- cbind(X, "intercept" = 1)
  N <- nrow(X)
  K <- ncol(X)
  shrink.slct <- 1:(K-1)            # shrink all coefficients except the intercept

  tXX <- crossprod(X)
  tXY <- crossprod(X, y)
  dXX <- diag(tXX)

  ###---------------- Prior mean for regression coefficients ----------------###
  b_pr <- rep(0, K)

  ###-------------- Prior variances for regression coefficients -------------###
  # Local scalings (HS)
  lam_b <- rep(1, K-1)
  nu_b  <- rep(1, K-1)
  # Global scalings (HS)
  tau_b  <- 1
  zeta_b <- 1

  # SSVS spike and slab scalings (semiautomatic)
  tau.l_b <- 0.01   # spike
  tau.u_b <- 100    # slab
  c_b <- as.numeric(var(y))/dXX     # proxy for the variance of the OLS coefficient 

  v_pr <- rep(1, K)
  v_pr[K] <- 1e4                    # intercept: essentially kept unshrunk
  v_pr.inv <- 1/v_pr

  ###------------------ Prior moments for error variance --------------------###
  s_pr <- 3    # Prior DoF
  S_pr <- 0.3  # Prior scaling
  s_po <- s_pr + N/2   # Posterior DoF

  ###----------------------------- MCMC setup -------------------------------###
  ntot     <- nburn + nsave
  save.set <- 1:nsave + nburn
  save.ind <- 0

  b_store     <- matrix(NA, nsave, K)
  sig2_store  <- matrix(NA, nsave, 1)
  kappa_store <- delta_store <- matrix(NA, nsave, K-1)
  m_eff_store <- matrix(NA, nsave, 1)

  # Initialization
  sig2_draw  <- as.numeric(var(y))
  delta_draw <- rep(NA, K-1)

  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  for (irep in 1:ntot){
    # Step 1: Sample regression coefficients
    V_po <- solve(tXX/sig2_draw + diag(v_pr.inv))
    b_po <- V_po %*% (tXY/sig2_draw + diag(v_pr.inv) %*% b_pr)
    b_draw <- b_po + t(chol(V_po)) %*% rnorm(K)

    # Step 2: Sample prior variances (all coefficients but the intercept)
    if (prior == "HS"){
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
      v_pr[v_pr < 1e-15] <- 1e-15
    }else if (prior == "SSVS"){
      for (kk in shrink.slct){
        # Bernoulli
        p0.kk <- dnorm(b_draw[kk] - b_pr[kk], 0, sqrt(tau.l_b*c_b[kk])) * 0.5
        p1.kk <- dnorm(b_draw[kk] - b_pr[kk], 0, sqrt(tau.u_b*c_b[kk])) * 0.5
        prob.kk <- p1.kk/(p0.kk + p1.kk)
        if (is.nan(prob.kk)) prob.kk <- 1
        # Assign value to delta
        if (prob.kk > runif(1)) delta.kk <- 1 else delta.kk <- 0
        v_pr[kk] <- c_b[kk]*((1 - delta.kk)*tau.l_b + delta.kk*tau.u_b)
        delta_draw[kk] <- delta.kk
      }
    }
    v_pr.inv <- 1/v_pr

    # Shrinkage factors and effective number of unshrunk coefficients
    kappa_b <- 1/(1 + dXX[shrink.slct]/sig2_draw*v_pr[shrink.slct])
    m_eff   <- sum(1 - kappa_b)

    # Step 3: Sample error variance
    eps_draw <- y - X %*% b_draw
    S_po <- S_pr + crossprod(eps_draw)/2
    sig2_draw <- 1/rgamma(1, s_po, S_po)

    # Storage
    if (irep %in% save.set){
      save.ind <- save.ind + 1
      b_store[save.ind, ]     <- b_draw
      sig2_store[save.ind, ]  <- sig2_draw
      kappa_store[save.ind, ] <- kappa_b
      delta_store[save.ind, ] <- delta_draw
      m_eff_store[save.ind, ] <- m_eff
    }
    setTxtProgressBar(pb, irep)
  }
  close(pb)

  list(b = b_store, sig2 = sig2_store, kappa = kappa_store,
       pip = delta_store, m_eff = m_eff_store)
}


#### --------------------------------------------------------------------- ####
#### --------------------    Applications    ----------------------------- ####
#### --------------------------------------------------------------------- ####
# Use the full balanced FRED-QD panel for h = 1 and h = 4 step ahead
# and a Bayesian linear regression with homoskedastic errors

fig.dir  <- "figs/"; dir.create(fig.dir, showWarnings = FALSE)
csv_path <- "data/2026-04-QD.csv"   # stored FRED-QD vintage
str.smp  <- 1973               # Sample start
end.smp  <- 2019 + 3/4         # Sample end (pre-COVID, as in the paper)
Y.slct   <- "GDPC1"            # Target variable
X.slct   <- NULL               # NULL = all fully observed series; or hand-pick,
                               # e.g. c("INDPRO", "UNRATE", "FEDFUNDS", "GS10", "S&P 500")
p        <- 1                  # No. of lags
res <- list()
for (hh in c(1, 4)){
  h <- hh
  source("fred_setup.R")
  res[[as.character(hh)]] <- list(HS   = mcmc.sampler(data$y, data$X, "HS"),
                                  SSVS = mcmc.sampler(data$y, data$X, "SSVS"),
                                  lbl  = colnames(data$X))
  rm(h)
}

#### --------------------------------------------------------------------- ####
#### ------ Real vs financial classification (by FRED-QD mnemonic) ------- ####
#### --------------------------------------------------------------------- ####

fin.mnem <- c("TB3MS","TB6MS","GS1","GS5","GS10","MORTGAGE30US","AAA","BAA",
              "BAA10YM","MORTG10YRx","TB6M3Mx","GS1TB3Mx","GS10TB3Mx",
              "CPF3MTB3Mx","TB3SMFFM","T5YFFM","AAAFFM","COMPAPFF",
              "M2REAL","BUSLOANSx","NONREVSLx","REALLNx","TOTALSLx",
              "DTCOLNVHFNM","DTCTHFNM","CONSPIx","INVEST","VIXCLSx",
              "NIKKEI225","NASDAQCOM","S&P 500","S&P div yield","S&P PE ratio",
              "EXSZUSx","EXJPUSx","EXUSUKx","EXCAUSx")

#### --------------------------------------------------------------------- ####
#### --------------- Fig (1): SSVS PIPs and HS shrinkage factor ---------- ####
#### --------------------------------------------------------------------- ####

pdf(paste0(fig.dir, "S02_GaR_selection.pdf"), width = 11, height = 5.2)
par(mfrow = c(1, 2), mar = c(4.2, 7.5, 2, 0.8), las = 1, mgp = c(2.6, 0.7, 0))
for (hh in c("1", "4")){
  pip  <- colMeans(res[[hh]]$SSVS$pip)
  onek <- 1 - colMeans(res[[hh]]$HS$kappa)
  lbl  <- res[[hh]]$lbl
  o    <- order(pip, decreasing = TRUE)[20:1]
  mnem <- sub("\\(t(-[0-9]+)?\\)$", "", lbl)      # strip the (t)/(t-1) suffix
  cols <- ifelse(mnem[o] %in% fin.mnem, rred, bblue)
  bp <- barplot(pip[o], horiz = TRUE, names.arg = lbl[o], col = adjustcolor(cols, 0.75),
                border = NA, xlim = c(0, 1.02), cex.names = 0.62,
                xlab = expression("SSVS PIP (bars), HS " * 1 - kappa[j] * " (dots)"),
                main = ifelse(hh == "1", "h = 1 (one-quarter ahead)", "h = 4 (one-year ahead)"))
  abline(v = 0.5, lty = 2, col = "grey55")
  points(onek[o], bp, pch = 21, bg = "white", col = "grey20", lwd = 1.1)
  box(col = "grey60")
}
dev.off()

#### --------------------------------------------------------------------- ####
#### -------- Fig (2): HS posterior median + 90% credible sets ----------- ####
#### --------------------------------------------------------------------- ####

pdf(paste0(fig.dir, "S02_GaR_hs_coefficients.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(4.2, 7.5, 2, 0.8), las = 1, mgp = c(2.6, 0.7, 0))
for (hh in c("1", "4")){
  b   <- res[[hh]]$HS$b[, -ncol(res[[hh]]$HS$b)]   # drop intercept
  lbl <- res[[hh]]$lbl
  bm  <- apply(b, 2, median); bl <- apply(b, 2, quantile, 0.05); bu <- apply(b, 2, quantile, 0.95)
  o   <- order(abs(bm), decreasing = TRUE)[15:1]
  mnem <- sub("\\(t(-[0-9]+)?\\)$", "", lbl)      # strip the (t)/(t-1) suffix
  cols <- ifelse(mnem[o] %in% fin.mnem, rred, bblue)
  m_eff <- mean(rowSums(1 - res[[hh]]$HS$kappa))
  plot(0, type = "n", xlim = range(c(bl[o], bu[o])), ylim = c(0.5, 15.5),
       axes = FALSE, xlab = "Posterior median and 90% credible set", ylab = "",
       main = ifelse(hh == "1", "h = 1 (one-quarter ahead)", "h = 4 (one-year ahead)"))
  abline(v = 0, col = "grey55")
  segments(bl[o], 1:15, bu[o], 1:15, col = adjustcolor(cols, 0.85), lwd = 2.4)
  points(bm[o], 1:15, pch = 19, col = cols, cex = 0.9)
  axis(1); axis(2, at = 1:15, labels = lbl[o], cex.axis = 0.62)
  box(col = "grey60")
  legend("bottomright", bty = "n", cex = 0.8, text.col = "grey25",
         legend = sprintf("m_eff ~ %.0f of K = %d", m_eff, length(lbl)))
}
dev.off()

cat("Figures written to", fig.dir, "\n")
