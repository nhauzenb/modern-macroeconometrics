###--------------------------------------------------------------------------###
###----------------------------- Data set-up --------------------------------###
###------------ FRED-QD (McCracken & Ng, 2020), stored vintage csv ----------###
###--------------------------------------------------------------------------###
# Options (set in the main script before source(); defaults below):
#   csv_path : stored FRED-QD vintage csv
#   str.smp  : sample start, e.g. 1973
#   end.smp  : sample end,   e.g. 2019 + 3/4
#   h        : forecast horizon in quarters
#   Y.slct   : target mnemonic; y(t) = annualized growth over the past h
#              quarters, regressed on predictors dated t-h (lagged perspective)
#   X.slct   : NULL = all fully observed (balanced) series;
#              or a vector of mnemonics, e.g. c("INDPRO", "UNRATE", "FEDFUNDS")
#   p        : number of predictor lags (1 = period-t information only)

if (!exists("h"))           stop("Set the forecast horizon 'h' (in quarters) before source('fred_setup.R').")
if (!exists("str.smp"))     str.smp     <- 1973
if (!exists("end.smp"))     end.smp     <- 2019 + 3/4
if (!exists("Y.slct"))      Y.slct      <- "GDPC1"
if (!exists("X.slct"))      X.slct      <- NULL
if (!exists("p"))           p           <- 1
if (!exists("standardize")) standardize <- TRUE

ppy <- 4  # Periods per year

# McCracken & Ng transformation function (tcodes 1-7)
transxf <- function(x, code) {
  n <- length(x); y <- rep(NA_real_, n); s <- 1e-6
  if      (code == 1) y <- x
  else if (code == 2) y[2:n] <- x[2:n] - x[1:(n - 1)]
  else if (code == 3) y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
  else if (code == 4 && min(x, na.rm = TRUE) > s) y <- log(x)
  else if (code == 5 && min(x, na.rm = TRUE) > s) {
    lx <- log(x); y[2:n] <- lx[2:n] - lx[1:(n - 1)]
  } else if (code == 6 && min(x, na.rm = TRUE) > s) {
    lx <- log(x); y[3:n] <- lx[3:n] - 2 * lx[2:(n - 1)] + lx[1:(n - 2)]
  } else if (code == 7) {
    y[3:n] <- (x[3:n] / x[2:(n - 1)] - 1) - (x[2:(n - 1)] / x[1:(n - 2)] - 1)
  }
  y
}

###--------------------------------------------------------------------------###
###------------- Read csv file and split metadata from data rows ------------###
###--------------------------------------------------------------------------###
raw   <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE)
key   <- tolower(trimws(raw[[1]]))
trow  <- which(grepl("transform", key))
tcode <- as.numeric(raw[trow, -1]); names(tcode) <- colnames(raw)[-1]

drow  <- which(!is.na(as.Date(raw[[1]], format = "%m/%d/%Y")))   # genuine data rows
fred  <- raw[drow, ]
dates <- as.Date(fred[[1]], format = "%m/%d/%Y")
Xlvl  <- as.matrix(sapply(fred[, -1], as.numeric))

###--------------------------------------------------------------------------###
###------------------------- Quarterly ts objects ---------------------------###
###--------------------------------------------------------------------------###
yq0  <- c(as.numeric(format(dates[1], "%Y")),
          ceiling(as.numeric(format(dates[1], "%m"))/3))
lvl  <- ts(Xlvl, start = yq0, frequency = ppy)                # levels
Yraw <- ts(sapply(seq_along(tcode), function(j) transxf(Xlvl[, j], tcode[j])),
           start = yq0, frequency = ppy)                      # transformed
colnames(Yraw) <- names(tcode)

###--------------------------------------------------------------------------###
###------ Target: Annualised growth over the past h quarters (lagged view) --###
###--------------------------------------------------------------------------###
# y(t) = (100*ppy/h) * (log Y_t - log Y_{t-h}), paired below with x(t-h)
lt  <- log(lvl[, Y.slct]); n <- length(lt); ann <- 100 * ppy / h
ygr <- ts(c(rep(NA, h), ann * (lt[(1 + h):n] - lt[1:(n - h)])),
          start = yq0, frequency = ppy)

# Sample window
Yraw <- window(Yraw, start = str.smp, end = end.smp)

###--------------------------------------------------------------------------###
###---------------------------- Create Y and X ------------------------------###
###--------------------------------------------------------------------------###
if (!is.null(X.slct)) {
  var.slct <- unique(c(Y.slct, X.slct))   # own growth is always included
  Yraw     <- Yraw[, var.slct]
} else {
  Yraw     <- Yraw[, colSums(is.na(Yraw)) == 0]   # balanced panel
  var.slct <- colnames(Yraw)
}
M <- length(var.slct)

# Lag structure: y(t) is paired with predictors dated t-h, ..., t-h-(p-1)
var.lbl <- paste0(var.slct, "(t-", h, ")")
if (p > 1) var.lbl <- c(var.lbl, paste0(rep(var.slct, p - 1),
                                        rep(paste0("(t-", h + 1:(p - 1), ")"), each = M)))
X <- ts(embed(as.matrix(Yraw), p),
        start = str.smp + (p - 1)/ppy, frequency = ppy)
colnames(X) <- var.lbl

# Align: predictors stop h quarters early (they enter with lag h)
X <- window(X,   end   = end.smp - h/ppy)
Y <- window(ygr, start = str.smp + (p - 1 + h)/ppy, end = end.smp)
d <- as.numeric(time(Y))   # dates refer to the target y(t)

# Standardize predictors (not the target)
if (standardize) X <- scale(X)

ok  <- is.finite(Y) & complete.cases(X)
data <- list(y = c(Y)[ok], X = X[ok, , drop = FALSE], dates = d[ok],
            target = Y.slct, h = h)
