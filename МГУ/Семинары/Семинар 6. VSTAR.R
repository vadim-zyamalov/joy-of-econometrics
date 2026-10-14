## -----------------------------------------------------------------------------
## https://github.com/vadim-zyamalov/vstar
##
## A1L <- .5 * diag(2)
## A2L <- .06 * diag(2)
## A1H <- matrix(c(.25, -.05, 0, -.45), ncol = 2)
## A2H <- matrix(c(.06, 0, 0, 0.06), ncol = 2)
## aA1L <- -.5 * diag(2)
## Omg <- diag(c(0.1, 0.1))
##
## B_2reg <- cbind(t(A1L), t(A1H))
## B_2reg_alt <- cbind(t(A1L), t(aA1L))
## B_3reg <- cbind(t(A1L), t(A1H), t(aA1L))
## B2_2reg <- cbind(t(cbind(A1L, A2L)), t(cbind(A1H, A2H)))
##
## -----------------------------------------------------------------------------

setwd("d:/YandexDisk/03. Преподавание/Эконометрика (МГУ)/Семинары/")
load("SEM6.Rdata", verbose = TRUE)

## ---- SETVAR -----------------------------------------------------------------
tsDyn::TVAR.LRtest(y.tar, lag = 1, thDelay = 1, mTh = 1)
tsDyn::TVAR(y.tar, lag = 1, include = "none", mTh = 1, model = "TAR") |> summary()

## ---- VSTAR ------------------------------------------------------------------
stvalues <- starvars::startingVLSTAR(y.2reg,
  p = 1, n.combi = 3,
  singlecgamma = TRUE, st = s, ncores = 6
)
res <- starvars::VLSTAR(
  y.2reg,
  p = 1,
  m = 2,
  st = s,
  starting = stvalues,
  ncores = 6,
  method = "NLS",
  singlecgamma = TRUE,
  n.iter = 100
)

stvalues <- starvars::startingVLSTAR(y.2reg.ar2,
  p = 2, n.combi = 3,
  singlecgamma = TRUE, st = s, ncores = 5
)
res <- starvars::VLSTAR(
  y.2reg.ar2,
  p = 2,
  m = 2,
  st = s,
  starting = stvalues,
  ncores = 5,
  method = "NLS",
  singlecgamma = TRUE
)
