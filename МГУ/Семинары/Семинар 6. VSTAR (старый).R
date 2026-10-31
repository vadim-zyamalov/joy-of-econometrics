## -----------------------------------------------------------------------------
## https://github.com/vadim-zyamalov/vstar
## -----------------------------------------------------------------------------

setwd("d:/YandexDisk/05. Преподавание/Эконометрика (МГУ)/")
rm(list = ls(all.names = TRUE))
cat("\014")

library(vstar)

load("SEM6.Rdata", verbose = TRUE)

## ---- SETVAR -----------------------------------------------------------------
tsDyn::TVAR.LRtest(y.tar, lag = 1, thDelay = 1, mTh = 1)
summary(tsDyn::TVAR(y.tar, lag = 1, include = "none", mTh = 1, model = "TAR"))

## ---- Линейная модель --------------------------------------------------------
model.lin <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = s,
  dataset = y.lin
)
round(linearity.test(model.lin, J = 1), 4)

## ---- 2 режима ---------------------------------------------------------------
model.2reg <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = s,
  dataset = y.2reg
)

round(linearity.test(model.2reg, J = 1), 4)

result.2reg <- vstar.grid(
  dataset = model.2reg,
  m = 2,
  gamma.limits = c(3, 5),
  points = 200,
  cores = 10
)
result.nls.2reg <- vstar.nls(result.2reg, tol = 1e-6, verbose = TRUE)
summary(result.nls.2reg)

round((coef(result.nls.2reg)$coef - B_2reg) / result.nls.2reg$sd, 4)

round(nonlinearity.test(result.nls.2reg, J = 3), 4)
round(
  serial.correlation.test(result.nls.2reg, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.2reg), 4)

## ---- 2 режима (SE) ----------------------------------------------------------
model.2reg.se <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = list("y1", 1),
  dataset = y.2reg.se
)

round(linearity.test(model.2reg.se, J = 1), 4)

result.2reg.se <- vstar.grid(
  dataset = model.2reg.se,
  m = 2,
  gamma.limits = c(3, 5),
  points = 200,
  cores = 10
)
result.nls.2reg.se <- vstar.nls(result.2reg.se, tol = 1e-6, verbose = TRUE)
summary(result.nls.2reg.se)

round((coef(result.nls.2reg.se)$coef - B_2reg) / result.nls.2reg.se$sd, 4)

round(nonlinearity.test(result.nls.2reg.se, J = 3), 4)
round(
  serial.correlation.test(result.nls.2reg.se, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.2reg.se), 4)

## ---- 2 режима (SE, AR2) ------------------------------------------------------
model.2reg.ar2 <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = list("y1", 1),
  dataset = y.2reg.ar2
)

round(linearity.test(model.2reg.ar2, J = 3), 4)

result.2reg.ar2 <- vstar.grid(
  dataset = model.2reg.ar2,
  m = 2,
  gamma.limits = c(3, 5),
  points = 200,
  cores = 10
)
result.nls.2reg.ar2 <- vstar.nls(result.2reg.ar2, tol = 1e-6, verbose = TRUE)
summary(result.nls.2reg.ar2)

round((coef(result.nls.2reg.ar2)$coef - B_2reg) / result.nls.2reg.ar2$sd, 4)

round(nonlinearity.test(result.nls.2reg.ar2, J = 3), 4)
round(
  serial.correlation.test(result.nls.2reg.ar2, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.2reg.ar2), 4)

## -------------------------
model.2reg.ar2 <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 2,
  trans = list("y1", 1),
  dataset = y.2reg.ar2
)

round(linearity.test(model.2reg.ar2, J = 1), 4)

result.2reg.ar2 <- vstar.grid(
  dataset = model.2reg.ar2,
  m = 2,
  gamma.limits = c(3, 5),
  points = 200,
  cores = 10
)
result.nls.2reg.ar2 <- vstar.nls(result.2reg.ar2, tol = 1e-6, verbose = TRUE)
summary(result.nls.2reg.ar2)

round((coef(result.nls.2reg.ar2)$coef - B2_2reg) / result.nls.2reg.ar2$sd, 4)

round(nonlinearity.test(result.nls.2reg.ar2, J = 3), 4)
round(
  serial.correlation.test(result.nls.2reg.ar2, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.2reg.ar2), 4)

## ---- Два режима (ошибка) ----------------------------------------------------
model.3reg <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = s,
  dataset = y.3reg
)

round(linearity.test(model.3reg, J = 1), 4)

result.3reg <- vstar.grid(
  dataset = model.3reg,
  m = 3,
  gamma.limits = c(3, 5),
  points = 10,
  cores = 10
)
result.nls.3reg <- vstar.nls(result.3reg, tol = 1e-6, verbose = TRUE)
summary(result.nls.3reg)

round((coef(result.nls.3reg)$coef - B_3reg) / result.nls.3reg$sd, 4)

round(nonlinearity.test(result.nls.3reg, J = 1), 4)
round(
  serial.correlation.test(result.nls.3reg, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.3reg), 4)

## ---- Автокорреляция ---------------------------------------------------------
model.acorr <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = s,
  dataset = y.acorr
)

round(linearity.test(model.acorr, J = 1), 4)

result.acorr <- vstar.grid(
  dataset = model.acorr,
  m = 2,
  gamma.limits = c(3, 5),
  points = 200,
  cores = 10
)
result.nls.acorr <- vstar.nls(result.acorr, tol = 1e-6, verbose = TRUE)
summary(result.nls.acorr)

round((coef(result.nls.acorr)$coef - B_2reg) / result.nls.acorr$sd, 4)

round(nonlinearity.test(result.nls.acorr, J = 3), 4)
round(
  serial.correlation.test(result.nls.acorr, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.acorr), 4)

# ---- Нестабильные коэффициенты -----------------------------------------------
model.chcoef <- vstar.prepare(
  endo = c("y1", "y2"),
  p = 1,
  trans = s,
  dataset = y.chcoef
)

round(linearity.test(model.chcoef, J = 1), 4)

result.chcoef <- vstar.grid(
  dataset = model.chcoef,
  m = 2,
  gamma.limits = c(3, 5),
  points = 200,
  cores = 10
)
result.nls.chcoef <- vstar.nls(result.chcoef, tol = 1e-6, verbose = TRUE)
summary(result.nls.chcoef)

round((coef(result.nls.chcoef)$coef - B_2reg) / result.nls.chcoef$sd, 4)

round(nonlinearity.test(result.nls.chcoef, J = 3), 4)
round(
  serial.correlation.test(result.nls.chcoef, J = 1, ortogonalize = TRUE),
  4
)
round(stability.test(result.nls.chcoef), 4)

## -----------------------------------------------------------------------------
G.func <- vstar:::get.G.function()
W <- as.matrix(cbind(1 - G.func(s, 4, median(s)), G.func(s, 4, median(s))))
fff <- sstvars::fitSTVAR(
  y.2reg,
  p = 1,
  M = 2,
  weight_function = "exogenous",
  weightfun_pars = W[2:950, ]
)
