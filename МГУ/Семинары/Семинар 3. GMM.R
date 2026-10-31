setwd("d:\\YandexDisk\\03. Преподавание\\Эконометрика (МГУ)\\")

library(zeallot)
library(ivreg)
library(gmm)

df <- read.csv("Лекции/data/schooling.csv")

lm(LWAGE ~ SCHOOLING + EXP + EXP2 + BLACK + SMSA + SOUTH, data = df) |> summary()

ivreg(LWAGE ~ BLACK + SMSA + SOUTH | SCHOOLING + EXP + EXP2 | AGE + AGE2 + NEARC, data = df) |> summary()

gmm(LWAGE ~ SCHOOLING + EXP + EXP2 + BLACK + SMSA + SOUTH, ~ BLACK + SMSA + SOUTH + AGE + AGE2 + NEARC, data = df, type = "iterative") |> summary()

## ---- Allen Wang -------------------------------------------------------------

N <- 2500
S <- 500

e <- rnorm(N)
y <- numeric(N)

for (t in 3:N) y[t] <- y[t - 1] - 0.21 * y[t - 2] + e[t] + .8 * e[t - 1]
y <- y[-(1:S)]

g <- function(pars, x) {
  N <- length(x)
  c(a1, a2, b) %<-% pars

  c1 <- numeric(N)
  c2 <- numeric(N)
  c3 <- numeric(N)
  c4 <- numeric(N)
  c5 <- numeric(N)

  for (t in 3:N) {
    c1[t] <- x[t] - a1 * x[t - 1] - a2 * x[t - 2]
    c2[t] <- (c1[t])^2 - (1 + b^2)
  }
  for (t in 6:N) {
    c3[t] <- c1[t] * c1[t - 1] - b
    c4[t] <- c1[t] * c1[t - 2]
    c5[t] <- c1[t] * c1[t - 3]
  }
  m <- cbind(c1, c2, c3, c4, c5)
  m[6:N, ]
}

G <- function(pars, x) {
  N <- length(x)
  c(a1, a2, b) %<-% pars
  c1 <- numeric(N)
  for (t in 3:N) c1[t] <- x[t] - a1 * x[t - 1] - a2 * x[t - 2]

  E1 <- numeric(5)
  E2 <- numeric(5)
  E3 <- c(0, -2 * b, -1, 0, 0)

  for (t in 6:N) {
    E1 <- E1 + c(
      -x[t - 1],
      -2 * x[t - 1] * c1[t],
      -x[t - 1] * c1[t - 1] - x[t - 2] * c1[t],
      -x[t - 1] * c1[t - 2] - x[t - 3] * c1[t],
      -x[t - 1] * c1[t - 3] - x[t - 4] * c1[t]
    )
    E2 <- E2 + c(
      -x[t - 2],
      -2 * x[t - 2] * c1[t],
      -x[t - 2] * c1[t - 1] - x[t - 3] * c1[t],
      -x[t - 2] * c1[t - 2] - x[t - 4] * c1[t],
      -x[t - 2] * c1[t - 3] - x[t - 5] * c1[t]
    )
  }
  E1 <- E1 / (N - 5)
  E2 <- E2 / (N - 5)
  cbind(E1, E2, E3)
}

arima(y, order = c(2, 0, 1), include.mean = FALSE)

gmmres <- gmm(g, x = y, t0 = c(a1 = 0, a2 = 0, b = 0), gradv = G, type = "iterative", method = "BFGS")
gmmres
gmmres$vcov |>
  diag() |>
  sqrt()
