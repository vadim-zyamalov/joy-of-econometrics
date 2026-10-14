setwd("d:/YandexDisk/05. Преподавание/Эконометрика (МГУ)/")
rm(list = ls(all.names = TRUE))
cat("\014")

library(vstar)

N <- 1000

s <- rnorm(N)

thr_2reg <- median(s)
thr_3reg <- quantile(s, prob = c(.33, .67))

print(thr_2reg)
print(thr_3reg)

G.func <- vstar:::get.G.function("L")

sm <- s %x% t(vstar:::unity(1))
gm <- vstar:::unity(N) %x% t(4)
cm <- vstar:::unity(N) %x% t(thr_2reg)

G.mat_2reg <- G.func(sm, gm, cm)

sm <- s %x% t(vstar:::unity(2))
gm <- vstar:::unity(N) %x% t(c(4, 4))
cm <- vstar:::unity(N) %x% t(thr_3reg)

G.mat_3reg <- G.func(sm, gm, cm)

A1L <- .5 * diag(2)
A2L <- .06 * diag(2)
A1H <- matrix(c(.25, -.05, 0, -.45), ncol = 2)
A2H <- matrix(c(.06, 0, 0, 0.06), ncol = 2)
aA1L <- -.5 * diag(2)
Omg <- diag(c(0.1, 0.1))

B_2reg <- cbind(t(A1L), t(A1H))
B_2reg_alt <- cbind(t(A1L), t(aA1L))
B_3reg <- cbind(t(A1L), t(A1H), t(aA1L))
B2_2reg <- cbind(t(cbind(A1L, A2L)), t(cbind(A1H, A2H)))

# ---- SETVAR ------------------------------------------------------------------
y.tar <- matrix(0, ncol = 2, nrow = 1001)
y.tar[1, ] <- rnorm(2)

for (k in 1:N) {
  y.lag <- y.tar[k, , drop = FALSE]
  d <- ifelse(y.lag[1, 1] >= 0, 1, 0)
  PSI <- t(t(c(d, 1 - d))) %x% diag(2)
  y.tar[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) + Omg %*% rnorm(2))
}
y.tar <- y.tar[52:(N + 1), 1:2, drop = FALSE]
colnames(y.tar) <- c("y1", "y2")

plot.ts(y.tar)

# ---- Linear ------------------------------------------------------------------
y.lin <- matrix(0, ncol = 2, nrow = 1001)
y.lin[1, ] <- rnorm(2)

for (k in 1:N) {
  y.lag <- y.lin[k, , drop = FALSE]
  y.lin[k + 1, ] <- t(A1L %*% t(y.lag) + rnorm(2))
}

y.lin <- y.lin[52:(N + 1), 1:2]
colnames(y.lin) <- c("y1", "y2")

plot.ts(y.lin)

# ---- 2 режима ----------------------------------------------------------------
y.2reg <- matrix(0, ncol = 2, nrow = 1001)
y.2reg[1, ] <- rnorm(2)

for (k in 1:N) {
  y.lag <- y.2reg[k, , drop = FALSE]
  PSI <- t(t(c(1, G.mat_2reg[k, ]))) %x% diag(2)
  y.2reg[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) + Omg %*% rnorm(2))
}
y.2reg <- y.2reg[52:(N + 1), 1:2, drop = FALSE]
colnames(y.2reg) <- c("y1", "y2")

plot.ts(y.2reg)

# ---- 2 режима (SE) -----------------------------------------------------------
y.2reg.se <- matrix(0, ncol = 2, nrow = 1001)
y.2reg.se[1, ] <- rnorm(2)

for (k in 1:N) {
  y.lag <- y.2reg.se[k, , drop = FALSE]
  PSI <- t(t(c(1, G.func(y.lag[1, 1], 4, 0)))) %x% diag(2)
  y.2reg.se[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) + Omg %*% rnorm(2))
}
y.2reg.se <- y.2reg.se[52:(N + 1), 1:2, drop = FALSE]
colnames(y.2reg.se) <- c("y1", "y2")

plot.ts(y.2reg.se)

# ---- 2 режима (SE, AR2) ------------------------------------------------------
y.2reg.ar2 <- matrix(0, ncol = 2, nrow = 1002)
y.2reg.ar2[1:2, ] <- rnorm(4)

for (k in 1:N) {
  y.lag <- cbind(y.2reg.ar2[k + 1, , drop = FALSE], y.2reg.ar2[k, , drop = FALSE])
  PSI <- t(t(c(1, G.func(y.lag[1, 1], 4, 0)))) %x% diag(2)
  y.2reg.ar2[k + 2, ] <- t(t(PSI) %*% t(B2_2reg) %*% t(y.lag) + Omg %*% rnorm(2))
}
y.2reg.ar2 <- y.2reg.ar2[53:(N + 2), 1:2, drop = FALSE]
colnames(y.2reg.ar2) <- c("y1", "y2")

plot.ts(y.2reg.ar2)

# ---- 3 режима ----------------------------------------------------------------
y.3reg <- matrix(0, ncol = 2, nrow = 1001)
y.3reg[1, ] <- rnorm(2)

for (k in 1:N) {
  y.lag <- y.3reg[k, , drop = FALSE]
  PSI <- t(t(c(1, G.mat_3reg[k, ]))) %x% diag(2)
  y.3reg[k + 1, ] <- t(t(PSI) %*% t(B_3reg) %*% t(y.lag) + Omg %*% rnorm(2))
}
y.3reg <- y.3reg[52:(N + 1), 1:2, drop = FALSE]
colnames(y.3reg) <- c("y1", "y2")

plot.ts(y.3reg)

## ----- Автокорреляция ошибок -------------------------------------------------
y.acorr <- matrix(0, ncol = 2, nrow = 1001)
y.acorr[1, ] <- rnorm(2)

e <- matrix(rnorm(2002), ncol = 2)

for (k in 1:N) {
  y.lag <- y.acorr[k, , drop = FALSE]
  e[k + 1, ] <- e[k + 1, ] + .8 * e[k, ]
  PSI <- t(t(c(1, G.mat_2reg[k, ]))) %x% diag(2)
  y.acorr[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) +
                          Omg %*% t(e[k + 1, , drop = FALSE]))
}
y.acorr <- y.acorr[52:(N + 1), 1:2, drop = FALSE]
colnames(y.acorr) <- c("y1", "y2")

plot.ts(y.acorr)

# ---- Нестабильные коэффициенты -----------------------------------------------
y.chcoef <- matrix(0, ncol = 2, nrow = 1001)
y.chcoef[1, ] <- rnorm(2)

for (k in 1:N) {
  y.lag <- y.chcoef[k, , drop = FALSE]
  PSI <- t(t(c(1, G.mat_2reg[k, ]))) %x% diag(2)
  if (k <= 525) {
    y.chcoef[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) +
                             Omg %*% rnorm(2))
  } else {
    y.chcoef[k + 1, ] <- t(t(PSI) %*% t(B_2reg_alt) %*% t(y.lag) +
                             Omg %*% rnorm(2))
  }
}
y.chcoef <- y.chcoef[52:(N + 1), 1:2, drop = FALSE]
colnames(y.chcoef) <- c("y1", "y2")

plot.ts(y.chcoef)

# ------------------------------------------------------------------------------
s <- s[51:N]

save(list = ls()[grep("^(B\\d?_.*|s|y\\.(?!lag).*)$", ls(), perl = TRUE)], file = "SEM4.Rdata")
