# PHASE 4-5 - Estimation on real data with HOLD-OUT, comparison between
# 3 models (M1 linear, M2 AR1, M3 AR1 + seasonality + Student-t),
# DIC, PPC, and FORECAST out-of-sample.


suppressPackageStartupMessages({
  library(rjags); library(coda); library(dplyr); library(lubridate)
})
load.module("dic")

# setwd() not needed: open progetto_fsl2.Rproj so the working directory is the project root
set.seed(20260607)

# ------------- dataset ---------------------------------------
dati <- readRDS("dataset_finale.rds") %>% arrange(Data)
cat("Final dataset:", nrow(dati), "days\n",
    "PM10 station: Roma Preneste (L.GO PERESTRELLO)\n",
    "Weather station: AL007 Boncompagni (Roma centro, 4.5 km da Preneste)\n")

# --- hold-out: the last 90 days ---
N_HOLDOUT <- 90
N_TOT     <- nrow(dati)
idx_tr    <- seq_len(N_TOT - N_HOLDOUT)
idx_te    <- (N_TOT - N_HOLDOUT + 1):N_TOT

dati_tr <- dati[idx_tr, ]
dati_te <- dati[idx_te, ]
cat("Training:", nrow(dati_tr), "days |",
    "Hold-out:", nrow(dati_te), "days (",
    as.character(min(dati_te$Data)), "->",
    as.character(max(dati_te$Data)), ")\n\n")

# ------------- features --------------------------------------
build_X <- function(dat, X_mean = NULL, X_sd = NULL) {
  X_raw <- as.matrix(dat[, c("Temp", "Umidita", "Vento", "Pioggia")])
  if (is.null(X_mean)) X_mean <- colMeans(X_raw)
  if (is.null(X_sd))   X_sd   <- apply(X_raw, 2, sd)
  list(X = scale(X_raw, center = X_mean, scale = X_sd),
       mean = X_mean, sd = X_sd)
}
prep <- build_X(dati_tr)
X_tr <- prep$X
Y_tr <- log(dati_tr$PM10)
N    <- length(Y_tr); K <- ncol(X_tr)

#seasonal terms (sin/cos of the day-of-year)
doy_tr  <- yday(dati_tr$Data)
sin1_tr <- sin(2 * pi * doy_tr / 365.25)
cos1_tr <- cos(2 * pi * doy_tr / 365.25)

#-----------------------------------------------
# fit JAGS + summary

fit_jags <- function(file, data_jags, inits, params,
                     n_burn = 3000, n_iter = 10000, thin = 2,
                     n_chains = 3, dic_iter = 10000) {
  m <- jags.model(file, data = data_jags, inits = inits,
                  n.chains = n_chains, quiet = TRUE)
  update(m, n_burn, progress.bar = "none")
  samp <- coda.samples(m, variable.names = params,
                       n.iter = n_iter, thin = thin,
                       progress.bar = "none")
  dic  <- dic.samples(m, n.iter = dic_iter, type = "pD",
                      progress.bar = "none")
  list(model = m, samples = samp, dic = dic)
}

# =================================================================
# M1: linear
# =================================================================
cat("M1 linear \n")
data1 <- list(N = N, K = K, Y = Y_tr, X = X_tr)
init1 <- function() list(alpha = 0, beta = rep(0, K), tau = 1)
fit1 <- fit_jags("model1_linear.bug", data1, init1,
                 c("alpha","beta","sigma"))

# =================================================================
# M2: AR(1) + weather
# =================================================================
cat("M2 AR(1) \n")
data2 <- list(N = N, K = K, Y = Y_tr, X = X_tr)
init2 <- function() list(alpha = 0, beta = rep(0, K),
                         phi = 0, tau = 1)
fit2 <- fit_jags("model2_ar1.bug", data2, init2,
                 c("alpha","beta","phi","sigma"))

# =================================================================
# M3: AR(1) + seasonality + Student-t
# =================================================================
cat("M3 AR(1) + seasonality + Student-t \n")
data3 <- list(N = N, K = K, Y = Y_tr, X = X_tr,
              sin1 = sin1_tr, cos1 = cos1_tr)
init3 <- function() list(alpha = 0, beta = rep(0, K),
                         gamma1 = 0, gamma2 = 0,
                         phi = 0, tau = 1, nu = 10)
fit3 <- fit_jags("model3_ar1_season_t.bug", data3, init3,
                 c("alpha","beta","gamma1","gamma2",
                   "phi","sigma","nu"))

# posterior summary + R-hat
riassumi <- function(fit, label) {
  cat("\n---", label, "---\n")
  s <- summary(fit$samples)
  cat("Mean / SD:\n");           print(round(s$statistics[,1:2], 3))
  cat("CI 95% (2.5%, 97.5%):\n");print(round(s$quantiles[, c(1, 5)], 3))
  gr <- gelman.diag(fit$samples, multivariate = FALSE)$psrf
  cat("R-hat (max):", round(max(gr[, 1]), 3), "\n")
  cat("DIC:\n");                  print(fit$dic)
}
riassumi(fit1, "M1 lineare")
riassumi(fit2, "M2 AR(1)")
riassumi(fit3, "M3 AR(1)+stag+t")

# traceplots 
saveplot <- function(samps, vars, file) {
  png(file, width = 1200, height = 200 * ceiling(length(vars)/2), res = 130)
  op <- par(mfrow = c(ceiling(length(vars)/2), 2), mar = c(3,3,2,1))
  for (v in vars) plot(samps[, v], main = v)
  par(op); dev.off()
}
saveplot(fit1$samples, c("alpha","beta[1]","beta[2]","beta[3]","beta[4]","sigma"),
         "trace_model1.png")
saveplot(fit2$samples, c("alpha","phi","beta[1]","beta[2]","beta[3]","beta[4]","sigma"),
         "trace_model2.png")
saveplot(fit3$samples, c("alpha","phi","beta[1]","beta[2]","beta[3]","beta[4]",
                         "gamma1","gamma2","sigma","nu"),
         "trace_model3.png")

#DIC table 
dic_val <- function(d) sum(d$deviance) + sum(d$penalty)
dic_tab <- data.frame(
  modello = c("M1 lineare", "M2 AR(1)", "M3 AR(1)+stag+t"),
  D_mean  = c(sum(fit1$dic$deviance), sum(fit2$dic$deviance), sum(fit3$dic$deviance)),
  pD      = c(sum(fit1$dic$penalty),  sum(fit2$dic$penalty),  sum(fit3$dic$penalty)),
  DIC     = c(dic_val(fit1$dic),      dic_val(fit2$dic),      dic_val(fit3$dic))
)
dic_tab$dDIC <- dic_tab$DIC - min(dic_tab$DIC)
cat("\nConfronto DIC \n")
print(round(dic_tab[, -1], 2), row.names = FALSE)
cat("(lower is better; dDIC = delta of best)\n")

# =================================================================
# POSTERIOR PREDICTIVE CHECK on training (M3)
# =================================================================
cat("\ PPC training (M3) \n")
ppc_samp <- coda.samples(fit3$model, "Y_rep",
                         n.iter = 2000, thin = 2,
                         progress.bar = "none")
rep_mat <- as.matrix(ppc_samp)
T_obs <- c(media = mean(exp(Y_tr)), sd = sd(exp(Y_tr)),
           massimo = max(exp(Y_tr)),
           p_oltre_50 = mean(exp(Y_tr) > 50))
T_rep <- t(apply(rep_mat, 1, function(r) {
  x <- exp(r); c(mean(x), sd(x), max(x), mean(x > 50))
}))
ppp <- colMeans(sweep(T_rep, 2, T_obs, ">="))
names(ppp) <- names(T_obs)
cat("Bayes p-value (ideal ~ 0.5):\n"); print(round(ppp, 3))

png("ppc_stats.png", width = 1400, height = 900, res = 130)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (j in seq_along(T_obs)) {
  hist(T_rep[, j], breaks = 40, col = "grey85", border = "white",
       main = paste0(names(T_obs)[j], "  (p = ", round(ppp[j], 3), ")"),
       xlab = "T(Y_rep)")
  abline(v = T_obs[j], col = "red", lwd = 2)
}
par(op); dev.off()

# =================================================================
#  out-of-sample FORECAST (M2 e M3) on the 90-day hold-out
# =================================================================
cat("\nForecast hold-out (90 days)\n")
X_te_full <- build_X(dati_te, X_mean = prep$mean, X_sd = prep$sd)$X
doy_te   <- yday(dati_te$Data)
sin1_te  <- sin(2 * pi * doy_te / 365.25)
cos1_te  <- cos(2 * pi * doy_te / 365.25)
Y_last_train <- Y_tr[length(Y_tr)]
N_te <- nrow(dati_te)

# extract posterior in matrix
post_mat <- function(samp) {
  mat <- as.matrix(samp)
  idx <- sample(seq_len(nrow(mat)), min(2000, nrow(mat)))
  mat[idx, , drop = FALSE]
}

#M2 forecast 
M2 <- post_mat(fit2$samples)
F2 <- matrix(NA_real_, nrow(M2), N_te)
for (it in seq_len(nrow(M2))) {
  a   <- M2[it, "alpha"]; phi <- M2[it, "phi"]
  b   <- M2[it, paste0("beta[", 1:K, "]")]
  sig <- M2[it, "sigma"]
  ylag <- Y_last_train
  for (t in seq_len(N_te)) {
    mu <- a + phi * ylag + sum(b * X_te_full[t, ])
    F2[it, t] <- rnorm(1, mu, sig)
    ylag <- F2[it, t]
  }
}

#M3 forecast
M3 <- post_mat(fit3$samples)
F3 <- matrix(NA_real_, nrow(M3), N_te)
for (it in seq_len(nrow(M3))) {
  a   <- M3[it, "alpha"]; phi <- M3[it, "phi"]
  b   <- M3[it, paste0("beta[", 1:K, "]")]
  sig <- M3[it, "sigma"]
  g1  <- M3[it, "gamma1"]; g2 <- M3[it, "gamma2"]
  nu  <- M3[it, "nu"]
  ylag <- Y_last_train
  for (t in seq_len(N_te)) {
    mu <- a + phi * ylag + sum(b * X_te_full[t, ]) +
          g1 * sin1_te[t] + g2 * cos1_te[t]
    # rt(1, df=nu) * sig + mu = location-scale t
    F3[it, t] <- mu + sig * rt(1, df = nu)
    ylag <- F3[it, t]
  }
}

#forecast matrix
Y_te <- log(dati_te$PM10)
metr <- function(F) {
  pm_pred <- exp(F)
  med <- apply(pm_pred, 2, median)
  lo  <- apply(pm_pred, 2, quantile, 0.025)
  hi  <- apply(pm_pred, 2, quantile, 0.975)
  oss <- dati_te$PM10
  c(
    RMSE     = sqrt(mean((med - oss)^2)),
    MAE      = mean(abs(med - oss)),
    coverage = mean(oss >= lo & oss <= hi),
    width_CI = mean(hi - lo)
  )
}
m2 <- metr(F2); m3 <- metr(F3)
cat("Forecast skill on the 90-day hold-out (PM10 ug/m3 scale):\n")
print(round(rbind(M2 = m2, M3 = m3), 2))
cat("(ideal coverage = 0.95; lower RMSE/MAE = better)\n")

# plot
med2 <- apply(exp(F2), 2, median)
lo2  <- apply(exp(F2), 2, quantile, 0.025); hi2 <- apply(exp(F2), 2, quantile, 0.975)
med3 <- apply(exp(F3), 2, median)
lo3  <- apply(exp(F3), 2, quantile, 0.025); hi3 <- apply(exp(F3), 2, quantile, 0.975)

png("forecast_holdout.png", width = 1400, height = 750, res = 130)
op <- par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(dati_te$Data, dati_te$PM10, pch = 19, col = "black",
     ylim = c(0, max(dati_te$PM10, hi3, hi2)),
     xlab = "Data", ylab = expression(PM[10] ~ "(µg/m³)"),
     main = "Forecast 90-day hold-out  (point = observated)")
polygon(c(dati_te$Data, rev(dati_te$Data)), c(lo2, rev(hi2)),
        col = adjustcolor("steelblue", 0.18), border = NA)
polygon(c(dati_te$Data, rev(dati_te$Data)), c(lo3, rev(hi3)),
        col = adjustcolor("tomato",   0.18), border = NA)
lines(dati_te$Data, med2, col = "steelblue", lwd = 1.6)
lines(dati_te$Data, med3, col = "tomato",   lwd = 1.6)
legend("topright",
       c("osservato","M2 AR(1) mediana","M3 AR1+stag+t mediana",
         "M2 CI 95%","M3 CI 95%"),
       col = c("black","steelblue","tomato",
               adjustcolor("steelblue", 0.5), adjustcolor("tomato", 0.5)),
       pch = c(19, NA, NA, 15, 15), lty = c(NA, 1, 1, NA, NA), bty = "n")
par(op); dev.off()

# =================================================================
#frequentist comparison on training set
# =================================================================
cat("\nfrequentist comparison on training set\n")
df_lm <- data.frame(Y = Y_tr, X_tr,
                    sin1 = sin1_tr, cos1 = cos1_tr)
colnames(df_lm)[2:5] <- c("Temp", "Umidita", "Vento", "Pioggia")
df_lm$Y_lag <- c(NA, head(Y_tr, -1))

fit_lm1 <- lm(Y ~ Temp + Umidita + Vento + Pioggia, data = df_lm)
fit_lm2 <- lm(Y ~ Y_lag + Temp + Umidita + Vento + Pioggia, data = df_lm[-1, ])
fit_lm3 <- lm(Y ~ Y_lag + Temp + Umidita + Vento + Pioggia + sin1 + cos1,
              data = df_lm[-1, ])
cat("\nlm M1\n"); print(round(coef(summary(fit_lm1)), 3))
cat("\nlm M2\n"); print(round(coef(summary(fit_lm2)), 3))
cat("\nlm M3\n"); print(round(coef(summary(fit_lm3)), 3))

#saved anything
saveRDS(list(fit1 = fit1, fit2 = fit2, fit3 = fit3,
             dic_tab = dic_tab, ppp = ppp,
             metr_M2 = m2, metr_M3 = m3,
             F2 = F2, F3 = F3,
             prep = prep, dati = dati, idx_tr = idx_tr, idx_te = idx_te,
             fit_lm1 = fit_lm1, fit_lm2 = fit_lm2, fit_lm3 = fit_lm3),
        "real_fit_results.rds")
cat("\n output: trace_model{1,2,3}.png, ppc_stats.png,",
    "forecast_holdout.png, real_fit_results.rds\n")
