# Residual diagnostics for Model 3 (AR1 + stag + t)

suppressPackageStartupMessages({
  library(rjags); library(coda); library(dplyr); library(lubridate)
})

# setwd() not needed: open progetto_fsl2.Rproj so the working directory is the project root

res   <- readRDS("real_fit_results.rds")
fit3  <- res$fit3
dati  <- res$dati
idx_tr <- res$idx_tr
dati_tr <- dati[idx_tr, ]
Y_tr  <- log(dati_tr$PM10)
N     <- length(Y_tr); K <- 4

#standardized predictors as in the fit
X_tr  <- scale(as.matrix(dati_tr[, c("Temp","Umidita","Vento","Pioggia")]),
               center = res$prep$mean, scale = res$prep$sd)
doy_tr  <- yday(dati_tr$Data)
sin1_tr <- sin(2 * pi * doy_tr / 365.25)
cos1_tr <- cos(2 * pi * doy_tr / 365.25)

# ---- precise estimates (posterior mean) ------------------------
mat   <- as.matrix(fit3$samples)
get   <- function(nm) mean(mat[, nm])
alpha <- get("alpha"); phi   <- get("phi")
beta  <- sapply(1:K, function(k) get(paste0("beta[", k, "]")))
gamma1 <- get("gamma1"); gamma2 <- get("gamma2")
sigma <- get("sigma");   nu    <- get("nu")
cat("M3 precise estimates:\n")
cat(sprintf("  alpha=%.3f  phi=%.3f  sigma=%.3f  nu=%.2f\n",
            alpha, phi, sigma, nu))
cat(sprintf("  beta=[%s]  gamma=(%.3f, %.3f)\n",
            paste(round(beta, 3), collapse = ", "), gamma1, gamma2))

# ----conditional fitted e standardized residuals ----------
mu_fit <- numeric(N)
mu_fit[1] <- alpha / (1 - phi) + sum(beta * X_tr[1, ]) +
             gamma1 * sin1_tr[1] + gamma2 * cos1_tr[1]
for (t in 2:N) {
  mu_fit[t] <- alpha + phi * Y_tr[t - 1] + sum(beta * X_tr[t, ]) +
               gamma1 * sin1_tr[t] + gamma2 * cos1_tr[t]
}
res_std <- (Y_tr - mu_fit) / sigma

cat("\nstandardized residuals: mean =", round(mean(res_std), 3),
    " SD =", round(sd(res_std), 3),
    "\n  (expected average~0, SD~sqrt(nu/(nu-2))) =",
    round(sqrt(nu/(nu - 2)), 3), "\n")

# ---- Ljung-Box: residual autocorrelation -------------------
lb1  <- Box.test(res_std, lag = 1,  type = "Ljung-Box")
lb10 <- Box.test(res_std, lag = 10, type = "Ljung-Box")
lb20 <- Box.test(res_std, lag = 20, type = "Ljung-Box")
cat("Ljung-Box p-values (H0: no autocorrelation):\n")
cat(sprintf("  lag 1:  p = %.3f\n", lb1$p.value))
cat(sprintf("  lag 10: p = %.3f\n", lb10$p.value))
cat(sprintf("  lag 20: p = %.3f\n", lb20$p.value))

# ---- diagnostic plot--------------------------------------
png("residuals_diag.png", width = 1500, height = 1000, res = 130)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# (1) Residuals vs data
plot(dati_tr$Data, res_std, pch = 19, cex = 0.4, col = "grey25",
     xlab = "Data", ylab = "residuo std",
     main = "Standardized residues over time")
abline(h = 0, col = "red")
abline(h = c(-3, 3), col = "red", lty = 2)

# (2) ACF
acf(res_std, lag.max = 40, main = "ACF residuals")
# (3) PACF
pacf(res_std, lag.max = 40, main = "PACF residuals")

# (4) QQ vs t_nu
qs <- ppoints(N)
q_emp <- sort(res_std)
q_th  <- qt(qs, df = nu)
plot(q_th, q_emp, pch = 19, cex = 0.4, col = "grey25",
     xlab = paste0("theoretical quantiles t_", round(nu, 1)),
     ylab = "empirical quantiles",
     main = "QQ-plot vs t-Student")
abline(0, 1, col = "red")

# (5) Residuals vs fitted
plot(mu_fit, res_std, pch = 19, cex = 0.4, col = "grey25",
     xlab = "expected value (log PM10)", ylab = "residual std",
     main = "Residuals vs fitted")
abline(h = 0, col = "red")
lines(lowess(mu_fit, res_std), col = "tomato", lwd = 2)

# (6) Histogram of residuals vs. theoretical t-density
hist(res_std, breaks = 60, freq = FALSE,
     col = "grey85", border = "white",
     main = paste0("Hist residuals vs t_", round(nu, 1), " theorical"),
     xlab = "residual std", xlim = c(-8, 8))
xx <- seq(-8, 8, by = 0.05)
lines(xx, dt(xx, df = nu), col = "tomato", lwd = 2)
lines(xx, dnorm(xx), col = "steelblue", lwd = 1.5, lty = 2)
legend("topright", c(paste0("t_", round(nu, 1)), "N(0,1)"),
       col = c("tomato","steelblue"), lty = c(1, 2), bty = "n")

par(op); dev.off()

# ---- summary------------------------------------
out <- c(
  media_residui     = mean(res_std),
  sd_residui        = sd(res_std),
  sd_attesa_t       = sqrt(nu/(nu - 2)),
  skewness          = mean((res_std - mean(res_std))^3) /
                       sd(res_std)^3,
  kurtosis          = mean((res_std - mean(res_std))^4) /
                       sd(res_std)^4,
  LB_lag1_p         = lb1$p.value,
  LB_lag10_p        = lb10$p.value,
  LB_lag20_p        = lb20$p.value,
  pct_oltre_3sd     = mean(abs(res_std) > 3),
  outlier_max_abs   = max(abs(res_std))
)
cat("\n--- diagnostic summary M3 ---\n")
print(round(out, 3))

saveRDS(list(res_std = res_std, mu_fit = mu_fit, summary = out),
        "residuals_M3.rds")
cat("\nsaved diagnostic in residuals_diag.png\n")

