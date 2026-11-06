# FASE 3 - Simulation-based recovery check

suppressPackageStartupMessages({
  library(rjags)
  library(coda)
})

# setwd() not needed: open progetto_fsl2.Rproj so the working directory is the project root
set.seed(20260607)

# ----------------------- true parameter ----------------------
N            <- 800
K            <- 4
true_alpha   <- 1.20
true_phi     <- 0.55
true_beta    <- c(-0.30, 0.15, -0.40, -0.10)   # Temp, Umid, Vento, Pioggia
true_sigma   <- 0.30
true_tau     <- 1 / true_sigma^2

# fake predictors (standardized~ N(0,1))
X <- matrix(rnorm(N * K), nrow = N, ncol = K)

# AR(1) simulation
Y <- numeric(N)
Y[1] <- rnorm(1,
              mean = true_alpha / (1 - true_phi) + sum(true_beta * X[1, ]),
              sd   = true_sigma)
for (t in 2:N) {
  mu_t <- true_alpha + true_phi * Y[t - 1] + sum(true_beta * X[t, ])
  Y[t] <- rnorm(1, mu_t, true_sigma)
}

cat("Simulated", N, "point.  Y range:",
    round(range(Y), 2), "\n\n")

# ------------------------- JAGS ------------------------------
data_jags <- list(N = N, K = K, Y = Y, X = X)
inits     <- function() list(alpha = 0, beta = rep(0, K),
                             phi = 0, tau = 1)

cat("Compilation JAGS (AR1 model)\n")
m <- jags.model("model2_ar1.bug",
                data     = data_jags,
                inits    = inits,
                n.chains = 3,
                quiet    = TRUE)

cat("Burn-in 2000\n")
update(m, 2000, progress.bar = "none")

cat("Sampling 8000\n")
samp <- coda.samples(m,
                     variable.names = c("alpha", "beta", "phi", "sigma"),
                     n.iter         = 8000,
                     thin           = 2,
                     progress.bar   = "none")

# --------- summary and comparison with real values --------------
s <- summary(samp)
post <- as.data.frame(s$statistics[, c("Mean", "SD")])
post$Q2.5  <- s$quantiles[, "2.5%"]
post$Q97.5 <- s$quantiles[, "97.5%"]
post$true  <- c(true_alpha, true_beta, true_phi, true_sigma)
post$copre <- (post$Q2.5 <= post$true) & (post$true <= post$Q97.5)
cat("\n Posterior estimates vs. true values \n")
print(round(post, 3))

# convergence diagnostics
gr <- gelman.diag(samp, multivariate = FALSE)
cat("\nGelman-Rubin (R-hat):\n")
print(gr$psrf)


# saving for control
saveRDS(list(samp = samp, true = post$true,
             names = rownames(post)),
        "sim_recovery_samples.rds")

png("sim_recovery_traceplot.png", width = 1200, height = 900, res = 130)
op <- par(mfrow = c(3, 2), mar = c(3, 3, 2, 1))
plot(samp[, "alpha"],   main = "alpha")
plot(samp[, "phi"],     main = "phi")
plot(samp[, "beta[1]"], main = "beta[1] (Temp)")
plot(samp[, "beta[2]"], main = "beta[2] (Umid)")
plot(samp[, "beta[3]"], main = "beta[3] (Vento)")
plot(samp[, "sigma"],   main = "sigma")
par(op)
dev.off()
cat("\nTraceplot saved in sim_recovery_traceplot.png\n")


