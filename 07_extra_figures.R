
#additional graphics for the presentation

suppressPackageStartupMessages({library(dplyr); library(lubridate)})
# setwd() not needed: open progetto_fsl2.Rproj so the working directory is the project root

dati <- readRDS("dataset_finale.rds") %>% arrange(Data)
cat("N =", nrow(dati), "days\n")

# 1) time series PM10 with line theshold 50 ug/m3 
png("pm10_timeseries.png", width = 1500, height = 700, res = 130)
par(mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
plot(dati$Data, dati$PM10, type = "h", col = "grey75",
     xlab = "Date", ylab = expression(PM[10] ~ "(µg/m³)"),
     main = "Daily PM10 at Rome–Preneste, 2021–2024",
     ylim = c(0, max(dati$PM10) * 1.1))

sm <- lowess(as.numeric(dati$Data), dati$PM10, f = 0.08)
lines(as.Date(sm$x, origin = "1970-01-01"), sm$y,
      col = "tomato", lwd = 2.5)
abline(h = 50, col = "steelblue", lwd = 1.5, lty = 2)
legend("topright", inset = 0.02,
       legend = c("daily observations",
                  "LOWESS smoother (annual cycle)",
                  "EU daily limit 50 µg/m³"),
       col = c("grey45", "tomato", "steelblue"),
       lty = c(1, 1, 2), lwd = c(1, 2.5, 1.5), bty = "n", cex = 0.9)
#peek annotations
i_max <- which.max(dati$PM10)
text(dati$Data[i_max], dati$PM10[i_max] + 5,
     paste0(round(dati$PM10[i_max], 0), " µg/m³"),
     col = "black", cex = 0.85, font = 2)
dev.off()
cat("Saved pm10_timeseries.png\n")

# 2) PM10 raw vs log(PM10) 
png("pm10_log_hist.png", width = 1500, height = 650, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
# raw
hist(dati$PM10, breaks = 60, col = "grey80", border = "white",
     xlab = expression(PM[10] ~ "(µg/m³)"), freq = FALSE,
     main = "Raw scale: strong right-skew",
     ylim = c(0, 0.07))
abline(v = mean(dati$PM10), col = "tomato", lwd = 2)
abline(v = median(dati$PM10), col = "steelblue", lwd = 2, lty = 2)
legend("topright",
       legend = c(paste0("mean = ", round(mean(dati$PM10), 1)),
                  paste0("median = ", round(median(dati$PM10), 1))),
       col = c("tomato", "steelblue"), lty = c(1, 2),
       lwd = 2, bty = "n", cex = 0.9)
# log
hist(log(dati$PM10), breaks = 60, col = "grey80", border = "white",
     xlab = expression(log(PM[10])), freq = FALSE,
     main = "Log scale: nearly symmetric",
     ylim = c(0, 1.0))
abline(v = mean(log(dati$PM10)), col = "tomato", lwd = 2)
abline(v = median(log(dati$PM10)), col = "steelblue", lwd = 2, lty = 2)
legend("topright",
       legend = c(paste0("mean = ", round(mean(log(dati$PM10)), 2)),
                  paste0("median = ", round(median(log(dati$PM10)), 2))),
       col = c("tomato", "steelblue"), lty = c(1, 2),
       lwd = 2, bty = "n", cex = 0.9)
dev.off()
cat("Saved pm10_log_hist.png\n")

# 3) bar chart DIC
res <- readRDS("real_fit_results.rds")
dic_vals <- c("M1 linear" = 1252.3,
              "M2 AR(1)" = 570.2,
              "M3 AR(1)+seas+t" = 364.5)
png("dic_barplot.png", width = 1200, height = 650, res = 130)
par(mar = c(4, 5, 2.5, 1), mgp = c(2.5, 0.7, 0))
cols <- c("#d9544d", "#e6a23c", "#1e7e3c")
bp <- barplot(dic_vals, horiz = TRUE,
              col = cols, border = NA, las = 1,
              xlab = "DIC (lower is better)",
              main = "Deviance Information Criterion",
              xlim = c(0, max(dic_vals) * 1.15))
text(dic_vals, bp, labels = sprintf("  %.1f", dic_vals),
     pos = 4, cex = 0.95, font = 2)
# delta arrows
delta12 <- dic_vals[1] - dic_vals[2]
delta23 <- dic_vals[2] - dic_vals[3]
mtext(paste0("Δ = ", round(delta12, 0)),
      at = (dic_vals[1] + dic_vals[2]) / 2, side = 3, line = -2,
      cex = 0.85, col = "#444")
mtext(paste0("Δ = ", round(delta23, 0)),
      at = (dic_vals[2] + dic_vals[3]) / 2, side = 3, line = -3,
      cex = 0.85, col = "#444")
dev.off()
cat("Saved dic_barplot.png\n")
