# Bayesian Modelling of Daily PM10 Concentrations in Rome (2021–2024)

Bayesian analysis of daily PM10 concentrations at the Roma‑Preneste air‑quality
station (EEA, *L.GO PERESTRELLO*, urban background), using daily weather
covariates from ARPA Lazio. Models are fitted in **JAGS** via `rjags`.

Three models are compared on real data with a hold‑out period:

- **M1** – linear regression on weather covariates
- **M2** – AR(1) model
- **M3** – AR(1) + seasonality + Student‑t errors (final model)

Comparison uses DIC, posterior predictive checks and out‑of‑sample forecasts.
Simulation‑based recovery checks validate the estimation procedure.

## Repository contents

| File | Description |
|------|-------------|
| `00_choose_station.R` | Picks the ARPA Lazio weather station closest to Roma‑Preneste |
| `01_data_prep.R` | Downloads PM10 (EEA, package `EEAaq`), reads ARPA weather data, builds `dataset_finale.rds/.csv` |
| `02_sim_recovery.R` | Simulation‑based parameter recovery check |
| `03_real_fit.R` | Fits M1, M2, M3 on real data: hold‑out, DIC, PPC, forecast |
| `04_sim_recovery_M3.R` | Simulation‑based recovery for the final model M3 |
| `05_residuals.R` | Residual diagnostics for M3 |
| `06_regen_figures.R` | Regenerates the final figures |
| `07_extra_figures.R` | Extra figures for the presentation |
| `report.pdf` | Project report |
| `presentation_final.html` | Final presentation |

## How to run

1. Install [JAGS](https://mcmc-jags.sourceforge.io/) and the R packages:
   ```r
   install.packages(c("rjags", "coda", "dplyr", "tidyr", "lubridate", "EEAaq"))
   ```
2. Download the ARPA Lazio daily weather files
   `dati-meteorologici_2021.csv` … `dati-meteorologici_2024.csv`
   and place them in the project folder.
3. Open `progetto_fsl2.Rproj` in RStudio and run the scripts in order
   (`00` → `07`).
