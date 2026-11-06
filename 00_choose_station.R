# Select the ARPA Lazio weather station closest to Rome-Preneste
# Coordinates from the ARPA Lazio "technical data sheets"

suppressPackageStartupMessages({library(dplyr)})

stz <- data.frame(
  cod_staz = sprintf("AL%03d", 1:9),
  nome     = c("CNR Tor Vergata (Roma)",
               "Latina",
               "Cavaliere (Roma)",
               "Castel di Guido (Roma)",
               "Istituto Jucci (Rieti)",
               "Aeroporto Mil. Frosinone",
               "Boncompagni (Roma centro)",
               "Aeroporto Mil. Viterbo",
               "Ceprano"),
  zona     = c("metropolitana","urbana","rurale","rurale","urbana",
               "urbana","urbana","rurale","rurale"),
  lat      = c(41.84, 41.48, 41.92, 41.89, 42.40,
               41.65, 41.90, 42.43, 41.55),
  lon      = c(12.64, 12.85, 12.65, 12.27, 12.85,
               13.32, 12.49, 12.07, 13.50)
)

# Roma-Preneste (L.GO PERESTRELLO, from EEA)
PRE_lat <- 41.88602
PRE_lon <- 12.54162

#Haversine distance (km)
hav <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to <- function(x) x * pi / 180
  dlat <- to(lat2 - lat1); dlon <- to(lon2 - lon1)
  a <- sin(dlat/2)^2 + cos(to(lat1)) * cos(to(lat2)) * sin(dlon/2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

stz$dist_km <- hav(PRE_lat, PRE_lon, stz$lat, stz$lon)
stz <- stz %>% arrange(dist_km)
cat("Distances from the control unit Roma-Preneste (41.886N, 12.542E):\n")
stz_print <- stz
stz_print$lat <- round(stz_print$lat, 2)
stz_print$lon <- round(stz_print$lon, 2)
stz_print$dist_km <- round(stz_print$dist_km, 1)
print(stz_print, row.names = FALSE)
cat("\nSelected weather station:", stz$cod_staz[1],
    "-", stz$nome[1],
    "(", round(stz$dist_km[1], 1), "km )\n")

saveRDS(stz, "stazioni_arpa.rds")

