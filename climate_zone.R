# ========================
# 0. LIBRARY LOADING
# ========================
library(readr)
library(dplyr)
library(ggplot2)
library(cluster)
library(RANN)
library(tibble)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
#=======================
# 1.GLOBAL RCC LOAD AND CALCULATE YANG CENTROIDS
#=======================

RCC <- read_csv("RCC.csv", show_col_types = FALSE)

centroides_yang5 <- RCC %>%
  group_by(tp) %>%
  summarise(
    AOD = mean(AOD, na.rm = TRUE),
    albedo = mean(albedo, na.rm = TRUE),
    cloud = mean(cloud, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(tp)

write_csv(centroides_yang5, "centroides_yang5.csv")

# ========================
# 2. FILTER ARGENTINA AND LOCAL CLUSTERING
# ========================

RCC_arg <- RCC %>%
  filter(lat >= -56, lat <= -20, lon >= -75, lon <= -53)

# ========================
# 3. NORMALIZAR VARIABLES Y CLUSTERING K=5
# ========================
set.seed(123)

variables_cluster <- dplyr::select(RCC_arg, AOD, albedo, cloud)
variables_norm <- scale(variables_cluster)
clust_k5 <- kmeans(variables_norm, centers = 5, nstart = 25)

RCC_arg <- RCC_arg %>%
  mutate(tp_local = clust_k5$cluster)

write_csv(RCC_arg, "RCC_argentina_tp_local.csv")

# ========================
# 4.AUXILIARY FUNCTIONS
# ========================
# --- Assign Yang classes to local clusters ---
asignar_tp_original_yang <- function(archivo_RCC_local, archivo_salida = NULL) {
  centroides_yang <- read_csv("centroides_yang5.csv", show_col_types = FALSE)
  RCC_local <- read_csv(archivo_RCC_local, show_col_types = FALSE)

  centroides_locales <- RCC_local %>%
    group_by(tp_local) %>%
    summarise(
      AOD = mean(AOD, na.rm = TRUE),
      albedo = mean(albedo, na.rm = TRUE),
      cloud = mean(cloud, na.rm = TRUE),
      .groups = "drop"
    )

  match_tp <- sapply(1:nrow(centroides_locales), function(i) {
    distancias <- apply(centroides_yang[, -1], 1, function(centro) {
      sum((as.numeric(centroides_locales[i, 2:4]) - centro)^2)
    })
    centroides_yang$tp[which.min(distancias)]
  })

  centroides_locales$tp_yang_match <- match_tp

  RCC_final <- RCC_local %>%
    left_join(centroides_locales[, c("tp_local", "tp_yang_match")], by = "tp_local")

  if (!is.null(archivo_salida)) {
    write_csv(RCC_final, archivo_salida)
  }

  return(RCC_final)
}

# --- Assign climate zone to stations ---
obtener_zona_climatica <- function(archivo_estaciones, archivo_RCC) {
  estaciones <- read_delim(archivo_estaciones, delim = ",", show_col_types = FALSE)
  RCC <- read_csv(archivo_RCC, show_col_types = FALSE)

  stopifnot(all(c("lon", "lat", "tp_local") %in% names(RCC)))
  stopifnot(all(c("longitud", "latitud") %in% names(estaciones)))

  nn <- nn2(data = RCC[, c("lon", "lat")],
            query = estaciones[, c("longitud", "latitud")],
            k = 1)
  idx <- nn$nn.idx[, 1]

  estaciones <- estaciones %>%
    mutate(
      tp = RCC$tp_local[idx],
      AOD = RCC$AOD[idx],
      albedo = RCC$albedo[idx],
      cloud = RCC$cloud[idx]
    )

  return(estaciones)
}

# ========================
# 5. EXECUTE FUNCTIONS
# ========================

RCC_arg_refinada <- asignar_tp_original_yang(
  archivo_RCC_local = "RCC_argentina_tp_local.csv",
  archivo_salida = "RCC_argentina_tp_yang_match.csv"
)

resultado <- obtener_zona_climatica(
  archivo_estaciones = "parametros_estaciones_h.csv",
  archivo_RCC = "RCC_argentina_tp_yang.csv"
)

write_csv(resultado, "parametros_estaciones_tp.csv")

# ========================
# 6. VISUALIZATION
# ========================

# Basemap of Argentina
argentina <- ne_countries(scale = "medium", country = "Argentina", returnclass = "sf")

ggplot() +
  geom_point(data = RCC_arg, aes(x = lon, y = lat, color = factor(tp_local)), size = 0.8) +
  geom_sf(data = argentina, fill = NA, color = "black", linewidth = 0.4) +
  coord_sf(xlim = c(-75, -53), ylim = c(-56, -21), expand = FALSE) +
  scale_color_brewer(palette = "Set1", name = "tp_local") +
  theme_minimal() +
  labs(title = "Climate zones on the outline of Argentina")
