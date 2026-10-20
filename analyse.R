# =============================================================================
#  Domaine vital du grand dauphin (Tursiops truncatus) dans le Parc naturel
#  marin du Cap Corse et de l'Agriate et chevauchement avec le trafic maritime
#
#  Analyse exploratoire à partir de données ouvertes
#  Auteure : Gwladys Helaine — octobre 2026
#
#  Étapes :
#   0. Paramètres et packages
#   1. Périmètre d'étude (Parc)
#   2. Bathymétrie et distance à la côte
#   3. Observations : téléchargement OBIS + GBIF, bancarisation, nettoyage
#   4. Domaine vital : noyau (KDE 50 % / 95 %) et polygone convexe (MCP 95 %)
#   5. Habitat : comparaison utilisation / disponibilité (profondeur, côte)
#   6. Pressions : trafic maritime (EMODnet) dans et hors du domaine vital
#   7. Cartes et figures
#   8. Exports (données propres, couches SIG, tableaux, résultats)
#
#  Lancer depuis le dossier du projet (ouvrir delfinu-exercice.Rproj ou
#  setwd() sur ce dossier), puis : source("analyse.R")
#
#  REPRODUCTIBILITÉ
#   - Logiciels : analyse réalisée sous R 4.5 (Linux) ; les versions exactes de
#     R et de chaque package sont enregistrées à chaque exécution dans
#     outputs/session_info.txt et outputs/versions_packages.csv.
#     Sous Linux, le package sf demande les bibliothèques système GDAL, GEOS et
#     PROJ (paquets libgdal-dev, libgeos-dev, libproj-dev, libudunits2-dev).
#   - Données : les observations OBIS/GBIF et la bathymétrie NOAA sont
#     téléchargées au premier lancement puis conservées dans data/raw/ ; les
#     lancements suivants réutilisent ces fichiers, ce qui garantit les mêmes
#     résultats. La date du téléchargement est enregistrée dans
#     outputs/sources_donnees.csv. Pour repartir de données à jour :
#     FORCER_TELECHARGEMENT <- TRUE (les résultats peuvent alors changer, les
#     portails étant alimentés en continu).
#     Le périmètre du Parc (PatriNat, base des espaces protégés, 07/2026) et la
#     densité de navires (EMODnet Human Activities, moyenne annuelle 2023) sont
#     à déposer à la main dans data/raw/parc/ et data/raw/trafic/.
#   - Hasard : le tirage des points de disponibilité est rendu reproductible
#     par set.seed(2026).
#   - Paramètres : tous les seuils et choix de méthode sont regroupés dans la
#     section 0 ; les modifier permet de tester la sensibilité des résultats.
# =============================================================================


# ---- 0. Paramètres et packages ---------------------------------------------
# Objectif : regrouper tous les choix de l'analyse (seuils de nettoyage,
#            résolutions, réglages des méthodes) et charger les outils.
# Pourquoi : un paramètre modifié ici se répercute partout ; on peut ainsi
#            tester la sensibilité des résultats sans toucher au reste du code.
# Sorties  : paramètres, packages chargés, dossiers du projet créés.

ANNEE_MIN          <- 2000   # observations plus anciennes écartées
INCERTITUDE_MAX_M  <- 5000   # précision de localisation minimale acceptée
TOLERANCE_TERRE_M  <- 20     # altitude max tolérée (grille bathy grossière près des côtes)
TAILLE_CELLULE_M   <- 1000   # agrégation : 1 obs. par cellule de 1 km et par jour
RESOLUTION_BATHY_MIN <- 0.25  # résolution de la bathymétrie NOAA, en minutes d'arc (0.25 = 15")
RESOLUTION_M       <- 250    # résolution des rasters reprojetés (500 si repli sur 1 minute)
N_MIN_KDE          <- 30     # en dessous, le domaine vital est très incertain
FACTEUR_H          <- 1      # multiplicateur de la largeur de bande href (1 = href)
N_DISPO_PAR_OBS    <- 10     # points de disponibilité tirés par observation
CRS_METRIQUE       <- "EPSG:2154"  # Lambert-93 (couvre la Corse)
FORCER_TELECHARGEMENT <- FALSE

# Emprise de téléchargement (WGS84) : un peu plus large que le Parc
BBOX <- c(xmin = 8.5, ymin = 42.5, xmax = 10.0, ymax = 43.5)

# Emprise approximative du Parc, utilisée SEULEMENT si la couche officielle
# n'a pas été déposée dans data/raw/parc/
BBOX_PARC_APPROX <- c(xmin = 8.8, ymin = 42.62, xmax = 9.75, ymax = 43.35)

set.seed(2026)

# robis (OBIS) est facultatif : s'il n'est pas installé, le script interroge
# directement l'API d'OBIS avec jsonlite
paquets <- c("sf", "terra", "dplyr", "readr", "ggplot2", "MASS",
             "rgbif", "marmap", "jsonlite")
manquants <- setdiff(paquets, rownames(installed.packages()))
if (length(manquants) > 0) {
  message("Installation des packages manquants : ", paste(manquants, collapse = ", "))
  tryCatch(install.packages(manquants),
           error = function(e) message("Installation impossible : ", conditionMessage(e)))
}

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
})
sf_use_s2(FALSE)  # calculs planaires en Lambert-93

for (d in c("data/raw/parc", "data/raw/trafic", "data/propre",
            "outputs/figures", "outputs/sig", "outputs/tableaux")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Petite fonction utilitaire : récupère une colonne si elle existe, sinon NA
colonne <- function(d, nom) if (nom %in% names(d)) d[[nom]] else rep(NA, nrow(d))

# Polygone d'une emprise, au format sf (WGS84)
emprise_vers_sf <- function(b) {
  st_as_sfc(st_bbox(c(xmin = b[["xmin"]], ymin = b[["ymin"]],
                      xmax = b[["xmax"]], ymax = b[["ymax"]]),
                    crs = st_crs(4326)))
}

resultats <- list(date_analyse = Sys.Date())


# ---- 1. Périmètre d'étude ----------------------------------------------------
# Objectif : lire le contour officiel du Parc (base PatriNat des espaces
#            protégés) et définir l'emprise de téléchargement des données.
# Méthode  : recherche automatique de la bonne couche et du bon polygone dans
#            le GeoPackage ; emprise élargie au Parc + 0,15° de marge.
# Sorties  : parc (polygone en Lambert-93), BBOX (emprise en WGS84).

fichiers_parc <- list.files("data/raw/parc", pattern = "\\.(shp|gpkg|geojson)$",
                            full.names = TRUE, ignore.case = TRUE)

if (length(fichiers_parc) > 0) {
  f_parc    <- fichiers_parc[1]
  motif     <- "cap ?corse|capicorsu"
  motif_pnm <- "parc naturel marin|parcu naturale marinu"

  # Texte de chaque ligne (tous les champs texte collés) pour y chercher les noms
  texte_lignes <- function(d) {
    d <- st_drop_geometry(d)
    car <- vapply(d, is.character, logical(1))
    if (!any(car)) return(rep("", nrow(d)))
    apply(d[, car, drop = FALSE], 1, paste, collapse = " ")
  }
  est_le_parc <- function(d) {
    t <- texte_lignes(d)
    grepl(motif, t, ignore.case = TRUE) & grepl(motif_pnm, t, ignore.case = TRUE)
  }

  # Un GeoPackage peut contenir plusieurs couches (une par type d'espace protégé) :
  # on prend celle des parcs naturels marins ("pnm"), sinon celle qui contient
  # une ligne "Parc naturel marin ... Cap Corse"
  couches <- st_layers(f_parc)$name
  couche  <- couches[grepl("pnm", couches, ignore.case = TRUE)][1]
  if (is.na(couche)) {
    for (cc in couches) {
      tmp <- tryCatch(st_read(f_parc, layer = cc, quiet = TRUE), error = function(e) NULL)
      if (inherits(tmp, "sf") && any(est_le_parc(tmp))) {
        couche <- cc
        break
      }
    }
  }
  if (is.na(couche)) {
    stop("Parc du Cap Corse introuvable dans ", f_parc, ". Couches disponibles : ",
         paste(couches, collapse = ", "))
  }
  parc <- st_read(f_parc, layer = couche, quiet = TRUE)
  message("Couche utilisée pour le Parc : ", couche)

  # On garde les lignes qui citent le Cap Corse ; dans une couche de parcs marins,
  # le nom seul suffit, sinon on exige aussi la mention "parc naturel marin"
  t_parc <- texte_lignes(parc)
  garde  <- grepl(motif, t_parc, ignore.case = TRUE)
  if (!grepl("pnm", couche, ignore.case = TRUE)) {
    garde <- garde & grepl(motif_pnm, t_parc, ignore.case = TRUE)
  }
  parc <- parc[garde, ]
  message("Polygone(s) retenu(s) pour le Parc :\n  ",
          paste(substr(t_parc[garde], 1, 150), collapse = "\n  "))
  if (nrow(parc) == 0) stop("Aucun polygone du Parc du Cap Corse trouvé dans la couche ", couche)
  if (nrow(parc) > 1) message(nrow(parc), " polygones retenus pour le Parc (fusionnés).")
  parc <- st_sf(geometry = st_union(st_make_valid(st_geometry(parc))))
  resultats$parc_officiel <- TRUE
  message("Périmètre du Parc lu depuis : ", fichiers_parc[1])
} else {
  parc <- st_sf(geometry = emprise_vers_sf(BBOX_PARC_APPROX))
  resultats$parc_officiel <- FALSE
  warning("Couche du Parc absente de data/raw/parc/ : emprise APPROXIMATIVE utilisée. ",
          "Voir le GUIDE pour télécharger le périmètre officiel.")
}
parc <- st_transform(parc, CRS_METRIQUE)

# L'emprise de téléchargement doit couvrir tout le Parc (+ 0,15° de marge)
bb_parc <- st_bbox(st_transform(parc, 4326))
BBOX <- c(xmin = min(BBOX[["xmin"]], bb_parc[["xmin"]] - 0.15),
          ymin = min(BBOX[["ymin"]], bb_parc[["ymin"]] - 0.15),
          xmax = max(BBOX[["xmax"]], bb_parc[["xmax"]] + 0.15),
          ymax = max(BBOX[["ymax"]], bb_parc[["ymax"]] + 0.15))
BBOX <- round(BBOX, 2)
resultats$surface_parc_km2 <- as.numeric(st_area(parc)) / 1e6
message("Surface du périmètre du Parc : ", round(resultats$surface_parc_km2), " km²")
message("Emprise de téléchargement : ", paste(names(BBOX), BBOX, collapse = ", "))

# Si l'emprise a changé depuis le dernier lancement, on retélécharge
fichier_emprise <- "data/raw/emprise_telechargement.rds"
if (!file.exists(fichier_emprise) ||
    !isTRUE(all.equal(readRDS(fichier_emprise), BBOX))) {
  message("Emprise modifiée : les données en cache seront retéléchargées.")
  FORCER_TELECHARGEMENT <- TRUE
}
if (FORCER_TELECHARGEMENT) {
  unlink(c(list.files("data/raw", pattern = "^bathymetrie_.*min\\.tif$", full.names = TRUE),
           "data/raw/observations_brutes.csv"))
}
saveRDS(BBOX, fichier_emprise)


# ---- 2. Bathymétrie et distance à la côte ------------------------------------
# Objectif : décrire le milieu (profondeur, distance à la côte) sur toute la
#            zone d'étude.
# Méthode  : bathymétrie NOAA ETOPO 2022 (via marmap), reprojetée en Lambert-93 ;
#            terre = altitude positive ; distance de chaque cellule de mer à la
#            terre la plus proche.
# Sorties  : bathy_l93, dist_cote (rasters), zone_etude (partie marine du Parc).

# Une bathymétrie personnelle (GEBCO, EMODnet...) peut être déposée sous
# data/raw/bathymetrie_perso.tif ; sinon on télécharge ETOPO 2022 via marmap.
fichier_perso <- "data/raw/bathymetrie_perso.tif"
telecharger_bathy <- function(res) {
  b <- marmap::getNOAA.bathy(lon1 = BBOX[["xmin"]], lon2 = BBOX[["xmax"]],
                             lat1 = BBOX[["ymin"]], lat2 = BBOX[["ymax"]],
                             resolution = res, keep = FALSE)
  xyz <- marmap::as.xyz(b)
  names(xyz) <- c("x", "y", "z")
  rast(xyz, type = "xyz", crs = "EPSG:4326")
}

if (file.exists(fichier_perso)) {
  fichier_bathy <- fichier_perso
  message("Bathymétrie personnelle utilisée : ", fichier_perso)
} else {
  res_bathy <- RESOLUTION_BATHY_MIN
  fichier_bathy <- sprintf("data/raw/bathymetrie_%smin.tif", res_bathy)
  if (!file.exists(fichier_bathy)) {
    message("Téléchargement de la bathymétrie (NOAA ETOPO 2022 via marmap, ",
            res_bathy, " minute d'arc)...")
    b <- tryCatch(telecharger_bathy(res_bathy), error = function(e) {
      message("Résolution ", res_bathy, "' indisponible (", conditionMessage(e),
              ") : repli sur 1 minute d'arc.")
      NULL
    })
    if (is.null(b)) {
      res_bathy <- 1
      fichier_bathy <- sprintf("data/raw/bathymetrie_%smin.tif", res_bathy)
      if (!file.exists(fichier_bathy)) b <- telecharger_bathy(res_bathy)
    }
    if (!is.null(b)) writeRaster(b, fichier_bathy, overwrite = TRUE)
  }
}
bathy <- rast(fichier_bathy)
# Résolution réelle (en mètres, dans le sens nord-sud) et grille de calcul adaptée
resultats$resolution_bathy_m <- if (is.lonlat(bathy)) res(bathy)[2] * 111320 else res(bathy)[2]
RESOLUTION_M <- if (resultats$resolution_bathy_m < 700) 250 else 500
message("Résolution de la bathymétrie : environ ", round(resultats$resolution_bathy_m), " m")
names(bathy) <- "profondeur"

# Convention : valeurs négatives en mer. Certaines sources donnent la profondeur
# en positif : on vérifie au centre du Parc et on inverse si besoin.
centre_parc <- st_transform(st_point_on_surface(parc), crs(bathy))
if (isTRUE(extract(bathy, vect(centre_parc))[1, 2] > 0)) {
  bathy <- -bathy
  names(bathy) <- "profondeur"
  message("Bathymétrie positive détectée : signe inversé.")
}

bathy_l93 <- project(bathy, CRS_METRIQUE, method = "bilinear", res = RESOLUTION_M)
names(bathy_l93) <- "profondeur"

terre    <- ifel(bathy_l93 >= 0, 1, NA)          # 1 = terre, NA = mer
mer      <- ifel(bathy_l93 <  0, 1, NA)          # 1 = mer
dist_cote <- mask(distance(terre), mer)          # distance (m) à la terre la plus proche
names(dist_cote) <- "dist_cote"

terre_sf <- st_as_sf(as.polygons(terre))
mer_sf   <- st_as_sf(as.polygons(mer))

# Zone d'étude = partie marine du Parc
zone_etude <- st_sf(geometry = st_union(st_intersection(st_geometry(parc),
                                                        st_union(st_geometry(mer_sf)))))
resultats$surface_zone_km2 <- as.numeric(st_area(zone_etude)) / 1e6


# ---- 3. Observations : téléchargement, bancarisation, nettoyage -------------
# Objectif : constituer un jeu d'observations fiable et traçable.
# Méthode  : téléchargement OBIS et GBIF, harmonisation des colonnes (Darwin
#            Core), puis filtres successifs notés dans un journal : absences,
#            spécimens, échouages, doublons entre portails, dates, précision,
#            points à terre, agrégation à 1 observation par km² et par jour
#            (pseudo-réplication), sélection dans le Parc.
# Sorties  : obs_parc (observations retenues), journal (nombre restant par étape).

wkt <- st_as_text(emprise_vers_sf(BBOX))
fichier_brut <- "data/raw/observations_brutes.csv"

# Interrogation directe de l'API d'OBIS (utilisée si le package robis est absent)
api_obis <- function(taille_page = 5000) {
  base  <- "https://api.obis.org/v3/occurrence"
  pages <- list()
  apres <- NULL
  repeat {
    url <- paste0(base,
                  "?scientificname=", URLencode("Tursiops truncatus", reserved = TRUE),
                  "&geometry=", URLencode(wkt, reserved = TRUE),
                  "&size=", taille_page,
                  if (!is.null(apres)) paste0("&after=", apres))
    reponse <- jsonlite::fromJSON(url, flatten = TRUE)
    res <- reponse$results
    if (is.null(res) || length(res) == 0 || NROW(res) == 0) break
    # Tout en texte pour pouvoir empiler les pages sans conflit de types
    res <- as.data.frame(lapply(res, function(x) {
      if (is.list(x)) vapply(x, function(v) paste(v, collapse = ";"), character(1))
      else as.character(x)
    }), stringsAsFactors = FALSE)
    pages[[length(pages) + 1]] <- res
    message("  OBIS : ", sum(vapply(pages, nrow, integer(1))), " enregistrements reçus")
    if (nrow(res) < taille_page) break
    apres <- res$id[nrow(res)]
  }
  bind_rows(pages)
}

telecharger_obis <- function() {
  if (requireNamespace("robis", quietly = TRUE)) {
    d <- robis::occurrence(scientificname = "Tursiops truncatus", geometry = wkt)
  } else {
    message("Package robis absent : interrogation directe de l'API OBIS.")
    d <- api_obis()
  }
  if (is.null(d) || nrow(d) == 0) return(NULL)
  tibble(source        = "OBIS",
         id_source     = as.character(colonne(d, "id")),
         longitude     = as.numeric(colonne(d, "decimalLongitude")),
         latitude      = as.numeric(colonne(d, "decimalLatitude")),
         date          = as.Date(substr(as.character(colonne(d, "eventDate")), 1, 10),
                                 format = "%Y-%m-%d"),
         incertitude_m = as.numeric(colonne(d, "coordinateUncertaintyInMeters")),
         type_obs      = as.character(colonne(d, "basisOfRecord")),
         jeu_donnees   = as.character(colonne(d, "datasetName")),
         effectif      = as.numeric(colonne(d, "individualCount")),
         statut        = as.character(colonne(d, "occurrenceStatus")),
         absence       = as.character(colonne(d, "absence")))
}

telecharger_gbif <- function() {
  r <- rgbif::occ_search(scientificName = "Tursiops truncatus", geometry = wkt,
                         hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                         limit = 10000)
  d <- r$data
  if (is.null(d) || nrow(d) == 0) return(NULL)

  # Le GBIF identifie les jeux de données par un code (datasetKey) :
  # on demande leur titre à l'API du GBIF
  cles   <- unique(na.omit(as.character(colonne(d, "datasetKey"))))
  titres <- vapply(cles, function(k) {
    tryCatch(as.character(jsonlite::fromJSON(paste0("https://api.gbif.org/v1/dataset/", k))$title),
             error = function(e) NA_character_)
  }, character(1))
  nom_jeu <- as.character(colonne(d, "datasetName"))
  nom_jeu <- ifelse(is.na(nom_jeu), unname(titres[as.character(colonne(d, "datasetKey"))]), nom_jeu)
  nom_jeu <- ifelse(is.na(nom_jeu), as.character(colonne(d, "institutionCode")), nom_jeu)

  tibble(source        = "GBIF",
         id_source     = as.character(colonne(d, "key")),
         longitude     = as.numeric(colonne(d, "decimalLongitude")),
         latitude      = as.numeric(colonne(d, "decimalLatitude")),
         date          = as.Date(substr(as.character(colonne(d, "eventDate")), 1, 10),
                                 format = "%Y-%m-%d"),
         incertitude_m = as.numeric(colonne(d, "coordinateUncertaintyInMeters")),
         type_obs      = as.character(colonne(d, "basisOfRecord")),
         jeu_donnees   = nom_jeu,
         effectif      = as.numeric(colonne(d, "individualCount")),
         statut        = as.character(colonne(d, "occurrenceStatus")),
         absence       = NA_character_)
}

if (!file.exists(fichier_brut) || FORCER_TELECHARGEMENT) {
  message("Téléchargement des observations OBIS et GBIF...")
  obs_obis <- tryCatch(telecharger_obis(), error = function(e) {
    message("OBIS indisponible : ", conditionMessage(e)); NULL })
  obs_gbif <- tryCatch(telecharger_gbif(), error = function(e) {
    message("GBIF indisponible : ", conditionMessage(e)); NULL })
  obs_brutes <- bind_rows(obs_obis, obs_gbif)
  if (nrow(obs_brutes) > 0) readr::write_csv(obs_brutes, fichier_brut)
} else {
  obs_brutes <- readr::read_csv(fichier_brut, show_col_types = FALSE,
                                col_types = readr::cols(.default = "c")) |>
    mutate(across(c(longitude, latitude, incertitude_m, effectif), as.numeric),
           date = as.Date(date))
}

# Données complémentaires éventuelles (partenaires, observations personnelles…)
# Format attendu : colonnes longitude, latitude, date (AAAA-MM-JJ) en WGS84
fichier_compl <- "data/raw/observations_complementaires.csv"
if (file.exists(fichier_compl)) {
  compl <- readr::read_csv(fichier_compl, show_col_types = FALSE) |>
    transmute(source = "Complémentaire", id_source = as.character(row_number()),
              longitude = as.numeric(longitude), latitude = as.numeric(latitude),
              date = as.Date(as.character(date)), incertitude_m = NA_real_,
              type_obs = "HUMAN_OBSERVATION", jeu_donnees = "Données complémentaires",
              effectif = NA_real_, statut = "PRESENT", absence = NA_character_)
  obs_brutes <- bind_rows(obs_brutes, compl)
}

if (is.null(obs_brutes) || nrow(obs_brutes) == 0) {
  stop("Aucune observation récupérée. Vérifier la connexion internet ",
       "ou déposer un fichier data/raw/observations_complementaires.csv.")
}

# -- Nettoyage, avec un journal de chaque étape (traçabilité)
journal <- data.frame(etape = character(), n_restant = integer())
noter <- function(d, etape) {
  n <- nrow(d)  # calculer d'abord : force les étapes précédentes du pipe à s'exécuter
  journal <<- rbind(journal, data.frame(etape = etape, n_restant = n))
  d
}

obs <- obs_brutes |>
  noter("Données brutes (OBIS + GBIF + complémentaires)") |>
  # Harmonisation du type d'observation (HumanObservation -> HUMAN_OBSERVATION)
  mutate(type_obs = toupper(gsub("([a-z])([A-Z])", "\\1_\\2", type_obs)),
         annee    = as.integer(format(date, "%Y"))) |>
  filter(!is.na(longitude), !is.na(latitude),
         between(longitude, -180, 180), between(latitude, -90, 90)) |>
  noter("Coordonnées présentes et valides") |>
  # Données d'absence ("cherché mais pas vu") : inutiles pour un domaine vital
  filter(toupper(coalesce(statut, "")) != "ABSENT",
         toupper(coalesce(absence, "")) != "TRUE") |>
  noter("Données d'absence exclues") |>
  filter(is.na(type_obs) | !grepl("SPECIMEN|FOSSIL", type_obs)) |>
  noter("Spécimens de collection exclus (musées)") |>
  # Les réseaux d'échouages publient des animaux morts sur la côte comme de
  # simples observations : ce ne sont pas des dauphins vivants en mer
  filter(is.na(jeu_donnees) | !grepl("chouage|stranding", jeu_donnees, ignore.case = TRUE)) |>
  noter("Données d'échouages exclues") |>
  # OBIS republie une partie des données GBIF : doublons inter-sources
  mutate(cle_doublon = paste(round(longitude, 3), round(latitude, 3), date)) |>
  distinct(cle_doublon, .keep_all = TRUE) |>
  noter("Doublons supprimés (même position à 100 m près et même date)") |>
  filter(!is.na(annee), annee >= ANNEE_MIN) |>
  noter(paste0("Date connue et année >= ", ANNEE_MIN)) |>
  filter(is.na(incertitude_m) | incertitude_m <= INCERTITUDE_MAX_M) |>
  noter(paste0("Précision de localisation <= ", INCERTITUDE_MAX_M / 1000, " km"))

obs_sf <- st_as_sf(obs, coords = c("longitude", "latitude"), crs = 4326,
                   remove = FALSE) |>
  st_transform(CRS_METRIQUE)

obs_sf$profondeur <- extract(bathy_l93, vect(obs_sf))[, 2]
obs_sf <- obs_sf |>
  filter(!is.na(profondeur), profondeur < TOLERANCE_TERRE_M) |>
  noter("Points à terre exclus")

# Plusieurs positions du même groupe le même jour = pseudo-réplication.
# On garde une observation par cellule de 1 km et par jour.
xy <- st_coordinates(obs_sf)
obs_sf <- obs_sf |>
  mutate(cellule = paste(floor(xy[, 1] / TAILLE_CELLULE_M),
                         floor(xy[, 2] / TAILLE_CELLULE_M), date)) |>
  distinct(cellule, .keep_all = TRUE) |>
  noter("Agrégation : 1 observation par km² et par jour")

obs_zone <- obs_sf  # toutes les observations propres de l'emprise (contexte)
obs_parc <- obs_sf[lengths(st_intersects(obs_sf, zone_etude)) > 0, ] |>
  noter("Observations situées dans le Parc")

n_parc <- nrow(obs_parc)
message(n_parc, " observations retenues dans le Parc.")
if (n_parc < 5) stop("Trop peu d'observations dans le Parc pour estimer un domaine vital.")
if (n_parc < N_MIN_KDE) {
  warning("Seulement ", n_parc, " observations : domaine vital très incertain, ",
          "à présenter comme une illustration de la méthode.")
}

resultats$journal       <- journal
resultats$n_parc        <- n_parc
resultats$n_zone        <- nrow(obs_zone)
resultats$n_brut        <- nrow(obs_brutes)
resultats$par_source    <- count(st_drop_geometry(obs_parc), source)
resultats$par_jeu       <- st_drop_geometry(obs_parc) |>
  count(jeu_donnees, sort = TRUE)
resultats$annees        <- range(obs_parc$annee)
resultats$par_annee     <- count(st_drop_geometry(obs_parc), annee)


# ---- 4. Domaine vital ----------------------------------------------------------
# Objectif : estimer la zone utilisée habituellement par les grands dauphins.
# Méthode  : estimation par noyau (KDE, largeur de bande de référence href),
#            isoplèthes 50 % (cœur) et 95 % ; polygone convexe minimum (MCP 95 %)
#            pour comparaison ; analyse de sensibilité avec une largeur de bande
#            réduite de moitié.
# Sorties  : domaine_vital (polygones), surfaces en km², kde_volume (raster).

xy_parc <- st_coordinates(obs_parc)

# 4a. Estimation par noyau (KDE) — largeur de bande de référence href
sigma <- sqrt(0.5 * (var(xy_parc[, 1]) + var(xy_parc[, 2])))
href  <- sigma * n_parc^(-1 / 6)
h     <- href * FACTEUR_H

emp <- st_bbox(st_buffer(zone_etude, 5000))
# MASS::kde2d divise h par 4 : on passe 4*h pour obtenir un écart-type = h
kd <- MASS::kde2d(xy_parc[, 1], xy_parc[, 2], h = 4 * h, n = 300,
                  lims = c(emp[["xmin"]], emp[["xmax"]], emp[["ymin"]], emp[["ymax"]]))

dx <- diff(kd$x[1:2]); dy <- diff(kd$y[1:2])
kde <- rast(nrows = length(kd$y), ncols = length(kd$x),
            xmin = min(kd$x) - dx / 2, xmax = max(kd$x) + dx / 2,
            ymin = min(kd$y) - dy / 2, ymax = max(kd$y) + dy / 2,
            crs = CRS_METRIQUE)
values(kde) <- as.vector(kd$z[, ncol(kd$z):1])   # lignes du raster = y décroissants
kde <- mask(kde, vect(zone_etude))                # densité limitée à la mer du Parc
names(kde) <- "densite"

# Isoplèthes : on cumule la densité des cellules les plus fortes jusqu'à p %
v  <- values(kde, mat = FALSE)
ok <- !is.na(v)
ordre <- order(v[ok], decreasing = TRUE)
cumul <- numeric(sum(ok))
cumul[ordre] <- cumsum(v[ok][ordre]) / sum(v[ok])
vol <- rep(NA_real_, length(v)); vol[ok] <- cumul
kde_volume <- setValues(kde, vol)
names(kde_volume) <- "volume"

isoplethe <- function(p) {
  pol <- st_as_sf(as.polygons(ifel(kde_volume <= p, 1, NA)))
  g <- st_union(st_geometry(pol))
  st_sf(niveau = paste0("KDE ", p * 100, " %"),
        surface_km2 = as.numeric(st_area(g)) / 1e6,
        geometry = g)
}
iso50 <- isoplethe(0.50)
iso95 <- isoplethe(0.95)
domaine_vital <- rbind(iso50, iso95)

# 4b. Polygone convexe minimum (MCP 95 %) — méthode simple de comparaison
centre <- colMeans(xy_parc)
d_centre <- sqrt((xy_parc[, 1] - centre[1])^2 + (xy_parc[, 2] - centre[2])^2)
mcp95 <- obs_parc[d_centre <= quantile(d_centre, 0.95), ] |>
  st_union() |> st_convex_hull() |>
  st_intersection(st_geometry(zone_etude))
mcp95 <- st_sf(niveau = "MCP 95 %", surface_km2 = as.numeric(st_area(mcp95)) / 1e6,
               geometry = st_union(mcp95))

resultats$href_km <- href / 1000
resultats$domaine_vital <- rbind(domaine_vital, mcp95) |>
  st_drop_geometry() |>
  mutate(part_du_parc_pct = 100 * surface_km2 / resultats$surface_zone_km2)

# 4c. Analyse de sensibilité : la largeur de bande de référence (href) tend à
# surlisser. On recalcule le domaine vital avec une largeur de bande réduite de
# moitié et on vérifie si le cœur reste au même endroit.
FACTEUR_SENSIBILITE <- 0.5

volume_kde <- function(hh) {
  k <- MASS::kde2d(xy_parc[, 1], xy_parc[, 2], h = 4 * hh, n = 300,
                   lims = c(emp[["xmin"]], emp[["xmax"]], emp[["ymin"]], emp[["ymax"]]))
  r <- rast(kde)                                   # même grille que le KDE de référence
  values(r) <- as.vector(k$z[, ncol(k$z):1])
  r <- mask(r, vect(zone_etude))
  v <- values(r, mat = FALSE); ok <- !is.na(v)
  o <- order(v[ok], decreasing = TRUE)
  cu <- numeric(sum(ok)); cu[o] <- cumsum(v[ok][o]) / sum(v[ok])
  vol <- rep(NA_real_, length(v)); vol[ok] <- cu
  setValues(r, vol)
}
contour_iso <- function(vol, p) st_union(st_geometry(st_as_sf(as.polygons(ifel(vol <= p, 1, NA)))))
aire_km2 <- function(g) if (length(g) == 0 || all(st_is_empty(g))) 0 else as.numeric(sum(st_area(g))) / 1e6

vol_sens   <- volume_kde(FACTEUR_SENSIBILITE * href)
coeur_sens <- contour_iso(vol_sens, 0.50)
dv95_sens  <- contour_iso(vol_sens, 0.95)
coeur_ref  <- st_geometry(iso50)
recouvrement <- aire_km2(suppressWarnings(st_intersection(coeur_sens, coeur_ref)))

resultats$sensibilite_kde <- data.frame(
  facteur_h = FACTEUR_SENSIBILITE,
  h_km = FACTEUR_SENSIBILITE * href / 1000,
  surface_coeur_km2 = aire_km2(coeur_sens),
  surface_dv95_km2  = aire_km2(dv95_sens),
  # Part du cœur « sensibilité » située dans le cœur de référence
  part_coeur_dans_reference_pct = 100 * recouvrement / aire_km2(coeur_sens)
)
message("Sensibilité (h x ", FACTEUR_SENSIBILITE, ") : cœur = ",
        round(resultats$sensibilite_kde$surface_coeur_km2), " km², dont ",
        round(resultats$sensibilite_kde$part_coeur_dans_reference_pct), " % dans le cœur de référence.")


# ---- 5. Habitat : utilisation vs disponibilité -----------------------------
# Objectif : tester si les zones fréquentées diffèrent du milieu disponible.
# Méthode  : comparaison des observations à des points tirés au hasard dans le
#            Parc (plan utilisation / disponibilité) ; tests de Wilcoxon ;
#            corrélation de Spearman puis régression logistique.
# Sorties  : hab (tableau utilisation / disponibilité), tests, rapport de cotes.

# Points de "disponibilité" tirés au hasard dans la mer du Parc
dispo <- st_sf(geometry = st_sample(zone_etude, size = N_DISPO_PAR_OBS * n_parc,
                                    type = "random"))

hab <- rbind(
  data.frame(presence = 1, type = "Observations",
             extract(c(bathy_l93, dist_cote), vect(obs_parc))[, -1]),
  data.frame(presence = 0, type = "Disponible (aléatoire)",
             extract(c(bathy_l93, dist_cote), vect(dispo))[, -1])
) |>
  filter(!is.na(profondeur), !is.na(dist_cote)) |>
  mutate(prof_m = -profondeur, dist_cote_km = dist_cote / 1000)

resultats$habitat_resume <- hab |>
  group_by(type) |>
  summarise(n = n(),
            prof_mediane_m = median(prof_m),
            prof_q25 = quantile(prof_m, 0.25), prof_q75 = quantile(prof_m, 0.75),
            dist_mediane_km = median(dist_cote_km),
            .groups = "drop")

resultats$test_prof <- wilcox.test(prof_m ~ presence, data = hab)
resultats$test_dist <- wilcox.test(dist_cote_km ~ presence, data = hab)

# Modèle de sélection d'habitat (utilisation/disponibilité) — régression logistique
# On vérifie d'abord la corrélation entre les deux variables
resultats$correlation_prof_dist <- cor(hab$prof_m, hab$dist_cote_km, method = "spearman")
if (abs(resultats$correlation_prof_dist) > 0.7) {
  formule <- presence ~ scale(log1p(prof_m))
  resultats$note_glm <- "Profondeur et distance à la côte trop corrélées (|rho| > 0.7) : seule la profondeur est gardée."
} else {
  formule <- presence ~ scale(log1p(prof_m)) + scale(dist_cote_km)
  resultats$note_glm <- "Profondeur (log) et distance à la côte, centrées-réduites."
}
glm_hab <- glm(formule, data = hab, family = binomial)
coefs <- summary(glm_hab)$coefficients
ic <- suppressMessages(confint.default(glm_hab))
resultats$glm <- data.frame(
  variable   = rownames(coefs),
  odds_ratio = exp(coefs[, 1]),
  ic95_bas   = exp(ic[, 1]),
  ic95_haut  = exp(ic[, 2]),
  p_value    = coefs[, 4],
  row.names  = NULL
)


# ---- 6. Pressions : trafic maritime -----------------------------------------
# Objectif : mesurer le chevauchement entre le domaine vital et le trafic.
# Méthode  : densité de navires EMODnet (AIS, h/km²/mois) comparée entre le cœur,
#            l'anneau 50-95 % et le reste du Parc (médianes) ; part des cellules
#            de fort trafic (3e quartile, puis 50 % et 10 % en sensibilité)
#            situées dans le domaine vital.
# Sorties  : resultats$trafic (tableau par zone et par catégorie de navires).

fichiers_trafic <- list.files("data/raw/trafic", pattern = "\\.tif$",
                              full.names = TRUE, ignore.case = TRUE)
# Ordre logique : tous navires (sert pour la carte et la note), plaisance, pêche, autres
nb <- basename(fichiers_trafic)
rang <- ifelse(grepl("tous", nb, ignore.case = TRUE), 1,
               ifelse(grepl("plaisance", nb, ignore.case = TRUE), 2,
                      ifelse(grepl("peche|pêche", nb, ignore.case = TRUE), 3, 4)))
fichiers_trafic <- fichiers_trafic[order(rang, nb)]
resultats$trafic_disponible <- length(fichiers_trafic) > 0

# Zones : cœur (KDE 50 %), domaine vital (50-95 %), reste du Parc
zones <- classify(kde_volume, rcl = matrix(c(0,    0.50, 1,
                                             0.50, 0.95, 2,
                                             0.95, 1.01, 3), ncol = 3, byrow = TRUE),
                  include.lowest = TRUE)
etiquettes_zones <- c("Cœur (KDE 50 %)", "Domaine vital (50-95 %)", "Reste du Parc")

if (resultats$trafic_disponible) {
  trafic_resume <- list()
  trafic_valeurs <- list()
  trafic_rasters <- list()

  for (f in fichiers_trafic) {
    nom <- tools::file_path_sans_ext(basename(f))
    r <- rast(f)[[1]]
    # Les rasters EMODnet couvrent toute l'Europe : on découpe d'abord dans leur
    # projection d'origine (plus rapide), puis on reprojette en Lambert-93
    r <- crop(r, project(vect(st_buffer(zone_etude, 20000)), crs(r)))
    r <- project(r, CRS_METRIQUE, method = "bilinear")
    r <- mask(crop(r, vect(st_buffer(zone_etude, 5000))), vect(zone_etude))
    names(r) <- "trafic"
    trafic_rasters[[nom]] <- r

    z <- resample(zones, r, method = "near")
    df <- data.frame(zone = values(z, mat = FALSE), trafic = values(r, mat = FALSE)) |>
      filter(!is.na(zone), !is.na(trafic)) |>
      mutate(zone = factor(etiquettes_zones[zone], levels = etiquettes_zones),
             couche = nom)
    trafic_valeurs[[nom]] <- df

    # Zones de fort trafic : cellules au-dessus du 3e quartile des valeurs du Parc
    # (seuil relatif, recalculé pour chaque catégorie de navires)
    seuil <- unname(quantile(df$trafic, 0.75))
    fort  <- df$trafic > seuil
    part_surface_dv <- mean(df$zone != "Reste du Parc")
    # Si la majorité des cellules ont une densité nulle (cas de la pêche équipée
    # de l'AIS), le seuil est proche de zéro et l'indicateur n'a pas de sens
    part_nulle <- mean(df$trafic == 0)
    indicateur_ok <- part_nulle < 0.5 && seuil > 0
    part_fort_dans_dv <- if (indicateur_ok) mean(df$zone[fort] != "Reste du Parc") else NA_real_
    # Sensibilité au choix du seuil : 50 % et 10 % des cellules les plus fréquentées
    part_fort_seuil <- function(prob) {
      if (!indicateur_ok) return(NA_real_)
      s <- unname(quantile(df$trafic, prob))
      mean(df$zone[df$trafic > s] != "Reste du Parc")
    }
    part_fort_q50 <- part_fort_seuil(0.50)
    part_fort_q90 <- part_fort_seuil(0.90)

    trafic_resume[[nom]] <- df |>
      group_by(zone) |>
      summarise(trafic_moyen = mean(trafic), trafic_median = median(trafic),
                .groups = "drop") |>
      mutate(couche = nom,
             seuil_fort_trafic = seuil,
             part_cellules_nulles_pct = 100 * part_nulle,
             part_fort_trafic_dans_dv_pct = 100 * part_fort_dans_dv,
             part_fort_q50_dans_dv_pct = 100 * part_fort_q50,
             part_fort_q90_dans_dv_pct = 100 * part_fort_q90,
             part_surface_dv_pct = 100 * part_surface_dv)
  }
  resultats$trafic <- bind_rows(trafic_resume)
  trafic_df <- bind_rows(trafic_valeurs)
} else {
  message("Aucun raster de trafic dans data/raw/trafic/ : partie 'pressions' ignorée.")
}


# ---- 7. Cartes et figures ----------------------------------------------------
# Objectif : produire des cartes et graphiques prêts à diffuser.
# Méthode  : ggplot2 ; habillage cartographique commun (échelle, nord
#            géographique, repères, sources et auteur en légende).
# Sorties  : outputs/figures/ (carte du domaine vital, habitat, trafic).

AUTEUR <- "G. Helaine"
mois_fr <- c("janvier", "février", "mars", "avril", "mai", "juin", "juillet",
             "août", "septembre", "octobre", "novembre", "décembre")
date_fr <- paste(mois_fr[as.integer(format(Sys.Date(), "%m"))], format(Sys.Date(), "%Y"))
signature <- paste0("Projection Lambert-93 (EPSG:2154). Réalisation : ", AUTEUR, ", ", date_fr, ".")
virgule <- function(x, d = 1) format(round(x, d), decimal.mark = ",", nsmall = d)

theme_carte <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_line(colour = "grey88", linewidth = 0.2),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
        axis.title = element_blank(),
        axis.text = element_text(size = 7, colour = "grey30"),
        legend.position = "right",
        legend.title = element_text(size = 8.5, face = "bold"),
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(0.45, "cm"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        plot.caption = element_text(colour = "grey40", size = 6.5, hjust = 0, lineheight = 1.1),
        plot.title.position = "plot",
        plot.caption.position = "plot")

lim <- st_bbox(st_buffer(zone_etude, 3000))
largeur <- lim[["xmax"]] - lim[["xmin"]]
hauteur <- lim[["ymax"]] - lim[["ymin"]]
emp_carte <- ext(lim[["xmin"]], lim[["xmax"]], lim[["ymin"]], lim[["ymax"]])
bathy_df <- as.data.frame(mask(crop(bathy_l93, emp_carte), crop(mer, emp_carte)), xy = TRUE)
mention_parc <- if (resultats$parc_officiel) {
  "périmètre du Parc : PatriNat, base des espaces protégés (07/2026)"
} else {
  "périmètre du Parc APPROXIMATIF (couche officielle non chargée)"
}

# Part du principal jeu de données (pour rappeler le biais d'effort sur la carte)
pct_jeu1 <- round(100 * resultats$par_jeu$n[1] / sum(resultats$par_jeu$n))

# -- Repères géographiques
lieux <- data.frame(
  nom   = c("Saint-Florent", "Bastia", "Cap\nCorse", "Désert\ndes Agriates"),
  type  = c("ville", "ville", "terre", "terre"),
  lon   = c(9.303, 9.450, 9.425, 9.17),
  lat   = c(42.681, 42.697, 42.93, 42.70),
  dx    = c(1500, 1800, 0, 0),
  dy    = c(-1800, 0, 0, 0),
  hjust = c(0, 0, 0.5, 0.5)
)
xy_lieux <- st_coordinates(st_transform(st_as_sf(lieux, coords = c("lon", "lat"), crs = 4326),
                                        CRS_METRIQUE))
lieux$x <- xy_lieux[, 1]
lieux$y <- xy_lieux[, 2]
lieux <- lieux[lieux$x > lim[["xmin"]] & lieux$x < lim[["xmax"]] &
               lieux$y > lim[["ymin"]] & lieux$y < lim[["ymax"]], ]

couches_lieux <- function() {
  list(
    geom_point(data = lieux[lieux$type == "ville", ], aes(x, y), shape = 22, size = 1.3,
               fill = "grey10", colour = "white", stroke = 0.3),
    geom_text(data = lieux[lieux$type == "ville", ], aes(x + dx, y + dy, label = nom, hjust = hjust),
              size = 2.5, colour = "grey10", vjust = 0.5, lineheight = 0.9),
    geom_text(data = lieux[lieux$type == "terre", ], aes(x + dx, y + dy, label = nom),
              size = 2.2, colour = "grey45", fontface = "italic", lineheight = 0.85)
  )
}

# -- Barre d'échelle et flèche du nord (orientée sur le nord géographique :
#    en Lambert-93, le nord vrai est décalé de plusieurs degrés en Corse)
centre_ll <- st_coordinates(st_transform(st_centroid(st_as_sfc(lim)), 4326))[1, 1:2]
p_nord <- st_coordinates(st_transform(st_sfc(st_point(centre_ll), st_point(centre_ll + c(0, 0.1)),
                                             crs = 4326), CRS_METRIQUE))
v_nord <- p_nord[2, 1:2] - p_nord[1, 1:2]
u_nord <- v_nord / sqrt(sum(v_nord^2))

habillage <- function() {
  x0 <- lim[["xmin"]] + 0.04 * largeur
  y0 <- lim[["ymin"]] + 0.06 * hauteur
  L  <- 20000
  h  <- 0.013 * hauteur
  xn <- lim[["xmax"]] - 0.05 * largeur
  yn <- lim[["ymax"]] - 0.17 * hauteur
  ln <- 0.09 * hauteur
  list(
    annotate("rect", xmin = x0 - 0.015 * largeur, xmax = x0 + L + 0.06 * largeur,
             ymin = y0 - 0.045 * hauteur, ymax = y0 + h + 0.02 * hauteur,
             fill = "white", alpha = 0.85, colour = NA),
    annotate("rect", xmin = x0, xmax = x0 + L / 2, ymin = y0, ymax = y0 + h,
             fill = "grey15", colour = "grey15", linewidth = 0.2),
    annotate("rect", xmin = x0 + L / 2, xmax = x0 + L, ymin = y0, ymax = y0 + h,
             fill = "white", colour = "grey15", linewidth = 0.2),
    annotate("text", x = c(x0, x0 + L / 2, x0 + L), y = y0 - 0.008 * hauteur,
             label = c("0", "10", "20 km"), vjust = 1, size = 2.3),
    annotate("segment", x = xn, y = yn, xend = xn + ln * u_nord[1], yend = yn + ln * u_nord[2],
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             linewidth = 1.8, colour = "white"),
    annotate("segment", x = xn, y = yn, xend = xn + ln * u_nord[1], yend = yn + ln * u_nord[2],
             arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
             linewidth = 0.6, colour = "grey10"),
    annotate("point", x = xn + 1.3 * ln * u_nord[1], y = yn + 1.3 * ln * u_nord[2],
             shape = 21, size = 4.2, fill = "white", colour = NA),
    annotate("text", x = xn + 1.3 * ln * u_nord[1], y = yn + 1.3 * ln * u_nord[2],
             label = "N", fontface = "bold", size = 3)
  )
}

cadre <- coord_sf(crs = st_crs(CRS_METRIQUE), xlim = unname(lim[c(1, 3)]),
                  ylim = unname(lim[c(2, 4)]), expand = FALSE)

titre_espece <- function(debut) {
  bquote(bold(.(debut) ~ "du grand dauphin (") * bolditalic("Tursiops truncatus") * bold(")"))
}

etiquettes_dv <- c("KDE 50 %" = "Cœur (isoplèthe 50 %)",
                   "KDE 95 %" = "Domaine vital (isoplèthe 95 %)")

# -- Carte 1 : domaine vital
obs_hors <- obs_zone[lengths(st_intersects(obs_zone, zone_etude)) == 0, ]
lab_ret  <- paste0("Dans le Parc, analysées (n = ", n_parc, ")")
lab_hors <- paste0("Hors du Parc, non analysées (n = ", nrow(obs_hors), ")")

carte_dv <- ggplot() +
  geom_raster(data = bathy_df, aes(x, y, fill = -profondeur)) +
  scale_fill_gradient(name = "Profondeur (m)", low = "#e3eef6", high = "#08306b",
                      trans = "sqrt", breaks = c(100, 500, 1000, 2000)) +
  geom_contour(data = bathy_df, aes(x, y, z = -profondeur, linetype = "Isobathe 200 m"),
               breaks = 200, colour = "grey55", linewidth = 0.3) +
  geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey70", linewidth = 0.15) +
  geom_sf(data = parc, aes(linetype = "Limite du Parc"), fill = NA, colour = "grey15",
          linewidth = 0.45, key_glyph = "path") +
  geom_sf(data = domaine_vital, aes(colour = niveau), fill = NA, linewidth = 0.8,
          key_glyph = "path") +
  geom_sf(data = obs_hors, aes(shape = lab_hors), size = 0.8, colour = "grey35",
          fill = "grey85", stroke = 0.3) +
  geom_sf(data = obs_parc, aes(shape = lab_ret), size = 1, colour = "black",
          fill = "white", stroke = 0.3) +
  scale_colour_manual(name = "Domaine vital (KDE)", values = c("KDE 50 %" = "#d7301f",
                                                               "KDE 95 %" = "#fc8d59"),
                      labels = etiquettes_dv) +
  scale_shape_manual(name = "Observations", values = setNames(c(21, 21), c(lab_ret, lab_hors)),
                     breaks = c(lab_ret, lab_hors)) +
  scale_linetype_manual(name = "Repères", values = c("Limite du Parc" = "dashed",
                                                     "Isobathe 200 m" = "solid"),
                        breaks = c("Limite du Parc", "Isobathe 200 m")) +
  guides(colour   = guide_legend(order = 1),
         shape    = guide_legend(order = 2, override.aes = list(fill = c("white", "grey85"),
                                                                colour = c("black", "grey35"),
                                                                size = 1.8)),
         linetype = guide_legend(order = 3, override.aes = list(colour = c("grey15", "grey55"))),
         fill     = guide_colourbar(order = 4, reverse = TRUE)) +
  couches_lieux() +
  habillage() +
  cadre +
  labs(title = titre_espece("Domaine vital estimé"),
       subtitle = paste0("Parc naturel marin du Cap Corse et de l'Agriate\nEstimation par noyau ",
                         "(largeur de bande href = ", virgule(h / 1000), " km) — ",
                         n_parc, " observations, ", resultats$annees[1], "–", resultats$annees[2]),
       caption = paste0("Sources : observations OBIS et GBIF (", pct_jeu1,
                        " % issues d'un même jeu de données opportunistes : la répartition reflète ",
                        "aussi l'effort d'observation) ;\nbathymétrie NOAA ETOPO 2022 (~", round(resultats$resolution_bathy_m, -1), " m) ; ",
                        mention_parc, ".\n", signature)) +
  theme_carte
ggsave("outputs/figures/carte_domaine_vital.png", carte_dv, width = 9, height = 6.4, dpi = 300,
       bg = "white")

# -- Figure : habitat utilisé vs disponible
fp_txt <- function(p) if (p < 0.001) "p < 0,001" else paste0("p = ", virgule(p, 3))
lab_obs <- paste0("Observations\n(n = ", sum(hab$presence == 1), ")")
lab_dis <- paste0("Milieu disponible\n(n = ", sum(hab$presence == 0), " points aléatoires)")
hab$type_lab <- factor(ifelse(hab$presence == 1, lab_obs, lab_dis), levels = c(lab_obs, lab_dis))
hab_long <- rbind(
  data.frame(type = hab$type_lab, variable = "Profondeur (m)", valeur = hab$prof_m),
  data.frame(type = hab$type_lab, variable = "Distance à la côte (km)", valeur = hab$dist_cote_km)
)
hab_long$variable <- factor(hab_long$variable, levels = c("Profondeur (m)", "Distance à la côte (km)"))

fig_hab <- ggplot(hab_long, aes(x = type, y = valeur, fill = type)) +
  geom_violin(alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ variable, scales = "free_y") +
  scale_y_sqrt() +
  scale_fill_manual(values = setNames(c("#d7301f", "#6baed6"), c(lab_obs, lab_dis)),
                    guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Profondeur et distance à la côte : observations et milieu disponible",
       subtitle = paste0("Tests de Wilcoxon — profondeur : ", fp_txt(resultats$test_prof$p.value),
                         " ; distance à la côte : ", fp_txt(resultats$test_dist$p.value)),
       caption = paste0("Milieu disponible : points tirés au hasard dans la partie marine du Parc (",
                        N_DISPO_PAR_OBS, " par observation). Axes en racine carrée.\n",
                        "Réalisation : ", AUTEUR, ", ", date_fr, ".")) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        plot.caption = element_text(colour = "grey40", size = 6.5, hjust = 0),
        plot.title.position = "plot", plot.caption.position = "plot",
        strip.text = element_text(face = "bold", size = 9.5))
ggsave("outputs/figures/habitat_profondeur_cote.png", fig_hab, width = 8, height = 4.4,
       dpi = 300, bg = "white")

# -- Trafic maritime
libelle_couche <- function(x) {
  vapply(x, function(n) {
    annee <- regmatches(n, regexpr("[0-9]{4}", n))
    type  <- sub("_?[0-9]{4}$", "", n)
    noms  <- c(tous_navires = "Tous navires", plaisance = "Plaisance", peche = "Pêche")
    type  <- if (type %in% names(noms)) noms[[type]] else gsub("_", " ", type)
    if (length(annee) == 1) paste0(type, " (", annee, ")") else type
  }, character(1), USE.NAMES = FALSE)
}

if (resultats$trafic_disponible) {
  nom1 <- names(trafic_rasters)[1]
  tr_df <- as.data.frame(trafic_rasters[[nom1]], xy = TRUE)

  carte_trafic <- ggplot() +
    geom_raster(data = tr_df, aes(x, y, fill = trafic)) +
    scale_fill_viridis_c(name = "Densité de navires\n(h / km² / mois)", trans = "log1p",
                         option = "magma", direction = -1,
                         breaks = c(0, 1, 10, 100, 1000)) +
    geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey70", linewidth = 0.15) +
    geom_sf(data = parc, aes(linetype = "Limite du Parc"), fill = NA, colour = "grey15",
            linewidth = 0.45, key_glyph = "path") +
    geom_sf(data = domaine_vital, aes(colour = niveau), fill = NA, linewidth = 0.8,
            key_glyph = "path") +
    scale_colour_manual(name = "Domaine vital (KDE)", values = c("KDE 50 %" = "#00a6a6",
                                                                 "KDE 95 %" = "#7fd4d4"),
                        labels = etiquettes_dv) +
    scale_linetype_manual(name = "Repères", values = c("Limite du Parc" = "dashed")) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2),
           fill = guide_colourbar(order = 3)) +
    couches_lieux() +
    habillage() +
    cadre +
    labs(title = titre_espece("Trafic maritime et domaine vital"),
         subtitle = paste0("Densité de navires équipés de l'AIS — ", libelle_couche(nom1),
                           ", moyenne annuelle, grille de 1 km"),
         caption = paste0("Sources : EMODnet Human Activities, densité de navires (AIS) ; ",
                          "observations OBIS et GBIF ; ", mention_parc, ".\n",
                          "L'AIS n'équipe pas la plupart des petites embarcations de plaisance ",
                          "et de pêche : le trafic côtier est sous-estimé.\n", signature)) +
    theme_carte
  ggsave("outputs/figures/carte_trafic.png", carte_trafic, width = 9, height = 6.4, dpi = 300,
         bg = "white")

  trafic_df$couche_lab <- factor(libelle_couche(trafic_df$couche),
                                 levels = libelle_couche(names(trafic_rasters)))
  fig_trafic <- ggplot(trafic_df, aes(x = zone, y = trafic, fill = zone)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
    facet_wrap(~ couche_lab, scales = "free_y") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000)) +
    scale_x_discrete(labels = c("Cœur\n(KDE 50 %)", "Domaine vital\n(50–95 %)", "Reste\ndu Parc")) +
    scale_fill_manual(values = c("#d7301f", "#fc8d59", "#bdbdbd"), guide = "none") +
    labs(x = NULL, y = "Densité de navires (h / km² / mois, échelle log)",
         title = "Intensité du trafic maritime selon la zone du domaine vital",
         subtitle = "Une valeur par cellule de 1 km² ; moyenne annuelle EMODnet (navires équipés de l'AIS)",
         caption = paste0("Le trafic des petites embarcations, souvent sans AIS, n'est pas représenté. ",
                          "Réalisation : ", AUTEUR, ", ", date_fr, ".")) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          plot.caption = element_text(colour = "grey40", size = 6.5, hjust = 0),
          plot.title.position = "plot", plot.caption.position = "plot",
          strip.text = element_text(face = "bold", size = 9.5),
          axis.text.x = element_text(size = 7.5))
  ggsave("outputs/figures/trafic_par_zone.png", fig_trafic, width = 8.5, height = 4.4,
         dpi = 300, bg = "white")
}


# ---- 8. Exports ----------------------------------------------------------------
# Objectif : rendre l'analyse vérifiable et réutilisable.
# Sorties  : données nettoyées, journal et dictionnaire (data/propre/) ;
#            couches SIG pour QGIS (outputs/sig/) ; tableaux de résultats
#            (outputs/tableaux/) ; résultats pour la note (outputs/resultats.rds) ;
#            versions des logiciels et dates des données sources (outputs/).

# 8a. Données propres + dictionnaire + journal de nettoyage
obs_export <- st_drop_geometry(obs_parc) |>
  select(source, id_source, longitude, latitude, date, annee, profondeur,
         incertitude_m, type_obs, jeu_donnees, effectif)
readr::write_csv(obs_export, "data/propre/observations_propres_parc.csv")
readr::write_csv(journal, "data/propre/journal_nettoyage.csv")

dictionnaire <- data.frame(
  champ = c("source", "id_source", "longitude", "latitude", "date", "annee",
            "profondeur", "incertitude_m", "type_obs", "jeu_donnees", "effectif"),
  type = c("texte", "texte", "réel", "réel", "date", "entier", "réel", "réel",
           "texte", "texte", "réel"),
  unite = c("", "", "degrés (WGS84)", "degrés (WGS84)", "AAAA-MM-JJ", "", "m (négatif en mer)",
            "m", "", "", "individus"),
  description = c("Portail d'origine (OBIS, GBIF, Complémentaire)",
                  "Identifiant de l'enregistrement dans la source",
                  "Longitude de l'observation", "Latitude de l'observation",
                  "Date de l'observation", "Année de l'observation",
                  "Profondeur extraite de la bathymétrie au point",
                  "Incertitude de localisation déclarée par la source",
                  "Type d'enregistrement (Darwin Core basisOfRecord)",
                  "Nom du jeu de données d'origine",
                  "Nombre d'individus déclaré (souvent non renseigné)")
)
readr::write_csv(dictionnaire, "data/propre/dictionnaire_donnees.csv")

# 8b. Couches SIG pour QGIS (un seul GeoPackage + rasters)
gpkg <- "outputs/sig/delfinu_exercice.gpkg"
# Supprimer l'ancien GeoPackage et ses fichiers annexes SQLite (-wal, -shm, -journal)
unlink(list.files("outputs/sig", pattern = "^delfinu_exercice\\.gpkg", full.names = TRUE))
st_write(zone_etude, gpkg, layer = "zone_etude", quiet = TRUE)
st_write(obs_parc[, c("source", "date", "annee", "profondeur", "jeu_donnees")],
         gpkg, layer = "observations_parc", quiet = TRUE)
st_write(obs_zone[, c("source", "date", "annee", "profondeur", "jeu_donnees")],
         gpkg, layer = "observations_emprise", quiet = TRUE)
st_write(domaine_vital, gpkg, layer = "domaine_vital_kde", quiet = TRUE)
st_write(mcp95, gpkg, layer = "domaine_vital_mcp95", quiet = TRUE)
writeRaster(kde_volume, "outputs/sig/kde_volume.tif", overwrite = TRUE)
writeRaster(bathy_l93, "outputs/sig/bathymetrie_l93.tif", overwrite = TRUE)
writeRaster(dist_cote, "outputs/sig/distance_cote_l93.tif", overwrite = TRUE)

# 8c. Tableaux
readr::write_csv(resultats$domaine_vital, "outputs/tableaux/surfaces_domaine_vital.csv")
readr::write_csv(resultats$sensibilite_kde, "outputs/tableaux/sensibilite_domaine_vital.csv")
readr::write_csv(resultats$habitat_resume, "outputs/tableaux/habitat_resume.csv")
readr::write_csv(resultats$glm, "outputs/tableaux/modele_habitat.csv")
if (resultats$trafic_disponible) {
  readr::write_csv(resultats$trafic, "outputs/tableaux/trafic_par_zone.csv")
}

# 8d. Résultats pour la note
saveRDS(resultats, "outputs/resultats.rds")
writeLines(capture.output(sessionInfo()), "outputs/session_info.txt")

# Versions exactes de R et des packages utilisés
versions <- data.frame(
  logiciel = c("R", paquets),
  version  = c(paste(R.version$major, R.version$minor, sep = "."),
               vapply(paquets, function(p) tryCatch(as.character(packageVersion(p)),
                                                    error = function(e) NA_character_),
                      character(1)))
)
readr::write_csv(versions, "outputs/versions_packages.csv")

# Traçabilité des données sources : fichier, date de téléchargement, taille
fichiers_sources <- list.files("data/raw", recursive = TRUE, full.names = TRUE)
sources <- data.frame(
  fichier = fichiers_sources,
  date_fichier = format(file.info(fichiers_sources)$mtime, "%Y-%m-%d %H:%M"),
  taille_Mo = round(file.info(fichiers_sources)$size / 1e6, 2)
)
readr::write_csv(sources, "outputs/sources_donnees.csv")

message("\nTerminé. Résultats dans outputs/ — lancer ensuite le rendu de note.qmd.")
print(resultats$journal)
print(resultats$domaine_vital)
