# =============================================================================
#  Grand dauphin (Tursiops truncatus) dans le Parc naturel marin d'Iroise :
#  zones de concentration des observations et chevauchement avec le trafic
#  maritime (hors pêche professionnelle)
#  Les paramètres propres au site sont regroupés dans le bloc « Site étudié ».
#
#  Analyse exploratoire à partir de données ouvertes
#  Auteure : Gwladys Helaine — octobre 2026
#  Code écrit avec l'aide d'un assistant d'IA (Claude, Anthropic) ; démarche,
#  contrôle des données et des résultats, choix de méthode et d'interprétation :
#  Gwladys Helaine.
#
#  Étapes :
#   0. Paramètres et packages
#   1. Périmètre d'étude (Parc)
#   2. Bathymétrie et distance à la côte
#   3. Observations : téléchargement OBIS + GBIF, bancarisation, nettoyage et
#      contrôle qualité (positions sur grille, doublons imprécis)
#   4. Zones de concentration : noyau (KDE 50 % / 95 %) sur l'ensemble du Parc,
#      puis par secteur (Molène, Sein) ; tests de sensibilité
#   5. Habitat : comparaison utilisation / disponibilité (profondeur, côte)
#   6. Pressions : trafic maritime (EMODnet) dans et hors des zones, avec
#      contrôle de l'effet de la proximité de la côte
#   7. Cartes et figures
#   8. Exports (données propres, couches SIG, tableaux, résultats)
#
#  Lancer depuis le dossier du projet (ouvrir delfinu-iroise.Rproj ou
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


# ---- Site étudié ------------------------------------------------------------
# Tout ce qui est propre au site est réuni ici : pour analyser un autre parc,
# il suffit de modifier ce bloc.

SITE <- list(
  nom_long = "Parc naturel marin d'Iroise",
  motif    = "iroise",                        # mot cherché dans la base PatriNat
  # Emprise de téléchargement (WGS84), élargie automatiquement si besoin
  bbox     = c(xmin = -5.7, ymin = 47.85, xmax = -4.2, ymax = 48.7),
  # Emprise approximative du Parc, utilisée seulement en l'absence de couche officielle
  bbox_parc_approx = c(xmin = -5.45, ymin = 47.95, xmax = -4.35, ymax = 48.6),
  # Isobathe repère et graduations de profondeur pour les cartes (mer d'Iroise
  # peu profonde : une isobathe de 200 m n'y serait pas pertinente)
  isobathe = 50,
  breaks_profondeur = c(10, 25, 50, 100),
  # Repères géographiques affichés sur les cartes
  lieux = data.frame(
    # type : "ville" (port, carré noir), "terre" (île ou cap), "mer" (passage)
    nom   = c("Île de Sein", "Le Conquet", "Douarnenez", "Camaret-sur-Mer",
              "Ouessant", "Archipel\nde Molène", "Presqu'île\nde Crozon", "Cap Sizun",
              "Raz de Sein", "Chenal du Four", "Phare\nd'Ar Men"),
    # type « phare » : repère ponctuel en mer (triangle), étiquette décalée au-dessus
    type  = c("ville", "ville", "ville", "ville", "terre", "terre", "terre", "terre",
              "mer", "mer", "phare"),
    # Coordonnées WGS84 : bourgs et ports pour les villes, centre approximatif pour
    # les îles, presqu'îles et passages ; Ar Men : 48°03'00" N, 4°59'54" O
    lon   = c(-4.851, -4.771, -4.329, -4.595, -5.075, -4.958, -4.45, -4.60, -4.775, -4.80, -4.9983),
    lat   = c(48.038, 48.360, 48.093, 48.278, 48.462, 48.397, 48.235, 48.035, 48.045, 48.45, 48.0500),
    dx    = c(0, 1500, -1500, -1500, 0, 0, 0, 0, 0, 0, 0),
    dy    = c(-1800, 0, 1500, 0, 0, 3500, 0, 2500, 0, 0, 2200),
    hjust = c(0.5, 0, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  ),
  # Secteurs analysés séparément (section 4d). Les suivis par photo-identification
  # décrivent deux groupes côtiers résidents : archipel de Molène (GECC, 2021) et
  # chaussée de Sein, jusqu'au phare d'Ar Men, où des dauphins de passage sont signalés (Buanic, OFB, 2021-2022). Les boîtes
  # ont été tracées d'après ces suivis, après une première estimation sur tout le
  # Parc qui montrait deux concentrations distinctes : elles servent à séparer
  # les deux groupes, ce ne sont pas des limites biologiques. Leur influence est
  # testée en section 4e (boîtes réduites et agrandies de 5 km).
  secteurs = data.frame(
    nom     = c("Molène", "Sein"),
    fichier = c("molene", "sein"),
    libelle = c("Archipel de Molène – Le Conquet", "Chaussée de Sein – Raz de Sein"),
    xmin = c(-5.20, -5.10), xmax = c(-4.60, -4.68),
    ymin = c(48.22, 47.95), ymax = c(48.52, 48.12)
  ),
  # Catégories de trafic exclues de l'analyse (la mission exclut la pêche
  # professionnelle) : motif cherché dans le nom des fichiers de trafic
  trafic_exclu = "peche|pêche"
)


# ---- 0. Paramètres et packages ---------------------------------------------
# Objectif : regrouper tous les choix de l'analyse (seuils de nettoyage,
#            résolutions, réglages des méthodes) et charger les outils.
# Pourquoi : un paramètre modifié ici se répercute partout ; on peut ainsi
#            tester la sensibilité des résultats sans toucher au reste du code.
# Sorties  : paramètres, packages chargés, dossiers du projet créés.

ANNEE_MIN          <- 2000   # observations plus anciennes écartées
INCERTITUDE_MAX_M  <- 5000   # précision de localisation minimale acceptée
TOLERANCE_TERRE_M  <- 20     # altitude max tolérée : grille bathy grossière près des côtes
                             # et estran découvert à marée basse (marnage de plusieurs
                             # mètres en Iroise)
PAS_GRILLE_DEG     <- 1 / 12 # 5 minutes d'arc : positions dont la longitude et la
                             # latitude tombent pile sur ce pas = centres de maille
                             # d'un atlas, pas des positions réelles
INCERTITUDE_FLOUTEE_M <- 5000 # précision annoncée à partir de laquelle une donnée est
                             # considérée comme peu précise (doublon possible d'un autre jeu)
DISTANCE_DOUBLON_FLOU_M <- 10000 # rayon de recherche d'une observation plus précise le même jour
TAILLE_CELLULE_M   <- 1000   # agrégation : 1 obs. par cellule de 1 km et par jour
RESOLUTION_BATHY_MIN <- 0.25  # résolution de la bathymétrie NOAA, en minutes d'arc (0.25 = 15")
RESOLUTION_M       <- 250    # résolution des rasters reprojetés (500 si repli sur 1 minute)
N_MIN_KDE          <- 30     # en dessous, les zones estimées sont très incertaines
FACTEUR_H          <- 1      # multiplicateur de la largeur de bande href (1 = href)
N_MIN_SECTEUR      <- 10     # minimum d'observations pour estimer les zones d'un secteur
N_DISPO_PAR_OBS    <- 10     # points de disponibilité tirés par observation
CRS_METRIQUE       <- "EPSG:2154"  # Lambert-93 (France métropolitaine)
FORCER_TELECHARGEMENT <- FALSE

# Emprises (définies dans le bloc « Site étudié »)
BBOX <- SITE$bbox
BBOX_PARC_APPROX <- SITE$bbox_parc_approx

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
  motif     <- SITE$motif
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
  # une ligne "Parc naturel marin ..." portant le nom du site
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
    stop(SITE$nom_long, " introuvable dans ", f_parc, ". Couches disponibles : ",
         paste(couches, collapse = ", "))
  }
  parc <- st_read(f_parc, layer = couche, quiet = TRUE)
  message("Couche utilisée pour le Parc : ", couche)

  # On garde les lignes qui citent le nom du site ; dans une couche de parcs marins,
  # le nom seul suffit, sinon on exige aussi la mention "parc naturel marin"
  t_parc <- texte_lignes(parc)
  garde  <- grepl(motif, t_parc, ignore.case = TRUE)
  if (!grepl("pnm", couche, ignore.case = TRUE)) {
    garde <- garde & grepl(motif_pnm, t_parc, ignore.case = TRUE)
  }
  parc <- parc[garde, ]
  message("Polygone(s) retenu(s) pour le Parc :\n  ",
          paste(substr(t_parc[garde], 1, 150), collapse = "\n  "))
  if (nrow(parc) == 0) stop("Aucun polygone du ", SITE$nom_long, " trouvé dans la couche ", couche)
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
#            positions ramenées sur une grille d'atlas, points à terre,
#            doublons imprécis entre jeux de données, agrégation à 1 observation par km² et par jour
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
  # Données d'absence ("cherché mais pas vu") : inutiles pour localiser l'espèce
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

# -- Positions ramenées sur une grille (centres de maille)
# Certains jeux anciens ne donnent pas la position de l'observation mais le centre
# d'une maille d'atlas (ex. 4°50' O, 48°25' N). Plusieurs dizaines d'observations
# se retrouvent alors empilées sur le même point, et la précision déclarée (50 m)
# ne correspond pas à la position réelle. Un point GPS réel tombe très rarement pile sur ce pas pour la
# longitude et la latitude : on peut donc les repérer automatiquement.
sur_grille <- function(v, pas) abs(v / pas - round(v / pas)) < 1e-3
obs$sur_grille <- sur_grille(obs$longitude, PAS_GRILLE_DEG) &
                  sur_grille(obs$latitude,  PAS_GRILLE_DEG)
resultats$exclus_grille <- obs |>
  filter(sur_grille) |>
  count(jeu_donnees, longitude, latitude, name = "n", sort = TRUE)
obs <- obs |>
  filter(!sur_grille) |>
  noter("Positions sur une grille de 5' exclues (centres de maille d'atlas)")

obs_sf <- st_as_sf(obs, coords = c("longitude", "latitude"), crs = 4326,
                   remove = FALSE) |>
  st_transform(CRS_METRIQUE)

obs_sf$profondeur <- extract(bathy_l93, vect(obs_sf))[, 2]
obs_sf <- obs_sf |>
  filter(!is.na(profondeur), profondeur < TOLERANCE_TERRE_M) |>
  noter("Points à terre exclus")

# -- Doublons imprécis entre jeux de données
# Certains inventaires nationaux republient des observations d'autres jeux avec
# une précision de 5 km : même date, position
# décalée de quelques km. Le filtre de doublons exact ne les voit pas. Une donnée
# de précision >= 5 km est retirée si une observation plus précise (ou de précision non déclarée) d'un autre jeu existe le
# même jour à moins de 10 km : on garde la version la plus précise.
idx_flou <- which(!is.na(obs_sf$incertitude_m) &
                  obs_sf$incertitude_m >= INCERTITUDE_FLOUTEE_M)
a_retirer <- integer(0)
if (length(idx_flou) > 0) {
  voisins <- st_is_within_distance(obs_sf[idx_flou, ], obs_sf,
                                   dist = DISTANCE_DOUBLON_FLOU_M)
  jumeau <- vapply(seq_along(idx_flou), function(k) {
    i <- idx_flou[k]
    j <- setdiff(voisins[[k]], i)
    any(obs_sf$date[j] == obs_sf$date[i] &
        coalesce(obs_sf$jeu_donnees[j], "") != coalesce(obs_sf$jeu_donnees[i], "") &
        coalesce(obs_sf$incertitude_m[j], 0) < obs_sf$incertitude_m[i])
  }, logical(1))
  a_retirer <- idx_flou[jumeau]
}
resultats$exclus_floutes <- count(st_drop_geometry(obs_sf[a_retirer, ]), jeu_donnees, name = "n")
if (length(a_retirer) > 0) obs_sf <- obs_sf[-a_retirer, ]
obs_sf <- noter(obs_sf, "Doublons imprécis exclus (même jour, < 10 km, autre jeu plus précis)")

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
if (n_parc < 5) stop("Trop peu d'observations dans le Parc pour estimer des zones de concentration.")
if (n_parc < N_MIN_KDE) {
  warning("Seulement ", n_parc, " observations : zones de concentration très incertaines, ",
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


# ---- 4. Zones de concentration des observations -------------------------------
# Objectif : estimer la zone utilisée habituellement par les grands dauphins.
# Méthode  : estimation par noyau (KDE, largeur de bande de référence href),
#            isoplèthes 50 % (zone centrale) et 95 % (zone étendue) ; polygone convexe minimum (MCP 95 %)
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
# surlisser. On recalcule les zones avec une largeur de bande réduite de
# moitié et on vérifie si la zone centrale reste au même endroit.
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
  # Part de la zone centrale « sensibilité » située dans la zone centrale de référence
  part_coeur_dans_reference_pct = 100 * recouvrement / aire_km2(coeur_sens)
)
message("Sensibilité (h x ", FACTEUR_SENSIBILITE, ") : zone centrale = ",
        round(resultats$sensibilite_kde$surface_coeur_km2), " km², dont ",
        round(resultats$sensibilite_kde$part_coeur_dans_reference_pct), " % dans la zone centrale de référence.")

# 4d. Zones de concentration par secteur (Molène, Sein)
# Pourquoi : estimé sur tout le Parc, le KDE mélange les deux groupes côtiers et
# des dauphins vus au large (dauphins de passage ou du large, que les données
# publiques ne permettent pas de distinguer). Une seule largeur de bande, calculée
# sur des points très dispersés, est alors trop grande pour chaque groupe.
# Méthode : pour chaque secteur défini dans le bloc « Site étudié », on ne garde
# que les observations du secteur, on recalcule href sur ces seules observations,
# puis le KDE (grille plus fine, ~200 m) limité à la mer du secteur. On teste
# aussi href / 2. On note la position du centre de la zone centrale pour la comparer aux
# zones décrites par les suivis.
# Sorties : resultats$secteurs (surfaces), resultats$secteurs_obs (effectifs,
#           années, jeux de données), secteurs_dv (polygones).
kde_secteur <- function(xy, masque, hh, n_grille = 200) {
  e <- st_bbox(masque)
  k <- MASS::kde2d(xy[, 1], xy[, 2], h = 4 * hh, n = n_grille,
                   lims = c(e[["xmin"]], e[["xmax"]], e[["ymin"]], e[["ymax"]]))
  ddx <- diff(k$x[1:2]); ddy <- diff(k$y[1:2])
  r <- rast(nrows = length(k$y), ncols = length(k$x),
            xmin = min(k$x) - ddx / 2, xmax = max(k$x) + ddx / 2,
            ymin = min(k$y) - ddy / 2, ymax = max(k$y) + ddy / 2, crs = CRS_METRIQUE)
  values(r) <- as.vector(k$z[, ncol(k$z):1])
  r <- mask(r, vect(masque))
  v <- values(r, mat = FALSE); ok <- !is.na(v)
  o <- order(v[ok], decreasing = TRUE)
  cu <- numeric(sum(ok)); cu[o] <- cumsum(v[ok][o]) / sum(v[ok])
  vol <- rep(NA_real_, length(v)); vol[ok] <- cu
  setValues(r, vol)
}
boite_secteur <- function(sec) {
  st_transform(st_as_sfc(st_bbox(c(xmin = sec$xmin, ymin = sec$ymin,
                                   xmax = sec$xmax, ymax = sec$ymax), crs = st_crs(4326))),
               CRS_METRIQUE)
}

secteurs_dv  <- list()   # polygones
secteurs_vol <- list()   # rasters de volume (href du secteur), pour le trafic
lignes_surf  <- list()
lignes_obs   <- list()
for (s in seq_len(nrow(SITE$secteurs))) {
  sec    <- SITE$secteurs[s, ]
  masque <- st_sf(geometry = st_intersection(st_geometry(zone_etude), boite_secteur(sec)))
  pts    <- obs_parc[lengths(st_intersects(obs_parc, masque)) > 0, ]
  n_s    <- nrow(pts)
  if (n_s < N_MIN_SECTEUR) {
    message("Secteur ", sec$nom, " : ", n_s, " observations, zones non estimées.")
    next
  }
  xy_s   <- st_coordinates(pts)
  href_s <- sqrt(0.5 * (var(xy_s[, 1]) + var(xy_s[, 2]))) * n_s^(-1 / 6)

  jeux <- count(st_drop_geometry(pts), jeu_donnees, sort = TRUE)
  lignes_obs[[sec$nom]] <- data.frame(
    secteur = sec$nom, n_obs = n_s,
    annee_min = min(pts$annee), annee_max = max(pts$annee),
    n_depuis_2014 = sum(pts$annee >= 2014),
    jeu_principal = jeux$jeu_donnees[1],
    part_jeu_principal_pct = 100 * jeux$n[1] / n_s,
    part_inaturalist_pct = 100 * mean(grepl("iNaturalist", pts$jeu_donnees)),
    href_km = href_s / 1000)

  for (f in c(1, FACTEUR_SENSIBILITE)) {
    vol <- kde_secteur(xy_s, masque, f * href_s)
    if (f == 1) secteurs_vol[[sec$nom]] <- vol
    for (p in c(0.50, 0.95)) {
      g <- contour_iso(vol, p)
      if (length(g) == 0 || all(st_is_empty(g))) next
      centre_coeur <- st_coordinates(st_transform(st_centroid(g), 4326))
      secteurs_dv[[length(secteurs_dv) + 1]] <- st_sf(
        data.frame(secteur = sec$nom, facteur_h = f, niveau = paste0("KDE ", p * 100, " %")),
        geometry = g)
      lignes_surf[[length(lignes_surf) + 1]] <- data.frame(
        secteur = sec$nom, facteur_h = f, h_km = f * href_s / 1000,
        niveau = paste0("KDE ", p * 100, " %"), surface_km2 = aire_km2(g),
        centre_lon = centre_coeur[1, 1], centre_lat = centre_coeur[1, 2])
    }
  }
}
secteurs_dv <- if (length(secteurs_dv) > 0) do.call(rbind, secteurs_dv) else NULL
resultats$secteurs     <- bind_rows(lignes_surf)
resultats$secteurs_obs <- bind_rows(lignes_obs)
message("Zones de concentration par secteur :")
print(resultats$secteurs)

# 4e. Sensibilité au découpage des secteurs
# Pourquoi : les limites des secteurs sont un choix de l'analyste. On vérifie
# que la zone centrale ne dépend pas de ce choix en refaisant le calcul avec des
# boîtes réduites puis agrandies de 5 km de chaque côté.
# Sortie : resultats$sensibilite_secteurs (surface de la zone centrale et part
#          de cette surface comprise dans la zone centrale de référence).
MARGES_SECTEUR_M <- c(-5000, 5000)
lignes_dec <- list()
geom_dec   <- list()   # contours et rectangles, pour la figure de sensibilité
for (s in seq_len(nrow(SITE$secteurs))) {
  sec <- SITE$secteurs[s, ]
  ref <- secteurs_dv[secteurs_dv$secteur == sec$nom & secteurs_dv$facteur_h == 1 &
                       secteurs_dv$niveau == "KDE 50 %", ]
  if (is.null(ref) || nrow(ref) == 0) next
  for (m in MARGES_SECTEUR_M) {
    boite_m <- st_buffer(boite_secteur(sec), m, joinStyle = "MITRE", mitreLimit = 5)
    masque  <- st_sf(geometry = st_intersection(st_geometry(zone_etude), boite_m))
    pts     <- obs_parc[lengths(st_intersects(obs_parc, masque)) > 0, ]
    if (nrow(pts) < N_MIN_SECTEUR) next
    xy <- st_coordinates(pts)
    hm <- sqrt(0.5 * (var(xy[, 1]) + var(xy[, 2]))) * nrow(pts)^(-1 / 6)
    g  <- contour_iso(kde_secteur(xy, masque, hm), 0.50)
    rec <- aire_km2(suppressWarnings(st_intersection(g, st_geometry(ref))))
    variante <- if (m < 0) "Rectangle réduit de 5 km" else "Rectangle agrandi de 5 km"
    geom_dec[[length(geom_dec) + 1]] <- st_sf(
      data.frame(secteur = sec$nom, variante = variante, objet = c("zone", "rectangle")),
      geometry = c(st_union(g), st_geometry(boite_m)))
    lignes_dec[[length(lignes_dec) + 1]] <- data.frame(
      secteur = sec$nom, marge_km = m / 1000, n_obs = nrow(pts), h_km = hm / 1000,
      surface_centrale_km2 = aire_km2(g),
      # part de la nouvelle zone centrale comprise dans celle de référence
      part_dans_reference_pct = 100 * rec / aire_km2(g),
      # part de la zone centrale de référence couverte par la nouvelle
      part_reference_couverte_pct = 100 * rec / aire_km2(st_geometry(ref)))
  }
}
resultats$sensibilite_secteurs <- bind_rows(lignes_dec)
# On ajoute l'analyse principale pour la figure
for (s in seq_len(nrow(SITE$secteurs))) {
  sec <- SITE$secteurs[s, ]
  ref <- secteurs_dv[secteurs_dv$secteur == sec$nom & secteurs_dv$facteur_h == 1 &
                       secteurs_dv$niveau == "KDE 50 %", ]
  if (is.null(ref) || nrow(ref) == 0) next
  geom_dec[[length(geom_dec) + 1]] <- st_sf(
    data.frame(secteur = sec$nom, variante = "Analyse principale", objet = c("zone", "rectangle")),
    geometry = c(st_union(st_geometry(ref)), boite_secteur(sec)))
}
geom_dec <- if (length(geom_dec) > 0) do.call(rbind, geom_dec) else NULL
message("Sensibilité au découpage des secteurs :")
print(resultats$sensibilite_secteurs)


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

# Modèle utilisation / disponibilité (association, pas mesure de sélection) — régression logistique
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
# Objectif : décrire le chevauchement entre les zones de concentration et le trafic.
# Méthode  : densité de navires EMODnet (AIS, h/km²/mois) comparée entre la zone centrale,
#            l'anneau 50-95 % et le reste du Parc (médianes) ; part des cellules
#            de fort trafic (3e quartile, puis 50 % et 10 % en sensibilité)
#            situées dans les zones de concentration. Les catégories exclues pour ce site
#            (pêche professionnelle) sont définies dans le bloc « Site étudié ».
# Sorties  : resultats$trafic (tableau par zone et par catégorie de navires).

fichiers_trafic <- list.files("data/raw/trafic", pattern = "\\.tif$",
                              full.names = TRUE, ignore.case = TRUE)
# Catégories exclues de l'analyse pour ce site (voir le bloc « Site étudié »)
if (nzchar(SITE$trafic_exclu)) {
  exclus <- grepl(SITE$trafic_exclu, basename(fichiers_trafic), ignore.case = TRUE)
  if (any(exclus)) message("Couches de trafic exclues : ",
                           paste(basename(fichiers_trafic)[exclus], collapse = ", "))
  fichiers_trafic <- fichiers_trafic[!exclus]
}
# Ordre logique : tous navires (sert pour la carte et la note), plaisance,
# voile, passagers, autres
nb <- basename(fichiers_trafic)
rang <- ifelse(grepl("tous", nb, ignore.case = TRUE), 1,
        ifelse(grepl("plaisance", nb, ignore.case = TRUE), 2,
        ifelse(grepl("voile", nb, ignore.case = TRUE), 3,
        ifelse(grepl("passager", nb, ignore.case = TRUE), 4, 5))))
fichiers_trafic <- fichiers_trafic[order(rang, nb)]
resultats$trafic_disponible <- length(fichiers_trafic) > 0

# Zones : zone centrale (KDE 50 %), zone étendue (50-95 %), reste du Parc.
# Elles sont définies à partir des zones estimées par secteur (section 4d), qui
# est le résultat principal : toutes les analyses de trafic portent ainsi sur
# les mêmes contours que ceux des cartes. Le KDE sur l'ensemble du Parc n'est
# utilisé qu'en repli, si aucun secteur n'a pu être estimé.
if (!is.null(secteurs_dv)) {
  dv_ref  <- secteurs_dv[secteurs_dv$facteur_h == 1, ]
  gabarit <- rast(kde_volume)
  r50 <- rasterize(vect(dv_ref[dv_ref$niveau == "KDE 50 %", ]), gabarit, field = 1)
  r95 <- rasterize(vect(dv_ref[dv_ref$niveau == "KDE 95 %", ]), gabarit, field = 1)
  zones <- ifel(!is.na(r50), 1, ifel(!is.na(r95), 2, 3))
  zones <- mask(zones, vect(zone_etude))
  resultats$zones_trafic <- "zones de concentration par secteur (href de chaque secteur)"
} else {
  zones <- classify(kde_volume, rcl = matrix(c(0,    0.50, 1,
                                               0.50, 0.95, 2,
                                               0.95, 1.01, 3), ncol = 3, byrow = TRUE),
                    include.lowest = TRUE)
  resultats$zones_trafic <- "zones de concentration sur l'ensemble du Parc"
}
etiquettes_zones <- c("Zone centrale (KDE 50 %)", "Zone étendue (KDE 50-95 %)", "Reste du Parc")

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

# 6c. Contrôle de l'effet « proximité de la côte »
# Pourquoi : le trafic est en général plus dense près des côtes et des ports. Les
# zones centrales étant côtières, une partie de l'écart avec le reste du Parc
# peut venir de cette seule proximité. On compare donc aussi les zones centrales
# à la partie du reste du Parc située à une distance de la côte comparable (pas
# plus loin que 90 % des cellules des zones centrales).
# Sortie : resultats$trafic_cotier
if (resultats$trafic_disponible) {
  lignes <- list()
  for (nom in names(trafic_rasters)) {
    r  <- trafic_rasters[[nom]]
    z  <- resample(zones, r, method = "near")
    dc <- resample(dist_cote, r, method = "bilinear")
    df <- data.frame(zone = values(z, mat = FALSE), trafic = values(r, mat = FALSE),
                     dist = values(dc, mat = FALSE)) |>
      filter(!is.na(zone), !is.na(trafic), !is.na(dist))
    dmax <- unname(quantile(df$dist[df$zone == 1], 0.9))
    cotier <- df$zone == 3 & df$dist <= dmax
    lignes[[nom]] <- data.frame(
      couche = nom, dist_max_km = dmax / 1000,
      median_centrale = median(df$trafic[df$zone == 1]),
      median_reste = median(df$trafic[df$zone == 3]),
      median_reste_cotier = median(df$trafic[cotier]),
      n_cellules_reste_cotier = sum(cotier))
  }
  resultats$trafic_cotier <- bind_rows(lignes)
  message("Contrôle de la proximité de la côte :")
  print(resultats$trafic_cotier)
}

# 6b. Trafic par secteur : zone centrale du secteur (KDE 50 %, href du secteur) comparé
# au reste de la mer du secteur
if (resultats$trafic_disponible && length(secteurs_vol) > 0) {
  lignes <- list()
  for (nom_s in names(secteurs_vol)) {
    for (nom_t in names(trafic_rasters)) {
      r   <- trafic_rasters[[nom_t]]
      vol <- resample(secteurs_vol[[nom_s]], r, method = "near")
      lignes[[paste(nom_s, nom_t)]] <- data.frame(vol = values(vol, mat = FALSE),
                                                 trafic = values(r, mat = FALSE)) |>
        filter(!is.na(vol), !is.na(trafic)) |>
        mutate(zone = ifelse(vol <= 0.50, "Zone centrale du secteur", "Reste du secteur")) |>
        group_by(zone) |>
        summarise(n_cellules = n(), trafic_median = median(trafic),
                  trafic_moyen = mean(trafic), .groups = "drop") |>
        mutate(secteur = nom_s, couche = nom_t, .before = 1)
    }
  }
  resultats$trafic_secteurs <- bind_rows(lignes)
}


# ---- 7. Cartes et figures ----------------------------------------------------
# Objectif : produire des cartes et graphiques prêts à diffuser.
# Méthode  : ggplot2 ; habillage cartographique commun (échelle, nord
#            géographique, repères, sources et auteur en légende).
# Sorties  : outputs/figures/ (carte des zones de concentration, habitat, trafic).

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
lieux <- SITE$lieux
xy_lieux <- st_coordinates(st_transform(st_as_sf(lieux, coords = c("lon", "lat"), crs = 4326),
                                        CRS_METRIQUE))
lieux$x <- xy_lieux[, 1]
lieux$y <- xy_lieux[, 2]
lieux <- lieux[lieux$x > lim[["xmin"]] & lieux$x < lim[["xmax"]] &
               lieux$y > lim[["ymin"]] & lieux$y < lim[["ymax"]], ]

# Les noms sont posés sur une étiquette blanche semi-transparente : ils restent
# lisibles sur le fond bleu foncé et au milieu des points d'observation.
couches_lieux <- function() {
  etiquette <- function(type, ...) {
    geom_label(data = lieux[lieux$type == type, ],
               aes(x + dx, y + dy, label = nom, hjust = hjust),
               fill = scales::alpha("white", 0.75), linewidth = 0,
               label.padding = unit(0.08, "lines"), label.r = unit(0.05, "lines"),
               vjust = 0.5, lineheight = 0.85, ...)
  }
  list(
    geom_point(data = lieux[lieux$type == "ville", ], aes(x, y), shape = 22, size = 1.3,
               fill = "grey10", colour = "white", stroke = 0.3),
    etiquette("ville", size = 2.4, colour = "grey10"),
    etiquette("terre", size = 2.2, colour = "grey30", fontface = "italic"),
    etiquette("mer",   size = 2.1, colour = "#08306b", fontface = "italic"),
    geom_point(data = lieux[lieux$type == "phare", ], aes(x, y), shape = 24, size = 1.6,
               fill = "white", colour = "grey10", stroke = 0.4),
    etiquette("phare", size = 2.1, colour = "#08306b", fontface = "italic")
  )
}

# -- Barre d'échelle et flèche du nord (orientée sur le nord géographique :
#    en Lambert-93, le nord vrai est décalé de quelques degrés loin du méridien central)
centre_ll <- st_coordinates(st_transform(st_centroid(st_as_sfc(lim)), 4326))[1, 1:2]
p_nord <- st_coordinates(st_transform(st_sfc(st_point(centre_ll), st_point(centre_ll + c(0, 0.1)),
                                             crs = 4326), CRS_METRIQUE))
v_nord <- p_nord[2, 1:2] - p_nord[1, 1:2]
u_nord <- v_nord / sqrt(sum(v_nord^2))

habillage <- function(lim = get("lim", envir = globalenv()), L = 20000) {
  largeur <- lim[["xmax"]] - lim[["xmin"]]
  hauteur <- lim[["ymax"]] - lim[["ymin"]]
  x0 <- lim[["xmin"]] + 0.04 * largeur
  y0 <- lim[["ymin"]] + 0.06 * hauteur
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
    # Une instruction par graduation : nécessaire pour les cartes à plusieurs panneaux
    annotate("text", x = x0, y = y0 - 0.008 * hauteur, label = "0", vjust = 1, size = 2.3),
    annotate("text", x = x0 + L / 2, y = y0 - 0.008 * hauteur, label = format(L / 2000), vjust = 1, size = 2.3),
    annotate("text", x = x0 + L, y = y0 - 0.008 * hauteur, label = paste(L / 1000, "km"), vjust = 1, size = 2.3),
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

titre_espece <- function(debut, fin = "") {
  bquote(bold(.(debut) ~ "de grand dauphin (") * bolditalic("Tursiops truncatus") * bold(")" * .(fin)))
}

# Vocabulaire des figures. Le KDE est calculé sur des observations
# opportunistes : il décrit la densité des observations, qui dépend à la fois de
# la présence des dauphins et de l'endroit où se trouvent les observateurs. Les
# figures parlent donc de « zones de concentration des observations » ; le lien
# avec le domaine vital des groupes est discuté dans la note.
TITRE_KDE <- "Zones de concentration des observations"
etiquettes_dv <- c("KDE 50 %" = "Zone centrale (isoplèthe 50 %)",
                   "KDE 95 %" = "Zone étendue (isoplèthe 95 %)")

# Nom court du principal jeu de données (cité dans les légendes)
# (les noms GBIF répètent souvent l'intitulé : « X - Y » ; on garde X)
nom_court <- function(x, n = 120) {
  x <- sub(" - .*$", "", x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
}
jeu1 <- nom_court(resultats$par_jeu$jeu_donnees[1])

# Observations du Parc situées dans un secteur ou non (le KDE par secteur
# n'utilise que les premières ; l'analyse d'habitat les utilise toutes)
if (!is.null(secteurs_dv)) {
  boites <- do.call(c, lapply(seq_len(nrow(SITE$secteurs)), function(i) boite_secteur(SITE$secteurs[i, ])))
  dans_secteur <- lengths(st_intersects(obs_parc, st_union(boites))) > 0
} else {
  boites <- NULL
  dans_secteur <- rep(TRUE, n_parc)
}
resultats$n_hors_secteur <- sum(!dans_secteur)

# Contours affichés sur les cartes générales : zones estimées par secteur
# (résultat principal). L'estimation sur l'ensemble du Parc n'est pas tracée,
# pour ne pas montrer deux domaines vitaux sur la même carte : ses surfaces
# restent dans outputs/tableaux/ et dans le GeoPackage (domaine_vital_kde).
# Mettre AFFICHER_CONTEXTE <- TRUE pour la tracer en tirets.
AFFICHER_CONTEXTE <- FALSE
if (!is.null(secteurs_dv)) {
  dv_carte <- secteurs_dv[secteurs_dv$facteur_h == 1, c("niveau")]
  if (AFFICHER_CONTEXTE) {
    contexte <- iso95[, "niveau"]
    contexte$niveau <- "Parc 95 %"
    dv_carte <- rbind(dv_carte, contexte)
  }
  etiquettes_dv <- c("KDE 50 %"  = "Zone centrale (isoplèthe 50 %)",
                     "KDE 95 %"  = "Zone étendue (isoplèthe 95 %)",
                     "Parc 95 %" = "Isoplèthe 95 %, estimation sur tout le Parc")
} else {
  dv_carte <- domaine_vital
}
traits_dv <- c("KDE 50 %" = "solid", "KDE 95 %" = "solid", "Parc 95 %" = "22")
# Niveaux réellement présents sur la carte (breaks et labels doivent avoir la même longueur)
niveaux_carte <- intersect(names(etiquettes_dv), unique(dv_carte$niveau))

# -- Carte 1 : zones de concentration des observations
obs_hors <- obs_zone[lengths(st_intersects(obs_zone, zone_etude)) == 0, ]
lab_ret  <- paste0("Dans un secteur : KDE et habitat (n = ", sum(dans_secteur), ")")
lab_hs   <- paste0("Dans le Parc, hors secteur : habitat seulement (n = ", sum(!dans_secteur), ")")
lab_hors <- paste0("Hors du Parc : non analysées (n = ", nrow(obs_hors), ")")
obs_parc$categorie <- ifelse(dans_secteur, lab_ret, lab_hs)
lab_isobathe <- paste0("Isobathe ", SITE$isobathe, " m")
# Légendes construites à partir des seules catégories présentes (breaks, labels
# et override.aes doivent avoir la même longueur)
cat_obs <- data.frame(lab = c(lab_ret, lab_hs, lab_hors), forme = c(21, 24, 21),
                      fond = c("white", "white", "grey85"), trait = c("black", "black", "grey35"))
cat_obs <- cat_obs[c(any(dans_secteur), any(!dans_secteur), nrow(obs_hors) > 0), ]
cat_rep <- data.frame(lab = c("Limite du Parc", "Secteur analysé", lab_isobathe),
                      type = c("solid", "longdash", "dotted"),
                      couleur = c("grey15", "grey30", "grey10"), epaisseur = c(0.45, 0.3, 0.5))
if (is.null(boites)) cat_rep <- cat_rep[cat_rep$lab != "Secteur analysé", ]


carte_dv <- ggplot() +
  geom_raster(data = bathy_df, aes(x, y, fill = -profondeur)) +
  scale_fill_gradient(name = "Profondeur (m)", low = "#e3eef6", high = "#08306b",
                      trans = "sqrt", breaks = SITE$breaks_profondeur) +
  # Isobathe en pointillés foncés : un trait gris clair fin disparaît sur le bleu
  geom_contour(data = bathy_df, aes(x, y, z = -profondeur, linetype = lab_isobathe),
               breaks = SITE$isobathe, colour = "grey10", linewidth = 0.4) +
  geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey70", linewidth = 0.15) +
  geom_sf(data = parc, aes(linetype = "Limite du Parc"), fill = NA, colour = "grey15",
          linewidth = 0.45, key_glyph = "path") +
  # Contexte (tout le Parc) en tirets, zones par secteur en trait plein
  # (couche ajoutée seulement si AFFICHER_CONTEXTE = TRUE : une couche vide fait
  # planter ggplot)
  (if (AFFICHER_CONTEXTE) geom_sf(data = dv_carte[dv_carte$niveau == "Parc 95 %", ],
                                  aes(colour = niveau), fill = NA, linewidth = 0.7,
                                  linetype = "22", key_glyph = "path")) +
  geom_sf(data = dv_carte[dv_carte$niveau != "Parc 95 %", ], aes(colour = niveau),
          fill = NA, linewidth = 0.8, key_glyph = "path") +
  (if (!is.null(boites)) geom_sf(data = st_sf(geometry = boites), aes(linetype = "Secteur analysé"),
                                 fill = NA, colour = "grey30", linewidth = 0.3, key_glyph = "path")) +
  geom_sf(data = obs_hors, aes(shape = lab_hors), size = 0.8, colour = "grey35",
          fill = "grey85", stroke = 0.3) +
  geom_sf(data = obs_parc, aes(shape = categorie), size = 1.1, colour = "black",
          fill = "white", stroke = 0.3) +
  scale_colour_manual(name = "Estimation par noyau (KDE)",
                      values = c("KDE 50 %" = "#d7301f", "KDE 95 %" = "#fc8d59",
                                 "Parc 95 %" = "#fdd0a2"),
                      breaks = niveaux_carte, labels = etiquettes_dv[niveaux_carte]) +
  # Formes différentes (cercle / triangle) : le remplissage est déjà utilisé
  # par la bathymétrie
  scale_shape_manual(name = "Observations", values = setNames(cat_obs$forme, cat_obs$lab),
                     breaks = cat_obs$lab) +
  scale_linetype_manual(name = "Repères", values = setNames(cat_rep$type, cat_rep$lab),
                        breaks = cat_rep$lab) +
  guides(colour   = guide_legend(order = 1, override.aes = list(
                      linetype = unname(traits_dv[niveaux_carte]))),
         shape    = guide_legend(order = 2, override.aes = list(
                      fill = cat_obs$fond, colour = cat_obs$trait, size = 1.8)),
         linetype = guide_legend(order = 3, override.aes = list(
                      colour = cat_rep$couleur, linewidth = cat_rep$epaisseur)),
         fill     = guide_colourbar(order = 4, reverse = TRUE)) +
  couches_lieux() +
  habillage() +
  cadre +
  labs(title = titre_espece(TITRE_KDE),
       subtitle = if (!is.null(secteurs_dv)) {
         paste0(SITE$nom_long, " — ", n_parc, " observations, ", resultats$annees[1], "–",
                resultats$annees[2], "\nEstimation par noyau dans chaque secteur, largeur de bande href (",
                paste0(virgule(resultats$secteurs_obs$href_km), " km ",
                       resultats$secteurs_obs$secteur, collapse = " ; "),
                ")", if (AFFICHER_CONTEXTE) paste0(" ; tirets : estimation sur tout le Parc (h = ",
                                                  virgule(h / 1000), " km)"))
       } else {
         paste0(SITE$nom_long, "\nEstimation par noyau ",
                "(largeur de bande href = ", virgule(h / 1000), " km) — ",
                n_parc, " observations, ", resultats$annees[1], "–", resultats$annees[2])
       },
       caption = paste0("Observations hétérogènes (portails OBIS et GBIF) : leur répartition reflète aussi l'effort d'observation. ",
                        pct_jeu1, " % proviennent d'un même jeu :\n« ", jeu1, " ».\nLes contours décrivent la densité des ",
                        "observations, pas directement le domaine vital des groupes.\n",
                        "Bathymétrie NOAA ETOPO 2022 (~", round(resultats$resolution_bathy_m, -1), " m) ; ",
                        mention_parc, ".\n", signature)) +
  theme_carte
ggsave("outputs/figures/carte_domaine_vital.png", carte_dv, width = 9.5, height = 6.6, dpi = 300,
       bg = "white")

# -- Figure : habitat utilisé vs disponible
fp_txt <- function(p) if (p < 0.001) "p < 0,001" else paste0("p = ", virgule(p, 3))
# « valeur observations vs valeur milieu disponible »
med_txt <- function(col, unite) {
  r <- resultats$habitat_resume
  paste0(virgule(r[[col]][r$type == "Observations"]), unite, " (observations) contre ",
         virgule(r[[col]][r$type != "Observations"]), unite, " (milieu disponible)")
}
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
  # Graduations rapprochées près de zéro : c'est là que se trouvent les médianes
  # des observations (échelle en racine carrée)
  scale_y_sqrt(breaks = c(0, 1, 2, 5, 10, 20, 30, 40, 50, 75, 100, 150)) +
  scale_fill_manual(values = setNames(c("#d7301f", "#6baed6"), c(lab_obs, lab_dis)),
                    guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Profondeur et distance à la côte : observations et milieu disponible",
       subtitle = paste0("Médianes — profondeur : ", med_txt("prof_mediane_m", " m"),
                         "\nMédianes — distance à la côte : ", med_txt("dist_mediane_km", " km"),
                         "\nTests de Wilcoxon : ", fp_txt(resultats$test_prof$p.value), " (profondeur) ; ",
                         fp_txt(resultats$test_dist$p.value), " (distance)"),
       caption = paste0("Milieu disponible : points tirés au hasard dans la partie marine du Parc (",
                        N_DISPO_PAR_OBS, " par observation). Axes en racine carrée.\nProfondeur et distance ",
                        "extraites d'une bathymétrie à ~", round(resultats$resolution_bathy_m, -1), " m de résolution.\n",
                        "Les p-values supposent des observations indépendantes : avec des observations proches ",
                        "dans l'espace et un effort hétérogène, elles sont optimistes.\n",
                        "Réalisation : ", AUTEUR, ", ", date_fr, ".")) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        plot.caption = element_text(colour = "grey40", size = 6.5, hjust = 0),
        plot.title.position = "plot", plot.caption.position = "plot",
        strip.text = element_text(face = "bold", size = 9.5))
ggsave("outputs/figures/habitat_profondeur_cote.png", fig_hab, width = 8, height = 5.3,
       dpi = 300, bg = "white")

# -- Trafic maritime
libelle_couche <- function(x) {
  vapply(x, function(n) {
    annee <- regmatches(n, regexpr("[0-9]{4}", n))
    type  <- sub("_?[0-9]{4}$", "", n)
    noms  <- c(tous_navires = "Tous navires", plaisance = "Plaisance", voile = "Voile",
               passagers = "Passagers", peche = "Pêche")
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
    # Mêmes contours que ceux utilisés pour les calculs de trafic (par secteur)
    geom_sf(data = dv_carte[dv_carte$niveau != "Parc 95 %", ], aes(colour = niveau),
            fill = NA, linewidth = 0.8, key_glyph = "path") +
    scale_colour_manual(name = "Estimation par noyau (KDE)", values = c("KDE 50 %" = "#00a6a6",
                                                                 "KDE 95 %" = "#7fd4d4"),
                        labels = etiquettes_dv[c("KDE 50 %", "KDE 95 %")]) +
    scale_linetype_manual(name = "Repères", values = c("Limite du Parc" = "solid")) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2),
           fill = guide_colourbar(order = 3)) +
    couches_lieux() +
    habillage() +
    cadre +
    labs(title = titre_espece(paste("Trafic maritime et", tolower(TITRE_KDE))),
         subtitle = paste0("Densité de navires équipés de l'AIS — ", libelle_couche(nom1),
                           ", moyenne annuelle, grille de 1 km ; zones calculées dans chaque secteur"),
         caption = paste0("Sources : EMODnet Human Activities, densité de navires (AIS) ; ",
                          "observations OBIS et GBIF ; ", mention_parc, ".\n",
                          "« Tous navires » inclut la pêche professionnelle, exclue des autres catégories analysées. ",
                          "L'AIS n'est obligatoire que pour une partie des navires\n",
                          "(navires à passagers, grands navires, pêche de 15 m et plus) : le trafic des petites ",
                          "embarcations de plaisance et de pêche est sous-estimé.\n", signature)) +
    theme_carte
  ggsave("outputs/figures/carte_trafic.png", carte_trafic, width = 9.5, height = 6.6, dpi = 300,
         bg = "white")

  # Carte par catégorie : plaisance, voile, navires à passagers (pêche exclue)
  # Objectif : voir quelles activités traversent les zones de concentration,
  # ce que la carte « tous navires » ne permet pas de distinguer.
  noms_cat <- names(trafic_rasters)[grepl("plaisance|voile|passager", names(trafic_rasters))]
  if (length(noms_cat) > 0) {
    cat_df <- do.call(rbind, lapply(noms_cat, function(n) {
      d <- as.data.frame(trafic_rasters[[n]], xy = TRUE)
      d$couche_lab <- libelle_couche(n)
      d
    }))
    cat_df$couche_lab <- factor(cat_df$couche_lab, levels = libelle_couche(noms_cat))

    carte_categories <- ggplot() +
      geom_raster(data = cat_df, aes(x, y, fill = trafic)) +
      scale_fill_viridis_c(name = "Densité de navires (h / km² / mois)", trans = "log1p",
                           option = "magma", direction = -1,
                           breaks = c(0, 1, 10, 100, 1000)) +
      geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey70", linewidth = 0.15) +
      geom_sf(data = parc, fill = NA, colour = "grey15", linewidth = 0.35) +
      geom_sf(data = dv_carte[dv_carte$niveau != "Parc 95 %", ], aes(colour = niveau),
              fill = NA, linewidth = 0.6, key_glyph = "path") +
      scale_colour_manual(name = NULL, values = c("KDE 50 %" = "#00a6a6", "KDE 95 %" = "#7fd4d4"),
                          labels = etiquettes_dv[c("KDE 50 %", "KDE 95 %")]) +
      facet_wrap(~ couche_lab, ncol = length(noms_cat)) +
      habillage() +
      cadre +
      labs(title = titre_espece(paste("Trafic maritime par catégorie et", tolower(TITRE_KDE))),
           subtitle = "Densité de navires équipés de l'AIS, moyenne annuelle 2023, grille de 1 km ; même échelle de couleurs pour les trois panneaux",
           caption = paste0("Sources : EMODnet Human Activities, densité de navires (AIS) ; observations OBIS et GBIF ; ",
                            mention_parc, ".\n",
                            "L'AIS n'est obligatoire que pour une partie des navires : la plaisance et la voile sont sous-estimées.\n",
                            signature)) +
      theme_carte +
      theme(legend.position = "bottom", legend.box = "vertical",
            legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.3, "cm"),
            axis.text = element_blank(), panel.grid = element_blank(),
            strip.text = element_text(face = "bold", size = 9.5))
    ggsave("outputs/figures/carte_trafic_categories.png", carte_categories, width = 11, height = 5.6,
           dpi = 300, bg = "white")
  }

  # Schéma de synthèse pour la première page de la note
  # Objectif : faire passer les deux messages principaux d'un coup d'œil.
  # Les couloirs de navigation sont les cellules les plus fréquentées du Parc
  # (5 % des cellules aux plus fortes densités EMODnet 2023), calculées pour la
  # plaisance et la voile réunies d'une part, les navires à passagers d'autre part ;
  # les cellules isolées sont retirées.
  top_df <- function(motif, p = 0.95) {
    n <- names(trafic_rasters)[grepl(motif, names(trafic_rasters))]
    if (length(n) == 0) return(NULL)
    r <- Reduce(`+`, trafic_rasters[n])
    d <- as.data.frame(r, xy = TRUE)
    names(d)[3] <- "v"
    d <- d[d$v >= quantile(d$v, p, na.rm = TRUE), ]
    # On retire les cellules isolées (sans voisine retenue) : elles ne forment pas un couloir
    taille <- res(r)[1]
    voisins <- as.matrix(dist(d[, c("x", "y")])) <= 1.5 * taille
    d[rowSums(voisins) > 1, ]
  }
  nautisme  <- top_df("plaisance|voile")
  passagers <- top_df("passager")
  couloirs <- rbind(
    if (!is.null(nautisme))  data.frame(nautisme[, c("x", "y")],  cat = "Plaisance et voile"),
    if (!is.null(passagers)) data.frame(passagers[, c("x", "y")], cat = "Navires à passagers")
  )
  couloirs$cat <- factor(couloirs$cat, levels = c("Plaisance et voile", "Navires à passagers"))
  couloirs <- couloirs[order(couloirs$cat), ]          # passagers dessinés au-dessus
  taille_cellule <- res(trafic_rasters[[1]])
  lieux_s <- lieux[lieux$nom %in% c("Ouessant", "Archipel\nde Molène", "Le Conquet", "Île de Sein"), ]
  zones_s <- dv_carte[dv_carte$niveau %in% c("KDE 50 %", "KDE 95 %"), ]

  schema_synthese <- ggplot() +
    geom_sf(data = parc, fill = "#f3f7fa", colour = NA) +
    geom_sf(data = zones_s[zones_s$niveau == "KDE 95 %", ], aes(fill = "Zone étendue des observations"), colour = NA) +
    geom_sf(data = zones_s[zones_s$niveau == "KDE 50 %", ], aes(fill = "Zone centrale des observations"), colour = NA) +
    geom_tile(data = couloirs, aes(x, y, fill = cat),
              width = taille_cellule[1], height = taille_cellule[2]) +
    scale_fill_manual(name = NULL,
                      values = c("Zone centrale des observations" = "#ef8a85",
                                 "Zone étendue des observations"  = "#fbd3d0",
                                 "Plaisance et voile"             = "#e08214",
                                 "Navires à passagers"            = "#54278f"),
                      breaks = c("Zone centrale des observations", "Zone étendue des observations",
                                 "Plaisance et voile", "Navires à passagers")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    geom_sf(data = zones_s[zones_s$niveau == "KDE 50 %", ], fill = NA, colour = "#b2182b", linewidth = 0.5) +
    geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey60", linewidth = 0.15) +
    geom_sf(data = parc, fill = NA, colour = "grey40", linewidth = 0.3) +
    geom_label(data = lieux_s, aes(x + dx, y + dy, label = nom, hjust = hjust),
               size = 2.7, fill = scales::alpha("white", 0.8), linewidth = 0,
               label.padding = unit(0.1, "lines"), lineheight = 0.85) +
    habillage() +
    cadre +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 9.5, colour = "#2f6690"),
          plot.subtitle = element_text(size = 8, colour = "grey25"),
          plot.title.position = "plot",
          legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 8.5), legend.key.size = unit(0.4, "cm"),
          plot.margin = margin(4, 4, 4, 4)) +
    labs(title = bquote(bold("Grand dauphin (") * bolditalic("Tursiops truncatus") * bold(") et trafic maritime")),
         subtitle = "Zones de concentration des observations et secteurs de trafic le plus dense")
  ggsave("outputs/figures/schema_synthese.png", schema_synthese, width = 4.6, height = 5.2, dpi = 300, bg = "white")

  trafic_df$couche_lab <- factor(libelle_couche(trafic_df$couche),
                                 levels = libelle_couche(names(trafic_rasters)))
  # Nombre de cellules de 1 km² par zone (identique pour toutes les catégories)
  n_zone <- table(trafic_df$zone[trafic_df$couche == nom1])
  lab_zones_n <- paste0(c("Zone centrale\n(KDE 50 %)", "Zone étendue\n(KDE 50–95 %)", "Reste\ndu Parc"),
                        "\nn = ", as.integer(n_zone[etiquettes_zones]), " cellules")
  fig_trafic <- ggplot(trafic_df, aes(x = zone, y = trafic, fill = zone)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
    facet_wrap(~ couche_lab, scales = "free_y") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000)) +
    scale_x_discrete(labels = lab_zones_n) +
    scale_fill_manual(values = c("#d7301f", "#fc8d59", "#bdbdbd"), guide = "none") +
    labs(x = NULL, y = "Densité de navires (h / km² / mois, échelle log)",
         title = "Intensité du trafic maritime selon les zones de concentration des observations",
         subtitle = paste0("Une valeur par cellule de 1 km² ; moyenne annuelle EMODnet (navires équipés de l'AIS)\n",
                           "Échelles verticales différentes selon la catégorie"),
         caption = paste0("Description, sans test statistique : les cellules voisines ne sont pas indépendantes. ",
                          "« Tous navires » inclut la pêche professionnelle.\n",
                          "Le trafic des petites embarcations, souvent sans AIS, n'est pas représenté. ",
                          "Les observations opportunistes étant concentrées dans les eaux côtières,\n",
                          "une partie du recouvrement avec le trafic peut venir de l'effort d'observation. ",
                          "Réalisation : ", AUTEUR, ", ", date_fr, ".")) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          plot.caption = element_text(colour = "grey40", size = 6.5, hjust = 0),
          plot.title.position = "plot", plot.caption.position = "plot",
          strip.text = element_text(face = "bold", size = 9.5),
          axis.text.x = element_text(size = 7.5))
  ggsave("outputs/figures/trafic_par_zone.png", fig_trafic, width = 8.5, height = 5.2,
         dpi = 300, bg = "white")
}

# -- Cartes par secteur (zoom Molène et Sein)
if (!is.null(secteurs_dv)) {
  for (s in seq_len(nrow(SITE$secteurs))) {
    sec <- SITE$secteurs[s, ]
    if (!sec$nom %in% secteurs_dv$secteur) next
    boite <- boite_secteur(sec)
    lim_s <- st_bbox(st_buffer(boite, 1500))
    info  <- resultats$secteurs_obs[resultats$secteurs_obs$secteur == sec$nom, ]
    dv_s  <- secteurs_dv[secteurs_dv$secteur == sec$nom, ]
    dv_s$trait <- paste0(ifelse(dv_s$niveau == "KDE 50 %", "Zone centrale, 50 %", "Zone étendue, 95 %"),
                         ifelse(dv_s$facteur_h == 1, " (href)", " (href / 2)"))
    dv_s <- dv_s[!(dv_s$niveau == "KDE 95 %" & dv_s$facteur_h != 1), ]
    pts_s <- obs_parc[lengths(st_intersects(obs_parc, boite)) > 0, ]
    # Bathymétrie recoupée sur l'emprise du secteur (et non sur celle de la carte
    # générale, qui peut être plus petite et laisser une bande blanche)
    emp_s   <- ext(lim_s[["xmin"]], lim_s[["xmax"]], lim_s[["ymin"]], lim_s[["ymax"]])
    bathy_s <- as.data.frame(mask(crop(bathy_l93, emp_s), crop(mer, emp_s)), xy = TRUE)
    styles <- c("Zone centrale, 50 % (href)" = "#d7301f", "Zone étendue, 95 % (href)" = "#fc8d59",
                "Zone centrale, 50 % (href / 2)" = "#67000d")
    traits <- c("Zone centrale, 50 % (href)" = "solid", "Zone étendue, 95 % (href)" = "solid",
                "Zone centrale, 50 % (href / 2)" = "22")
    alerte_n <- if (info$n_obs < N_MIN_KDE) {
      paste0("\nEffectif faible (moins de ", N_MIN_KDE, " observations) : estimation indicative")
    } else ""
    carte_s <- ggplot() +
      geom_raster(data = bathy_s, aes(x, y, fill = -profondeur)) +
      scale_fill_gradient(name = "Profondeur (m)", low = "#e3eef6", high = "#08306b",
                          trans = "sqrt", breaks = SITE$breaks_profondeur) +
      geom_contour(data = bathy_s, aes(x, y, z = -profondeur), breaks = SITE$isobathe,
                   colour = "grey10", linewidth = 0.35, linetype = "dotted") +
      geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey70", linewidth = 0.15) +
      geom_sf(data = parc, fill = NA, colour = "grey15", linewidth = 0.45) +
      geom_sf(data = boite, fill = NA, colour = "grey40", linewidth = 0.3, linetype = "longdash") +
      geom_sf(data = dv_s, aes(colour = trait, linetype = trait), fill = NA, linewidth = 0.8,
              key_glyph = "path") +
      geom_sf(data = pts_s, size = 1.2, shape = 21, colour = "black", fill = "white", stroke = 0.3) +
      scale_colour_manual(name = "Estimation par noyau (KDE)", values = styles, breaks = names(styles)) +
      scale_linetype_manual(name = "Estimation par noyau (KDE)", values = traits, breaks = names(traits)) +
      guides(fill = guide_colourbar(reverse = TRUE)) +
      couches_lieux() +
      habillage(lim = lim_s, L = 10000) +
      coord_sf(crs = st_crs(CRS_METRIQUE), xlim = unname(lim_s[c(1, 3)]),
               ylim = unname(lim_s[c(2, 4)]), expand = FALSE) +
      labs(title = titre_espece(TITRE_KDE),
           subtitle = paste0("Secteur ", sec$libelle, "\n",
                             "Estimation sur les seules observations du secteur (n = ",
                             info$n_obs, ", ", info$annee_min, "–", info$annee_max,
                             ") ; href du secteur = ", virgule(info$href_km), " km", alerte_n),
           caption = paste0("Points blancs : observations du secteur. Tirets longs : secteur analysé, tracé autour ",
                            "du groupe côtier décrit par les suivis\npar photo-identification (sensibilité testée à ± 5 km). Pointillés : isobathe ",
                            SITE$isobathe, " m. Les contours décrivent la densité des observations, ",
                            "qui dépend aussi de l'effort d'observation.\n",
                            "Jeu de données principal : « ", nom_court(info$jeu_principal, 70), " » (",
                            round(info$part_jeu_principal_pct), " % des observations)",
                            if (!grepl("iNaturalist", info$jeu_principal))
                              paste0(" ; iNaturalist : ", round(info$part_inaturalist_pct), " %"),
                            ". Sources : OBIS, GBIF ;\n",
                            "bathymétrie NOAA ETOPO 2022 ; ", mention_parc, ".\n", signature)) +
      theme_carte
    ggsave(paste0("outputs/figures/carte_secteur_", sec$fichier, ".png"), carte_s,
           width = 8.5, height = 6.4, dpi = 300, bg = "white")
  }
}


# -- Figure : sensibilité au découpage des secteurs
# Une carte par secteur : les trois rectangles (réduit, principal, agrandi) en
# tirets, et la zone centrale (isoplèthe 50 %) obtenue avec chacun, en trait plein.
if (!is.null(geom_dec)) {
  couleurs_var <- c("Rectangle réduit de 5 km" = "#2c7fb8",
                    "Analyse principale"       = "#d7301f",
                    "Rectangle agrandi de 5 km" = "#7a0177")
  for (nom_s in unique(geom_dec$secteur)) {
    sec   <- SITE$secteurs[SITE$secteurs$nom == nom_s, ]
    gd    <- geom_dec[geom_dec$secteur == nom_s, ]
    gd$variante <- factor(gd$variante, levels = names(couleurs_var))
    lim_d <- st_bbox(st_buffer(gd[gd$objet == "rectangle", ], 1500))
    tab   <- resultats$sensibilite_secteurs[resultats$sensibilite_secteurs$secteur == nom_s, ]
    ref_s <- resultats$secteurs[resultats$secteurs$secteur == nom_s &
                                  resultats$secteurs$facteur_h == 1 &
                                  resultats$secteurs$niveau == "KDE 50 %", ]
    resume <- paste0("Zone centrale : ", round(tab$surface_centrale_km2[tab$marge_km < 0]), " km² (réduit) ; ",
                     round(ref_s$surface_km2), " km² (principal) ; ",
                     round(tab$surface_centrale_km2[tab$marge_km > 0]), " km² (agrandi)")
    fig_d <- ggplot() +
      geom_sf(data = terre_sf, fill = "#e8e4dc", colour = "grey70", linewidth = 0.15) +
      geom_sf(data = parc, fill = NA, colour = "grey15", linewidth = 0.4) +
      geom_sf(data = gd[gd$objet == "rectangle", ], aes(colour = variante), fill = NA,
              linewidth = 0.4, linetype = "longdash", show.legend = FALSE) +
      geom_sf(data = gd[gd$objet == "zone", ], aes(colour = variante), fill = NA,
              linewidth = 0.9, key_glyph = "path") +
      geom_sf(data = obs_parc, size = 0.9, shape = 21, colour = "black", fill = "white", stroke = 0.3) +
      scale_colour_manual(name = "Zone centrale (isoplèthe 50 %)
obtenue avec :", values = couleurs_var,
                          breaks = names(couleurs_var), drop = FALSE) +
      couches_lieux() +
      habillage(lim = lim_d, L = 10000) +
      coord_sf(crs = st_crs(CRS_METRIQUE), xlim = unname(lim_d[c(1, 3)]),
               ylim = unname(lim_d[c(2, 4)]), expand = FALSE) +
      labs(title = paste0("Sensibilité au découpage du secteur ", sec$libelle),
           subtitle = paste0("Même calcul avec le rectangle réduit, principal ou agrandi de 5 km de chaque côté\n", resume),
           caption = paste0("Tirets longs : les trois rectangles, chacun de la couleur de la zone qu'il produit. ",
                            "Points blancs : observations du Parc.\n",
                            "Lecture : si les zones sont emboîtées, leur position ne dépend pas du découpage ; ",
                            "seule leur étendue varie, avec la largeur de bande.\n", signature)) +
      theme_carte
    ggsave(paste0("outputs/figures/sensibilite_secteur_", sec$fichier, ".png"), fig_d,
           width = 8.5, height = 6.4, dpi = 300, bg = "white")
  }
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
# Détail des observations écartées par les deux contrôles qualité spécifiques
readr::write_csv(resultats$exclus_grille,  "data/propre/exclus_positions_grille.csv")
readr::write_csv(resultats$exclus_floutes, "data/propre/exclus_doublons_floutes.csv")

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
gpkg <- "outputs/sig/delfinu_iroise.gpkg"
# Supprimer l'ancien GeoPackage et ses fichiers annexes SQLite (-wal, -shm, -journal)
unlink(list.files("outputs/sig", pattern = "^delfinu_iroise\\.gpkg", full.names = TRUE))
st_write(zone_etude, gpkg, layer = "zone_etude", quiet = TRUE)
st_write(obs_parc[, c("source", "date", "annee", "profondeur", "jeu_donnees")],
         gpkg, layer = "observations_parc", quiet = TRUE)
st_write(obs_zone[, c("source", "date", "annee", "profondeur", "jeu_donnees")],
         gpkg, layer = "observations_emprise", quiet = TRUE)
st_write(domaine_vital, gpkg, layer = "domaine_vital_kde", quiet = TRUE)
st_write(mcp95, gpkg, layer = "domaine_vital_mcp95", quiet = TRUE)
if (!is.null(secteurs_dv)) st_write(secteurs_dv, gpkg, layer = "domaine_vital_secteurs", quiet = TRUE)
writeRaster(kde_volume, "outputs/sig/kde_volume.tif", overwrite = TRUE)
writeRaster(bathy_l93, "outputs/sig/bathymetrie_l93.tif", overwrite = TRUE)
writeRaster(dist_cote, "outputs/sig/distance_cote_l93.tif", overwrite = TRUE)

# 8c. Tableaux
readr::write_csv(resultats$domaine_vital, "outputs/tableaux/surfaces_domaine_vital.csv")
readr::write_csv(resultats$sensibilite_kde, "outputs/tableaux/sensibilite_domaine_vital.csv")
readr::write_csv(resultats$habitat_resume, "outputs/tableaux/habitat_resume.csv")
readr::write_csv(resultats$glm, "outputs/tableaux/modele_habitat.csv")
readr::write_csv(resultats$secteurs, "outputs/tableaux/domaine_vital_secteurs.csv")
readr::write_csv(resultats$secteurs_obs, "outputs/tableaux/observations_secteurs.csv")
readr::write_csv(resultats$sensibilite_secteurs, "outputs/tableaux/sensibilite_secteurs.csv")
if (!is.null(resultats$trafic_secteurs)) {
  readr::write_csv(resultats$trafic_secteurs, "outputs/tableaux/trafic_secteurs.csv")
}
if (!is.null(resultats$trafic_cotier)) {
  readr::write_csv(resultats$trafic_cotier, "outputs/tableaux/trafic_controle_cote.csv")
}
if (resultats$trafic_disponible) {
  readr::write_csv(resultats$trafic, "outputs/tableaux/trafic_par_zone.csv")
}

# 8d. Résultats pour la note et la présentation
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

message("\nTerminé. Résultats dans outputs/ — lancer ensuite postgis.R (base de données) puis le rendu de note.qmd.")
print(resultats$journal)
message("Positions sur grille exclues, par jeu et par point :")
print(resultats$exclus_grille)
message("Doublons imprécis exclus, par jeu :")
print(resultats$exclus_floutes)
print(resultats$domaine_vital)
message("Par secteur :")
print(resultats$secteurs_obs)
print(resultats$secteurs)
if (!is.null(resultats$trafic_secteurs)) print(resultats$trafic_secteurs)
