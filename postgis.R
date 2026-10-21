# =============================================================================
#  Bancarisation des résultats dans une base PostgreSQL / PostGIS
#
#  À lancer APRÈS analyse.R : source("postgis.R")
#
#  Prérequis (une seule fois, voir le GUIDE) :
#   - PostgreSQL et PostGIS installés, base « delfinu » créée avec l'extension
#     postgis (CREATE EXTENSION postgis;)
#   - packages R : DBI, RPostgres
#
#  Ce script :
#   1. se connecte à la base « delfinu » ;
#   2. crée un schéma « iroise » et y écrit les tables spatiales (observations,
#      périmètre du Parc, zones de concentration des observations, dans les
#      tables « domaine_vital* », nom hérité de la première version) avec leurs
#      index spatiaux ;
#   3. ajoute une contrainte de cohérence (année) ;
#   4. exécute les requêtes d'exemple de requetes.sql et affiche les résultats.
# =============================================================================

paquets <- c("DBI", "RPostgres", "sf", "readr")
manquants <- setdiff(paquets, rownames(installed.packages()))
if (length(manquants) > 0) install.packages(manquants)

library(DBI)
library(sf)

BASE   <- "delfinu"
SCHEMA <- "iroise"
GPKG   <- "outputs/sig/delfinu_iroise.gpkg"

if (!file.exists(GPKG)) stop("Lancer d'abord analyse.R (fichier ", GPKG, " absent).")

# ---- 1. Connexion -------------------------------------------------------------
# Sous Linux, la connexion locale utilise le compte de l'utilisateur courant
con <- dbConnect(RPostgres::Postgres(), dbname = BASE)
# (pas de on.exit() ici : dans un script lancé avec source(), il fermerait la
#  connexion immédiatement ; on la ferme à la fin du script)

if (nrow(dbGetQuery(con, "SELECT 1 FROM pg_extension WHERE extname = 'postgis'")) == 0) {
  stop("L'extension PostGIS n'est pas activée dans la base ", BASE,
       " : exécuter CREATE EXTENSION postgis; dans psql.")
}
message("Connecté à la base ", BASE, " — ",
        dbGetQuery(con, "SELECT postgis_full_version() AS v")$v |> substr(1, 40))

# ---- 2. Écriture des tables -----------------------------------------------------
dbExecute(con, paste0("CREATE SCHEMA IF NOT EXISTS ", SCHEMA))

ecrire <- function(x, table) {
  nom <- Id(schema = SCHEMA, table = table)
  st_write(x, con, nom, delete_layer = TRUE, quiet = TRUE)
  # Index spatial : accélère les requêtes géographiques (ST_Within, ST_DWithin...)
  geom <- attr(x, "sf_column")
  dbExecute(con, sprintf("CREATE INDEX IF NOT EXISTS %s_geom_idx ON %s.%s USING GIST (%s)",
                         table, SCHEMA, table, geom))
  message("Table ", SCHEMA, ".", table, " : ", nrow(x), " ligne(s)")
}

obs     <- st_read(GPKG, layer = "observations_emprise", quiet = TRUE)
parc    <- st_read(GPKG, layer = "zone_etude", quiet = TRUE)
dv      <- st_read(GPKG, layer = "domaine_vital_kde", quiet = TRUE)
dv_sec  <- st_read(GPKG, layer = "domaine_vital_secteurs", quiet = TRUE)
names(obs)    <- tolower(names(obs))
names(dv)     <- tolower(names(dv))
names(dv_sec) <- tolower(names(dv_sec))

ecrire(obs,  "observations")
ecrire(parc, "parc")
ecrire(dv,   "domaine_vital")            # estimation sur tout le Parc (contexte)
ecrire(dv_sec, "domaine_vital_secteurs") # résultat principal : Molène, Sein

# ---- 3. Contraintes de cohérence ------------------------------------------------
# Une base protège les données : une valeur impossible est refusée à la saisie
dbExecute(con, sprintf("ALTER TABLE %s.observations DROP CONSTRAINT IF EXISTS annee_valide", SCHEMA))
dbExecute(con, sprintf("ALTER TABLE %s.observations ADD CONSTRAINT annee_valide
                        CHECK (annee BETWEEN 1900 AND EXTRACT(YEAR FROM CURRENT_DATE))", SCHEMA))
dbExecute(con, sprintf("COMMENT ON TABLE %s.observations IS
  'Observations nettoyées de Tursiops truncatus (OBIS, GBIF) dans l''emprise d''étude, Lambert-93'", SCHEMA))

# ---- 4. Requêtes d'exemple ------------------------------------------------------
# Les requêtes sont dans requetes.sql, séparées par des lignes « -- @nom »
texte <- readLines("requetes.sql", encoding = "UTF-8")
debuts <- grep("^-- @", texte)
fins   <- c(debuts[-1] - 1, length(texte))
# Les résultats d'une requête ne sont pas stockés dans la base : ils sont
# calculés à la demande. On les enregistre donc ici, un CSV par requête, plus
# un journal texte qui reprend toutes les requêtes et leurs résultats.
dossier_req <- "outputs/tableaux/requetes"
dir.create(dossier_req, recursive = TRUE, showWarnings = FALSE)
journal_req <- file.path(dossier_req, "journal_requetes.txt")
writeLines(c(paste("Requêtes exécutées sur la base", BASE, "le", format(Sys.time(), "%d/%m/%Y %H:%M")),
             ""), journal_req, useBytes = TRUE)

for (i in seq_along(debuts)) {
  bloc <- texte[debuts[i]:fins[i]]
  nom  <- sub("^-- @", "", bloc[1])
  sql  <- paste(bloc[!grepl("^\\s*--", bloc)], collapse = "\n")
  if (!nzchar(trimws(sql))) next
  message("\n== ", nom)
  resultat <- dbGetQuery(con, sql)
  print(resultat)
  # Un fichier CSV par requête (nom simplifié : sans accents ni espaces)
  fichier <- sprintf("%02d_%s.csv", i, gsub("^_|_$", "", gsub("[^a-z0-9]+", "_",
                     tolower(iconv(nom, "UTF-8", "ASCII//TRANSLIT")))))
  readr::write_csv(resultat, file.path(dossier_req, fichier))
  # Journal : intitulé, requête SQL, résultat
  cat("== ", nom, "\n", sql, "\n\n", paste(capture.output(print(resultat)), collapse = "\n"),
      "\n\n", sep = "", file = journal_req, append = TRUE)
}
message("\nRésultats des requêtes enregistrés dans ", dossier_req)

dbDisconnect(con)
message("\nBase prête. Dans QGIS : Couche > Ajouter une couche > PostGIS, ",
        "nouvelle connexion vers la base « ", BASE, " ».")
