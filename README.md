# Domaine vital du grand dauphin dans le PNM du Cap Corse et de l'Agriate

Analyse exploratoire, à partir de données ouvertes, du domaine vital du grand dauphin (*Tursiops truncatus*) dans le Parc naturel marin du Cap Corse et de l'Agriate et de son chevauchement avec le trafic maritime. Exercice préparatoire au projet Delfinu.

**Auteur :** Gwladys Helaine — octobre 2026

## Contenu

| Fichier / dossier | Rôle |
|---|---|
| `analyse.R` | Script complet : téléchargement, nettoyage, domaine vital, habitat, pressions, exports |
| `note.qmd` | Note de synthèse (HTML / PDF) générée à partir des résultats |
| `presentation.qmd` | Présentation (reveal.js) |
| `data/propre/` | Observations nettoyées, journal de nettoyage, dictionnaire de données |
| `outputs/figures/` | Cartes et graphiques |
| `outputs/sig/` | Couches SIG (GeoPackage, rasters Lambert-93) |
| `outputs/tableaux/` | Résultats chiffrés |

## Reproduire

1. Déposer le périmètre du Parc (INPN) dans `data/raw/parc/` et, en option, les rasters de densité de navires EMODnet dans `data/raw/trafic/`.
2. Lancer `source("analyse.R")` depuis le dossier du projet.
3. Rendre `note.qmd` et `presentation.qmd` avec Quarto.

## Sources de données

- Observations : OBIS, GBIF
- Bathymétrie : NOAA ETOPO (via le package marmap)
- Trafic maritime : EMODnet Human Activities, densité de navires
- Périmètre : INPN, espaces protégés

## Méthodes

Estimation du domaine vital par noyau (isoplèthes 50 % et 95 %, largeur de bande href) et polygone convexe minimum ; comparaison utilisation / disponibilité de l'habitat (tests de Wilcoxon, régression logistique) ; comparaison de la densité de navires entre zones du domaine vital. Outils : R (sf, terra, dplyr, ggplot2), QGIS, Quarto.

Les limites (biais d'effort d'observation, absence de photo-identification, résolution) sont détaillées dans la note.
