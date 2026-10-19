# Domaine vital du grand dauphin dans le PNM du Cap Corse et de l'Agriate

Analyse exploratoire, à partir de données ouvertes, du domaine vital du grand dauphin (*Tursiops truncatus*) dans le Parc naturel marin du Cap Corse et de l'Agriate et de son chevauchement avec le trafic maritime. Exercice préparatoire au projet Delfinu.

**Auteure :** Gwladys Helaine — octobre 2026

**À lire en premier :** [`note.pdf`](note.pdf), la note de synthèse (méthodes, résultats, limites, conclusion).

## Principaux résultats

- 1 329 enregistrements bruts (OBIS, GBIF) → **161 observations** retenues dans le Parc après nettoyage (échouages, doublons, localisations imprécises, points à terre).
- Cœur du domaine vital (noyau 50 %) : **136 km²**, dans le golfe de Saint-Florent ; domaine vital (noyau 95 %) : **1 165 km²**.
- Observations côtières : profondeur médiane de 60 m contre 1 962 m pour le milieu disponible.
- Robustesse : avec une largeur de bande réduite de moitié, 98 % du cœur recalculé reste dans le cœur de référence.
- Densité médiane de navires (AIS) **2,0 fois** plus élevée dans le cœur du domaine vital que dans le reste du Parc, et **3,3 fois** pour la plaisance ; la surreprésentation des zones de fort trafic dans le domaine vital se maintient quel que soit le seuil testé.
- Limites principales : 78 % des observations proviennent d'une seule base opportuniste (la répartition reflète aussi l'effort d'observation), et 83 % datent de 2001 à 2010 alors que le trafic date de 2023. Les résultats illustrent une méthode plus qu'ils ne décrivent la situation actuelle.

## Contenu du dépôt

| Fichier | Rôle |
|---|---|
| `analyse.R` | Script complet : téléchargement, nettoyage, domaine vital, habitat, pressions, figures, exports |
| `note.qmd` / `note.pdf` | Note de synthèse (source Quarto et version PDF) |
| `presentation.qmd` | Présentation (reveal.js) |
| `carte_domaine_vital.png`, `carte_trafic.png` | Cartes |
| `habitat_profondeur_cote.png`, `trafic_par_zone.png` | Graphiques |
| `observations_propres_parc.csv` | Observations nettoyées (161) |
| `journal_nettoyage.csv` | Nombre d'observations restant après chaque règle de nettoyage |
| `dictionnaire_donnees.csv` | Description de chaque colonne des données nettoyées |
| `surfaces_domaine_vital.csv`, `habitat_resume.csv`, `modele_habitat.csv`, `trafic_par_zone.csv` | Résultats chiffrés |
| `sensibilite_domaine_vital.csv` | Analyse de sensibilité du domaine vital (largeur de bande réduite de moitié) |
| `delfinu-exercice.Rproj` | Projet RStudio |

## Reproduire l'analyse

1. Ouvrir `delfinu-exercice.Rproj` dans RStudio.
2. Déposer le périmètre du Parc (`sig_metrop.gpkg`, base des espaces protégés de PatriNat) dans `data/raw/parc/` et, en option, les rasters de densité de navires EMODnet dans `data/raw/trafic/`.
3. Lancer `source("analyse.R")` : le script télécharge les observations et la bathymétrie, puis recrée tous les dossiers et résultats (`data/propre/`, `outputs/`).
4. Rendre `note.qmd` et `presentation.qmd` avec Quarto.

## Sources de données

- Observations : OBIS, GBIF (dont base OGREVA de la DREAL Corse, Observatoire Pelagis, GIS3M, Tethys Research Institute)
- Bathymétrie : NOAA ETOPO 2022 (via le package marmap)
- Trafic maritime : EMODnet Human Activities, densité de navires (AIS), moyenne annuelle 2023
- Périmètre : PatriNat, base des espaces protégés (07/2026)

## Méthodes

Estimation du domaine vital par noyau (isoplèthes 50 % et 95 %, largeur de bande href) et polygone convexe minimum, avec analyse de sensibilité à la largeur de bande ; comparaison utilisation / disponibilité de l'habitat (tests de Wilcoxon, régression logistique) ; comparaison de la densité de navires entre zones du domaine vital, avec test de sensibilité au seuil de fort trafic. Outils : R (sf, terra, dplyr, ggplot2), QGIS, Quarto.
