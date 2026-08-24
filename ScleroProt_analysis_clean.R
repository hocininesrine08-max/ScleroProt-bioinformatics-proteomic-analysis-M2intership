# ==============================================================================
# ScléroProt - Analyse protéomique des compartiments EV et FREE
# ==============================================================================
# Analyse comparative des groupes Healthy controls, Raynaud, VEDOSS et SSc.
# Les étapes comprennent le prétraitement, l'imputation, l'ACP, les analyses
# différentielles et la visualisation des signatures protéiques.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

packages <- c(
  "readxl", "ggplot2", "ggpubr", "gridExtra", "ggrepel",
  "dplyr", "tidyr", "grid", "imputeLCMD", "FactoMineR",
  "factoextra", "VennDiagram", "pheatmap"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Packages manquants : ",
    paste(missing_packages, collapse = ", "),
    ". Installer les packages avant l'exécution du script."
  )
}

invisible(lapply(packages, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# 2. Paramètres généraux et importation
# ------------------------------------------------------------------------------

data_dir <- "data"

# Adapter les noms de fichiers si les matrices EV et FREE sont séparées.
file_ev <- file.path(data_dir, "Matrix-VES-prot (1).xlsx")
file_free <- file.path(data_dir, "Matrix-VES-prot (1).xlsx")

mon_data_ev <- read_excel(file_ev)
mon_data_free <- read_excel(file_free)

# Index de début des variables protéiques dans les matrices d'origine.
index_depart_commun <- 8

# ------------------------------------------------------------------------------
# 3. Fonctions de préparation et de visualisation
# ------------------------------------------------------------------------------

palette_sexe <- c(
  "Femme" = "#FFC0CB",
  "Homme" = "#B0E0E6"
)

contour_sexe <- c(
  "Femme" = "#C71585",
  "Homme" = "#4682B4"
)

palette_pid <- c(
  "PID0" = "#FDFD96",
  "PID1" = "#C3B1E1",
  "0" = "grey90"
)

contour_pid <- c(
  "PID0" = "#9B870C",
  "PID1" = "#512888",
  "0" = "grey40"
)

palette_detail <- c(
  "Healthy controls" = "#A3D9C9",
  "VEDOSS" = "#FFD1A9",
  "Raynaud" = "#E6A8D7",
  "Idiopathic Raynaud" = "#E6A8D7",
  "Raynaud connective" = "#D8BFD8",
  "SSc Limited cutaneous" = "#F5A3A3",
  "SSc diffuse cutaneous" = "#E6B8B8"
)

contour_detail <- c(
  "Healthy controls" = "#2D6A4F",
  "VEDOSS" = "#B25E00",
  "Raynaud" = "#800080",
  "Idiopathic Raynaud" = "#800080",
  "Raynaud connective" = "#4B0082",
  "SSc Limited cutaneous" = "#9A1F1F",
  "SSc diffuse cutaneous" = "#5C1111"
)

palette_pie <- c(
  "Données Présentes" = "#A3D9C9",
  "Données Manquantes (NA)" = "#F5A3A3"
)

contour_pie <- c(
  "Données Présentes" = "#2D6A4F",
  "Données Manquantes (NA)" = "#9A1F1F"
)

generer_barplot <- function(data, colonne, titre,
                            palette_remplissage, palette_contour) {
  df_counts <- data %>%
    group_by(.data[[colonne]]) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Prop = Count / sum(Count),
      Label = paste0(Count, " (", round(Prop * 100, 1), "%)")
    )

  ggplot(
    df_counts,
    aes(
      x = .data[[colonne]],
      y = Count,
      fill = .data[[colonne]],
      color = .data[[colonne]]
    )
  ) +
    geom_col(width = 0.7, linewidth = 0.8) +
    geom_text(
      aes(label = Label),
      vjust = -0.4,
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(values = palette_remplissage, na.value = "grey90") +
    scale_color_manual(values = palette_contour, na.value = "grey40") +
    theme_minimal() +
    labs(title = titre, x = "", y = "Nombre de patients") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none"
    )
}

nettoyer_clinique <- function(data) {
  if ("sexe" %in% colnames(data)) {
    data$sexe_clean <- ifelse(
      tolower(data$sexe) %in% c("f", "femme", "female"),
      "Femme",
      "Homme"
    )
  } else if ("sex" %in% colnames(data)) {
    data$sexe_clean <- ifelse(
      tolower(data$sex) %in% c("f", "femme", "female"),
      "Femme",
      "Homme"
    )
  } else {
    data$sexe_clean <- NA_character_
  }

  data$Groupe_Detail <- case_when(
    tolower(data$identification) == "healthy controls" ~ "Healthy controls",
    tolower(data$identification) == "vedoss" ~ "VEDOSS",
    tolower(data$identification) == "raynaud" ~ "Raynaud",
    tolower(data$identification) == "idiopathic raynaud" ~ "Idiopathic Raynaud",
    tolower(data$identification) == "raynaud connective" ~ "Raynaud connective",
    tolower(data$identification) == "ssc diffuse cutaneous" ~ "SSc diffuse cutaneous",
    tolower(data$identification) == "ssc limited cutaneous" ~ "SSc Limited cutaneous",
    TRUE ~ as.character(data$identification)
  )

  if (!"severite" %in% colnames(data)) {
    data$severite <- "0"
  }

  data
}

convertir_en_numerique <- function(matrice) {
  df <- as.data.frame(matrice)
  for (col in colnames(df)) {
    df[[col]] <- as.numeric(as.character(df[[col]]))
  }
  as.matrix(df)
}

dessiner_camembert <- function(df_pie, titre) {
  ggplot(df_pie, aes(x = "", y = Nombre, fill = Statut, color = Statut)) +
    geom_col(width = 1, linewidth = 0.8) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = palette_pie) +
    scale_color_manual(values = contour_pie) +
    theme_void() +
    labs(title = titre) +
    geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      fontface = "bold",
      size = 4.5
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# ------------------------------------------------------------------------------
# 4. Harmonisation clinique et descriptif des cohortes
# ------------------------------------------------------------------------------

mon_data_ev <- nettoyer_clinique(mon_data_ev)
mon_data_free <- nettoyer_clinique(mon_data_free)

g1 <- generer_barplot(
  mon_data_ev, "sexe_clean", "EV : Répartition du sexe",
  palette_sexe, contour_sexe
)
g2 <- generer_barplot(
  mon_data_ev, "Groupe_Detail", "EV : Groupes cliniques",
  palette_detail, contour_detail
)
g3 <- generer_barplot(
  mon_data_ev, "severite", "EV : Atteinte pulmonaire",
  palette_pid, contour_pid
)

grid.arrange(g1, g2, g3, ncol = 3)

g4 <- generer_barplot(
  mon_data_free, "sexe_clean", "FREE : Répartition du sexe",
  palette_sexe, contour_sexe
)
g5 <- generer_barplot(
  mon_data_free, "Groupe_Detail", "FREE : Groupes cliniques",
  palette_detail, contour_detail
)
g6 <- generer_barplot(
  mon_data_free, "severite", "FREE : Atteinte pulmonaire",
  palette_pid, contour_pid
)

grid.arrange(g4, g5, g6, ncol = 3)

# ------------------------------------------------------------------------------
# 5. Filtrage des protéines selon les valeurs manquantes
# ------------------------------------------------------------------------------

matrice_ev_brute <- mon_data_ev[
  , index_depart_commun:(ncol(mon_data_ev) - 2)
]

matrice_free_brute <- mon_data_free[
  , index_depart_commun:(ncol(mon_data_free) - 2)
]

taux_na_ev <- colMeans(is.na(matrice_ev_brute))
taux_na_free <- colMeans(is.na(matrice_free_brute))

df_pie_ev <- data.frame(
  Statut = c("Données Présentes", "Données Manquantes (NA)"),
  Nombre = c(
    sum(!is.na(matrice_ev_brute)),
    sum(is.na(matrice_ev_brute))
  )
) %>%
  mutate(Label = paste0(round(Nombre / sum(Nombre) * 100, 1), "%"))

df_pie_free <- data.frame(
  Statut = c("Données Présentes", "Données Manquantes (NA)"),
  Nombre = c(
    sum(!is.na(matrice_free_brute)),
    sum(is.na(matrice_free_brute))
  )
) %>%
  mutate(Label = paste0(round(Nombre / sum(Nombre) * 100, 1), "%"))

grid.arrange(
  dessiner_camembert(df_pie_ev, "Global EV"),
  dessiner_camembert(df_pie_free, "Global FREE"),
  ncol = 2
)

# Seuil de conservation : au maximum 70 % de valeurs manquantes.
ev_solides <- names(taux_na_ev[taux_na_ev <= 0.70])
free_solides <- names(taux_na_free[taux_na_free <= 0.70])

matrice_ev_filtree <- matrice_ev_brute[, ev_solides, drop = FALSE]
matrice_free_filtree <- matrice_free_brute[, free_solides, drop = FALSE]

# ------------------------------------------------------------------------------
# 6. Réintégration des protéines d'intérêt
# ------------------------------------------------------------------------------

prot_a_sauver_ev <- c("IGHV1-58", "AMY1A;AMY1B;AMY1C", "FCER1G", "NAGLU")
prot_a_sauver_free <- c("PRPS1", "CALML5")

matrice_ev_repêchee <- mon_data_ev[
  , unique(c(colnames(matrice_ev_filtree), prot_a_sauver_ev)),
  drop = FALSE
]

matrice_free_repêchee <- mon_data_free[
  , unique(c(colnames(matrice_free_filtree), prot_a_sauver_free)),
  drop = FALSE
]

mat_ev_safe <- convertir_en_numerique(matrice_ev_repêchee)
mat_free_safe <- convertir_en_numerique(matrice_free_repêchee)

# ------------------------------------------------------------------------------
# 7. Imputation des valeurs manquantes par QRILC
# ------------------------------------------------------------------------------

matrice_ev_imputee <- impute.QRILC(mat_ev_safe)[[1]]
matrice_free_imputee <- impute.QRILC(mat_free_safe)[[1]]

clinique_ev <- mon_data_ev[, c(1:7, (ncol(mon_data_ev) - 1):ncol(mon_data_ev))]
clinique_free <- mon_data_free[, c(1:7, (ncol(mon_data_free) - 1):ncol(mon_data_free))]

clinique_ev$Statut_Simple <- ifelse(
  tolower(clinique_ev$identification) == "healthy controls",
  "Healthy controls",
  ifelse(tolower(clinique_ev$identification) == "vedoss", "VEDOSS", "SSc")
)

clinique_free$Statut_Simple <- ifelse(
  tolower(clinique_free$identification) == "healthy controls",
  "Healthy controls",
  ifelse(tolower(clinique_free$identification) == "vedoss", "VEDOSS", "SSc")
)

mon_data_ev_final <- cbind(clinique_ev, as.data.frame(matrice_ev_imputee))
mon_data_free_final <- cbind(clinique_free, as.data.frame(matrice_free_imputee))

# ------------------------------------------------------------------------------
# 8. ACP et identification des patients potentiellement aberrants
# ------------------------------------------------------------------------------

generer_acp_outliers <- function(data_finale, nom_matrice) {
  cols_a_exclure <- c(
    "sexe", "sex", "severite",
    "evolution-SSc", "connectivite"
  )

  mat_prot <- data_finale %>%
    select(where(is.numeric)) %>%
    select(-any_of(cols_a_exclure))

  res_pca <- PCA(mat_prot, scale.unit = TRUE, graph = FALSE)

  coords <- as.data.frame(res_pca$ind$coord[, 1:2])
  colnames(coords) <- c("PC1", "PC2")

  coords$identification <- data_finale$identification
  coords$Statut_Simple <- data_finale$Statut_Simple
  coords$Distance_Centre <- sqrt(coords$PC1^2 + coords$PC2^2)

  seuil_outlier <- mean(coords$Distance_Centre) +
    2.5 * sd(coords$Distance_Centre)

  coords$Est_Outlier <- coords$Distance_Centre > seuil_outlier
  coords$Label_Outlier <- ifelse(
    coords$Est_Outlier,
    coords$identification,
    ""
  )

  p_acp <- ggplot(
    coords,
    aes(x = PC1, y = PC2, color = Statut_Simple)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(aes(size = Distance_Centre), alpha = 0.7) +
    geom_point(
      data = filter(coords, Est_Outlier),
      shape = 1,
      size = 5,
      color = "black",
      stroke = 1.2
    ) +
    geom_text_repel(
      aes(label = Label_Outlier),
      color = "black",
      fontface = "bold",
      size = 3,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c(
        "SSc" = "#FF4B4B",
        "Healthy controls" = "#4B9CD3",
        "VEDOSS" = "#50B848"
      )
    ) +
    scale_size_continuous(range = c(1.5, 4.5)) +
    theme_minimal() +
    labs(
      title = paste0("ACP - Outliers (", nom_matrice, ")"),
      x = paste0("PC1 (", round(res_pca$eig[1, 2], 1), "%)"),
      y = paste0("PC2 (", round(res_pca$eig[2, 2], 1), "%)")
    )

  print(p_acp)
  coords
}

outliers_ev <- generer_acp_outliers(mon_data_ev_final, "Matrice EV")
outliers_free <- generer_acp_outliers(mon_data_free_final, "Matrice FREE")

# Liste des patients à exclure après inspection de l'ACP.
patients_a_exclure <- character(0)

mon_data_ev_clean <- mon_data_ev_final %>%
  filter(!identification %in% patients_a_exclure)

mon_data_free_clean <- mon_data_free_final %>%
  filter(!identification %in% patients_a_exclure)

# ------------------------------------------------------------------------------
# 9. Analyse différentielle SSc vs Healthy controls
# ------------------------------------------------------------------------------

analyser_ssc_vs_hc <- function(data_finale, nom_matrice) {
  df_sub <- data_finale %>%
    filter(Statut_Simple %in% c("SSc", "Healthy controls"))

  idx_prot <- (
    which(colnames(df_sub) == "Statut_Simple") + 1
  ):ncol(df_sub)

  noms_prot <- colnames(df_sub)[idx_prot]

  resultats <- lapply(noms_prot, function(prot) {
    valeurs_ssc <- as.numeric(
      df_sub[[prot]][df_sub$Statut_Simple == "SSc"]
    )
    valeurs_hc <- as.numeric(
      df_sub[[prot]][df_sub$Statut_Simple == "Healthy controls"]
    )

    if (sum(!is.na(valeurs_ssc)) < 2 ||
        sum(!is.na(valeurs_hc)) < 2) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      Log2FC = mean(valeurs_ssc, na.rm = TRUE) -
        mean(valeurs_hc, na.rm = TRUE),
      PValue = t.test(valeurs_ssc, valeurs_hc)$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$PValue, method = "BH")
  resultats$LogP <- -log10(resultats$PValue)

  resultats$Statut <- "Non significatif"
  resultats$Statut[
    resultats$PValue < 0.05 & resultats$Log2FC > 0.5
  ] <- "Up-régulé dans SSc"
  resultats$Statut[
    resultats$PValue < 0.05 & resultats$Log2FC < -0.5
  ] <- "Down-régulé dans SSc"

  top_labels <- resultats %>%
    filter(PValue < 0.05, abs(Log2FC) > 0.5) %>%
    arrange(PValue) %>%
    head(10)

  p_volcano <- ggplot(
    resultats,
    aes(x = Log2FC, y = LogP, fill = Statut, color = Statut)
  ) +
    geom_point(shape = 21, alpha = 0.85, size = 2.5, stroke = 0.8) +
    geom_text_repel(
      data = top_labels,
      aes(label = Proteine),
      color = "black",
      fontface = "bold",
      size = 3
    ) +
    scale_fill_manual(
      values = c(
        "Non significatif" = "grey80",
        "Up-régulé dans SSc" = "#F5A3A3",
        "Down-régulé dans SSc" = "#A3D9C9"
      )
    ) +
    scale_color_manual(
      values = c(
        "Non significatif" = "grey50",
        "Up-régulé dans SSc" = "#9A1F1F",
        "Down-régulé dans SSc" = "#2D6A4F"
      )
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey40"
    ) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot - SSc vs HC (", nom_matrice, ")"),
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  print(p_volcano)
  resultats
}

stats_ev_final <- analyser_ssc_vs_hc(mon_data_ev_clean, "Matrice EV")
stats_free_final <- analyser_ssc_vs_hc(mon_data_free_clean, "Matrice FREE")

# ------------------------------------------------------------------------------
# 10. Top 15 des protéines surexprimées dans les SSc
# ------------------------------------------------------------------------------

extraire_top_surexprimees <- function(stats_df, nom_matrice, n = 15) {
  stats_df %>%
    filter(PValue < 0.05, Log2FC > 0) %>%
    arrange(FDR) %>%
    head(n) %>%
    select(Proteine, Log2FC, PValue, FDR) %>%
    mutate(Matrice = nom_matrice)
}

top_up_ev <- extraire_top_surexprimees(stats_ev_final, "EV")
top_up_free <- extraire_top_surexprimees(stats_free_final, "FREE")

print(top_up_ev)
print(top_up_free)

generer_barplot_surexpression <- function(stats_df, nom_matrice, n = 15) {
  df_plot <- stats_df %>%
    filter(PValue < 0.05, Log2FC > 0) %>%
    arrange(desc(Log2FC)) %>%
    head(n)

  if (nrow(df_plot) == 0) {
    return(NULL)
  }

  df_plot$Proteine <- factor(
    df_plot$Proteine,
    levels = rev(df_plot$Proteine)
  )

  ggplot(df_plot, aes(x = Log2FC, y = Proteine, fill = Log2FC)) +
    geom_col(
      color = "#9A1F1F",
      linewidth = 0.6,
      width = 0.7,
      alpha = 0.85
    ) +
    scale_fill_gradient(
      low = "#FFE3E3",
      high = "#F5A3A3"
    ) +
    geom_text(
      aes(
        label = paste0(
          "FDR = ",
          formatC(FDR, format = "e", digits = 1)
        )
      ),
      hjust = 0,
      size = 2.8,
      fontface = "italic"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.35))
    ) +
    theme_minimal() +
    labs(
      title = paste0(
        "Top ", n,
        " protéines surexprimées (SSc vs HC) - ",
        nom_matrice
      ),
      x = "Log2 Fold Change",
      y = ""
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

p_up_ev <- generer_barplot_surexpression(stats_ev_final, "EV")
p_up_free <- generer_barplot_surexpression(stats_free_final, "FREE")

if (!is.null(p_up_ev) && !is.null(p_up_free)) {
  grid.arrange(p_up_ev, p_up_free, ncol = 2)
}

# ------------------------------------------------------------------------------
# 11. Analyse SSc diffuse cutanée vs SSc limitée cutanée
# ------------------------------------------------------------------------------

analyser_diffuse_vs_limitee <- function(data_clean, nom_matrice) {
  df_sub <- data_clean %>%
    filter(
      Groupe_Detail %in% c(
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      )
    )

  cols_exclure <- c(
    "sexe", "sex", "severite",
    "evolution-SSc", "connectivite"
  )

  mat_prot <- df_sub %>%
    select(where(is.numeric)) %>%
    select(-any_of(cols_exclure))

  groupes <- df_sub$Groupe_Detail

  resultats <- lapply(colnames(mat_prot), function(prot) {
    valeurs_diffuse <- mat_prot[[prot]][
      groupes == "SSc diffuse cutaneous"
    ]
    valeurs_limitee <- mat_prot[[prot]][
      groupes == "SSc Limited cutaneous"
    ]

    if (sum(!is.na(valeurs_diffuse)) < 2 ||
        sum(!is.na(valeurs_limitee)) < 2) {
      return(NULL)
    }

    test <- wilcox.test(
      valeurs_diffuse,
      valeurs_limitee,
      exact = FALSE
    )

    data.frame(
      Proteine = prot,
      Log2FC = log2(
        mean(valeurs_diffuse, na.rm = TRUE) /
          mean(valeurs_limitee, na.rm = TRUE)
      ),
      PValue = test$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$PValue, method = "BH")

  resultats$Statut <- case_when(
    resultats$PValue < 0.05 & resultats$Log2FC > 0.5 ~
      "Up dans Diffuse",
    resultats$PValue < 0.05 & resultats$Log2FC < -0.5 ~
      "Down dans Diffuse",
    TRUE ~ "Non significatif"
  )

  p_volcano <- ggplot(
    resultats,
    aes(x = Log2FC, y = -log10(PValue), color = Statut)
  ) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_color_manual(
      values = c(
        "Up dans Diffuse" = "#FF4B4B",
        "Down dans Diffuse" = "#4B9CD3",
        "Non significatif" = "grey70"
      )
    ) +
    geom_text_repel(
      data = filter(
        resultats,
        PValue < 0.01,
        abs(Log2FC) > 0.5
      ),
      aes(label = Proteine),
      color = "black",
      size = 3,
      fontface = "bold",
      max.overlaps = 15
    ) +
    theme_minimal() +
    labs(
      title = paste0(
        "SSc diffuse vs SSc limitée - ",
        nom_matrice
      ),
      x = "Log2 Fold Change (dcSSc / lcSSc)",
      y = "-log10(p-value)"
    )

  print(p_volcano)
  resultats
}

stats_diff_ev <- analyser_diffuse_vs_limitee(
  mon_data_ev_clean, "EV"
)
stats_diff_free <- analyser_diffuse_vs_limitee(
  mon_data_free_clean, "FREE"
)

# ------------------------------------------------------------------------------
# 12. Analyse PID0 vs PID1 chez les patients SSc
# ------------------------------------------------------------------------------

analyser_pid0_vs_pid1 <- function(data_finale, nom_matrice) {
  df_sub <- data_finale %>%
    filter(
      Statut_Simple == "SSc",
      !is.na(severite),
      severite %in% c("PID0", "PID1")
    ) %>%
    mutate(
      PID_Clean = ifelse(
        severite == "PID1",
        "PID1 (Avec)",
        "PID0 (Sans)"
      )
    )

  idx_prot <- (
    which(colnames(df_sub) == "Statut_Simple") + 1
  ):(ncol(df_sub) - 1)

  noms_prot <- colnames(df_sub)[idx_prot]

  resultats <- lapply(noms_prot, function(prot) {
    valeurs_pid1 <- as.numeric(
      df_sub[[prot]][df_sub$severite == "PID1"]
    )
    valeurs_pid0 <- as.numeric(
      df_sub[[prot]][df_sub$severite == "PID0"]
    )

    if (sum(!is.na(valeurs_pid1)) < 2 ||
        sum(!is.na(valeurs_pid0)) < 2) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      Log2FC = mean(valeurs_pid1, na.rm = TRUE) -
        mean(valeurs_pid0, na.rm = TRUE),
      PValue = t.test(valeurs_pid1, valeurs_pid0)$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$PValue, method = "BH")
  resultats$LogP <- -log10(resultats$PValue)

  resultats$Statut <- case_when(
    resultats$PValue < 0.05 & resultats$Log2FC > 0.3 ~
      "Surexprimé dans PID1",
    resultats$PValue < 0.05 & resultats$Log2FC < -0.3 ~
      "Surexprimé dans PID0",
    TRUE ~ "Non significatif"
  )

  top_labels <- resultats %>%
    filter(PValue < 0.05, abs(Log2FC) > 0.3) %>%
    arrange(PValue) %>%
    head(10)

  p_volcano <- ggplot(
    resultats,
    aes(x = Log2FC, y = LogP, fill = Statut, color = Statut)
  ) +
    geom_point(shape = 21, alpha = 0.9, size = 2.5) +
    geom_text_repel(
      data = top_labels,
      aes(label = Proteine),
      color = "black",
      fontface = "bold",
      size = 3
    ) +
    scale_fill_manual(
      values = c(
        "Non significatif" = "grey80",
        "Surexprimé dans PID1" = "#FF0000",
        "Surexprimé dans PID0" = "#0000FF"
      )
    ) +
    scale_color_manual(
      values = c(
        "Non significatif" = "grey50",
        "Surexprimé dans PID1" = "#8B0000",
        "Surexprimé dans PID0" = "#00008B"
      )
    ) +
    geom_vline(
      xintercept = c(-0.3, 0.3),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey40"
    ) +
    theme_minimal() +
    labs(
      title = paste0("PID1 vs PID0 - ", nom_matrice),
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    )

  print(p_volcano)
  resultats
}

stats_pid_ev <- analyser_pid0_vs_pid1(mon_data_ev_clean, "EV")
stats_pid_free <- analyser_pid0_vs_pid1(mon_data_free_clean, "FREE")

# ------------------------------------------------------------------------------
# 13. Analyse Healthy controls vs VEDOSS
# ------------------------------------------------------------------------------

analyser_vedoss_vs_hc <- function(data_clean, nom_matrice) {
  df_sub <- data_clean %>%
    filter(
      Groupe_Detail %in% c(
        "VEDOSS",
        "Healthy controls"
      )
    )

  noms_prot <- df_sub %>%
    select(where(is.numeric)) %>%
    colnames()

  resultats <- lapply(noms_prot, function(prot) {
    val_v <- as.numeric(
      df_sub[[prot]][df_sub$Groupe_Detail == "VEDOSS"]
    )
    val_h <- as.numeric(
      df_sub[[prot]][df_sub$Groupe_Detail == "Healthy controls"]
    )

    if (sum(!is.na(val_v)) < 2 ||
        sum(!is.na(val_h)) < 2) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      Log2FC = mean(val_v, na.rm = TRUE) -
        mean(val_h, na.rm = TRUE),
      PValue = t.test(val_v, val_h)$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$PValue, method = "BH")
  resultats$LogP <- -log10(resultats$PValue)

  resultats$Statut <- case_when(
    resultats$PValue < 0.05 & resultats$Log2FC > 0.3 ~
      "Surexprimé (VEDOSS)",
    resultats$PValue < 0.05 & resultats$Log2FC < -0.3 ~
      "Sous-exprimé (VEDOSS)",
    TRUE ~ "Non significatif"
  )

  top_labels <- resultats %>%
    filter(PValue < 0.05, abs(Log2FC) > 0.3) %>%
    arrange(PValue) %>%
    head(10)

  p <- ggplot(
    resultats,
    aes(x = Log2FC, y = LogP, fill = Statut, color = Statut)
  ) +
    geom_point(shape = 21, alpha = 0.8, size = 3) +
    geom_text_repel(
      data = top_labels,
      aes(label = Proteine),
      color = "black",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "Surexprimé (VEDOSS)" = "#FF0000",
        "Sous-exprimé (VEDOSS)" = "#0000FF",
        "Non significatif" = "grey85"
      )
    ) +
    scale_color_manual(
      values = c(
        "Surexprimé (VEDOSS)" = "#8B0000",
        "Sous-exprimé (VEDOSS)" = "#00008B",
        "Non significatif" = "grey50"
      )
    ) +
    theme_minimal() +
    labs(
      title = paste0("VEDOSS vs HC - ", nom_matrice),
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    )

  print(p)
  resultats
}

stats_vedoss_hc_ev <- analyser_vedoss_vs_hc(
  mon_data_ev_clean, "EV"
)
stats_vedoss_hc_free <- analyser_vedoss_vs_hc(
  mon_data_free_clean, "FREE"
)

# ------------------------------------------------------------------------------
# 14. Analyse VEDOSS vs SSc
# ------------------------------------------------------------------------------

analyser_vedoss_vs_ssc <- function(data_clean, nom_matrice) {
  df_sub <- data_clean %>%
    filter(
      Groupe_Detail %in% c(
        "VEDOSS",
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      )
    ) %>%
    mutate(
      Groupe_Comparaison = ifelse(
        Groupe_Detail == "VEDOSS",
        "VEDOSS",
        "SSc"
      )
    )

  noms_prot <- df_sub %>%
    select(where(is.numeric)) %>%
    colnames()

  resultats <- lapply(noms_prot, function(prot) {
    val_v <- as.numeric(
      df_sub[[prot]][df_sub$Groupe_Comparaison == "VEDOSS"]
    )
    val_s <- as.numeric(
      df_sub[[prot]][df_sub$Groupe_Comparaison == "SSc"]
    )

    if (sum(!is.na(val_v)) < 2 ||
        sum(!is.na(val_s)) < 2) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      Log2FC = mean(val_v, na.rm = TRUE) -
        mean(val_s, na.rm = TRUE),
      PValue = t.test(val_v, val_s)$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$PValue, method = "BH")
  resultats$LogP <- -log10(resultats$PValue)

  resultats$Statut <- case_when(
    resultats$PValue < 0.05 & resultats$Log2FC > 0.3 ~
      "Surexprimé (VEDOSS)",
    resultats$PValue < 0.05 & resultats$Log2FC < -0.3 ~
      "Sous-exprimé (VEDOSS)",
    TRUE ~ "Non significatif"
  )

  top_labels <- resultats %>%
    filter(PValue < 0.05, abs(Log2FC) > 0.3) %>%
    arrange(PValue) %>%
    head(10)

  p <- ggplot(
    resultats,
    aes(x = Log2FC, y = LogP, fill = Statut, color = Statut)
  ) +
    geom_point(shape = 21, alpha = 0.8, size = 3) +
    geom_text_repel(
      data = top_labels,
      aes(label = Proteine),
      color = "black",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "Surexprimé (VEDOSS)" = "#FF0000",
        "Sous-exprimé (VEDOSS)" = "#0000FF",
        "Non significatif" = "grey85"
      )
    ) +
    scale_color_manual(
      values = c(
        "Surexprimé (VEDOSS)" = "#8B0000",
        "Sous-exprimé (VEDOSS)" = "#00008B",
        "Non significatif" = "grey50"
      )
    ) +
    theme_minimal() +
    labs(
      title = paste0("VEDOSS vs SSc - ", nom_matrice),
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    )

  print(p)
  resultats
}

stats_vedoss_ssc_ev <- analyser_vedoss_vs_ssc(
  mon_data_ev_clean, "EV"
)
stats_vedoss_ssc_free <- analyser_vedoss_vs_ssc(
  mon_data_free_clean, "FREE"
)

# ------------------------------------------------------------------------------
# 15. Analyse à trois groupes : Raynaud, VEDOSS et SSc
# ------------------------------------------------------------------------------

generer_analyse_3groupes <- function(data, nom_matrice) {
  df_multi <- data %>%
    filter(
      Groupe_Detail %in% c(
        "Idiopathic Raynaud",
        "Raynaud connective",
        "VEDOSS",
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      )
    ) %>%
    mutate(
      Groupe = case_when(
        Groupe_Detail %in% c(
          "Idiopathic Raynaud",
          "Raynaud connective"
        ) ~ "Raynaud",
        Groupe_Detail == "VEDOSS" ~ "VEDOSS",
        TRUE ~ "SSc"
      ),
      Groupe = factor(
        Groupe,
        levels = c("Raynaud", "VEDOSS", "SSc")
      )
    )

  noms_prot <- df_multi %>%
    select(where(is.numeric)) %>%
    colnames()

  res_stats <- lapply(noms_prot, function(prot) {
    df_temp <- data.frame(
      Groupe = df_multi$Groupe,
      Valeur = as.numeric(df_multi[[prot]])
    ) %>%
      filter(!is.na(Valeur))

    if (
      length(unique(df_temp$Groupe)) < 3 ||
      any(table(df_temp$Groupe) < 2)
    ) {
      return(NULL)
    }

    fit <- try(
      aov(Valeur ~ Groupe, data = df_temp),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      PValue = summary(fit)[[1]][["Pr(>F)"]][1]
    )
  })

  stats_3gr <- bind_rows(res_stats)
  stats_3gr$FDR <- p.adjust(stats_3gr$PValue, method = "BH")

  top10 <- stats_3gr %>%
    arrange(FDR) %>%
    head(10) %>%
    pull(Proteine)

  df_long <- df_multi %>%
    select(Groupe, all_of(top10)) %>%
    pivot_longer(
      -Groupe,
      names_to = "Proteine",
      values_to = "Expression"
    ) %>%
    filter(!is.na(Expression))

  p <- ggplot(
    df_long,
    aes(x = Groupe, y = Expression, fill = Groupe)
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.7,
      outlier.shape = NA
    ) +
    geom_jitter(width = 0.1, alpha = 0.3) +
    stat_compare_means(
      comparisons = list(
        c("Raynaud", "VEDOSS"),
        c("VEDOSS", "SSc"),
        c("Raynaud", "SSc")
      ),
      method = "t.test",
      label = "p.signif"
    ) +
    facet_wrap(~Proteine, scales = "free_y", ncol = 5) +
    scale_fill_manual(
      values = c(
        "Raynaud" = "#FFA500",
        "VEDOSS" = "#FDFD96",
        "SSc" = "#E6B8B8"
      )
    ) +
    theme_minimal() +
    labs(
      title = paste0("Top 10 protéines - ", nom_matrice),
      x = "",
      y = "Expression"
    )

  list(plot = p, stats = stats_3gr)
}

res_3groupes_ev <- generer_analyse_3groupes(
  mon_data_ev_clean, "EV"
)
res_3groupes_free <- generer_analyse_3groupes(
  mon_data_free_clean, "FREE"
)

print(res_3groupes_ev$plot)
print(res_3groupes_free$plot)

# ------------------------------------------------------------------------------
# 16. Analyse VEDOSS vs Raynaud
# ------------------------------------------------------------------------------

analyser_vedoss_vs_raynaud <- function(data, nom_matrice) {
  df_comp <- data %>%
    filter(
      Groupe_Detail %in% c(
        "VEDOSS",
        "Idiopathic Raynaud",
        "Raynaud connective"
      )
    ) %>%
    mutate(
      Groupe = factor(
        ifelse(
          Groupe_Detail == "VEDOSS",
          "VEDOSS",
          "Raynaud"
        ),
        levels = c("Raynaud", "VEDOSS")
      )
    )

  noms_prot <- df_comp %>%
    select(where(is.numeric)) %>%
    colnames()

  resultats <- lapply(noms_prot, function(prot) {
    val_v <- df_comp[[prot]][df_comp$Groupe == "VEDOSS"]
    val_r <- df_comp[[prot]][df_comp$Groupe == "Raynaud"]

    val_v <- val_v[!is.na(val_v)]
    val_r <- val_r[!is.na(val_r)]

    if (length(val_v) < 2 || length(val_r) < 2) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      Log2FC = mean(val_v) - mean(val_r),
      PValue = t.test(val_v, val_r)$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$PValue, method = "BH")

  top10 <- resultats %>%
    arrange(PValue) %>%
    head(10) %>%
    pull(Proteine)

  df_long <- df_comp %>%
    select(Groupe, all_of(top10)) %>%
    pivot_longer(
      -Groupe,
      names_to = "Proteine",
      values_to = "Expression"
    )

  p_box <- ggplot(
    df_long,
    aes(x = Groupe, y = Expression, fill = Groupe)
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.7,
      outlier.shape = NA
    ) +
    geom_jitter(width = 0.1, alpha = 0.3) +
    stat_compare_means(
      method = "t.test",
      label = "p.signif"
    ) +
    facet_wrap(~Proteine, scales = "free_y", ncol = 5) +
    scale_fill_manual(
      values = c(
        "VEDOSS" = "#FDFD96",
        "Raynaud" = "#FFA500"
      )
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = paste0(
        "Top 10 protéines discriminantes - VEDOSS vs Raynaud (",
        nom_matrice, ")"
      ),
      x = "",
      y = "Expression"
    )

  resultats$Statut <- case_when(
    resultats$PValue < 0.05 & resultats$Log2FC > 0.5 ~ "Up",
    resultats$PValue < 0.05 & resultats$Log2FC < -0.5 ~ "Down",
    TRUE ~ "NS"
  )

  p_volcano <- ggplot(
    resultats,
    aes(x = Log2FC, y = -log10(PValue), color = Statut)
  ) +
    geom_point(alpha = 0.6) +
    geom_text_repel(
      data = filter(resultats, Statut != "NS"),
      aes(label = Proteine),
      size = 3
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c(
        "Up" = "red",
        "Down" = "#007BFF",
        "NS" = "grey"
      )
    ) +
    theme_minimal() +
    labs(
      title = paste0("VEDOSS vs Raynaud - ", nom_matrice),
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    )

  print(p_box)
  print(p_volcano)

  resultats
}

stats_vedoss_raynaud_ev <- analyser_vedoss_vs_raynaud(
  mon_data_ev_clean, "EV"
)
stats_vedoss_raynaud_free <- analyser_vedoss_vs_raynaud(
  mon_data_free_clean, "FREE"
)

# ------------------------------------------------------------------------------
# 17. Analyse ANOVA : HC, VEDOSS et SSc
# ------------------------------------------------------------------------------

generer_figure_top5_stats <- function(data, nom_matrice, n = 5) {
  df_select <- data %>%
    filter(
      Groupe_Detail %in% c(
        "Healthy controls",
        "VEDOSS",
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      )
    ) %>%
    mutate(
      Groupe_Simplifie = factor(
        case_when(
          Groupe_Detail == "Healthy controls" ~ "HC",
          Groupe_Detail == "VEDOSS" ~ "VEDOSS",
          TRUE ~ "SSc"
        ),
        levels = c("HC", "VEDOSS", "SSc")
      )
    ) %>%
    select(Groupe_Simplifie, where(is.numeric))

  noms_prot <- setdiff(
    colnames(df_select),
    "Groupe_Simplifie"
  )

  pvals_anova <- sapply(noms_prot, function(prot) {
    df_temp <- df_select[, c("Groupe_Simplifie", prot)]
    df_temp <- df_temp[complete.cases(df_temp), ]

    if (length(unique(df_temp$Groupe_Simplifie)) < 2) {
      return(1)
    }

    summary(
      aov(df_temp[[prot]] ~ df_temp$Groupe_Simplifie)
    )[[1]][["Pr(>F)"]][1]
  })

  top_noms <- names(sort(pvals_anova))[seq_len(
    min(n, length(pvals_anova))
  )]

  df_long <- df_select %>%
    select(Groupe_Simplifie, all_of(top_noms)) %>%
    pivot_longer(
      -Groupe_Simplifie,
      names_to = "Proteine",
      values_to = "Expression"
    ) %>%
    filter(!is.na(Expression))

  ggplot(
    df_long,
    aes(
      x = Groupe_Simplifie,
      y = Expression,
      fill = Groupe_Simplifie
    )
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.7,
      outlier.shape = NA
    ) +
    geom_jitter(width = 0.1, alpha = 0.3) +
    stat_compare_means(
      method = "t.test",
      comparisons = list(
        c("HC", "VEDOSS"),
        c("VEDOSS", "SSc"),
        c("HC", "SSc")
      ),
      label = "p.signif"
    ) +
    facet_wrap(~Proteine, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = c(
        "HC" = "#A3D9C9",
        "VEDOSS" = "#FDFD96",
        "SSc" = "#E6B8B8"
      )
    ) +
    theme_minimal() +
    labs(
      title = paste0("Top ", n, " protéines - ", nom_matrice),
      x = "",
      y = "Expression"
    ) +
    theme(legend.position = "none")
}

print(generer_figure_top5_stats(mon_data_ev_clean, "EV"))
print(generer_figure_top5_stats(mon_data_free_clean, "FREE"))

# ------------------------------------------------------------------------------
# 18. Comparaison des protéines significatives entre EV et FREE
# ------------------------------------------------------------------------------

sig_ev <- stats_ev_final %>%
  filter(FDR < 0.05) %>%
  pull(Proteine) %>%
  trimws()

sig_free <- stats_free_final %>%
  filter(FDR < 0.05) %>%
  pull(Proteine) %>%
  trimws()

proteines_communes <- intersect(sig_ev, sig_free)

tableau_commun <- data.frame(
  Proteine = proteines_communes,
  FDR_EV = stats_ev_final$FDR[
    match(proteines_communes, stats_ev_final$Proteine)
  ],
  FDR_FREE = stats_free_final$FDR[
    match(proteines_communes, stats_free_final$Proteine)
  ]
)

print(tableau_commun)

tableau_long <- tableau_commun %>%
  pivot_longer(
    cols = c(FDR_EV, FDR_FREE),
    names_to = "Compartiment",
    values_to = "FDR"
  )

dotplot_commun <- ggplot(
  tableau_long,
  aes(
    x = Compartiment,
    y = reorder(Proteine, FDR),
    size = -log10(FDR),
    color = -log10(FDR)
  )
) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Protéines communes significatives : EV vs FREE",
    subtitle = "Comparaison SSc vs Healthy controls",
    x = "",
    y = "Protéines",
    size = "-log10(FDR)",
    color = "-log10(FDR)"
  )

print(dotplot_commun)

venn_plot <- venn.diagram(
  x = list(EV = sig_ev, FREE = sig_free),
  category.names = c("EV significatif", "FREE significatif"),
  filename = NULL,
  fill = c("#69b3a2", "#404080"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

grid.newpage()
grid.draw(venn_plot)

# Top 25 protéines significatives dans au moins un compartiment.
tableau_comparatif <- full_join(
  stats_ev_final %>%
    select(Proteine, FDR_EV = FDR),
  stats_free_final %>%
    select(Proteine, FDR_FREE = FDR),
  by = "Proteine"
) %>%
  filter(FDR_EV < 0.05 | FDR_FREE < 0.05)

top25_prot <- tableau_comparatif %>%
  mutate(
    importance = rowMeans(
      cbind(-log10(FDR_EV), -log10(FDR_FREE)),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(importance)) %>%
  head(25) %>%
  pull(Proteine)

proteines_a_garder <- unique(
  c(top25_prot, "SFTPB", "PRDX2")
)

tableau_fdr_final <- data.frame(
  Proteine = proteines_a_garder
) %>%
  left_join(
    stats_ev_final %>%
      select(Proteine, FDR_EV = FDR),
    by = "Proteine"
  ) %>%
  left_join(
    stats_free_final %>%
      select(Proteine, FDR_FREE = FDR),
    by = "Proteine"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

print(tableau_fdr_final)

# ------------------------------------------------------------------------------
# 19. Protéines d'intérêt ciblées
# ------------------------------------------------------------------------------

prot_ev_cibles <- c(
  "BTD", "F2RL3", "F8", "ITIH1", "ITIH2",
  "LBP", "NAP1L1", "PON1", "PPP2R5C", "VWF"
)

prot_free_cibles <- c(
  "APCS", "GPX3", "GSN", "ITIH1", "ITIH2",
  "LBP", "LDHA", "PRG4", "SFTPB", "VCL"
)

get_fdr_specifique <- function(data, liste_prot) {
  df_clean <- data %>%
    filter(
      Groupe_Detail %in% c(
        "Idiopathic Raynaud",
        "Raynaud connective",
        "VEDOSS",
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      )
    ) %>%
    mutate(
      Groupe = case_when(
        Groupe_Detail %in% c(
          "Idiopathic Raynaud",
          "Raynaud connective"
        ) ~ "Raynaud",
        Groupe_Detail == "VEDOSS" ~ "VEDOSS",
        TRUE ~ "SSc"
      )
    )

  res_list <- lapply(liste_prot, function(prot) {
    if (!prot %in% colnames(df_clean)) {
      return(NULL)
    }

    p_rv <- t.test(
      df_clean[[prot]][df_clean$Groupe == "Raynaud"],
      df_clean[[prot]][df_clean$Groupe == "VEDOSS"]
    )$p.value

    p_vs <- t.test(
      df_clean[[prot]][df_clean$Groupe == "VEDOSS"],
      df_clean[[prot]][df_clean$Groupe == "SSc"]
    )$p.value

    p_sr <- t.test(
      df_clean[[prot]][df_clean$Groupe == "SSc"],
      df_clean[[prot]][df_clean$Groupe == "Raynaud"]
    )$p.value

    data.frame(
      Proteine = prot,
      P_Raynaud_VEDOSS = p_rv,
      P_VEDOSS_SSc = p_vs,
      P_SSc_Raynaud = p_sr
    )
  })

  res <- bind_rows(res_list)

  if (nrow(res) == 0) {
    return(res)
  }

  res$FDR_Raynaud_VEDOSS <- p.adjust(
    res$P_Raynaud_VEDOSS,
    method = "BH"
  )
  res$FDR_VEDOSS_SSc <- p.adjust(
    res$P_VEDOSS_SSc,
    method = "BH"
  )
  res$FDR_SSc_Raynaud <- p.adjust(
    res$P_SSc_Raynaud,
    method = "BH"
  )

  res
}

resultats_cibles_ev <- get_fdr_specifique(
  mon_data_ev_clean, prot_ev_cibles
)
resultats_cibles_free <- get_fdr_specifique(
  mon_data_free_clean, prot_free_cibles
)

print(resultats_cibles_ev)
print(resultats_cibles_free)

# ------------------------------------------------------------------------------
# 20. Analyse ciblée de FN1
# ------------------------------------------------------------------------------

if ("FN1" %in% colnames(mon_data_ev_clean)) {
  df_fn1 <- mon_data_ev_clean %>%
    filter(!is.na(Groupe_Detail), !is.na(FN1))

  groupes_fn1 <- c(
    "Healthy controls",
    "VEDOSS",
    "SSc diffuse cutaneous",
    "SSc Limited cutaneous"
  )

  df_fn1 <- df_fn1 %>%
    filter(Groupe_Detail %in% groupes_fn1)

  comparaisons_fn1 <- list(
    c("Healthy controls", "VEDOSS"),
    c("VEDOSS", "SSc diffuse cutaneous"),
    c("VEDOSS", "SSc Limited cutaneous")
  )

  p_fn1 <- sapply(comparaisons_fn1, function(pair) {
    t.test(
      df_fn1$FN1[df_fn1$Groupe_Detail == pair[1]],
      df_fn1$FN1[df_fn1$Groupe_Detail == pair[2]]
    )$p.value
  })

  fdr_fn1 <- p.adjust(p_fn1, method = "BH")

  p_fn1_ev <- ggplot(
    df_fn1,
    aes(x = Groupe_Detail, y = FN1, fill = Groupe_Detail)
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.7,
      outlier.shape = NA
    ) +
    geom_jitter(width = 0.1, alpha = 0.3) +
    stat_compare_means(
      comparisons = comparaisons_fn1,
      method = "t.test",
      label = "p.signif"
    ) +
    scale_fill_manual(values = palette_detail) +
    theme_minimal() +
    labs(
      title = "Fibronectine (FN1) dans les EV",
      x = "",
      y = "Expression (Log2)"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "bold"
      ),
      legend.position = "none"
    )

  print(p_fn1_ev)

  print(
    data.frame(
      Comparaison = c(
        "Healthy controls vs VEDOSS",
        "VEDOSS vs SSc diffuse cutaneous",
        "VEDOSS vs SSc Limited cutaneous"
      ),
      P_Value = p_fn1,
      FDR = fdr_fn1
    )
  )
}

# ------------------------------------------------------------------------------
# 21. Évolution de la SSc
# ------------------------------------------------------------------------------

analyser_evolution_ssc <- function(data) {
  df_ssc <- data %>%
    filter(
      Groupe_Detail %in% c(
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      ),
      !is.na(`evolution-SSc`)
    )

  prot_cols <- df_ssc %>%
    select(where(is.numeric)) %>%
    colnames()

  resultats <- lapply(prot_cols, function(prot) {
    test <- try(
      t.test(
        df_ssc[[prot]] ~ df_ssc$`evolution-SSc`
      ),
      silent = TRUE
    )

    if (inherits(test, "try-error")) {
      return(NULL)
    }

    data.frame(
      Proteine = prot,
      P_Value = test$p.value
    )
  })

  resultats <- bind_rows(resultats)
  resultats$FDR <- p.adjust(resultats$P_Value, method = "BH")

  resultats %>%
    arrange(FDR)
}

discrim_ev_evolution <- analyser_evolution_ssc(
  mon_data_ev_clean
)
discrim_free_evolution <- analyser_evolution_ssc(
  mon_data_free_clean
)

print(
  head(
    discrim_ev_evolution %>%
      filter(FDR < 0.05),
    10
  )
)

print(
  head(
    discrim_free_evolution %>%
      filter(FDR < 0.05),
    10
  )
)

plot_evolution <- function(data, stats_res, nom_matrice, n = 5) {
  top_prot <- head(stats_res$Proteine, n)

  df_long <- data %>%
    filter(
      Groupe_Detail %in% c(
        "SSc diffuse cutaneous",
        "SSc Limited cutaneous"
      ),
      !is.na(`evolution-SSc`)
    ) %>%
    select(`evolution-SSc`, all_of(top_prot)) %>%
    pivot_longer(
      -`evolution-SSc`,
      names_to = "Proteine",
      values_to = "Expression"
    )

  ggplot(
    df_long,
    aes(
      x = `evolution-SSc`,
      y = Expression,
      fill = `evolution-SSc`
    )
  ) +
    geom_boxplot(alpha = 0.7) +
    stat_compare_means(
      method = "t.test",
      label = "p.signif"
    ) +
    facet_wrap(~Proteine, scales = "free_y") +
    theme_minimal() +
    labs(
      title = paste0(
        "Protéines discriminantes selon l'évolution de la SSc - ",
        nom_matrice
      ),
      x = "Durée de la maladie",
      y = "Expression"
    )
}

print(
  plot_evolution(
    mon_data_ev_clean,
    discrim_ev_evolution,
    "EV"
  )
)

print(
  plot_evolution(
    mon_data_free_clean,
    discrim_free_evolution,
    "FREE"
  )
)

# ------------------------------------------------------------------------------
# 22. Comorbidités / connectivités
# ------------------------------------------------------------------------------

comorbidites <- c("lupus", "sharp", "dermatomyosite")

tableau_croise <- mon_data_free_clean %>%
  filter(connectivite %in% comorbidites) %>%
  group_by(Groupe_Detail, connectivite) %>%
  summarise(
    Nombre = n(),
    .groups = "drop"
  )

print(tableau_croise)

ggplot(
  tableau_croise,
  aes(x = Groupe_Detail, y = Nombre, fill = connectivite)
) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(
    title = "Association entre type de SSc et comorbidités",
    x = "Type de SSc",
    y = "Nombre de patients",
    fill = "Comorbidité"
  )

association_prot_maladie <- mon_data_free_clean %>%
  filter(connectivite %in% comorbidites) %>%
  group_by(connectivite) %>%
  summarise(
    across(
      c(APOA2, APOC1, F5, PNP, RNASE4),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

print(
  pivot_longer(
    association_prot_maladie,
    -connectivite,
    names_to = "Proteine",
    values_to = "Expression_Moyenne"
  )
)

data_heatmap <- association_prot_maladie %>%
  as.data.frame()

rownames(data_heatmap) <- data_heatmap$connectivite
data_heatmap$connectivite <- NULL

pheatmap(
  data_heatmap,
  scale = "column",
  main = "Signature protéique par comorbidité",
  color = colorRampPalette(
    c("blue", "white", "red")
  )(50)
)

ggplot(
  mon_data_free_clean %>%
    filter(connectivite %in% comorbidites),
  aes(x = Groupe_Detail, fill = connectivite)
) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(
    title = "Proportion des comorbidités selon le type de SSc",
    x = "Type de SSc",
    y = "Proportion",
    fill = "Comorbidité"
  )

# ------------------------------------------------------------------------------
# Fin du script
# ------------------------------------------------------------------------------
