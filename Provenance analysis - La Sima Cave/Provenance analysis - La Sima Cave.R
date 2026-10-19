# ==============================================================================
# Provenance analysis - La Sima Cave
# XRF data (oxides, wt%) from geological sources (CZ) and artefacts (SIM)
# ==============================================================================

# Packages
packages <- c("ggplot2", "dplyr", "tidyr", "compositions",
              "ggrepel", "RColorBrewer", "pheatmap", "psych")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}

# Output folder
output_dir <- if (dir.exists(file.path(Sys.getenv("USERPROFILE"), "Desktop"))) {
  file.path(Sys.getenv("USERPROFILE"), "Desktop")
} else if (dir.exists(file.path(Sys.getenv("USERPROFILE"), "Escritorio"))) {
  file.path(Sys.getenv("USERPROFILE"), "Escritorio")
} else if (dir.exists("~/Desktop")) {
  "~/Desktop"
} else {
  getwd()
}
cat("Files will be saved to:", output_dir, "\n")


# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------

data <- data.frame(
  Sample = c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37),
  Na2O   = c(0.0774, 0.417, 0.192, 0.272, 0, 0.0464, 0.105, 0.297, 0, 0.0466, 0.0365, 0, 0),
  MgO    = c(0.0561, 0.407, 0.373, 0.133, 0.0633, 0.0214, 0.199, 0.146, 0.072, 0.245, 0.051, 0.127, 0.102),
  Al2O3  = c(0.58, 10.7, 5.91, 1.69, 0.485, 0.335, 2.93, 1.07, 0.566, 4.28, 0.658, 1.23, 0.919),
  SiO2   = c(79.7, 44.1, 55, 69.4, 78, 67.4, 73.4, 63.5, 57.8, 56.7, 81.1, 60, 69),
  P2O5   = c(0.0297, 0.141, 0.154, 0.0666, 0.0863, 0.909, 0.0153, 0.283, 0.0478, 0.113, 0.044, 0.168, 0.0831),
  SO3    = c(0.0287, 0.101, 0.0516, 0.0434, 0.0235, 0.0224, 0.017, 0.122, 0.0156, 0.0241, 0.0303, 0.029, 0.0286),
  Cl     = c(0.495, 0.302, 0.0558, 0.142, 0.022, 0.0364, 0.106, 0.147, 0.0167, 0.0214, 0, 0, 0.0388),
  K2O    = c(0.0383, 0.73, 1.02, 0.216, 0.149, 0.0588, 0.126, 0.185, 0.0746, 0.248, 0.0793, 0.211, 0.113),
  CaO    = c(0.0328, 0.113, 0.1, 0.0607, 1.46, 0.523, 0.0399, 1.97, 0.745, 0.531, 0.134, 0.706, 1.41),
  TiO2   = c(0, 0.0861, 0.0696, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  MnO    = c(0, 0.0103, 0.0136, 0.00462, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Fe2O3  = c(0.0478, 1.11, 0.937, 0.137, 0.1, 0.0388, 0.0558, 0.561, 0.0663, 0.56, 0.104, 0.137, 0.106),
  NiO    = c(0, 0.0012, 0, 0.00285, 0.00263, 0, 0, 0, 0, 0, 0, 0, 0),
  CuO    = c(0, 0, 0, 0.00287, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  ZnO    = c(0, 0.00326, 0.00208, 0.00131, 0, 0, 0, 0.0221, 0, 0.00419, 0, 2.93e-05, 0),
  As2O3  = c(0, 0, 0, 0, 0, 0.00103, 0, 0, 0, 0, 0, 0, 0),
  SrO    = c(0, 0.0011, 0, 0, 0, 0.00194, 0.00145, 0.00244, 0, 0.00189, 0, 0, NA),
  BaO    = c(0, 0.0313, 0, 0, 0, 0.0181, 0, 0.00317, 0, 0, 0, 0, 0),
  CO2    = c(19.3, 41.8, 36.1, 27.8, 19.6, 30.6, 23, 31.7, 40.6, 37.2, 17.8, 37.4, 28.2),
  Nomenclature = c("CZ2-1", "CZ1-1", "CZ1-3",
                   "SIM.22-619", "SIM.22-744", "SIM.22-863",
                   "CZ3-2", "SIM-21.126", "SIM-21.134",
                   "SIM-22.362", "SIM-22.590", "SIM-22.719", "SIM-22.753")
)


# ------------------------------------------------------------------------------
# 2. Preparation
# ------------------------------------------------------------------------------

elements <- c("Na2O", "MgO", "Al2O3", "SiO2", "P2O5", "SO3", "Cl",
              "K2O", "CaO", "TiO2", "MnO", "Fe2O3", "NiO", "CuO",
              "ZnO", "As2O3", "SrO", "BaO")

df <- data[, c("Sample", "Nomenclature", elements)]
df <- df[complete.cases(df[elements]), ]
df$Type <- ifelse(grepl("^CZ", df$Nomenclature), "Source", "Artefact")

# Replace zeros with half the minimum positive value per element
for (el in elements) {
  v <- df[[el]]
  min_pos <- min(v[v > 0], na.rm = TRUE)
  if (is.finite(min_pos)) df[[el]][df[[el]] == 0] <- min_pos / 2
}


# ------------------------------------------------------------------------------
# 3. Boxplots by element
# ------------------------------------------------------------------------------

df_long <- df %>%
  pivot_longer(cols = all_of(elements), names_to = "Element", values_to = "Concentration")

p_box <- ggplot(df_long, aes(x = Type, y = Concentration, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
              size = 1.2, alpha = 0.6, aes(color = Type)) +
  scale_fill_manual(values = c("Source" = "#2E86AB", "Artefact" = "#A23B72")) +
  scale_color_manual(values = c("Source" = "#2E86AB", "Artefact" = "#A23B72")) +
  facet_wrap(~Element, scales = "free_y", ncol = 6) +
  labs(title = "Elemental concentrations in sources and artefacts",
       x = "", y = "Concentration (wt%)", fill = "Type", color = "Type") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.position = "bottom")

ggsave(file.path(output_dir, "Boxplots_elements.png"),
       p_box, width = 14, height = 10, dpi = 300)


# ------------------------------------------------------------------------------
# 4. Compositional transformation (CLR)
# ------------------------------------------------------------------------------

mat <- as.matrix(df[, elements])
rownames(mat) <- df$Nomenclature
clr_data <- scale(clr(acomp(mat)))


# ------------------------------------------------------------------------------
# 5. PCA with all elements
# ------------------------------------------------------------------------------

pca_all <- prcomp(clr_data, center = FALSE, scale. = FALSE)
var_all <- round(100 * summary(pca_all)$importance[2, 1:5], 1)

R <- cor(clr_data)
cat("\n--- Bartlett's test of sphericity ---\n")
print(psych::cortest.bartlett(R, n = nrow(clr_data)))
cat("\n--- KMO ---\n")
print(psych::KMO(R))


# ------------------------------------------------------------------------------
# 6. PCA without P and Ca
# ------------------------------------------------------------------------------

elements_noPCa <- setdiff(elements, c("P2O5", "CaO"))
mat_noPCa <- as.matrix(df[, elements_noPCa])
rownames(mat_noPCa) <- df$Nomenclature
clr_noPCa <- scale(clr(acomp(mat_noPCa)))
pca_noPCa <- prcomp(clr_noPCa, center = FALSE, scale. = FALSE)
var_noPCa <- round(100 * summary(pca_noPCa)$importance[2, 1:5], 1)


# ------------------------------------------------------------------------------
# 7. Assign each artefact to the nearest source
# ------------------------------------------------------------------------------

assign_sources <- function(pca_res, n_pcs = 5) {
  coords <- as.data.frame(pca_res$x[, 1:n_pcs])
  coords$Sample <- rownames(coords)
  coords$Type <- df$Type[match(coords$Sample, df$Nomenclature)]
  
  sources <- coords[coords$Type == "Source", , drop = FALSE]
  artefacts <- coords[coords$Type == "Artefact", , drop = FALSE]
  
  d <- matrix(NA, nrow = nrow(artefacts), ncol = nrow(sources))
  rownames(d) <- artefacts$Sample
  colnames(d) <- sources$Sample
  
  for (i in seq_len(nrow(artefacts))) {
    for (j in seq_len(nrow(sources))) {
      d[i, j] <- sqrt(sum((artefacts[i, 1:n_pcs] - sources[j, 1:n_pcs])^2))
    }
  }
  
  data.frame(
    Artefact = rownames(d),
    Source = colnames(d)[apply(d, 1, which.min)],
    Distance = round(apply(d, 1, min), 3),
    row.names = NULL
  )
}

assign_all <- assign_sources(pca_all)
assign_noPCa <- assign_sources(pca_noPCa)

cat("\n--- Assignment with all elements ---\n"); print(assign_all)
cat("\n--- Assignment without P and Ca ---\n"); print(assign_noPCa)

write.csv(assign_all,   file.path(output_dir, "Assignment_all.csv"),  row.names = FALSE)
write.csv(assign_noPCa, file.path(output_dir, "Assignment_noPCa.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# 8. PCA biplot with all elements
# ------------------------------------------------------------------------------

coords <- as.data.frame(pca_all$x[, 1:2])
coords$Sample <- rownames(coords)
coords$Type <- df$Type[match(coords$Sample, df$Nomenclature)]

loadings <- as.data.frame(pca_all$rotation[, 1:2])
loadings$Element <- rownames(loadings)

p_biplot <- ggplot(coords, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(color = Type), size = 4) +
  geom_text_repel(aes(label = Sample, color = Type), size = 3.5,
                  box.padding = 0.5, max.overlaps = 20) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * 4, yend = PC2 * 4),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgreen", linewidth = 0.8) +
  geom_text(data = loadings,
            aes(x = PC1 * 4.5, y = PC2 * 4.5, label = Element),
            color = "darkgreen", size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("Source" = "#2E86AB", "Artefact" = "#A23B72")) +
  labs(title = "PCA biplot (all elements)",
       x = paste0("PC1 (", var_all[1], "%)"),
       y = paste0("PC2 (", var_all[2], "%)"),
       color = "Type") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "PCA_biplot.png"),
       p_biplot, width = 12, height = 9, dpi = 300)


# ------------------------------------------------------------------------------
# 9. PCA biplot without P and Ca
# ------------------------------------------------------------------------------

coords2 <- as.data.frame(pca_noPCa$x[, 1:2])
coords2$Sample <- rownames(coords2)
coords2$Type <- df$Type[match(coords2$Sample, df$Nomenclature)]

loadings2 <- as.data.frame(pca_noPCa$rotation[, 1:2])
loadings2$Element <- rownames(loadings2)

p_biplot2 <- ggplot(coords2, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(color = Type), size = 4) +
  geom_text_repel(aes(label = Sample, color = Type), size = 3.5,
                  box.padding = 0.5, max.overlaps = 20) +
  geom_segment(data = loadings2,
               aes(x = 0, y = 0, xend = PC1 * 4, yend = PC2 * 4),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgreen", linewidth = 0.8) +
  geom_text(data = loadings2,
            aes(x = PC1 * 4.5, y = PC2 * 4.5, label = Element),
            color = "darkgreen", size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("Source" = "#2E86AB", "Artefact" = "#A23B72")) +
  labs(title = "PCA biplot (without P and Ca)",
       x = paste0("PC1 (", var_noPCa[1], "%)"),
       y = paste0("PC2 (", var_noPCa[2], "%)"),
       color = "Type") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "PCA_biplot_noPCa.png"),
       p_biplot2, width = 12, height = 9, dpi = 300)


# ------------------------------------------------------------------------------
# 10. Assignment map (all elements)
# ------------------------------------------------------------------------------

source_names <- unique(df$Nomenclature[df$Type == "Source"])
source_colors <- setNames(brewer.pal(n = length(source_names), name = "Set1"), source_names)

arrows <- merge(coords, assign_all, by.x = "Sample", by.y = "Artefact")
source_coords <- coords[coords$Sample %in% unique(arrows$Source),
                        c("PC1", "PC2", "Sample")]
names(source_coords)[1:2] <- c("PC1_end", "PC2_end")
arrows <- merge(arrows, source_coords, by.x = "Source", by.y = "Sample")

p_map <- ggplot(coords, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(data = subset(coords, Type == "Source"),
             aes(color = Sample), size = 5, shape = 17) +
  geom_point(data = subset(coords, Type == "Artefact"),
             color = "gray40", size = 3.5, shape = 16) +
  geom_segment(data = arrows,
               aes(x = PC1, y = PC2, xend = PC1_end, yend = PC2_end,
                   color = Source),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1, alpha = 0.8) +
  geom_text_repel(data = subset(coords, Type == "Artefact"),
                  aes(label = Sample), size = 3.2,
                  box.padding = 0.5, segment.color = "grey50", max.overlaps = 20) +
  geom_text_repel(data = subset(coords, Type == "Source"),
                  aes(label = Sample, color = Sample), size = 4.5,
                  fontface = "bold", box.padding = 0.8, segment.color = NA) +
  scale_color_manual(values = source_colors) +
  labs(title = "Assignment of artefacts to sources",
       x = paste0("PC1 (", var_all[1], "%)"),
       y = paste0("PC2 (", var_all[2], "%)"),
       color = "Source") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "Assignment_all.png"),
       p_map, width = 14, height = 10, dpi = 300)


# ------------------------------------------------------------------------------
# 11. Assignment map (without P and Ca)
# ------------------------------------------------------------------------------

arrows2 <- merge(coords2, assign_noPCa, by.x = "Sample", by.y = "Artefact")
source_coords2 <- coords2[coords2$Sample %in% unique(arrows2$Source),
                          c("PC1", "PC2", "Sample")]
names(source_coords2)[1:2] <- c("PC1_end", "PC2_end")
arrows2 <- merge(arrows2, source_coords2, by.x = "Source", by.y = "Sample")

p_map2 <- ggplot(coords2, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(data = subset(coords2, Type == "Source"),
             aes(color = Sample), size = 5, shape = 17) +
  geom_point(data = subset(coords2, Type == "Artefact"),
             color = "gray40", size = 3.5, shape = 16) +
  geom_segment(data = arrows2,
               aes(x = PC1, y = PC2, xend = PC1_end, yend = PC2_end,
                   color = Source),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1, alpha = 0.8) +
  geom_text_repel(data = subset(coords2, Type == "Artefact"),
                  aes(label = Sample), size = 3.2,
                  box.padding = 0.5, segment.color = "grey50", max.overlaps = 20) +
  geom_text_repel(data = subset(coords2, Type == "Source"),
                  aes(label = Sample, color = Sample), size = 4.5,
                  fontface = "bold", box.padding = 0.8, segment.color = NA) +
  scale_color_manual(values = source_colors) +
  labs(title = "Assignment of artefacts to sources (without P and Ca)",
       x = paste0("PC1 (", var_noPCa[1], "%)"),
       y = paste0("PC2 (", var_noPCa[2], "%)"),
       color = "Source") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "Assignment_noPCa.png"),
       p_map2, width = 14, height = 10, dpi = 300)


# ------------------------------------------------------------------------------
# 12. Scree plot
# ------------------------------------------------------------------------------

scree <- data.frame(
  Component = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  Variance = summary(pca_all)$importance[2, 1:10] * 100
)

p_scree <- ggplot(scree, aes(Component, Variance, group = 1)) +
  geom_line(color = "#2E86AB", linewidth = 1) +
  geom_point(color = "#2E86AB", size = 3) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")),
            vjust = -1, size = 3.5) +
  labs(title = "Variance explained by component",
       x = "Component", y = "Variance explained (%)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "Scree_plot.png"), p_scree, width = 10, height = 6, dpi = 300)


# ------------------------------------------------------------------------------
# 13. Heatmap
# ------------------------------------------------------------------------------

pheatmap(as.matrix(clr_data),
         annotation_row = data.frame(Type = df$Type,
                                     row.names = df$Nomenclature),
         main = "Heatmap of CLR-transformed data",
         fontsize_row = 10,
         filename = file.path(output_dir, "Heatmap_CLR.png"),
         width = 10, height = 8)


# ------------------------------------------------------------------------------
# 14. Summary table
# ------------------------------------------------------------------------------

summary_table <- data.frame(
  Sample = df$Nomenclature,
  Type = df$Type,
  SiO2 = df$SiO2, Al2O3 = df$Al2O3, Fe2O3 = df$Fe2O3,
  K2O = df$K2O, CaO = df$CaO, P2O5 = df$P2O5
)
write.csv(summary_table, file.path(output_dir, "Summary_table.csv"), row.names = FALSE)

cat("\nDone. All files are in:", output_dir, "\n")