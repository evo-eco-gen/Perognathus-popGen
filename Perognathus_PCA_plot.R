# ───────────────────────────────────────────────────────────────────────────────
# COMPLETE R SCRIPT: PCA SCATTER WITH CUSTOM COLORS, SHAPES, AND CENTROIDS (NO LABELS)
# ───────────────────────────────────────────────────────────────────────────────

# 0. PACKAGES
library(dplyr)
library(ggplot2)

# 1. PARAMETERS
prefix <- "Perognathus.124.lon.ino.amp.noFlavus.AUTOSOMES.1pc.pruned"
nPCs   <- 20

# 2. READ PCA COORDINATES & EIGENVALUES
PCA <- read.table(paste0(prefix, ".eigs"),
                  stringsAsFactors = FALSE)
colnames(PCA) <- c("ID",
                   paste0("PC", 1:nPCs),
                   "case.control")
# keep only ID + PCs
PCA <- PCA[, 1:(nPCs + 1)]

# clean up ID field (keep only the individual)
PCA$ID <- sapply(strsplit(PCA$ID, ":"), `[`, 2)

# 3. READ POPULATION ASSIGNMENTS
clst <- read.table("lia.2.clst.txt", stringsAsFactors = FALSE)
# clst: V1 = ID, V2 = population code (underscores)
PCA$sample <- clst$V2[match(PCA$ID, clst$V1)]

# drop the old ID column & rownames
PCA$ID <- NULL
rownames(PCA) <- NULL

# 4. CREATE LEGEND LABEL (spaces instead of underscores)
PCA <- PCA %>%
  mutate(pop_legend = gsub("_", " ", sample))

# 5. COMPUTE POPULATION CENTROIDS (PC1 vs PC2)
centroids <- PCA %>%
  group_by(pop_legend) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )

# 6. DEFINE YOUR CUSTOM COLOR & SHAPE PALETTES
pop_colors <- c(
  "inornatus E"         = "#DA3424",
  "inornatus W"         = "#F3A740",
  "amplus NW AZ"        = "#00A0DF",
  "amplus ColPlat"      = "#EA2E13",
  "amplus S AZ"         = "#AE95BF",
  "amplus Salt Gila"    = "#8F438F",
  "longimembris UpColR" = "#DA3424",
  "longimembris AZ hap5"= "#F4D690",
  "longimembris AZ hap6"= "#FAAA0A",
  "longimembris VirRiv" = "#DC9512",
  "longimembris GrBas"  = "#3CA0DB",
  "longimembris Mojave" = "#1A6AD2",
  "longimembris SW CA"  = "#84CAD9",
  "longimembris SE CA"  = "#C1E204",
  "longimembris SW AZ"  = "#FCEB09"
)

pop_shapes <- c(
  "inornatus E"         = 15,   # filled square
  "inornatus W"         = 16,   # filled circle
  "amplus NW AZ"        = 17,   # filled triangle up
  "amplus ColPlat"      = 18,   # filled diamond
  "amplus S AZ"         = 19,   # large filled circle
  "amplus Salt Gila"    = 20,   # small filled circle
  "longimembris UpColR" = 21,   # filled circle w/ border
  "longimembris AZ hap5"= 22,   # filled square w/ border
  "longimembris AZ hap6"= 23,   # filled diamond w/ border
  "longimembris VirRiv" = 24,   # filled triangle up w/ border
  "longimembris GrBas"  = 25,   # filled triangle down w/ border
  "longimembris Mojave" = 7,    # filled hexagon
  "longimembris SW CA"  = 8,    # filled star
  "longimembris SE CA"  = 6,    # filled upside hexagon
  "longimembris SW AZ"  = 14    # filled asterisk/cross
)

# 7. DRAW PCA WITH ggplot2
p <- ggplot(PCA,
            aes(x = PC1, y = PC2,
                color = pop_legend,
                shape = pop_legend)) +
  # samples as large, opaque points
  geom_point(size = 7, alpha = 1) +
  # centroids on top (same color+shape, bigger; no legend entry)
  geom_point(data = centroids,
             aes(x = PC1, y = PC2,
                 color = pop_legend,
                 shape = pop_legend),
             size = 9, stroke = 2,
             show.legend = FALSE) +
  # apply your custom palettes
  scale_color_manual(values = pop_colors) +
  scale_shape_manual(values = pop_shapes) +
  # classic theme + axis/legend formatting
  theme_classic() +
  labs(color = "Population",
       shape = "Population") +
  theme(
    axis.title     = element_text(size = 14),
    axis.text      = element_text(size = 12),
    legend.text    = element_text(size = 12),
    legend.key.size= unit(0.8, "lines")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  )

# 8. SAVE & PRINT
ggsave("PCA_plot.png", p, width = 7, height = 5, dpi = 300)
print(p)
