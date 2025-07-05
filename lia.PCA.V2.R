
### Running PCA
Comparing ADMIXTURE and PCA results often helps give insight and confirmation regarding population structure in a sample. To run PCA, a standard package that is well-suited for SNP data is the `smartpca` package maintained by Nick Patterson and Alkes Price (at http://data.broadinstitute.org/alkesgroup/EIGENSOFT/). To run it, we first set-up a basic `smartpca` parameter file from the command-line of a `bash` shell:

```{r, engine = 'bash', eval = FALSE}
PREFIX=H938_Euro.LDprune
echo genotypename: out/$PREFIX.bed > out/$PREFIX.par
echo snpname: out/H938_Euro.LDprune.bim >> out/$PREFIX.par
echo indivname: out/H938_Euro.LDprune.PCA.fam >> out/$PREFIX.par
echo snpweightoutname: out/H938_Euro.LDprune.snpeigs \
>> out/$PREFIX.par
echo evecoutname: out/H938_Euro.LDprune.eigs >> out/$PREFIX.par
echo evaloutname: out/H938_Euro.LDprune.eval >> out/$PREFIX.par
echo phylipoutname: out/H938_Euro.LDprune.fst >> out/$PREFIX.par
echo numoutevec: 20 >> out/$PREFIX.par
echo numoutlieriter: 0 >> out/$PREFIX.par
echo outlieroutname: out/H938_Euro.LDprune.out >> out/$PREFIX.par
echo altnormstyle: NO >> out/$PREFIX.par
echo missingmode: NO >> out/$PREFIX.par
echo nsnpldregress: 0 >> out/$PREFIX.par
echo noxdata: YES >> out/$PREFIX.par
echo nomalexhet: YES >> out/$PREFIX.par
```   
This input parameter file runs `smartpca` in its most basic mode (i.e. no automatic outlier removal or adjustments for LD - features which you might want to explore later). 

As a minor issue, `smartpca` ignores individuals in the `.fam` file if they are marked as missing in the phenotypes column. This `awk` command provides a new `.fam` file that will automatically include all individuals.
    
```{r, engine = 'bash', eval = FALSE}
awk '{print $1,$2,$3,$4,$5,1}' out/H938_Euro.LDprune.fam \
> out/H938_Euro.LDprune.PCA.fam
```
Now run `smartpca` with the following command. 
    
```{r, engine = 'bash', eval = FALSE}
smartpca -p ./out/H938_Euro.LDprune.par
```   
You will find the output files in the `out` sub-directory as specified in the parameter file.

### Plotting PCA results with PCAviz
The PCAviz package can be found at https://github.com/NovembreLab/PCAviz. It provides a simple interface for quickly creating plots from PCA results. It encodes several of our favored best practices for plotting PCA (such as using abbreviations for point characters and plotting median positions of each labelled group).  To install the package use:
```{r, eval=FALSE,echo=TRUE}
install.packages("devtools")
devtools::install_github("NovembreLab/PCAviz",
                         build_vignettes = TRUE)
```

The following command in `R` generates plots showing each individual sample's position in the PCA space and the median position of each labelled group in PCA space:

```{r, eval = TRUE, echo=TRUE, results='hide',warning=FALSE,message=FALSE, fig.height=10, fig.width=10, out.width = "100%", fig.cap="Pairwise plots of PC scores generated using the PCAviz package."}
library(PCAviz)
library(cowplot)
prefix <- "Perognathus.124.lon.ino.amp.noFlavus.AUTOSOMES.1pc.pruned"
nPCs <- 20

# Read in individual coordinates on PCs and eignvalues
PCA <- read.table(paste(prefix, ".eigs", sep = ""))
names(PCA) <- c("ID", paste("PC", (1:nPCs), sep = ""), 
                  "case.control")
PCA <- PCA[, 1:(nPCs + 1)] # Remove case/control column
eig.val <- sqrt(unlist(read.table(
  paste(prefix, ".eval", sep = "")))[1:nPCs])
sum.eig <- sum(unlist(read.table(
  paste(prefix, ".eval", sep = ""))))

# Read in snp weightings matrix
snpeigs <- read.table(paste(prefix, ".snpeigs", sep = ""))
names(snpeigs) <- c("ID", "chr", "pos", 
                    paste("PC", (1:nPCs), sep = ""))
snpeigs$chr <- factor(snpeigs$chr)
rownames(snpeigs) <- snpeigs$ID
snpeigs <- snpeigs[, -1]

# Note smartpca pushes the plink family and individual
# ids together so we need to extract out the ids afresh
tmp <- unlist(sapply(as.character(PCA$ID), strsplit, ":"))
ids <- tmp[seq(2, length(tmp), by = 2)]
PCA$ID <- ids

# Read in the group/cluster labels
clst <- read.table("lia.2.clst.txt")
# Order them to match the ids of PCA object
clst_unord <- clst$V2[match(ids, clst$V1)] 
PCA <- as.data.frame(PCA)
PCA <- cbind(PCA, clst_unord)
names(PCA)[ncol(PCA)] <- "sample"

################
### OPTION 1: PCAviz object if we desire to plot with sample names
# Build the PCAviz object
PCA$ID <- NULL
rownames(PCA) <- NULL

lia <- pcaviz(dat = PCA, sdev = eig.val, var = sum.eig)
lia <- pcaviz_abbreviate_var(lia, "sample")

# Make PCA plots
geom.point.summary.params <- list(
  shape = 16, stroke = 1, size = 5,
  alpha = .7, show.legend = TRUE  # Show the legend!
)

# An example plot block
plot1 <- plot(
	lia,
	coords = paste0("PC", c(1, 2)),
	color = "sample",
	geom.point.summary.params = geom.point.summary.params,
	summary.fun = "mean",
	geom.text.summary.params = list(size = 4, fontface = "bold"),
	scale.pc.axes = 0.6
)
################

### OPTION 2: Simple plotting with ggplot2 to ensure control over graphic parameters

plot1 <-ggplot(PCA, aes_string(x = "PC1", y = "PC2", color = "sample")) +
  geom_point(size = 4, alpha = 0.7) +
  theme_classic() +
  labs(color = "Population")


# Make PCA plots
geom.point.summary.params <- list(
  shape = 16, stroke = 1, size = 5,
  alpha = .7, show.legend = TRUE  # Show the legend!
)

################
### OPTION 3: Simple plotting with ggplot2 to ensure control over graphic parameters,
# extended to show centroid labels over populations for easier interpretation


##########
ITALICISE NAMES:
# ---- 1. Prepare legend and centroid labels ----

# Function to build plotmath label for centroid
make_plotmath_label <- function(x) {
  parts <- unlist(strsplit(x, "_"))
  if (length(parts) > 1) {
    first <- paste0("italic(", parts[1], ")")
    rest <- paste(parts[-1], collapse=" ")
    return(paste0(first, '~"', rest, '"'))
  } else {
    return(paste0("italic(", parts[1], ")"))
  }
}

# Add legend label (plain text, spaces) and centroid label (plotmath)
PCA <- PCA %>%
  mutate(
    pop_legend = gsub("_", " ", sample),
    pop_label = sapply(sample, make_plotmath_label)
  )

# ---- 2. Calculate centroids ----
centroids <- PCA %>%
  group_by(sample, pop_legend, pop_label) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 3. Plot ----
p <- ggplot(PCA, aes(x = PC1, y = PC2, color = pop_legend)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text(
    data = centroids,
    aes(label = pop_label),
    color = "black", fontface = "bold", size = 5, parse = TRUE
  ) +
  theme_classic() +
  labs(color = "Population") +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size=4)))

print(p)

#################
ADD THE COLOUR SCHEME
pop_colors <- c("inornatus E" = "#DA3424",
"inornatus W" = "#F3A740",
"amplus NW AZ" = "#00A0DF",
"amplus ColPlat" = "#EA2E13",
"amplus S AZ" = "#AE95BF",
"amplus Salt Gila" = "#8F438F",
"longimembris UpColR" = "#DA3424",
"longimembris AZ hap5" = "#F4D690",
"longimembris AZ hap6" = "#FAAA0A",
"longimembris VirRiv" = "#DC9512",
"longimembris GrBas" = "#3CA0DB",
"longimembris Mojave" = "#1A6AD2",
"longimembris SW CA" = "#84CAD9",
"longimembris SE CA" = "#C1E204",
"longimembris SW AZ" = "#FCEB09")

p <- ggplot(PCA, aes(x = PC1, y = PC2, color = pop_legend)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text(
    data = centroids,
    aes(label = pop_label),
    color = "black", fontface = "bold", size = 5, parse = TRUE
  ) +
  theme_classic() +
  labs(color = "Population") +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values = pop_colors)

print(p) + scale_color_manual(values = pop_colors)

####################
USE SHAPES

pop_shapes <- c(
  "inornatus E" = 15,           # Filled square
  "inornatus W" = 16,           # Filled circle
  "amplus NW AZ" = 17,          # Filled triangle up
  "amplus ColPlat" = 18,        # Filled diamond
  "amplus S AZ" = 19,           # Large filled circle
  "amplus Salt Gila" = 7,       # Filled hexagon
  "longimembris UpColR" = 8,    # Filled star
  "longimembris AZ hap5" = 20,  # Filled small circle
  "longimembris AZ hap6" = 21,  # Filled circle w/ border
  "longimembris VirRiv" = 22,   # Filled square w/ border
  "longimembris GrBas" = 23,    # Filled diamond w/ border
  "longimembris Mojave" = 24,   # Filled triangle up w/ border
  "longimembris SW CA" = 25,    # Filled triangle down w/ border
  "longimembris SE CA" = 6,     # Filled upside hexagon
  "longimembris SW AZ" = 14,     # Filled asterisk/cross
  "extra pop 1" = 1,   # Filled circle (default)
  "extra pop 2" = 2,   # Filled triangle
  "extra pop 3" = 5,   # Filled diamond
  "extra pop 4" = 9,   # Filled circle with plus
  "extra pop 5" = 10,  # Filled triangle with plus
  "extra pop 6" = 11   # Filled diamond with plus
)


p <- ggplot(PCA, aes(x = PC1, y = PC2, color = pop_legend, shape = pop_legend, fill = pop_legend)) +
  geom_point(size = 6, alpha = 1) +
  geom_text(
    data = centroids,
    aes(label = pop_label),
    color = "black", fontface = "bold", size = 5, parse = TRUE
  ) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values = pop_colors) +
  scale_shape_manual(values = pop_shapes) +
  theme_classic() +
  labs(color = "Population", shape = "Population", fill = "Population") +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4))
  )

print(p)

##########
NO CENTROIDS

p <- ggplot(PCA, aes(x = PC1, y = PC2, color = pop_legend, shape = pop_legend, fill = pop_legend)) +
  geom_point(size = 6, alpha = 1) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values = pop_colors) +
  scale_shape_manual(values = pop_shapes) +
  theme_classic() +
  labs(color = "Population", shape = "Population", fill = "Population") +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4))
  )

print(p)









# No text labels per point   ADJUST PAST PLOT1
###############
plot1 <- plot(
	lia,
	coords = paste0("PC", c(1, 2)),
	color = "sample",
	geom.point.summary.params = geom.point.summary.params,
	summary.fun = "mean",
	geom.text.summary.params = list(size = 4, fontface = "bold"),
	scale.pc.axes = 0.6
)
################
 
 plot2 <- plot(
  lia,
  coords = paste0("PC", c(2, 3)),
  color = "sample",
  geom.point.summary.params = geom.point.summary.params,
  geom.text.summary.params = NULL,
  geom.label.summary.params = NULL,
  scale.pc.axes = 0.6
)
plot3 <- plot(
  lia,
  coords = paste0("PC", c(4, 5)),
  color = "sample",
  geom.point.summary.params = geom.point.summary.params,
  geom.text.summary.params = NULL,
  geom.label.summary.params = NULL,
  scale.pc.axes = 0.6
)
plot4 <- plot(
  lia,
  coords = paste0("PC", c(5, 6)),
  color = "sample",
  geom.point.summary.params = geom.point.summary.params,
  geom.text.summary.params = NULL,
  geom.label.summary.params = NULL,
  scale.pc.axes = 0.6
)

plot_grid(plot1, plot2, plot3, plot4)
```

First one may notice several populations are separated with PC1 and PC2, with the more isolated populations being those that were most distinguished from the others by `ADMIXTURE`. PC4 distinguishes a subset of Orcadian individuals and PC5 distinguishes two Adygei individuals.  PC6 corresponds to the cryptic structure observed within Sardinians in the `ADMIXTURE` analysis. 

As an alternative visualization, it can be helpful to see the distribution of PC coordinates per population for each labelled group in the data:
```{r, eval = TRUE,fig.height=8,fig.width=8, out.width = "100%",fig.cap="Violin plots of PC values generated using the PCAviz package."}
pcaviz_violin(hgdp, pc.dims = paste0("PC", c(1:3)), 
              plot.grid.params = list(nrow = 3))
```

As mentioned above in the section on LD, it is useful to inspect the PC loadings to ensure that they broadly represent variation across the genome, rather than one or a small number of genomic regions [@Duforet-Frebourg16]. SNPs that are selected in the same direction as genome-wide structure can show high loadings, but what is particularly pathological is if the only SNPs that show high loadings are all concentrated in a single region of the genome, as might occur if the PCA is explaining genomic structure (such as an inversion) rather than population structure.

```{r, eval = TRUE, fig.height=6, fig.width=8, out.width = "100%", fig.cap="PC loading plots generated using the PCAviz package."}
for (i in 1:5) {
  plotname <- paste("plot", i, sep = "")
  plot <- pcaviz_loadingsplot(hgdp,
    pc.dim = paste0("PC", i),
    min.rank = 0.8, gap = 200, color = "chr",
    geom.point.params = list(show.legend = FALSE)
  ) +
    xlab("SNPs") + ylab(paste0("PC", i, " loading"))
  assign(plotname, plot)
}
# grep common legend
plot <- pcaviz_loadingsplot(hgdp,
  pc.dim = paste0("PC", 1),
  min.rank = 0.8, gap = 200, color = "chr"
) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", 
        legend.justification = "center")
plot_legend <- get_legend(plot)
# plot loadings
prow <- plot_grid(plot1, plot2, plot3, plot4, plot5,
                  nrow = 5, align = "vh")
plot_grid(prow, plot_legend, ncol = 1, rel_heights = c(1, .2))
```

The proportion of total variance explained by each PC is a useful metric for understanding structure in a sample and for evaluating how many PCs one might want to include in downstream analyses. This can be computed as  $\lambda_i / \sum_{k} \lambda_k$, with $\lambda_i$ being eigenvalues in decreasing order, and is plotted below:
```{r, eval = TRUE, echo=TRUE,warnings=FALSE,fig.width=3, fig.height=3,fig.cap="Proportion of variance explained by each PC.  Plot generated using the PCAviz package."}
screeplot(hgdp, type = "pve") + 
  ylim(0, 0.018) +
  ylab('Proportion of Variance Explained') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.line = element_line(size = 1, linetype = "solid"))
```

The results show that the top PCs only explain a small fraction of the variance (<1.5%) and that after about $K=6$ the variance explained per PC becomes relatively constant; roughly in line with the visual inspection of the `admixture` results that revealed $K=6$ may be reasonable for this dataset. 


#  Discussion

Our protocol above is relatively straightforward and presents the most basic implementation of these analyses. Each analysis sofware (`ADMIXTURE` and `smartpca`) and each visualization package (`pong` and `PCAviz`) contain numerous other options that may be suitable for specific analyses and we encourage the readers to spend time in the manuals of each.  Nonetheless, what we have presented is a useful start and a standard pipeline that we use in our research.

Two broad perspectives we find helpful useful to keep in mind are: 1) How the admixture model and PCA framework are related to each other indirectly as different forms of sparse factor analysis [@Engelhardt10]; 2) How the PCA framework in particular can be considered as a form of efficient data compression.  Both of these perspectives can be helpful in interpreting the outputs of the methods and for appreciating how these approaches best serve as helpful visual exploratory tools for analyzing structure in genetic data.  These methods are ultimately relatively simple statistical tools being used to summarize complex realities.  They are part of the toolkit for analysis, and often are extremely useful for framing specific models of population structure that can be further investigated using more detailed and explicit approaches (such as those based on coalescent or diffusion theory, Chapters 7 on MSMC and 8 on CoalHMM). 

\newpage
# References


https://github.com/NovembreLab/PCAviz/blob/master/R/pcaviz.R