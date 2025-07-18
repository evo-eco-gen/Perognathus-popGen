#!/usr/bin/env Rscript
# Usage: Rscript parallel_bcf_diversity_stats.R pop_assignments.txt scaffold_list.txt bcf_dir reference.fasta.fai output_prefix

library(PopGenome)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: Rscript parallel_bcf_diversity_stats.R popmap.txt scaffolds.txt bcf_dir ref.fai output_prefix")

popmap <- read.table(args[1], col.names=c("sample","pop"), stringsAsFactors=FALSE)
pop_counts <- table(popmap$pop)
multi_pops <- names(pop_counts)[pop_counts > 1]
if (any(pop_counts == 1)) cat("Removing singleton populations:", paste(names(pop_counts)[pop_counts==1], collapse = ", "), "\n")
popmap <- subset(popmap, pop %in% multi_pops)
pop_list <- split(popmap$sample, popmap$pop)
if (length(pop_list) == 0) stop("ERROR: No multi-sample populations remain after filtering.")

chrs <- scan(args[2], what=character(), sep="\n")
fai <- read.table(args[4], header=FALSE, stringsAsFactors=FALSE)
names(fai) <- c("scaffold","length","rest")
scaf_lengths <- setNames(fai$length, fai$scaffold)

genome_list <- lapply(chrs, function(chr) {
  bcf <- file.path(args[3], paste0(chr, ".vcf.gz"))
  if (!chr %in% names(scaf_lengths)) {
    cat("Skipping", chr, "- not in .fai\n")
    return(NULL)
  }
  len <- scaf_lengths[chr]
  cat("Processing", chr, "length", len, "\n")
  g <- tryCatch({
    readVCF(bcf, numcols=5000, include.unknown=TRUE,
            tid=chr, frompos=1, topos=len, parallel=TRUE)
  }, error = function(e) {
    cat("readVCF error for", chr, ":", e$message, "\n")
    return(NULL)
  })
  if (is.null(g)) return(NULL)

  g <- tryCatch({
    g <- set.populations(g, pop_list, diploid=TRUE)
    # ONLY nucleotide diversity (pi): skip haplotype
    g <- diversity.stats(g, pi = TRUE, subsites = FALSE, keep.site.info = FALSE)
    g <- neutrality.stats(g)
    g <- F_ST.stats(g, mode = "nucleotide")
    g <- detail.stats(g)
    g
  }, error = function(e) {
    cat("Stats error for", chr, ":", e$message, "\n")
    return(NULL)
  })

  g
})
genome_list <- Filter(Negate(is.null), genome_list)
if (length(genome_list) == 0) stop("ERROR: No scaffolds processed successfully.")

genome_all <- concatenate.classes(genome_list)
genome_all <- diversity.stats(genome_all, pi = TRUE, subsites = FALSE, keep.site.info = FALSE)
genome_all <- neutrality.stats(genome_all)
genome_all <- F_ST.stats(genome_all)
genome_all <- detail.stats(genome_all)

res <- list(
  pi         = genome_all@nuc.diversity.within,
  theta      = genome_all@theta_Watterson,
  TajimaD    = genome_all@Tajima.D,
  FuLiD      = genome_all@Fu.Li.D,
  FuLiF      = genome_all@Fu.Li.F,
  FayWuH     = genome_all@Fay.Wu.H,
  S          = genome_all@n.segregating.sites,
  LD         = mean(genome_all@Kelly.Z_nS, na.rm=TRUE),
  Fst_matrix = get.F_ST(genome_all)
)

sink(paste0(args[5], "_stats.txt"))
print(res[-length(res)])
cat("\nLD:", res$LD, "\n\nFst matrix:\n")
print(res$Fst_matrix)
sink()
write.csv(res$Fst_matrix, paste0(args[5], "_Fst_matrix.csv"), quote = FALSE)

cat("Done. Results written with prefix:", args[5], "\n")