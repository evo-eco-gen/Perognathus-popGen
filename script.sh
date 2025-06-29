### PEROGNATHUS PCA
### Running a PCA analysis of Perognathus silky mice data using EIGENSOFT.

### a pure x86_64 env under Rosetta to avoid compatibility issues with MAc ARM2 architecture
CONDA_SUBDIR=osx-64 mamba create -n perognathus   -c conda-forge -c bioconda   --strict-channel-priority   plink openblas
conda activate perognathus
file $(which plink)

# Install gsl and smartpca in lockstep
mamba install -c conda-forge -c bioconda eigensoft gsl --force-reinstall

smartpca version: 18140, 8.0.0
plink  v1.90b6.18

# 1. Convert the VCF into bed/bim format.
plink \
  --vcf Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.vcf.gz \
  --biallelic-only strict \
  --geno 0.9742 \
  --indep-pairwise 50 5 0.1 \
  --make-bed \
  --allow-extra-chr \
  --out Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc


# 2. Assign unique numeric chromosome codes to every scaffold.
# get a list of all non‐numeric “chromosomes” in  .bim:
cut -f1 Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.bim | grep -vE '^[1-9][0-9]?$' | sort -u > scaffolds.txt

# assign them all to chr 20 (valid autosome code):
awk 'BEGIN{OFS="\t"} {print $1, 20}' scaffolds.txt > scaffold2chr.txt

# 3. Rewrite the BIM “chromosome” column.
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR { map[$1]=$2; next }
     {
       # if this variant’s “chrom” is in our scaffold list, change it to 20
       if($1 in map) $1 = map[$1]
       print
     }' scaffold2chr.txt Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.bim \
  > Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.bim

# copy over the matching .bed/.fam so PLINK still sees a consistent trio:
cp Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.bed Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.bed
cp Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.fam Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.fam

# 4. Ensure unique SNP IDs in column 2 (chr_pos).
awk 'BEGIN{FS=OFS="\t"} { $2 = $1 "_" $4; print }' \
  Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.bim \
> Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.fixid.bim

# swap in:
mv Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.fixid.bim Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.bim


# 5. Sort variants by (chromosome, position).
# 5a) Make transpose files from uni‐chrom dataset:
plink \
  --bfile Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr \
  --recode transpose --allow-extra-chr \
  --out Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr

# 5b) Sort the resulting .tped by chr (col 1) then bp (col 4):
sort -k1,1n -k4,4n Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.tped > Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.tped

# 5c) Copy the matching .tfam:
cp Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.tfam Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.tfam

# 5d) Rebuild a sorted .bed/.bim/.fam from the sorted tped/tfam:
plink \
  --tped Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.tped \
  --tfam Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.tfam \
  --make-bed --allow-extra-chr \
  --out Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted


# 6. Fix the .fam "phenotype" column 6 from -9 to 1
awk 'BEGIN{OFS="\t"} { $6 = 1; print }' \
  Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.fam \
> Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.fixfam.fam

mv Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.fixfam.fam Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted.fam 


# 7. LD-prune
plink --allow-extra-chr \
  --bfile Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted \
  --indep-pairwise 50 5 0.1 \
  --out Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.prune

# 650000 of 779247 variants removed.

plink \
  --bfile Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.uniChr.sorted \
  --extract Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.prune.in \
  --make-bed --allow-extra-chr \
  --out Perognathus.125.flavu.lon.ino.amp.AUTOSOMES.1pc.pruned


# Run Eigensoft:
smartpca -p ./out/flia.smartpca.txt

