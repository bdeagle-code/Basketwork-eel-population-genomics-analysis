### Diastobranchus October 2026 DArT code for CSIRO DAP.R
#   The code runs with files available at https://data.csiro.au/collection/csiro:78574 
#   in the folder "R code and related files" 
#
#   This version of the code was used in analysis for the initial submission of a paper.
#   It has been put on GitHub to provide any updated code 
#
#############################################
#   D.capensis SNP data set analysis
Folder= "C:/Documents/00 People Papers Projects/Diastobranchus/RAD-seq 2026"
setwd (Folder);  dir()

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("SNPRelate")

# dartRverse_install("dartR.popgen", repository = "github", branch = "dev")
library(adegenet)
library(vcfR)
library(SNPRelate)
library(HardyWeinberg)
library (dartR.base)
library (StAMPP)
library(dartR.popgen)
library(vegan)
library(geosphere)
library(hierfstat)
library(rnaturalearth)
library (sf)
library(ggplot2)
library(maps)

#############################################
#############################################
###   Read in the SNP genotype data
vcf <- read.vcfR("NXGSQCAGRF25040091-3_variants.vcf", verbose = FALSE )
head(vcf)
gt <- extract.gt(vcf, element = "GT")
colnames(gt)

# Write out the file names and match with other metadata
# write.csv(colnames(gt), file="names2026.csv")
Dc_GenL = vcfR2genlight(vcf)
Dc_GenL

# Sample names
Dc_GenL@ind.names

##########################################################################################
##########################################################################################
###   Read in sample metadata
RAD_info1= read.table("EelGenoInfo2026.csv",sep=",",header=T, comment.char = "")
dim(RAD_info1)
head(RAD_info1)
table(RAD_info1$Location)

# Adjust names in genlight object
indNames(Dc_GenL) <- RAD_info1$Short_name
pop(Dc_GenL) <- RAD_info1$Location

###########################################################
######################################################################################################################
#Code to concatenate sequence for "phylogenetic" analysis     (MEGA analysis)
# Extract genotypes as a matrix, convert to nucleotides (A or T) and subset data so there are 3 eels per population
Dc_GenL
Dc_GenL@gen
geno_mat <- as.matrix(Dc_GenL)

head(geno_mat[, 1:10])
dim(geno_mat)

allele_freq <- colMeans(geno_mat, na.rm = TRUE) / 2
head(allele_freq)
hist(allele_freq, breaks =50)

hist(geno_mat, main = "Genotype distribution", xlab = "Genotype (0/1/2)")

indNames(Dc_GenL)
pop(Dc_GenL)

GenoMat_ACGT = t(geno_mat)
colnames (GenoMat_ACGT) = indNames(Dc_GenL)
head(GenoMat_ACGT[, 1:30])

# Just a quick look at the replicate samples in the data
colnames(GenoMat_ACGT[,c(23,24,25,26)])
Indexx = apply(GenoMat_ACGT[,c(23,24,25,26)],1,sum)>0
GenoMat_ACGT[Indexx,c(23,24,25,26)]

GenoMat_ACGT=replace(GenoMat_ACGT, GenoMat_ACGT== 0 ,"AA")
GenoMat_ACGT=replace(GenoMat_ACGT, GenoMat_ACGT== 1 ,"AT")
GenoMat_ACGT=replace(GenoMat_ACGT, GenoMat_ACGT== 2 ,"TT")
GenoMat_ACGT[is.na(GenoMat_ACGT)] = "NN"
GenoMat_ACGT; colnames (GenoMat_ACGT)
# Pick three from each population
ThreeIdx =  which(colnames (GenoMat_ACGT) %in% c("Pat_01", "GAB_01", "GAB_02", "IO_01", "IO_02", "IO_03",
                                                 "Pat_02", "Pat_03", "NE_Tas", "GAB_03", "SE_Tas",
                                                 "Mac_01", "Mac_02", "Vic_01", "Vic_02", "Vic_03",
                                                 "NZS_01", "NZS_02", "NZS_03", "WI_01", "WI_02",
                                                 "WI_03", "Kk_01", "Kk_02", "Kk_03", "SI_01",
                                                 "SI_02", "SI_03", "N_Mac_01", "N_Mac_02"))

ConCatseq <- function(x, outfile = "SNP_sequences.fasta") {
  # concatenate SNPs for each individual
  seqs <- apply(x, 2, paste0, collapse = "")
  # create FASTA format
  fasta <- paste0(">", names(seqs), "\n", seqs)
  # write FASTA file
  writeLines(fasta, outfile)
  # return sequences invisibly
  invisible(fasta)
}
# run function
ConCatseq(GenoMat_ACGT[,ThreeIdx])

#############################################
# Create Maps - look at where the samples are from
# get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(world) +
  geom_sf(fill = "grey95") +
  geom_point(data = RAD_info1,
             aes(x = Long,
                 y = Lat,
                 color = Site_col),
             size = 3) +
  scale_color_identity() +
  coord_sf(xlim = c(50, 180), ylim = c(-57, -30)) +
  theme_minimal()

# Map Aus/NZ
ggplot(world) +
  geom_sf(fill = "grey95") +
  geom_point(data = RAD_info1,
             aes(x = Long,
                 y = Lat,
                 color = Site_col),
             size = 4) +
  scale_color_identity() +
  coord_sf(xlim = c(130, 180), ylim = c(-57, -30)) +
  theme_minimal()

ggplot(world) +
  geom_sf(fill = "grey95") +
  geom_point(data = RAD_info1,
             aes(x = Long,
                 y = Lat,
                 fill = Site_col),
             shape = 21,
             color = "black",
             stroke = 0.6,
             size = 4) +
  scale_fill_identity() +
  coord_sf(xlim = c(130, 180), ylim = c(-57, -30)) +
  theme_minimal()

#############################################
#### Back to SNP data in the genlight object

Dc_GenL <- gl.compliance.check(Dc_GenL)
gl.set.verbosity(3) #how much information dartR prints to the console
nLoc(Dc_GenL)
locNames(Dc_GenL)
nInd(Dc_GenL)
popNames(Dc_GenL)
pop(Dc_GenL)
table(pop(Dc_GenL))

# ??gl.report
gl.report.allna(Dc_GenL)
gl.report.callrate(Dc_GenL)
# gl.report.rdepth(Dc_GenL)  # Fatal Error: Read depth not included among the locus metrics

# Filter based on call rate for loci 
Dc_GenL <- gl.filter.callrate(Dc_GenL, method = "loc", threshold = 0.95)
#Summary of filtered dataset
# Call Rate for loci > 0.95 
# Original No. of loci : 195211 
# No. of loci retained: 84960 

# Filter based on call rate for indiv
# Dc_GenL <- gl.filter.callrate(Dc_GenL, method = "ind", threshold = 0.95)
# No. of individuals retained: 93 

# DON'T Filter based on MAF - some population have low numbers
# Dc_GenL <- gl.filter.maf(Dc_GenL, threshold = 0.02)
# Summary of filtered dataset
# MAF for loci > 0.02 
# Initial number of loci: 84960 
# Number of loci deleted: 34038 
# Final number of loci: 50922 

Dc_GenL <- gl.filter.monomorphs(Dc_GenL) # no monomorphs

# Dc_GenL <- gl.filter.ld(Dc_GenL, threshold = 0.2) # No measure of LD to filter on

Dc_GenL <- gl.filter.hwe(Dc_GenL, alpha = 0.001)
# population Vic has less than 5 individuals... skipped
# Loci examined: 50922 
# Deleted 1063 loci with significant departure from HWE at alpha = 0.001 applied locus by locus
# Loci retained: 49859 
# Adjustment of p-values for multiple comparisons vary with sample size



# Could remove snps in LD - but this only picks up a small number 
# convert genlight → GDS
# snpgdsCreateGeno("temp.gds",
#                 genmat = as.matrix(Dc_GenL),
#                 sample.id = indNames(Dc_GenL),
#                 snp.id = locNames(Dc_GenL),
#                 snpfirstdim = FALSE)

#genofile <- snpgdsOpen("temp.gds")

# LD pruning
#snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2)

# extract pruned SNPs
#snpset.id <- unlist(snpset)

# subset original genlight
#Dc_GenL_pruned <- Dc_GenL[, snpset.id]

###########################################
#Remove two replicate individuals
indNames(Dc_GenL) 
gl_rm <- gl.drop.ind(Dc_GenL,c("GAB_05b","GAB_06b"))
idx = which(RAD_info1$Short_name %in% c("GAB_05b", "GAB_06b"))

###########################################
#  Read in colours for plotting
# Okabe–Ito palette
col_rm_new = RAD_info1$Site_col[-idx]
 c("black",
"#0072B2", # blue         Patience
"#56B4E9", # sky blue     GAB
"#E69F00", # orange       NZ south
"#F0E442", # yellow       White Island
"#D55E00", # vermillion   Kaikoura
"#009E73"  # bluish green Stewart Island
)

###########################################
# Run PCA  ****Select the number of axes******
pca = glPca(gl_rm, nf=6)

plot(pca$scores[,1], pca$scores[,2], xlab = "PC1", ylab = "PC2", pch = 19, col= col_rm_new)
identify(pca$scores[,1], pca$scores[,2], label = indNames(gl_rm) )

#gg plot
pca_12=pca
df <- data.frame(
  PC1 = pca_12$scores[,1],
  PC2 = pca_12$scores[,2],
  pop = pop(gl_rm)
)

df$pop <- pop(gl_rm)
var <- 100 * pca_12$eig / sum(pca_12$eig)

ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(colour = col_rm_new), size = 3) +
  scale_colour_identity() +
  
  stat_ellipse(
    aes(group = pop, colour = col_rm_new),
    type = "norm",
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  xlab(paste0("PC1 (", round(var[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var[2], 1), "%)")) +
  theme_minimal()



# Look at pairwise plots of PCs
pairs(pca_12$scores[,1:4], pch = 19)

# Look at loading
plot(abs(pca_12$loadings[,1]),
     ylab = "PC1 loading",
     xlab = "SNP index",
     pch = 16)

# Could remove snps in LD - but this only picks up a small number 
# convert genlight → GDS
# snpgdsCreateGeno("temp.gds",
#                 genmat = as.matrix(Dc_GenL),
#                 sample.id = indNames(Dc_GenL),
#                 snp.id = locNames(Dc_GenL),
#                 snpfirstdim = FALSE)

#genofile <- snpgdsOpen("temp.gds")

# LD pruning
#snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2)

# extract pruned SNPs
#snpset.id <- unlist(snpset)

# subset original genlight
#Dc_GenL_pruned <- Dc_GenL[, snpset.id]


#############################################
#############################################
# Try with just bigger pops and equal numbers 
gl_rm  # dataset with the replicates removed 

######
######  12 fish from 4 Pops
gl_12 <- gl.keep.ind(Dc_GenL,indNames(Dc_GenL)[RAD_info1$Pop_12==1 ])
gl_12 <- gl.filter.monomorphs(gl_12 ) # Monomorphic loci: 9562 
#Sanity check
cbind(RAD_info1$Location[RAD_info1$Pop_12==1 ] ,indNames(Dc_GenL)[RAD_info1$Pop_12==1 ])


######  8 fish from 6 Pops
gl_8 = gl.keep.ind(Dc_GenL,indNames(Dc_GenL)[RAD_info1$Pop_8==1 ])
gl_8 <- gl.filter.monomorphs(gl_8) 
#Sanity check
cbind(RAD_info1$Location[RAD_info1$Pop_8==1 ] ,indNames(Dc_GenL)[RAD_info1$Pop_8==1 ])

######  8 fish from 6 Pops - plus Macca X2 and Indian Ocean
#gl_8p = gl.keep.ind(Dc_GenL,indNames(Dc_GenL)[RAD_info1$Pop_8p==1 ])
#col_8p = RAD_info1$Site_col[RAD_info1$Pop_8p==1 ]


CURRENT=  gl_8   #gl_rm 

CURRENT@pop

# Do PCA on "CURRENT" 
pca <- glPca(CURRENT, nf = 10)

var <- pca$eig / sum(pca$eig) * 100

df <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  pop = as.factor(pop(CURRENT))
)

####### Do PCA separate plot for each of the data sets to make simple
### SNP_12_eels
df$pop <- pop(CURRENT)

pop_cols <- c(
  "GAB" = "#56B4E9",
  "NZ_south" = "#E69F00",
  "Patience" = "#0072B2",
  "Stewart_Island" = "#009E73"
)

p <- ggplot(df, aes(PC1, PC2)) +
  
  geom_point(aes(colour = pop), size = 3, alpha = 0.8) +
  
  stat_ellipse(
    aes(group = pop, colour = pop),
    type = "norm",
    linewidth = 1
  ) +
  
  scale_colour_manual(values = pop_cols) +
  
  xlab(paste0("PC1 (", round(var[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(var[2], 2), "%)")) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_line(colour = "grey90"),
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )
p

#Add sample names
p + geom_text(aes(label = rownames(df)), vjust = -0.5, size = 3)

###########################################
### SNP_8_eels
df$pop <- pop(CURRENT)

pop_cols <- c(
  "GAB" = "#56B4E9",
  "NZ_south" = "#E69F00",
  "Patience" = "#0072B2",
  "Stewart_Island" = "#009E73",
  "White_Island" = "#F0E442",
  "Kaikoura" = "#D55E00"
)

pop_shapes <- c(
  "GAB" = 16,
  "NZ_south" = 17,
  "Patience" = 15,
  "Stewart_Island" = 18,
  "White_Island" = 19,
  "Kaikoura" = 3
)

p <- ggplot(df, aes(PC1, PC2)) +
  
  geom_point(
    aes(colour = pop, shape = pop),
    size = 3,
    alpha = 0.8
  ) +
  
  stat_ellipse(
    aes(group = pop, colour = pop),
    type = "norm",
    linewidth = 1
  ) +
  
  scale_colour_manual(values = pop_cols) +
  scale_shape_manual(values = pop_shapes) +
  
  xlab(paste0("PC1 (", round(var[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(var[2], 2), "%)")) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

p

#Add sample names
p + geom_text(aes(label = rownames(df)), vjust = -0.5, size = 3)


###########################################
### All_eels
df$pop <- pop(CURRENT)
levels(as.factor(df$pop))

pop_cols <- c(
  "GAB" = "#56B4E9",
  "NZ_south" = "#E69F00",
  "Patience" = "#0072B2",
  "Stewart_Island" = "#009E73",
  "White_Island" = "#F0E442",
  "Kaikoura" = "#D55E00",
  "SE_Tas"  ="black",
  "Macca" ="black",
  "N_Macca" ="black",
  "NE_Tas" ="black",
  "Vic" ="black",
  "IndianO" ="black"
)

p <- ggplot(df, aes(PC1, PC2)) +
  
  geom_point(aes(colour = pop), size = 3, alpha = 0.8) +
  
  stat_ellipse(
    aes(group = pop, colour = pop),
    type = "norm",
    linewidth = 1
  ) +
  
  scale_colour_manual(values = pop_cols) +
  
  xlab(paste0("PC1 (", round(var[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(var[2], 2), "%)")) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_line(colour = "grey90"),
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )
p

#Add sample names
p + geom_text(aes(label = rownames(df)), vjust = -0.5, size = 3)

###########################################
# Scree plot
var <- pca$eig / sum(pca$eig) * 100

df_var <- data.frame(
  PC = 1:length(var),
  variance = var
)

ggplot(df_var[1:6,], aes(x = PC, y = variance)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("% variance explained") +
  xlab("Principal Component")

###########################################
###########################################
gl_rm  # data set with the replicates removed 

nInd(gl_rm)
nLoc(gl_rm)
popNames(gl_rm)

######  12 fish from 4 Pops
gl_12 

######  8 fish from 6 Pops
gl_8 

###########################################
##########################
# AMOVA
gl.amova(gl_rm)
gl.amova(gl_12)
gl.amova(gl_8)


#  FST (dartR)
#  gl.fst.pop(CURRENT, method = "weir.cockerham")

class(gl_rm) #[1] "dartR"
fst_perm <- gl.fst.pop(gl_rm, nboots = 1000)
fst_perm

fst_perm12 <- gl.fst.pop(gl_12, nboots = 1000)
fst_perm12

fst_perm8 <- gl.fst.pop(gl_8, nboots = 1000)
fst_perm8



#############################################
##############
# Check IBD
class(gl_8)
nInd(gl_8)
nLoc(gl_8)
popNames(gl_8)


#############################################
######  8 fish from 6 Pops - plus Patience 2023 eels
# gl_8_P_2023 = gl.keep.ind(Dc_GenL,indNames(Dc_GenL)[RAD_info1$Pop_8_Patience2==1 ])


# fst_mat <- gl.fst.pop(gl_8_P_2023)

fst_mat <- gl.fst.pop(gl_8)
gen_dist <- as.dist(fst_mat)
rownames(fst_mat)

coords <- data.frame(
  pop = c("Patience", "GAB", "NZ_south", "White_Island", "Kaikoura", "Stewart_Island"),
  lon = c(147.3773333, 131.7821667, 170.0315, 177.255305, 173.853583, 167.0433 ),
  lat = c(-44.12066667, -34.37766667, -48.4559, -37.50715, -42.575778, -46.5314)
 )

# Includes Patience twice
#coords <- data.frame(
#  pop = c("Patience", "GAB", "Patience2023", "NZ_south", "White_Island", "Kaikoura", "Stewart_Island"),
#  lon = c(147.3773333, 131.7821667, 147.3773333, 170.0315, 177.255305, 173.853583, 167.0433 ),
#  lat = c(-44.12066667, -34.37766667, -44.12066667, -48.4559, -37.50715, -42.575778, -46.5314)
# )

coords <- coords[match(rownames(fst_mat), coords$pop), ]

geo_dist <- distm(coords[, c("lon","lat")], fun = distHaversine)

# Convert to km and distance object
geo_dist <- as.dist(geo_dist / 1000)
str(geo_dist)

time_dist = geo_dist
time_dist[1:15] = c(0.1,9,6,5,9,9,6,5,9,15,14,0.1,1,15,14)
# time_dist[1:21] = c(0.1,7,9,6,5,9,7,9,6,5,9,1,14,13,1,15,14,0.1,1,15,14)



mantel_res <- mantel(gen_dist, time_dist, method = "pearson", permutations = 9999)

print(mantel_res)

plot(time_dist, gen_dist,
     xlab = "Geographic distance (km)",
     ylab = "Genetic distance (FST)",
     pch = 19)

abline(lm(as.numeric(gen_dist) ~ as.numeric(time_dist)), col = "blue", lwd = 2)


df_ibd <- data.frame(
  geo_dist = geo_dist,
  gen_dist = gen_dist
)

####
### A nicer Plot of the MAntel test result 
p <- ggplot(df_ibd, aes(x = time_dist, y = gen_dist)) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
 
   annotate(
    "text",
    x = 15,
    y = 0,
    label = "Mantel r = 0.483\np = 0.056",
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  labs(
    x = "Difference in collection date (years)",
    y = expression(Genetic~distance~(F[ST]))
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.title = element_text()
  )
p


#############################################
#############################################
# DAPC

grp <- pop(CURRENT)

dapc1 <- dapc(gl_12, grp)

scatter(dapc1)
