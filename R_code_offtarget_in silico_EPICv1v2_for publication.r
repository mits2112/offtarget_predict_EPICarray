# packages 
library(minfi)
library(ewastools) 
library(dplyr)
library(IlluminaHumanMethylationEPICv2manifest) # EPIC v2
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38) # EPIC v2
library(IlluminaHumanMethylationEPICmanifest) # EPIC v1
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38) # EPIC v1 
library(limma)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)

##########################################
### Processing of EPIC v1 and EPIC v2 Arrays ####
##########################################

#----EPICv2 -----
# Read IDAT files and create EPICv2 extended RGset
path <- "/Users/hjj/Desktop/EPICv2_Junjie/data/analysis"
targets <- read.metharray.sheet(path)
baseDir <- system.file(path, package = "minfiData")
RGset <- read.metharray.exp(targets = targets, verbose = TRUE, extended = TRUE)
RGset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#test
betas <- getBeta(RGset)
dim(betas)
betas[1:5, 1:5]


# QC
## Control metrics
baseDir <- "/Users/hjj/Desktop/EPICv2_Junjie/data/analysis"   
idat_dir <- file.path(baseDir, "IDAT")                        
sample_info <- read.metharray.sheet(baseDir)
sample_info$Basename <- file.path(
idat_dir,
paste0(sample_info$Slide, "_", sample_info$Array)
)
ewas_meth <- read_idats(sample_info$Basename)
ctrls <- control_metrics(ewas_meth)
sample_failure(ctrls) 
## methylated/unmethylated intensity check (minfi)
mset_raw <- preprocessRaw(RGset) 
qc <- getQC(mset_raw)
badSampleCutoff <- 10.5 
meds <- (qc$mMed + qc$uMed)/2
whichBad <- which((meds < badSampleCutoff)) 
controlqc_minfi <- plotQC(qc, badSampleCutoff = badSampleCutoff)
## Sex prediction 
check_sexV2 <- function (raw) 
{
    if (!all(c("M", "U", "manifest") %in% names(raw))) 
        stop("Invalid argument")
    with(raw, {
        chrX = manifest[chr == "chrX", index] ## manifest chr name updated
        chrY = manifest[chr == "chrY", index]## manifest chr name updated
        auto = manifest[!chr %in% c("chrX", "chrY"), index]
        chrX = M[chrX, , drop = FALSE] + U[chrX, , drop = FALSE]
        chrY = M[chrY, , drop = FALSE] + U[chrY, , drop = FALSE]
        auto = M[auto, , drop = FALSE] + U[auto, , drop = FALSE]
        chrX = colMeans(chrX, na.rm = TRUE)
        chrY = colMeans(chrY, na.rm = TRUE)
        auto = colMeans(auto, na.rm = TRUE)
        chrX = chrX/auto
        chrY = chrY/auto
        return(list(X = chrX, Y = chrY))
   })
}
predicted_sex <- check_sexV2(ewas_meth)
ewastools_predicted_sex <- predict_sex(predicted_sex$X,predicted_sex$Y)


# Preprocessing

## create noob object (minfi)
noob <- preprocessNoob(RGset)

## noob and functional normalization (minfi)
nPC=2 #default or choose the top PCs that accounted for 90% or more of the variance.
MSet.fnorm <- preprocessFunnorm(RGset, nPCs = nPC, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)

## choosing the replicate with lowest detection p-value as the representative beta value
unique_probes <- readRDS("/Users/hjj/Downloads/EPIC_v2_unique_probes.rds")
dup_probes <- readRDS("/Users/hjj/Downloads/EPIC_v2_dup_probes.rds")
detp <- minfi::detectionP(RGset)
rep_probes <- intersect(rownames(detp), dup_probes$TargetID) %>% sort()
length(rep_probes)
replicate_probes_detP <- detp[rep_probes, ]
input_beta_all <- getBeta(MSet.fnorm)
input_beta <- input_beta_all
select_probes <- intersect(rownames(input_beta), unique_probes$TargetID)
length(select_probes) 
unique_probes <- unique_probes[select_probes, ]
input_beta <- input_beta[select_probes, ]
dim(input_beta)
input_beta_unique <- input_beta
betas <- input_beta_all
dup_probes <- dup_probes %>% filter(TargetID %in% rownames(betas))
replicate_probes_detP <- replicate_probes_detP[intersect(rownames(betas), rownames(replicate_probes_detP)),]
betas_dup <- betas[rownames(replicate_probes_detP),]
x_detp <- replicate_probes_detP%>% as.data.frame()
x_detp$TargetID <- rownames(x_detp)
x_detp <- left_join(x_detp, dup_probes)
x_detp_matrix <- NULL
x_probe_select <- NULL
for (i in 1:ncol(replicate_probes_detP)){
    print(i)
    x_tmp <- x_detp[,c(i, ncol(x_detp)-1,ncol(x_detp))]
    colnames(x_tmp)[1] <- "pval"
    x_tmp <- cbind(x_tmp, data.frame(beta = betas_dup[, colnames(x_detp)[i]]))
    tmp <- slice_min(group_by(x_tmp, NAME), pval, n = 1, with_ties = FALSE)
    tmp <- as.data.frame(tmp)
    probe_select <- tmp %>% select(NAME, TargetID,pval)
    colnames(probe_select)[2:3] <- c(colnames(x_detp)[i], paste0(colnames(x_detp)[i], "_detP"))
    if(i ==1){
        x_probe_select <- probe_select
    }else{
        x_probe_select <- left_join(x_probe_select,probe_select)
    }
    y <- tmp[, 4] %>% as.matrix()
    rownames(y) <- tmp$TargetID
    colnames(y) <- colnames(x_detp)[i]
    x_detp_matrix <- cbind(x_detp_matrix, y)
}
dim(x_detp_matrix)
EPICv2_beta_matrix_by_detp <- rbind(input_beta_unique, x_detp_matrix)
beta_v2 <- EPICv2_beta_matrix_by_detp
dim(beta_v2)     
head(beta_v2)

##M-value calculation
unique_probes <- readRDS("/Users/hjj/Downloads/EPIC_v2_unique_probes.rds")
dup_probes <- readRDS("/Users/hjj/Downloads/EPIC_v2_dup_probes.rds")
input_M_all <- getM(MSet.fnorm)   
detp <- minfi::detectionP(RGset)
select_probes <- intersect(rownames(input_M_all), unique_probes$TargetID)
M_unique <- input_M_all[select_probes, ]
dup_probes <- dup_probes %>% filter(TargetID %in% rownames(input_M_all))
replicate_probes_detP <- detp[dup_probes$TargetID, ]
M_dup <- input_M_all[rownames(replicate_probes_detP), ]  
x_detp <- replicate_probes_detP %>% as.data.frame()
x_detp$TargetID <- rownames(x_detp)
x_detp <- left_join(x_detp, dup_probes)
x_detp_matrix_M <- NULL
x_probe_select <- NULL
for (i in 1:ncol(replicate_probes_detP)) {
    print(i)
    x_tmp <- x_detp[, c(i, ncol(x_detp) - 1, ncol(x_detp))]
    colnames(x_tmp)[1] <- "pval"
    x_tmp <- cbind(x_tmp, data.frame(Mval = M_dup[, colnames(x_detp)[i]]))
    tmp <- slice_min(group_by(x_tmp, NAME), pval, n = 1, with_ties = FALSE)
    tmp <- as.data.frame(tmp)
    probe_select <- tmp %>% select(NAME, TargetID, pval)
    colnames(probe_select)[2:3] <- c(colnames(x_detp)[i], paste0(colnames(x_detp)[i], "_detP"))
    if (i == 1) {
        x_probe_select <- probe_select
    } else {
        x_probe_select <- left_join(x_probe_select, probe_select)
    }
    y <- tmp[, 4] %>% as.matrix()
    rownames(y) <- tmp$TargetID
    colnames(y) <- colnames(x_detp)[i]
    x_detp_matrix_M <- cbind(x_detp_matrix_M, y)
}
EPICv2_M_matrix <- rbind(M_unique, x_detp_matrix_M)
M_v2 <- EPICv2_M_matrix
dim(M_v2)   
head(M_v2)

#----EPICv1 -----
# Read IDAT files and create RGset
path <- "/Users/hjj/Desktop/EPICv1_Honda/Junjie_analysis" 
targets <- read.metharray.sheet(path)  
baseDir <- system.file(path, package = "minfiData")
RGset <- read.metharray.exp(targets = targets, verbose = TRUE, extended = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

#test
betas <- getBeta(RGset)
dim(betas)
betas[1:5, 1:5]


# QC

## Control metrics
baseDir <- "/Users/hjj/Desktop/EPICv1_Honda/Junjie_analysis"   
idat_dir <- file.path(baseDir, "IDAT")                        
sample_info <- read.metharray.sheet(baseDir)
sample_info$Basename <- file.path(
idat_dir,
paste0(sample_info$Slide, "_", sample_info$Array)
)
ewas_meth <- read_idats(sample_info$Basename)
ctrls <- control_metrics(ewas_meth)
sample_failure(ctrls) 

## methylated/unmethylated intensity check (minfi)
mset_raw <- preprocessRaw(RGset) 
qc <- getQC(mset_raw)
badSampleCutoff <- 10.5 #default value
meds <- (qc$mMed + qc$uMed)/2
whichBad <- which((meds < badSampleCutoff)) 
controlqc_minfi <- plotQC(qc, badSampleCutoff = badSampleCutoff)

## Sex prediction (ewastools)
check_sexV2 <- function (raw) 
{
    if (!all(c("M", "U", "manifest") %in% names(raw))) 
        stop("Invalid argument")
    with(raw, {
        chrX = manifest[chr == "chrX", index] ## manifest chr name updated
        chrY = manifest[chr == "chrY", index]## manifest chr name updated
        auto = manifest[!chr %in% c("chrX", "chrY"), index]
        chrX = M[chrX, , drop = FALSE] + U[chrX, , drop = FALSE]
        chrY = M[chrY, , drop = FALSE] + U[chrY, , drop = FALSE]
        auto = M[auto, , drop = FALSE] + U[auto, , drop = FALSE]
        chrX = colMeans(chrX, na.rm = TRUE)
        chrY = colMeans(chrY, na.rm = TRUE)
        auto = colMeans(auto, na.rm = TRUE)
        chrX = chrX/auto
        chrY = chrY/auto
        return(list(X = chrX, Y = chrY))
    })
}
predicted_sex <- check_sexV2(ewas_meth)
ewastools_predicted_sex <- predict_sex(predicted_sex$X,predicted_sex$Y)

# Preprocessing

## create noob object (minfi)
noob <- preprocessNoob(RGset)

## noob and functional normalization (minfi)
nPC=2 #default or choose the top PCs that accounted for 90% or more of the variance.
MSet.fnorm <- preprocessFunnorm(RGset, nPCs = nPC, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)

## EPICv1_beta_matrix 
beta_v1 <- getBeta(MSet.fnorm)
dim(beta_v1)    
##EPICv1_M_matrix
M_v1 <- getM(MSet.fnorm)
dim(M_v1)

#Common probes of EPICv1 and EPICv2
probe_v2_cg <- sub("_(BC|TC)[0-9]+$", "", rownames(beta_v2))
common_cg <- intersect(rownames(beta_v1), probe_v2_cg)
beta_v1_common <- beta_v1[common_cg, ]
idx_v2 <- probe_v2_cg %in% common_cg
beta_v2_common <- beta_v2[idx_v2, ]
head(beta_v1_common)
head(beta_v2_common)
dim(beta_v1_common)
dim(beta_v2_common)
probe_v2_cg <- sub("_(BC|TC)[0-9]+$", "", rownames(M_v2))
common_cg <- intersect(rownames(M_v1), probe_v2_cg)
M_v1_common <- M_v1[common_cg, ]
M_v2_common <- M_v2[probe_v2_cg %in% common_cg, ]
head(M_v1_common)
head(M_v2_common)
dim(M_v1_common)   
dim(M_v2_common)  

###############################
### EPIC v2 Off-Target Prediction ###
###############################
probe.features.epicv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
probe.features <- probe.features.epicv2
head(probe.features)
nrow(probe.features)   
FilteredProbes <- probe.features[
  rownames(probe.features) %in% rownames(beta_v2_common),
  c("chr", "pos", "strand", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
#Rename Filtered Probes
colnames(FilteredProbes) <- c("CHR","MAPINFO","Strand","gene","feature","cgi")
nrow(FilteredProbes)   
FilteredProbes$CHR <- as.factor(FilteredProbes$CHR)
levels(FilteredProbes$CHR)

gr_FilteredProbes<- GRanges(seqnames = FilteredProbes$CHR,                            ranges=IRanges(start=FilteredProbes$MAPINFO,width=1),
                            strand  = c("*"),
                            genes   = FilteredProbes$gene,
                            feature = FilteredProbes$feature,
                            cgi     = FilteredProbes$cgi,
                            probe_ID = as.factor(rownames(FilteredProbes)) )
head(gr_FilteredProbes)

#Nearby Probe Identification (EPIC v2) 
#offtarget_folder <- "/Users/hjj/Desktop/EPICv2_Junjie/myog_gRNA_offtarget_bedfiles_hg38/"
offtarget_folder <- "/Users/hjj/Desktop/EPICv2_Junjie/gRNA_offtarget_bedfiles_hg38/"

limit_of_distance <- 500

#sgRNA_names <- c("myog")
sgRNA_names <- c("pj3","pj6","pj7","pj8")

offtargets_8mer_probes  <- as.factor(vector())
offtargets_12mer_probes <- as.factor(vector())
offtargets_20mer_probes <- as.factor(vector())

for (k in c(8,12,20)) {

  print(paste0(k, "mer_offtargets"))
  
  ov_probes_all <- as.vector(vector())
  ov_genes_all  <- as.vector(vector())
  
  offtargets_all <- as.factor(vector())
  offtargets_probes <- as.factor(vector())
  
  for ( i in sgRNA_names) {

    print("***************************")

    print(paste0("offtarget genes for gRNA ",i," ; ",k,"mer_offtarget"))
    
    offtarget_input <- paste0(offtarget_folder,i,"_",k,"mer_crisprdir_result.bed")
    
    offtargets <- read.delim(file = offtarget_input, header=F, sep="\t",skip = 1)
    
    if (length(offtargets) == 1 ) {
      rm(offtargets)
      print(paste0("total site number  = 0"))
      print(paste0("total probe number = 0"))
      print(paste0("total gene number  = 0"))
    }
    else {
      
      print(paste0("total site number  = ", nrow(offtargets)))    
      
      colnames(offtargets)[1:3] <- c("CHR","start","end")
      offtargets <- offtargets[offtargets$CHR %in% levels(FilteredProbes$CHR),]
      gr_offtargets <- GRanges(seqname = offtargets$CHR,
                               ranges  = IRanges(start = offtargets$start, end = offtargets$end),
                               strand  = "*" )

      gr_offtargets <- resize(gr_offtargets, width = 2*limit_of_distance, fix="center")
  
      ov_probes <- findOverlaps(gr_offtargets,gr_FilteredProbes)
      ov_probes <- subjectHits(ov_probes)
      ov_probes <- unique(ov_probes)
      ov_probes <- ov_probes[!is.na(ov_probes)]
      ov_probes <- ov_probes[!(ov_probes == "")]
      print(paste0("# off-target probes; ",length(ov_probes)))
  
genes_raw <- mcols(gr_FilteredProbes[ov_probes])$genes
genes_raw <- genes_raw[!is.na(genes_raw)]
genes_raw <- genes_raw[genes_raw != ""]
genes_one_per_probe <- sapply(strsplit(genes_raw, ";"), `[`, 1)
genes_one_per_probe <- genes_one_per_probe[!is.na(genes_one_per_probe)]
genes_one_per_probe <- genes_one_per_probe[genes_one_per_probe != ""]
ov_genes <- unique(genes_one_per_probe)
print(paste0("# off-target genes; ",length(ov_genes)))
print(head(ov_genes))
  
      write.table (ov_genes, file 
= paste0 ("EPICv2_offtargets_",i,"_",k,"mer_range_dis_",limit_of_distance,".txt"), 
                   sep = "\t",
                   quote = FALSE, 
                   row.names = FALSE,
                   col.names = FALSE )

      ov_probes_all <- c(ov_probes_all,ov_probes)
      ov_genes_all  <- c(ov_genes_all,ov_genes)
  
    }
  }
  
  print("------- total --------")
  
  print(paste0("# all off-target probes with redundancy; ",length(ov_probes_all)))
  ov_probes_all <- unique(ov_probes_all)
  ov_probes_all <- ov_probes_all[!is.na(ov_probes_all)]
  ov_probes_all <- ov_probes_all[!(ov_probes_all == "")]
  print(paste0("# all off-target probes; ",length(ov_probes_all)))
  offtargets_probes <- mcols(gr_FilteredProbes[ov_probes_all])$probe_ID
  
  ov_probes_all_genes <- mcols(gr_FilteredProbes[unique(ov_probes_all)])$genes
  ov_probes_all_genes <- unique(ov_probes_all_genes)
  ov_probes_all_genes <- ov_probes_all_genes[!is.na(ov_probes_all_genes)]
  ov_probes_all_genes <- ov_probes_all_genes[!(ov_probes_all_genes == "")]
  
  ov_genes_all <- unique(ov_genes_all)
  ov_genes_all <- ov_genes_all[!is.na(ov_genes_all)]
  ov_genes_all <- ov_genes_all[!(ov_genes_all == "")]
  print(paste0("# all off-target genes; ",length(ov_genes_all)))
  offtargets_genes  <- ov_genes_all

  assign(paste0("offtargets_",k,"mer_probes"), offtargets_probes)
  assign(paste0("offtargets_",k,"mer_genes"),  offtargets_genes)

  print("*********************")
  
  }

print(paste0("# offtargets_8mer_probes; ", length(offtargets_8mer_probes)))
print(paste0("# offtargets_12mer_probes; ", length(offtargets_12mer_probes)))
print(paste0("# offtargets_20mer_probes; ", length(offtargets_20mer_probes)))
print("*********************")

#####################
### EPIC v2 dot plot ###
#####################

M_use <- M_v2_common
beta_use <- beta_v2_common

# Positive vs Control  
control_samples  <- c("208567300043_R01C01", "208567300043_R02C01", "208567300043_R03C01")
positive_samples <- c("208567300043_R06C01", "208567300043_R07C01","208567300043_R08C01")

AllAve_Positive <- rowMeans(beta_use[, positive_samples], na.rm = TRUE)
AllAve_Control <- rowMeans(beta_use[, control_samples], na.rm = TRUE)

AllProbes <- data.frame(
  Control = AllAve_Control,
  Positive = AllAve_Positive)

offtargets_8mer_probes_v2 <- offtargets_8mer_probes
offtargets_12mer_probes_v2 <- offtargets_12mer_probes

idx_ns    <- !(rownames(beta_use) %in% c(offtargets_8mer_probes_v2,
                                        offtargets_12mer_probes_v2))
idx_8mer  <-  rownames(beta_use) %in% offtargets_8mer_probes_v2
idx_12mer <-  rownames(beta_use) %in% offtargets_12mer_probes_v2

#plot
png("offtarget_Positive_vs_Control_plot.png", width = 1000, height = 1000)

plot(
  x = AllProbes$Control,
  y = AllProbes$Positive,
  type = "n",
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp  = 1,
  xlab = "Mean methylation (Control)",
  ylab = "Mean methylation (Positive)",
  main = "Scatter plot of mean methylation"
)

points(
  AllProbes$Control[idx_ns],
  AllProbes$Positive[idx_ns],
  pch = 16,
  cex = 0.6,
  col = adjustcolor("grey50", alpha.f = 0.2)
)

points(
  AllProbes$Control[idx_8mer],
  AllProbes$Positive[idx_8mer],
  pch = 16,
  cex = 0.6,
  col = adjustcolor("blue", alpha.f = 1)
)

points(
  AllProbes$Control[idx_12mer],
  AllProbes$Positive[idx_12mer],
  pch = 16,
  cex = 2.0,
  col = adjustcolor("red", alpha.f = 1)
)

abline(a = 0, b = 1, col = "black", lwd = 1, lty = 2)
abline(a =  0.2, b = 1, col = "black", lwd = 1, lty = 2)
abline(a = -0.2, b = 1, col = "black", lwd = 1, lty = 2)

fit <- lm(Positive ~ Control, data = AllProbes)
R2  <- summary(fit)$r.squared
abline(fit, col = "black", lwd = 2)

text(
  x = 0.90, y = 0,
  labels = paste0("R² = ", round(R2, 4)),
  col = "black", cex = 1.2, adj = 0
)

dev.off()

##Count the percentage of 12mer and 8mer offtargets (EPIC v2)
deltaBeta_PvC <- AllProbes$Positive - AllProbes$Control   

cutoff <- 0.2

idx_hypo  <- deltaBeta_PvC < -cutoff
idx_mid   <- abs(deltaBeta_PvC) < cutoff
idx_hyper <- deltaBeta_PvC >  cutoff
idx_all <- idx_hypo | idx_mid | idx_hyper

is_8mer  <- rownames(beta_use) %in% offtargets_8mer_probes_v2
is_12mer <- rownames(beta_use) %in% offtargets_12mer_probes_v2
is_other <- !(is_8mer | is_12mer)

#8mer percentages
pct_8mer_hypo <- 100 * sum(idx_hypo  & is_8mer) / sum(idx_hypo)
pct_8mer_mid  <- 100 * sum(idx_mid   & is_8mer) / sum(idx_mid)
pct_8mer_hyper<- 100 * sum(idx_hyper & is_8mer) / sum(idx_hyper)
pct_8mer_total <- 100 * sum(idx_all & is_8mer) / sum(idx_all)

#12mer percentages
pct_12mer_hypo <- 100 * sum(idx_hypo  & is_12mer) / sum(idx_hypo)
pct_12mer_mid  <- 100 * sum(idx_mid   & is_12mer) / sum(idx_mid)
pct_12mer_hyper<- 100 * sum(idx_hyper & is_12mer) / sum(idx_hyper)
pct_12mer_total <- 100 * sum(idx_all & is_12mer) / sum(idx_all)

#other
pct_other_hypo <- 100 * sum(idx_hypo & is_other) / sum(idx_hypo)
pct_other_mid <- 100 * sum(idx_mid & is_other) / sum(idx_mid)
pct_other_hyper <- 100 * sum(idx_hyper & is_other) / sum(idx_hyper)
pct_other_total <- 100 * sum(idx_all & is_other) / sum(idx_all)

offtarget_pct <- data.frame(
  DeltaBeta_Group = c("12mer", "8mer", "other"),
 `Δβ> 0.2` = c(pct_12mer_hyper, pct_8mer_hyper, pct_other_hyper),
 `|Δβ| < 0.2`  = c(pct_12mer_mid, pct_8mer_mid, pct_other_mid), 
`Δβ< -0.2` = c(pct_12mer_hypo, pct_8mer_hypo, pct_other_hypo),
`total` = c(pct_12mer_total, pct_8mer_total, pct_other_total),
check.names = FALSE
)

offtarget_pct


##DMP(limma): significant probe/gene number of offtargets 
#Define groups (8 samples)
M_use <- M_v2_common
beta_use <- beta_v2_common
groupB <- factor(
  c("Control","Control","Control",
    "Negative","Negative",
    "Positive","Positive","Positive"),
  levels = c("Control","Negative","Positive")
)

designB <- model.matrix(~0 + groupB)
colnames(designB) <- levels(groupB)

contrastB <- makeContrasts(
  PvC = Positive - Control,
  NvC = Negative - Control,
  PvN = Positive - Negative,
  levels = designB
)

#limma: DMPs (based on M-values)
fit <- lmFit(M_use, designB)
fit <- contrasts.fit(fit, contrastB)
fit <- eBayes(fit)

# DMP
dmp_PvC <- topTable(fit, coef = "PvC", number = Inf, adjust.method = "BH")

#extract the adj.P.Val vector
adjP_PvC <- dmp_PvC[rownames(beta_use), "adj.P.Val"]

# Align adj.P.Val to probe order 
adjP_PvC_full <- adjP_PvC
names(adjP_PvC_full) <- rownames(beta_use)
# Keep only probes with valid statistical results
valid_idx <- !is.na(adjP_PvC_full)

#Define three Δβ intervals
cutoff <- 0.2
sig_idx <- valid_idx & adjP_PvC_full < 0.01

idx_hypo_sig  <- sig_idx & deltaBeta_PvC < -cutoff
idx_hyper_sig <- sig_idx & deltaBeta_PvC >  cutoff

#Define off-target probe sets
is_8mer  <- rownames(beta_use) %in% offtargets_8mer_probes_v2
is_12mer <- rownames(beta_use) %in% offtargets_12mer_probes_v2
is_other <- !(is_8mer | is_12mer)

sum(idx_hyper_sig & is_12mer)
sum(idx_hyper_sig & is_8mer)

anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

count_genes_from_probes <- function(probe_ids, anno) {
genes <- anno[probe_ids, "UCSC_RefGene_Name"]
genes <- genes[!is.na(genes) & genes != ""]
gene_vec <- unique(unlist(strsplit(genes, ";")))
length(gene_vec)
}

n_gene_hyper_8mer <- count_genes_from_probes(
rownames(beta_use)[idx_hyper_sig & is_8mer], anno)
n_gene_hyper_12mer <- count_genes_from_probes(
rownames(beta_use)[idx_hyper_sig & is_12mer], anno)

#Hyper_8mer_genes in csv
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
probes_hyper_8mer <- rownames(beta_use)[idx_hyper_sig & is_8mer]
genes_hyper_8mer <- anno[probes_hyper_8mer, "UCSC_RefGene_Name"]
genes_hyper_8mer <- genes_hyper_8mer[!is.na(genes_hyper_8mer) & genes_hyper_8mer != ""]
gene_list_hyper_8mer <- sort(
unique(unlist(strsplit(genes_hyper_8mer, ";")))
)

length(gene_list_hyper_8mer)  
head(gene_list_hyper_8mer)

write.csv(
data.frame(Gene = gene_list_hyper_8mer),
file = "Hyper_8mer_genes_EPICv2.csv",
row.names = FALSE
)
 
###############################
### EPIC v1 Off-Target Prediction ###
###############################
probe.features.epicv1 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
probe.features <- probe.features.epicv1
head(probe.features)
nrow(probe.features)    
FilteredProbes <- probe.features[
  rownames(probe.features) %in% rownames(beta_v1_common),
  c("chr", "pos", "strand", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")
]
colnames(FilteredProbes) <- c("CHR","MAPINFO","Strand","gene","feature","cgi")
nrow(FilteredProbes)  
FilteredProbes$CHR <- as.factor(FilteredProbes$CHR)
levels(FilteredProbes$CHR)

gr_FilteredProbes<- GRanges(seqnames = FilteredProbes$CHR,                            ranges=IRanges(start=FilteredProbes$MAPINFO,width=1),
                            strand  = c("*"),
                            genes   = FilteredProbes$gene,
                            feature = FilteredProbes$feature,
                            cgi     = FilteredProbes$cgi,
                            probe_ID = as.factor(rownames(FilteredProbes)) )
head(gr_FilteredProbes)

# Nearby Probe Identification (EPIC v1) 
#offtarget_folder <- "/Users/hjj/Desktop/EPICv1_Honda/myog_offtarget_bedfiles_hg19/"
offtarget_folder <- "/Users/hjj/Desktop/EPICv1_Honda/gRNA_offtarget_bedfiles_hg19/"

limit_of_distance <- 500

#sgRNA_names <- c("myog")
sgRNA_names <- c("pj3","pj6","pj7","pj8")

offtargets_8mer_probes  <- as.factor(vector())
offtargets_12mer_probes <- as.factor(vector())
offtargets_20mer_probes <- as.factor(vector())

for (k in c(8,12,20)) {

  print(paste0(k, "mer_offtargets"))
  
  ov_probes_all <- as.vector(vector())
  ov_genes_all  <- as.vector(vector())
  
  offtargets_all <- as.factor(vector())
  offtargets_probes <- as.factor(vector())
  
  for ( i in sgRNA_names) {

    print("***************************")

    print(paste0("offtarget genes for gRNA ",i," ; ",k,"mer_offtarget"))
    
    offtarget_input <- paste0(offtarget_folder,i,"_",k,"mer_crisprdir_result.bed")
    
    offtargets <- read.delim(file = offtarget_input, header=F, sep="\t",skip = 1)
    
    if (length(offtargets) == 1 ) {
      rm(offtargets)
      print(paste0("total site number  = 0"))
      print(paste0("total probe number = 0"))
      print(paste0("total gene number  = 0"))
    }
    else {
      
      print(paste0("total site number  = ", nrow(offtargets)))    
      
      colnames(offtargets)[1:3] <- c("CHR","start","end")
      offtargets <- offtargets[offtargets$CHR %in% levels(FilteredProbes$CHR),]
      gr_offtargets <- GRanges(seqname = offtargets$CHR,
                               ranges  = IRanges(start = offtargets$start, end = offtargets$end),
                               strand  = "*" )

      gr_offtargets <- resize(gr_offtargets, width = 2*limit_of_distance, fix="center")
  
      ov_probes <- findOverlaps(gr_offtargets,gr_FilteredProbes)
      ov_probes <- subjectHits(ov_probes)
      ov_probes <- unique(ov_probes)
      ov_probes <- ov_probes[!is.na(ov_probes)]
      ov_probes <- ov_probes[!(ov_probes == "")]
      print(paste0("# off-target probes; ",length(ov_probes)))
  
genes_raw <- mcols(gr_FilteredProbes[ov_probes])$genes
genes_raw <- genes_raw[!is.na(genes_raw)]
genes_raw <- genes_raw[genes_raw != ""]
genes_split <- unlist(strsplit(genes_raw, ";"))
genes_split <- genes_split[!is.na(genes_split)]
genes_split <- genes_split[genes_split != ""]
genes_unique <- unique(genes_split)
ov_genes <- genes_unique
print(paste0("# off-target genes; ",length(ov_genes)))
print(head(ov_genes))
  
write.table (ov_genes,file = 
paste0 ("EPICv1_offtargets_",i,"_",k,"mer_range_dis_",limit_of_distance,".txt"), 
                   sep = "\t",
                   quote = FALSE, 
                   row.names = FALSE,
                   col.names = FALSE )

       ov_probes_all <- c(ov_probes_all,ov_probes)
       ov_genes_all  <- c(ov_genes_all,ov_genes)
  
    }
  }
  
  print("------- total --------")
  
  print(paste0("# all off-target probes with redundancy; ",length(ov_probes_all)))
  ov_probes_all <- unique(ov_probes_all)
  ov_probes_all <- ov_probes_all[!is.na(ov_probes_all)]
  ov_probes_all <- ov_probes_all[!(ov_probes_all == "")]
  print(paste0("# all off-target probes; ",length(ov_probes_all)))
  offtargets_probes <- mcols(gr_FilteredProbes[ov_probes_all])$probe_ID
  
  ov_probes_all_genes <- mcols(gr_FilteredProbes[unique(ov_probes_all)])$genes
  ov_probes_all_genes <- unique(ov_probes_all_genes)
  ov_probes_all_genes <- ov_probes_all_genes[!is.na(ov_probes_all_genes)]
  ov_probes_all_genes <- ov_probes_all_genes[!(ov_probes_all_genes == "")]
  #print(paste0("# all off-target genes; ",length(ov_probes_all_genes)))

  #print(paste0("# all off-target genes; ",length(unique(ov_genes_all))))
  ov_genes_all <- unique(ov_genes_all)
  ov_genes_all <- ov_genes_all[!is.na(ov_genes_all)]
  ov_genes_all <- ov_genes_all[!(ov_genes_all == "")]
  print(paste0("# all off-target genes; ",length(ov_genes_all)))
  offtargets_genes  <- ov_genes_all

  assign(paste0("offtargets_",k,"mer_probes"), offtargets_probes)
  assign(paste0("offtargets_",k,"mer_genes"),  offtargets_genes)
  
  print("*********************")
  
  }

print(paste0("# offtargets_8mer_probes; ", length(offtargets_8mer_probes)))
print(paste0("# offtargets_12mer_probes; ", length(offtargets_12mer_probes)))
print(paste0("# offtargets_20mer_probes; ", length(offtargets_20mer_probes)))

print("*********************")


#####################
### EPIC v1 dot plot ###
#####################
M_use <- M_v1_common
beta_use <- beta_v1_common

## KAL vs NN
kal_samples  <- c("205096370071_R02C01", "205096370071_R04C01", "205096370071_R05C01")
nn_samples <- c("205096370071_R01C01", "205096370071_R03C01")

AllAve_KAL  <- rowMeans(beta_use[, kal_samples],  na.rm = TRUE)
AllAve_NN <- rowMeans(beta_use[, nn_samples], na.rm = TRUE)

AllProbes <- data.frame(
  KAL  = AllAve_KAL,
  NN = AllAve_NN)

offtargets_8mer_probes_v1 <- offtargets_8mer_probes
offtargets_12mer_probes_v1 <- offtargets_12mer_probes 

idx_ns    <- !(rownames(beta_use) %in% c(offtargets_8mer_probes_v1,
                                        offtargets_12mer_probes_v1))
idx_8mer  <-  rownames(beta_use) %in% offtargets_8mer_probes_v1
idx_12mer <-  rownames(beta_use) %in% offtargets_12mer_probes_v1

#plot
png("offtarget_KAL_vs_NN_plot.png", width = 1000, height = 1000)

plot(
  x = AllProbes$NN,
  y = AllProbes$KAL,
  type = "n",
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp  = 1,
  xlab = "Mean methylation (NN)",
  ylab = "Mean methylation (KAL)",
  main = "Scatter plot of mean methylation"
)

points(
  AllProbes$NN[idx_ns],
  AllProbes$KAL[idx_ns],
  pch = 16,
  cex = 0.6,
  col = adjustcolor("grey50", alpha.f = 0.2)
)

points(
  AllProbes$NN[idx_8mer],
  AllProbes$KAL[idx_8mer],
  pch = 16,
  cex = 0.6,
  col = adjustcolor("blue", alpha.f = 1)
)

points(
  AllProbes$NN[idx_12mer],
  AllProbes$KAL[idx_12mer],
  pch = 16,
  cex = 2.0,
  col = adjustcolor("red", alpha.f = 1)
)

abline(a = 0, b = 1, col = "black", lwd = 1, lty = 2)
abline(a =  0.2, b = 1, col = "black", lwd = 1, lty = 2)
abline(a = -0.2, b = 1, col = "black", lwd = 1, lty = 2)

fit <- lm(KAL ~ NN, data = AllProbes)
R2  <- summary(fit)$r.squared
abline(fit, col = "black", lwd = 2)

text(
  x = 0.90, y = 0,
  labels = paste0("R² = ", round(R2, 4)),
  col = "black", cex = 1.2, adj = 0
)

dev.off()


##Count the percentage of 12mer and 8mer offtargets (EPIC v1)

deltaBeta_KvN <- AllProbes$KAL - AllProbes$NN   # Δβ 

cutoff <- 0.2

idx_hypo  <- deltaBeta_KvN < -cutoff
idx_mid   <- abs(deltaBeta_KvN) < cutoff
idx_hyper <- deltaBeta_KvN >  cutoff
idx_all <- idx_hypo | idx_mid | idx_hyper

is_8mer  <- rownames(beta_use) %in% offtargets_8mer_probes_v1
is_12mer <- rownames(beta_use) %in% offtargets_12mer_probes_v1
is_other <- !(is_8mer | is_12mer)

#Calculate percentages
#8mer percentages
pct_8mer_hypo <- 100 * sum(idx_hypo  & is_8mer) / sum(idx_hypo)
pct_8mer_mid  <- 100 * sum(idx_mid   & is_8mer) / sum(idx_mid)
pct_8mer_hyper<- 100 * sum(idx_hyper & is_8mer) / sum(idx_hyper)
pct_8mer_total <- 100 * sum(idx_all & is_8mer) / sum(idx_all)

#12mer percentages
pct_12mer_hypo <- 100 * sum(idx_hypo  & is_12mer) / sum(idx_hypo)
pct_12mer_mid  <- 100 * sum(idx_mid   & is_12mer) / sum(idx_mid)
pct_12mer_hyper<- 100 * sum(idx_hyper & is_12mer) / sum(idx_hyper)
pct_12mer_total <- 100 * sum(idx_all & is_12mer) / sum(idx_all)

#other
pct_other_hypo <- 100 * sum(idx_hypo & is_other) / sum(idx_hypo)
pct_other_mid <- 100 * sum(idx_mid & is_other) / sum(idx_mid)
pct_other_hyper <- 100 * sum(idx_hyper & is_other) / sum(idx_hyper)
pct_other_total <- 100 * sum(idx_all & is_other) / sum(idx_all)

offtarget_pct <- data.frame(
  DeltaBeta_Group = c("12mer", "8mer", "other"),
 `Δβ> 0.2` = c(pct_12mer_hyper, pct_8mer_hyper, pct_other_hyper),
 `|Δβ| < 0.2`  = c(pct_12mer_mid, pct_8mer_mid, pct_other_mid), 
`Δβ< -0.2` = c(pct_12mer_hypo, pct_8mer_hypo, pct_other_hypo),
`total` = c(pct_12mer_total, pct_8mer_total, pct_other_total),
check.names = FALSE
)

offtarget_pct

#DMP(limma): significant probe/gene number  
M_use <- M_v1_common
beta_use <- beta_v1_common
groupB <- factor(
  c("Ctrl","Ctrl",
"NN","NN",
    "KAL","KAL","KAL"),
  levels = c("Ctrl","NN","KAL")
)

designB <- model.matrix(~0 + groupB)
colnames(designB) <- levels(groupB)

contrastB <- makeContrasts(
  KAL_vs_NN = KAL - NN,
  KAL_vs_Ctrl = KAL - Ctrl,
  NN_vs_Ctrl = NN - Ctrl,
  levels = designB
)

fit <- lmFit(M_use, designB)
fit <- contrasts.fit(fit, contrastB)
fit <- eBayes(fit)

# DMP
dmp_KvN <- topTable(fit, coef = "KAL_vs_NN", number = Inf, adjust.method = "BH")

adjP_KvN <- dmp_KvN[rownames(beta_use), "adj.P.Val"]

adjP_KvN_full <- adjP_KvN
names(adjP_KvN_full) <- rownames(beta_use)
valid_idx <- !is.na(adjP_KvN_full)

cutoff <- 0.2
sig_idx <- valid_idx & adjP_KvN_full < 0.01

idx_hypo_sig  <- sig_idx & deltaBeta_KvN < -cutoff
idx_hyper_sig <- sig_idx & deltaBeta_KvN >  cutoff

is_8mer  <- rownames(beta_use) %in% offtargets_8mer_probes_v1
is_12mer <- rownames(beta_use) %in% offtargets_12mer_probes_v1
is_other <- !(is_8mer | is_12mer)

sum(idx_hyper_sig & is_12mer)
#[1] 11
sum(idx_hyper_sig & is_8mer)
#[1] 754

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
count_genes_from_probes <- function(probe_ids, anno) {
genes <- anno[probe_ids, "UCSC_RefGene_Name"]
genes <- genes[!is.na(genes) & genes != ""]
gene_vec <- unique(unlist(strsplit(genes, ";")))
length(gene_vec)
}
n_gene_hyper_8mer <- count_genes_from_probes(
rownames(beta_use)[idx_hyper_sig & is_8mer], anno)
n_gene_hyper_12mer <- count_genes_from_probes(
rownames(beta_use)[idx_hyper_sig & is_12mer], anno)
n_gene_hyper_8mer 
n_gene_hyper_12mer

#Hyper_8mer_genes in csv
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
probes_hyper_8mer <- rownames(beta_use)[idx_hyper_sig & is_8mer]
genes_hyper_8mer <- anno[probes_hyper_8mer, "UCSC_RefGene_Name"]
genes_hyper_8mer <- genes_hyper_8mer[!is.na(genes_hyper_8mer) & genes_hyper_8mer != ""]
gene_list_hyper_8mer <- sort(
unique(unlist(strsplit(genes_hyper_8mer, ";")))
)

length(gene_list_hyper_8mer)  
head(gene_list_hyper_8mer)

write.csv(
data.frame(Gene = gene_list_hyper_8mer),
file = "Hyper_8mer_genes_EPICv1.csv",
row.names = FALSE
)

##################################################
### Comparison of EPIC v1 and EPIC v2 Off-Target Sites ###
##################################################
#Comparing the 8-mer and 12-mer off-target percentage of deltabeta>0.2
# Positive vs Control (EPICv2)
control_samples  <- c("208567300043_R01C01", "208567300043_R02C01", "208567300043_R03C01")
positive_samples <- c("208567300043_R06C01", "208567300043_R07C01","208567300043_R08C01")

AllAve_Positive <- rowMeans(beta_v2_common[, positive_samples], na.rm = TRUE)
AllAve_Control <- rowMeans(beta_v2_common[, control_samples], na.rm = TRUE)

AllProbes <- data.frame(
  Control = AllAve_Control,
  Positive = AllAve_Positive)

deltaBeta_PvC <- AllProbes$Positive - AllProbes$Control   # Δβ

cutoff <- 0.2

idx_hypo  <- deltaBeta_PvC < -cutoff
idx_mid   <- abs(deltaBeta_PvC) < cutoff
idx_hyper <- deltaBeta_PvC >  cutoff
idx_all <- idx_hypo | idx_mid | idx_hyper

is_8mer  <- rownames(beta_v2_common) %in% offtargets_8mer_probes_v2
is_12mer <- rownames(beta_v2_common) %in% offtargets_12mer_probes_v2
is_other <- !(is_8mer | is_12mer)

##8mer percentages
pct_8mer_hypo <- 100 * sum(idx_hypo  & is_8mer) / sum(is_8mer)
pct_8mer_mid  <- 100 * sum(idx_mid   & is_8mer) / sum(is_8mer)
pct_8mer_hyper<- 100 * sum(idx_hyper & is_8mer) / sum(is_8mer)
pct_8mer_total <- 100 * sum(idx_all & is_8mer) / sum(is_8mer)

n_8mer_hyper_v2 <- sum(idx_hyper & is_8mer)
n_8mer_total_v2 <- sum(idx_all & is_8mer)

##12mer percentages
pct_12mer_hypo <- 100 * sum(idx_hypo  & is_12mer) / sum(is_12mer)
pct_12mer_mid  <- 100 * sum(idx_mid   & is_12mer) / sum(is_12mer)
pct_12mer_hyper<- 100 * sum(idx_hyper & is_12mer) / sum(is_12mer)
pct_12mer_total <- 100 * sum(idx_all & is_12mer) / sum(is_12mer)

n_12mer_hyper_v2 <- sum(idx_hyper & is_12mer)
n_12mer_total_v2 <- sum(idx_all & is_12mer)

##other
pct_other_hypo <- 100 * sum(idx_hypo & is_other) / sum(is_other)
pct_other_mid <- 100 * sum(idx_mid & is_other) / sum(is_other)
pct_other_hyper <- 100 * sum(idx_hyper & is_other) / sum(is_other)
pct_other_total <- 100 * sum(idx_all & is_other) / sum(is_other)

n_other_hyper_v2 <- sum(idx_hyper & is_other)
n_other_total_v2 <- sum(idx_all & is_other)

#total
pct_hypo <- 100 * sum(idx_hypo) / sum(idx_all)
pct_mid <- 100 * sum(idx_mid) / sum(idx_all)
pct_hyper <- 100 * sum(idx_hyper) / sum(idx_all)
pct_total <- 100 * sum(idx_all) / sum(idx_all)

offtarget_pct_v2 <- data.frame(
  DeltaBeta_Group = c("12mer", "8mer", "other", "total"),
 `Δβ> 0.2` = c(pct_12mer_hyper, pct_8mer_hyper, pct_other_hyper, pct_hyper),
 `|Δβ| < 0.2`  = c(pct_12mer_mid, pct_8mer_mid, pct_other_mid, pct_mid), 
`Δβ< -0.2` = c(pct_12mer_hypo, pct_8mer_hypo, pct_other_hypo, pct_hypo),
`total` = c(pct_12mer_total, pct_8mer_total, pct_other_total, pct_total),
check.names = FALSE
)
offtarget_pct_v2

#KAL vs NN (EPICv1)
kal_samples  <- c("205096370071_R02C01", "205096370071_R04C01", "205096370071_R05C01")
nn_samples <- c("205096370071_R01C01", "205096370071_R03C01")

AllAve_KAL  <- rowMeans(beta_v1_common[, kal_samples],  na.rm = TRUE)
AllAve_NN <- rowMeans(beta_v1_common[, nn_samples], na.rm = TRUE)

AllProbes <- data.frame(
  KAL  = AllAve_KAL,
  NN = AllAve_NN)

deltaBeta_KvN <- AllProbes$KAL - AllProbes$NN   # Δβ 

cutoff <- 0.2

idx_hypo  <- deltaBeta_KvN < -cutoff
idx_mid   <- abs(deltaBeta_KvN) < cutoff
idx_hyper <- deltaBeta_KvN >  cutoff
idx_all <- idx_hypo | idx_mid | idx_hyper

is_8mer  <- rownames(beta_v1_common) %in% offtargets_8mer_probes_v1
is_12mer <- rownames(beta_v1_common) %in% offtargets_12mer_probes_v1
is_other <- !(is_8mer | is_12mer)

##8mer percentages
pct_8mer_hypo <- 100 * sum(idx_hypo  & is_8mer) / sum(is_8mer)
pct_8mer_mid  <- 100 * sum(idx_mid   & is_8mer) / sum(is_8mer)
pct_8mer_hyper<- 100 * sum(idx_hyper & is_8mer) / sum(is_8mer)
pct_8mer_total <- 100 * sum(idx_all & is_8mer) / sum(is_8mer)

n_8mer_hyper_v1 <- sum(idx_hyper & is_8mer)
n_8mer_total_v1 <- sum(idx_all & is_8mer)

##12mer percentages
pct_12mer_hypo <- 100 * sum(idx_hypo  & is_12mer) / sum(is_12mer)
pct_12mer_mid  <- 100 * sum(idx_mid   & is_12mer) / sum(is_12mer)
pct_12mer_hyper<- 100 * sum(idx_hyper & is_12mer) / sum(is_12mer)
pct_12mer_total <- 100 * sum(idx_all & is_12mer) / sum(is_12mer)

n_12mer_hyper_v1 <- sum(idx_hyper & is_12mer)
n_12mer_total_v1 <- sum(idx_all & is_12mer)

##other
pct_other_hypo <- 100 * sum(idx_hypo & is_other) / sum(is_other)
pct_other_mid <- 100 * sum(idx_mid & is_other) / sum(is_other)
pct_other_hyper <- 100 * sum(idx_hyper & is_other) / sum(is_other)
pct_other_total <- 100 * sum(idx_all & is_other) / sum(is_other)

n_other_hyper_v1 <- sum(idx_hyper & is_other)
n_other_total_v1 <- sum(idx_all & is_other)

##total
pct_hypo <- 100 * sum(idx_hypo) / sum(idx_all)
pct_mid <- 100 * sum(idx_mid) / sum(idx_all)
pct_hyper <- 100 * sum(idx_hyper) / sum(idx_all)
pct_total <- 100 * sum(idx_all) / sum(idx_all)

offtarget_pct_v1 <- data.frame(
  DeltaBeta_Group = c("12mer", "8mer", "other", "total"),
 `Δβ> 0.2` = c(pct_12mer_hyper, pct_8mer_hyper, pct_other_hyper, pct_hyper),
 `|Δβ| < 0.2`  = c(pct_12mer_mid, pct_8mer_mid, pct_other_mid, pct_mid), 
`Δβ< -0.2` = c(pct_12mer_hypo, pct_8mer_hypo, pct_other_hypo, pct_hypo),
`total` = c(pct_12mer_total, pct_8mer_total, pct_other_total, pct_total),
check.names = FALSE
)
offtarget_pct_v1


#Fisher’s exact test (statistic analysis of hyper-probe percentage of EPICv1 vs EPICv2)
##8mer-offtargets
mat <- matrix(c(
  n_8mer_hyper_v1,
  n_8mer_total_v1 - n_8mer_hyper_v1,
  n_8mer_hyper_v2,
  n_8mer_total_v2 - n_8mer_hyper_v2
), nrow = 2, byrow = TRUE)

fisher.test(mat)

## 12mer-offtargets
mat <- matrix(c(
  n_12mer_hyper_v1,
  n_12mer_total_v1 - n_12mer_hyper_v1,
  n_12mer_hyper_v2,
  n_12mer_total_v2 - n_12mer_hyper_v2
), nrow = 2, byrow = TRUE)

fisher.test(mat)

## other-offtargets
mat <- matrix(c(
  n_other_hyper_v1,
  n_other_total_v1 - n_other_hyper_v1,
  n_other_hyper_v2,
  n_other_total_v2 - n_other_hyper_v2
), nrow = 2, byrow = TRUE)

fisher.test(mat)


epicv1 <- data.frame(
  Platform = "EPICv1",
  Type = c("12mer", "8mer", "other"),
  up   = c(16.417910, 6.803995, 3.518918),   # Δβ > 0.2
  mid  = c(83.58209, 93.19600, 96.46855),   # |Δβ| < 0.2
  down = c(0.00000000, 0.00000000, 0.01253286) # Δβ < -0.2
)

epicv2 <- data.frame(
  Platform = "EPICv2",
  Type = c("12mer", "8mer", "other"),
  up   = c(4.4776119, 1.5604102, 0.2472780),   
  mid  = c(95.52239, 98.42176, 99.74526),
  down = c(0.000000000, 0.017833259, 0.007463401)
)

df <- rbind(epicv1, epicv2)

df_long <- pivot_longer(
  df,
  cols = c("up", "mid", "down"),
  names_to = "Delta",
  values_to = "Percent"
)
df_up <- subset(df_long, Delta == "up")

signif_data <- data.frame(
  Type = c("12mer", "8mer", "other"),
  p_value = c(0.04484, 2.2e-16, 2.2e-16),
  signif = c("*", "***", "***"),  # * p < 0.05, *** p < 0.001
  y_position = c(18, 8, 5)
)

##bar plot
p <- ggplot(df_up, aes(x = Type, y = Percent, fill = Platform)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6,
    alpha = 0.6   
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percent)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4
  ) +
  geom_segment(
    data = signif_data,
    aes(x = as.numeric(factor(Type)) - 0.3, 
        xend = as.numeric(factor(Type)) + 0.3,
        y = y_position, 
        yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_text(
    data = signif_data,
    aes(x = as.numeric(factor(Type)), 
        y = y_position + 0.15, 
        label = signif),
    inherit.aes = FALSE,
    size = 8,
    vjust = 0
  ) +
  scale_fill_manual(values = c(
    "EPICv1" = "#377EB8",
    "EPICv2" = "#E41A1C"
  )) +
  theme_classic(base_size = 14) +
  labs(
    y = "Δβ > 0.2 (%)",
    x = "",
    fill = "",
    title = "Off-target effects (Δβ > 0.2)"
  ) +
  coord_cartesian(ylim = c(0, max(df_up$Percent) * 1.25))

print(p)
ggsave("offtarget_bar_plot.png", p, width = 6, height = 5, dpi = 300)


##forest plot (show CI and OR from fisher’s exact test)
df_or <- data.frame(
Type = c("8mer", "12mer", "Other"),
OR   = c(4.60, 4.15, 14.72),
CI_low  = c(3.89, 1.03, 14.02),
CI_high = c(5.47, 24.32, 15.45)
)

p <- ggplot(df_or, aes(x = OR, y = Type)) +
    geom_errorbar(
        aes(xmin = CI_low, xmax = CI_high),
        orientation = "y",
        width = 0.2,
        linewidth = 1
    ) +
    geom_point(size = 3, color = "#DC2626") +
    geom_text(
        aes(
            label = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)
        ),
        vjust = 2,   
        size = 4
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_x_log10() +
    theme_classic(base_size = 14) +
    labs(
        x = "Odds Ratio (log scale)",
        y = "",
        title = "Off-target effect size (EPICv1 vs EPICv2)"
)
print(p)
ggsave("offtarget_forest_plot.png", p, width = 6, height = 5, dpi = 300)



###########################
### Off-Target Site Heatmap ###
###########################

#Heatmap of offtarget_8mer_genes
#EPIC v2
x <- read.table(file="/Users/hjj/Desktop/RNA-seq/60_HTSeq_STAR_TPM_protein_coding.txt", header=T, sep='\t')
x <- x[ , !(colnames(x)=="gene_id")]
head(x)

Hyper_8mer_genes <- read.delim("/Users/hjj/Desktop/EPICv1v2_analysis_251015/Hyper_8mer_genes_EPICv2.csv", header = TRUE,sep = ",")
head(Hyper_8mer_genes)

genelists <- unlist(Hyper_8mer_genes)
head(genelists)
length(genelists)

TPM_min = 10
height=300

nrow(x)  
# Keep only target genes
xSelect <- x[ x$gene_name %in% genelists, ]
nrow(xSelect) 

keep_cols <- c("gene_name",
grep("^CRoffOnly", colnames(xSelect), value = TRUE),
grep("^CRoffwgRNAs", colnames(xSelect), value = TRUE),
grep("^NoTf", colnames(xSelect), value = TRUE))
xSelect <- xSelect[, keep_cols]
xSelect <- xSelect[rowMeans(xSelect[, grep("^NoTf", colnames(xSelect))]) >= TPM_min,]
xmat <- as.matrix(xSelect[, -1])
rownames(xmat) <- xSelect$gene_name
head(xmat)

col_groups <- sub("_.*", "", colnames(xmat))
table(col_groups)

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(xmat)

mat_colors <- list(
group = c(
NoTf         = "blue",
CRoffwgRNAs    = "red"))

xmat_rel_log2 <-  xmat
## CRoffwgRNAs vs NoTf
xmat_rel_log2[,3] <-  ifelse(xmat[,5] == 0, NA, log2(xmat[,3]/xmat[,5]))
xmat_rel_log2[,4] <-  ifelse(xmat[,6] == 0, NA, log2(xmat[,4]/xmat[,6]))
## NoTf vs NoTf（baseline=0）
xmat_rel_log2[,5] <-  ifelse(xmat[,5] == 0, NA, log2(xmat[,5]/xmat[,5]))
xmat_rel_log2[,6] <-  ifelse(xmat[,6] == 0, NA, log2(xmat[,6]/xmat[,6]))
head(xmat_rel_log2,5)

xmat_rel_log2_col_ordered <- xmat_rel_log2[, c(3, 5, 4, 6)]

mat_rel_breaks <- seq(-5, 5, length.out = 1000)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
cellheight = height/nrow(xmat_rel_log2_col_ordered)

##heatmap
my_heatmap <- pheatmap(
mat               = xmat_rel_log2_col_ordered,
color             = cividis(1000),
breaks            = mat_rel_breaks,
border_color      = NA,
cluster_cols      = FALSE,
cluster_rows      = FALSE,
treeheight_row    = 0,
show_colnames     = TRUE,     
show_rownames     = TRUE,
fontsize_col       = 10,
fontsize_row      =10,
annotation_col    = mat_col,
annotation_colors = mat_colors,
drop_levels       = TRUE,
cellwidth         = 20,
cellheight        = cellheight,
main              = paste0("log2(TPM / NoTf TPM), NoTf TPM ≥ ", TPM_min),
gaps_cols         = c(
length(grep("^CRoffOnly", colnames(xmat_rel_log2_col_ordered))),
length(grep("^CRoffOnly|^CRoffwgRNAs", colnames(xmat_rel_log2_col_ordered)))
)
)

my_heatmap

save_pheatmap_png <- function(x, filename, width=2400, height=3000, res = 300) {
png(filename, width = width, height = height, res = res)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()
}

save_pheatmap_png(my_heatmap, paste0("offtargets_8mer_genes",
"_TPM_",TPM_min,
"_rel_log2_heatmap_EPICv2.png"))

write.table (rownames(xmat_rel_log2_col_ordered), 
             file = paste0("offtarget_gene_list_for_heatmap_EPICv2.txt"), 
             sep = "\t",
             quote = FALSE,
             row.names = FALSE,
             col.names = FALSE)




#Heatmap of offtarget_8mer_genes
#EPIC v1
x <- read.table(file="/Users/hjj/Desktop/RNA-seq/60_HTSeq_STAR_TPM_protein_coding.txt", header=T, sep='\t')
x <- x[ , !(colnames(x)=="gene_id")]
head(x)

Hyper_8mer_genes <- read.delim("/Users/hjj/Desktop/EPICv1v2_analysis_251015/Hyper_8mer_genes_EPICv1.csv", header = TRUE,sep = ",")
head(Hyper_8mer_genes)

genelists <- unlist(Hyper_8mer_genes)
head(genelists)
length(genelists)

TPM_min = 10
height=500

nrow(x)  
xSelect <- x[ x$gene_name %in% genelists, ]
nrow(xSelect) 

keep_cols <- c(
  "gene_name",
  grep("^KAL", colnames(xSelect), value = TRUE),
  grep("^NN",   colnames(xSelect), value = TRUE)
)
xSelect <- xSelect[, keep_cols]


xSelect <- xSelect[
  rowMeans(xSelect[, grep("^NN", colnames(xSelect))]) >= TPM_min,]
xmat <- as.matrix(xSelect[, -1])
rownames(xmat) <- xSelect$gene_name
col_groups <- sub("_.*", "", colnames(xmat))
table(col_groups)
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(xmat)
mat_colors <- list(group = c(NN  = "blue", KAL = "red"))

xmat_rel_log2 <-  xmat
## KAL vs NN
xmat_rel_log2[,1] <-  ifelse(xmat[,3] == 0, NA, log2(xmat[,1]/xmat[,3]))
xmat_rel_log2[,2] <-  ifelse(xmat[,4] == 0, NA, log2(xmat[,2]/xmat[,4]))
## NN vs NN（baseline=0）
xmat_rel_log2[,3] <-  ifelse(xmat[,3] == 0, NA, log2(xmat[,3]/xmat[,3]))
xmat_rel_log2[,4] <-  ifelse(xmat[,4] == 0, NA, log2(xmat[,4]/xmat[,4]))
head(xmat_rel_log2,3)

xmat_rel_log2_col_ordered <- xmat_rel_log2[, c(1, 3, 2, 4)]
mat_rel_breaks <- seq(-5, 5, length.out = 1000)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
cellheight = height / nrow(xmat_rel_log2_col_ordered)

KAL_cols <- grep("^KAL", colnames(xmat_rel_log2_col_ordered))

##heatmap
my_heatmap <- pheatmap(
  mat               = xmat_rel_log2_col_ordered,
  color             = cividis(1000),
  breaks            = seq(-5, 5, length.out = 1000),
  border_color      = NA,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  treeheight_row    = 0,
  show_colnames   = TRUE,
  show_rownames   = TRUE,
  fontsize_col       = 10,
fontsize_row      =3,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  cellwidth         = 20,
  cellheight        = cellheight,
  main              = paste0(
    "log2(KAL TPM / NN TPM), NN TPM ≥ ", TPM_min),
gaps_cols         = length(KAL_cols)   
)

my_heatmap

save_pheatmap_png <- function(x, filename, width=5000, height=6000, res = 300) {
png(filename, width = width, height = height, res = res)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()
}

save_pheatmap_png(
  my_heatmap,
  paste0(
    "offtargets_8mer_genes", "_TPM_",TPM_min,
    "_rel_log2_heatmap_EPICv1.png"))

write.table (rownames(xmat_rel_log2_col_ordered), 
             file = paste0("offtarget_gene_list_for_heatmap_EPICv1.txt"), 
             sep = "\t",
             quote = FALSE,
             row.names = FALSE,
             col.names = FALSE)


