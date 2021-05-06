#setwd("/usr4/bf527/sktpndt/Individual_Project/project3_Biologist")
library(dplyr)

# Heatmap clustering ------------------------------------------------------
# Read in CSVs
ahr <- read.csv("CMC_AHR_deseq_norm_counts.csv")
cxr <- read.csv("CORN__CARPXR_deseq_norm_counts.csv")
cyt <- read.csv("SALINE_Cytotoxic_deseq_norm_counts.csv")

# label them separately (easy to ID later)
colnames(ahr)[2:7] <- paste(colnames(ahr)[2:7],"_AHR", sep="")
colnames(cxr)[2:7] <- paste(colnames(cxr)[2:7],"_CXR", sep="")
colnames(cyt)[2:7] <- paste(colnames(cyt)[2:7],"_CYT", sep="")

# Merging them together
merged_df <- 
  dplyr::left_join(ahr, cxr, by=c("X"="X")) %>% 
  dplyr::left_join(y=cyt, by=c("X"="X"))

# Applying filters
# Filter by counts
## Finding appropriate threshold for average counts
plot(log(rowMeans(merged_df[2:7])[order(rowMeans(merged_df[2:19]))]),
     ylab="Log of Mean count",
     main="Distribution of Log(RowMean)")
# Graph starts picking up around 4 on the y axis
# Using threshold of 4
merged_df_filt <- merged_df[log(rowMeans(merged_df[2:19]))>4,]

# Filtering by coefficient of variation (CV)
cv <- function(x){
  result <- sd(x)/mean(x)
  return(result)
}

# Finding an appropriate CV threshold
cvs <- apply(merged_df_filt[2:19], 1, cv)
l <- 250
vals <- seq(min(cvs), max(cvs), length.out = 250)
lengths <- rep(0,l) # vector of lengths
for(i in 1:length(vals)){
  test <- which(cvs<vals[i])
  lengths[i] <- length(test)
}
plot(vals, lengths,
     xlab="CV Thresholds",
     ylab="Number of genes remaining",
     main="Effect of CV Threshold on size of data") 
# Looks like 0.5 is a good cutoff

# Returning rows with CV less than 0.5
cv_tf <- function(x,thresh){
  cvr <- cv(x)
  result <- cvr<thresh
  return(result)
}
merged_df_filt <- merged_df_filt[apply(merged_df_filt[2:19], 1, cv_tf, 0.5),]
# 8386 genes remaining

# Heatmap
heatmap(data.matrix(merged_df_filt[,2:19]))


# Getting DE Genes for GSEA analysis --------------------------------------

# Read csv
cmc <- read.csv("resCMC_AHR_deseq_results.csv")
corn <- read.csv("resCORN_CARPXR_deseq_results.csv")
saline <- read.csv("resSALINE_Cytotoxic_deseq_results.csv")

# Filtering for significance + log2FC
# Following thresholds outlined by the paper
# pvalue < 0.05, abs(log(FC))>1.5
cmc_filt <- cmc %>% 
  dplyr::filter(pvalue<0.05) %>% 
  dplyr::filter(abs(log2FoldChange)>1.5) %>% 
  dplyr::pull(X)
corn_filt <- corn %>% 
  dplyr::filter(pvalue<0.05) %>% 
  dplyr::filter(abs(log2FoldChange)>1.5) %>% 
  dplyr::pull(X)
saline_filt <- saline %>% 
  dplyr::filter(pvalue<0.05) %>% 
  dplyr::filter(abs(log2FoldChange)>1.5) %>% 
  dplyr::pull(X)


# Write to files
write(cmc_filt, file="cmc_res.txt")
write(corn_filt, file="corn_res.txt")
write(saline_filt, file="saline_res.txt")


