#setwd("/usr4/bf527/sktpndt/Individual_Project/project2_Analyst")
library(dplyr)

# Loading in the data
gene_exp <- read.table("gene_exp.diff", header = TRUE)

# Sorting by q-value
gene_exp <- gene_exp[order(gene_exp$q_value),]

# Getting out top 10 differentially expressed genes (report)
t10_diff <- gene_exp[1:10,] %>% 
  dplyr::select(gene, value_1, value_2, log2.fold_change., p_value, q_value)

# Histogram of log2(FC) unfiltered vs filtered
hist(gene_exp$log2.fold_change.,
     main = "Distribution of Fold Change values",
     xlab = "Log2(FC)",
     breaks=20)

# Filtering for only significant genes
gene_exp_filt <- gene_exp %>% 
  dplyr::filter(significant == "yes")

hist(gene_exp_filt$log2.fold_change.,
     main = "Distribution of Fold Change values",
     xlab = "Log2(FC)",
     breaks=20)

# Filtering for positive and negative fold changes
gene_exp_pos <- gene_exp_filt %>% 
  dplyr::filter(log2.fold_change.>0) %>% 
  dplyr::pull(gene) # 642 genes

gene_exp_neg <- gene_exp_filt %>% 
  dplyr::filter(log2.fold_change.<0) %>% 
  dplyr::pull(gene) # 995 genes

# Writing to txt files
write(gene_exp_pos, file="gene_exp_pos.txt", sep=",")
write(gene_exp_neg, file="gene_exp_neg.txt", sep=",")


## Finding the best threshold
# l <- 250
# vals <- seq(0, 10, length.out = l)
# lengths <- rep(0,l) # vector of lengths
# for(i in 1:length(vals)){
#   test <- which(abs(gene_exp$log2.fold_change.)>=vals[i])
#   lengths[i] <- length(test)
# }
# plot(vals, lengths,
#      xlab = "Log2(Fold Change)",
#      ylab = "Number of genes remaining",
#      main = "Effect of Log2(Fold Change) Thresholds on Number of Genes")
# abline(v=3, col="red", lty=3)
# abline(h=lengths[min(which(vals>=3))], col="red", lty=3)
# text(y=lengths[min(which(vals>=3))]+1000, x=3,
#      label="log2(FC)=3",
#      cex=0.8,
#      pos=4)