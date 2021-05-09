setwd("/usr4/bf527/sktpndt/Individual_Project/project2_Analyst")
library(dplyr)

# Loading in the data
gene_exp <- read.table("gene_exp.diff", header = TRUE)

# Sorting by q-value
gene_exp <- gene_exp[order(gene_exp$q_value),]

# Getting out top 10 differentially expressed genes (report)
t10_diff <- gene_exp[1:10,] %>% 
  dplyr::select(gene, value_1, value_2, log2.fold_change., p_value, q_value)

# Write to csv file
#write.csv(t10_diff, "t10_diff.csv")

# Histogram of log2(FC) unfiltered vs filtered
par(mfrow=c(1,2))
hist(gene_exp$log2.fold_change.,
     main = "Unfiltered data",
     font.main=1,
     xlab = "Log2(FC)",
     breaks=20)

# Filtering for only significant genes
gene_exp_filt <- gene_exp %>%
  dplyr::filter(significant == "yes")

hist(gene_exp_filt$log2.fold_change.,
     main = "Filtered data",
     font.main=1,
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
# write(gene_exp_pos, file="gene_exp_pos.txt", sep=",")
# write(gene_exp_neg, file="gene_exp_neg.txt", sep=",")
