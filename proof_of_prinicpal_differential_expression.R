# Load necessary libraries
library(edgeR)
library(limma)
library(tidyverse)

# Set working directory to the location of your count files
setwd('/Users/scampione/Projects/Buck_Institute/Ribo_Profiling/PMID31358845')


######################################################
# RNA-Seq Samples
######################################################

# List directories containing RNA count files
rna_subdirs <- list.files(pattern = "_rna$", full.names = TRUE)
rna_count_files <- c()

# Collect RNA count files
# for (dir in rna_subdirs) {
#   rna_count_files <- c(rna_count_files, 
#                        list.files(path = dir, 
#                                   pattern = "_gene_counts\\.txt$", 
#                                   full.names = TRUE))
# }


rna_count_files <- c("./SRR6761663_rna/SRR6761663_gene_counts.txt",
                      "./SRR6761664_rna/SRR6761664_gene_counts.txt",
                      "./SRR6761665_rna/SRR6761665_gene_counts.txt",
                      "./SRR6761666_rna/SRR6761666_gene_counts.txt")


# Adjust subdirectories format if necessary
rna_subdirs <- gsub("./", "", rna_subdirs)

# Prepare sample information
rna_sample_info <- data.frame(
  sra = gsub("_rna", "", rna_subdirs),
  samples = c('RNA_stress_rep1', 'RNA_stress_rep2', 
              'RNA_normal_rep1', 'RNA_normal_rep2'),
  condition = c('stress', 'stress', 
                'wt', 'wt'))

# Define column names for the count data
col_names <- c("gene_id", "count")

# Read in the count data
rna_count_list <- lapply(rna_count_files, function(f) {
  read_tsv(f, 
           col_names = col_names, 
           col_types = cols(
             gene_id = col_character(),
             count = col_double()))
})

# Create count matrix and sample information for edgeR
rna_genes <- rna_count_list[[1]] %>% pull(gene_id)
rna_count_matrix <- do.call(cbind, lapply(rna_count_list, function(df) df %>% pull(count)))
rownames(rna_count_matrix) <- rna_genes # row names as gene names
colnames(rna_count_matrix) <- rna_sample_info$samples

# output <- data.frame(rna_count_matrix)
# output$gene_id <- rna_genes
# prlgs <- c("RPL1A","RPL1B","RPL2A","RPL2B","RPL12A","RPL12B","RPL18A","RPL18B","RPL19A","RPL19B",
#            "RPL20A","RPL20B","RPL23A","RPL23B","RPL35A","RPL35B","RPL40A","RPL40B","RPL41A","RPL41B",
#            "RPL42A","RPL42B","RPL43A","RPL43B","RPS4A","RPS4B","RPS6A","RPS6B","RPS8A","RPS8B","RPS11A",
#            # "RPS11B","RPS16A","RPS16B","RPS18A","RPS18B","RPS23A","RPS23B","RPS24A","RPS24B","RPS30A",
#            "RPS30B")
# paralogs_rna_counts <- output[output$gene_id %in% prlgs,]
# paralogs_rna_counts
# write_csv(paralogs_rna_counts, 'paralogs_rna_counts.csv')

# Create DGEList object
rna_y <- DGEList(counts = rna_count_matrix, group = rna_sample_info$condition)

# Filter to remove lowly expressed genes
rna_keep <- filterByExpr(rna_y)
rna_y <- rna_y[rna_keep,]

# TMM normalization
rna_y <- calcNormFactors(rna_y, method="TMM")

# Calculate log2CPM with prior.count = 0.5 to avoid log of zero
rna_log2cpm <- cpm(rna_y, log = TRUE, prior.count = 0.5)

# Proceed with voom transformation
rna_v <- voom(rna_y, design = model.matrix(~0 + rna_sample_info$condition), plot = TRUE)

# Design matrix
rna_design <- model.matrix(~0 + rna_sample_info$condition)
colnames(rna_design) <- levels(factor(rna_sample_info$condition))

# Differential expression analysis with Limma
rna_vfit <- lmFit(rna_v, rna_design)

contrasts <- makeContrasts(StressVsWT = stress - wt, levels = rna_design)

# Fit model with contrasts
rna_fit <- contrasts.fit(rna_vfit, contrasts)

# Apply empirical Bayes moderation
rna_fit <- eBayes(rna_fit)

# Get results
rna_results <- topTable(rna_fit, sort.by = "none", n = Inf)

# Check the results
head(rna_results)

# Calculate the standard deviation of log2FC for thresholding
rna_std_log2FC <- sd(rna_results$logFC)

# Apply filtering based on adjusted p-value and |log2FC| > 1 SD(log2FC)
rna_sig_genes <- rna_results[with(rna_results, adj.P.Val < 0.05 & abs(logFC) > rna_std_log2FC),]
count(rna_sig_genes)

rna_sig_genes_up <- rna_sig_genes[with(rna_sig_genes, logFC > 0),]
rna_sig_genes_down <- rna_sig_genes[with(rna_sig_genes, logFC < 0),]
count(rna_sig_genes_up)
count(rna_sig_genes_down)

rna_sig_genes_up$genes <- rownames(rna_sig_genes_up)
rna_sig_genes_down$genes <- rownames(rna_sig_genes_down)

# write_csv(rna_sig_genes_up, "rna_sig_genes_up.csv")
# write_csv(rna_sig_genes_down, "rna_sig_genes_down.csv")




# RNA-seq samples
RNA_stress_rep1 <- rna_log2cpm[, "RNA_stress_rep1"]
RNA_stress_rep2 <- rna_log2cpm[, "RNA_stress_rep2"]
RNA_normal_rep1 <- rna_log2cpm[, "RNA_normal_rep1"]
RNA_normal_rep2 <- rna_log2cpm[, "RNA_normal_rep2"]

# Calculate the density for each sample
density_stress_rep1 <- density(RNA_stress_rep1)
density_stress_rep2 <- density(RNA_stress_rep2)
density_normal_rep1 <- density(RNA_normal_rep1)
density_normal_rep2 <- density(RNA_normal_rep2)

# Plot the first density to set up the plot
plot(density_stress_rep1, main="Density Plot for Samples", xlab="Log-cpm", ylab="Density", 
     xlim=range(c(density_stress_rep1$x, 
                  density_stress_rep2$x, 
                  density_normal_rep1$x, 
                  density_normal_rep2$x)), 
     ylim=c(0, max(c(density_stress_rep1$y, 
                     density_stress_rep2$y, 
                     density_normal_rep1$y, 
                     density_normal_rep2$y))), 
     col="red")

# Add the other densities to the same plot
lines(density_stress_rep2, col="cornflowerblue")
lines(density_normal_rep1, col="darkgreen")
lines(density_normal_rep2, col="orange")

# Add a legend
legend("topright", legend=c("RNA_S1", "RNA_S2", "RNA_N1", "RNA_N2"), 
       col=c("red", "blue", "green", "orange"), lty=1)







######################################################
# Ribo-Seq Samples
######################################################


# List directories containing RNA count files
ribo_subdirs <- list.files(pattern = "_ribo$", full.names = TRUE)
ribo_count_files <- c()


# Collect Ribo. Profiling count files



# ribo_count_files <- c("./SRR6761667_ribo/quality_threshold_1_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761667_full.txt",
#                       "./SRR6761668_ribo/quality_threshold_1_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761668_full.txt",
#                       "./SRR6761669_ribo/quality_threshold_1_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761669_full.txt",
#                       "./SRR6761670_ribo/quality_threshold_1_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761670_full.txt")


ribo_count_files <- c("./SRR6761667_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761667_full.txt",
                       "./SRR6761668_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761668_full.txt",
                       "./SRR6761669_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761669_full.txt",
                       "./SRR6761670_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761670_full.txt")

# ribo_count_files <- c("./SRR6761667_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761667_subset.txt",
#                       "./SRR6761668_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761668_subset.txt",
#                       "./SRR6761669_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761669_subset.txt",
#                       "./SRR6761670_ribo/yes_stranded_Blevins_Tavella_gtf_counts_SRR6761670_subset.txt")

# Adjust subdirectories format if necessary
ribo_subdirs <- gsub("./", "", ribo_subdirs)

# Prepare sample information
ribo_sample_info <- data.frame(
  sra = gsub("_ribo", "", ribo_subdirs),
  samples = c('RP_stress_rep1', 'RP_stress_rep2', 
              'RP_normal_rep1', 'RP_normal_rep2'),
  condition = c('stress', 'stress', 
                'wt', 'wt'))

# Define column names for the count data
col_names <- c("gene_id", "count")

# Read in the count data
ribo_count_list <- lapply(ribo_count_files, function(f) {
  read_tsv(f, 
           col_names = col_names, 
           col_types = cols(
             gene_id = col_character(),
             count = col_double()))
})

# Create count matrix and sample information for edgeR
ribo_genes <- ribo_count_list[[1]] %>% pull(gene_id)
ribo_count_matrix <- do.call(cbind, lapply(ribo_count_list, function(df) df %>% pull(count)))
rownames(ribo_count_matrix) <- ribo_genes # row names as gene names
colnames(ribo_count_matrix) <- ribo_sample_info$samples

# output <- data.frame(ribo_count_matrix)
# output$gene_id <- ribo_genes
# paralogs_ribo_counts <- output[output$gene_id %in% prlgs,]
# paralogs_ribo_counts
# write_csv(paralogs_ribo_counts, 'paralogs_ribo_counts.csv')


# Create DGEList object
ribo_y <- DGEList(counts = ribo_count_matrix, group = ribo_sample_info$condition)

# Filter to remove lowly expressed genes
ribo_keep <- filterByExpr(ribo_y)
ribo_y <- ribo_y[ribo_keep,]

# TMM normalization
ribo_y <- calcNormFactors(ribo_y, method = "TMM")

# Calculate log2CPM with prior.count = 0.5 to avoid log of zero
ribo_log2cpm <- cpm(ribo_y, log = TRUE, prior.count = 0.5)

# Proceed with voom transformation
ribo_v <- voom(ribo_y, design = model.matrix(~0 + ribo_sample_info$condition), plot = TRUE)


# Design matrix
ribo_design <- model.matrix(~0 + ribo_sample_info$condition)
colnames(ribo_design) <- levels(factor(ribo_sample_info$condition))



# Differential expression analysis with Limma
ribo_vfit <- lmFit(ribo_v, ribo_design)


contrasts <- makeContrasts(StressVsWT = stress - wt, levels = ribo_design)

# Fit model with contrasts
ribo_fit <- contrasts.fit(ribo_vfit, contrasts)

# Apply empirical Bayes moderation
ribo_fit <- eBayes(ribo_fit)

# Get results
ribo_results <- topTable(ribo_fit, sort.by = "none", n = Inf)
ribo_results <- ribo_results[ribo_results$logFC > -5,]
# Check the results
# head(ribo_results)


# Calculate the standard deviation of log2FC for thresholding
ribo_std_log2FC <- sd(ribo_results$logFC)

# Apply filtering based on adjusted p-value and |log2FC| > 1 SD(log2FC)
ribo_sig_genes <- ribo_results[with(ribo_results, adj.P.Val < 0.05 & abs(logFC) > ribo_std_log2FC),]
count(ribo_sig_genes)

ribo_sig_genes_up <- ribo_sig_genes[with(ribo_sig_genes, logFC > 0),]
ribo_sig_genes_down <- ribo_sig_genes[with(ribo_sig_genes, logFC < 0),]
count(ribo_sig_genes_up)
count(ribo_sig_genes_down)

ribo_sig_genes_up$genes <- rownames(ribo_sig_genes_up)
ribo_sig_genes_down$genes <- rownames(ribo_sig_genes_down)


# write_csv(ribo_sig_genes_up, "ribo_sig_genes_up.csv")
# write_csv(ribo_sig_genes_down, "ribo_sig_genes_down.csv")







# Figure 4C

# Merge the results based on gene identifiers
rna_results$gene_id <- rownames(rna_results)
ribo_results$gene_id <- rownames(ribo_results)

merged_results <- merge(rna_results, ribo_results, by="gene_id", suffixes=c("_rna", "_ribo"))


# merged_results$category <- with(merged_results, ifelse(
#   adj.P.Val_rna <  0.05  &  adj.P.Val_ribo <  0.05 & abs(logFC_rna) > sd(logFC_rna) & abs(logFC_ribo) > sd(logFC_ribo) & logFC_rna * logFC_ribo > 0, 'homodirectional',
#   ifelse(
#     (adj.P.Val_rna <  0.05 & abs(logFC_rna) > sd(logFC_rna) & abs(logFC_ribo) > sd(logFC_ribo) & logFC_rna * logFC_ribo < 0) , 'opposite_changes',
#     ifelse(
#       adj.P.Val_ribo <  0.05 & adj.P.Val_rna >= 0.05 & abs(logFC_ribo) > sd(logFC_ribo), 'translatome',
#       ifelse(
#         adj.P.Val_rna <  0.05 & adj.P.Val_ribo >=  0.05 & abs(logFC_rna) > sd(logFC_rna), 'transcriptome',
#         'no_change'
#       )
#     )
#   )
# ))


merged_results$category <- with(merged_results, ifelse(
  abs(logFC_rna) > sd(logFC_rna) & abs(logFC_ribo) > sd(logFC_ribo) & logFC_rna * logFC_ribo > 0, 'homodirectional',
  ifelse(
    (abs(logFC_rna) > sd(logFC_rna) & abs(logFC_ribo) > sd(logFC_ribo) & logFC_rna * logFC_ribo < 0) , 'opposite_changes',
    ifelse(
      abs(logFC_ribo) > sd(logFC_ribo), 'translatome',
      ifelse(
        abs(logFC_rna) > sd(logFC_rna), 'transcriptome',
        'no_change'
      )
    )
  )
))



count_same_dir <- nrow(merged_results[merged_results$category == 'homodirectional',])
count_same_dir

translatome_count_up <- nrow(merged_results[merged_results$category == 'translatome' & merged_results$logFC_ribo > 0, ])
translatome_count_up
translatome_count_down <- nrow(merged_results[merged_results$category == 'translatome' & merged_results$logFC_ribo < 0, ])
translatome_count_down
transcriptome_count_up <- nrow(merged_results[merged_results$category == 'transcriptome' & merged_results$logFC_rna > 0, ])
transcriptome_count_up
transcriptome_count_down <- nrow(merged_results[merged_results$category == 'transcriptome' & merged_results$logFC_rna < 0, ])
transcriptome_count_down

count_translatome_up_transcriptome_down <- nrow(merged_results[merged_results$category == 'opposite_changes' & merged_results$logFC_ribo > 0, ])
count_translatome_up_transcriptome_down

count_translatome_down_transcriptome_up <- nrow(merged_results[merged_results$category == 'opposite_changes' & merged_results$logFC_ribo < 0, ])
count_translatome_down_transcriptome_up

count_translatome_up_transcriptome_up <- nrow(merged_results[merged_results$category == 'homodirectional' & merged_results$logFC_ribo > 0, ])
count_translatome_up_transcriptome_up

count_translatome_down_transcriptome_down <- nrow(merged_results[merged_results$category == 'homodirectional' & merged_results$logFC_ribo < 0, ])
count_translatome_down_transcriptome_down


count(merged_results[merged_results$category == "translatome",])
count(merged_results[merged_results$category == "transcriptome",])

plot(merged_results$logFC_rna, merged_results$logFC_ribo, pch=20, col='white',
     xlab="log2 FC (transcriptome)", ylab="log2 FC (translatome)", 
     main="Correlation between log2 FC gene expression values",
     xlim=c(-4.75,4.75),
     ylim=c(-4.75,4.75))

# Add grid
abline(h = seq(-6, 6, by = 1), col = "lightgray", lty = "dotted")
abline(v = seq(-4, 4, by = 1), col = "lightgray", lty = "dotted")


# Add color
points(merged_results$logFC_rna[merged_results$category == 'no_change'], 
       merged_results$logFC_ribo[merged_results$category == 'no_change'], 
       pch=20, col='grey')

points(merged_results$logFC_rna[merged_results$category == 'homodirectional'], 
       merged_results$logFC_ribo[merged_results$category == 'homodirectional'], 
       pch=20, col='cornflowerblue')

points(merged_results$logFC_rna[merged_results$category == 'opposite_changes'], 
       merged_results$logFC_ribo[merged_results$category == 'opposite_changes'], 
       pch=20, col='red')

points(merged_results$logFC_rna[merged_results$category == 'translatome'], 
       merged_results$logFC_ribo[merged_results$category == 'translatome'], 
       pch=20, col='gold')

points(merged_results$logFC_rna[merged_results$category == 'transcriptome'], 
       merged_results$logFC_ribo[merged_results$category == 'transcriptome'], 
       pch=20, col='darkgreen')


# Add a legend
legend("topright", legend=c("homodirectional", "opposite changes", "translatome", "transcriptome", "no change"), 
       col=c('cornflowerblue', 'red', 'gold', 'darkgreen', 'grey'), pch=20)







# Fig 4A: Distribution of gene expression fold change (FC) values 

ribo_dens <- density(merged_results$logFC_ribo)
rna_dens <- density(merged_results$logFC_rna)

# Define colors with transparency using the alpha argument
col_ribo <- rgb(1, 0.5, 0.5, alpha=0.5) # Pink with transparency
col_rna <- rgb(173/255, 216/255, 230/255, alpha=0.5) # Light blue with transparency


# Plot the first density
plot(ribo_dens, xlim=c(-5,5), main="Distribution of gene expression fold change (FC) values", xlab="logFC", ylab="Density")
# Add the second density
lines(rna_dens)

# Use polygon() to fill the area under the curve with the first color
polygon(ribo_dens$x, ribo_dens$y, col=col_ribo, border=NA)

# Use polygon() to fill the area under the second curve with the second color
polygon(rna_dens$x, rna_dens$y, col=col_rna, border=NA)

# Add a legend with semi-transparent colors
legend("topright", legend=c("log2 FC(RP)", "log2 FC(RNA)"), 
       fill=c(col_ribo, col_rna), border=NA)









# Counts Table
setwd("/Users/scampione/Projects/Buck_Institute/Ribo_Profiling/PMID31358845/all_count_files")

# count_files <- list.files(pattern = ".txt$", full.names = TRUE)

count_files <- c("./SRR6761663_gene_counts.txt", 
                 "./SRR6761664_gene_counts.txt",
                 "./SRR6761665_gene_counts.txt",
                 "./SRR6761666_gene_counts.txt",
                 "./quality_threshold_0_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761667_full.txt",
                 "./quality_threshold_0_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761668_full.txt",
                 "./quality_threshold_0_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761669_full.txt",
                 "./quality_threshold_0_yes_stranded_Blevins_Tavella_gtf_counts_SRR6761670_full.txt")


sample_info <- data.frame(
  sra = c('SRR6761663', 'SRR6761664', 'SRR6761665', 'SRR6761666',
          'SRR6761667', 'SRR6761668', 'SRR6761669', 'SRR6761670'),
  samples = c('RNA_stress_rep1', 'RNA_stress_rep2', 
              'RNA_normal_rep1', 'RNA_normal_rep2',
              'RP_stress_rep1', 'RP_stress_rep2', 
              'RP_normal_rep1', 'RP_normal_rep2'),
  condition = c('stress', 'stress', 
                'wt', 'wt',
                'stress', 'stress', 
                'wt', 'wt'),
  type = c('rna','rna','rna','rna',
           'ribo','ribo','ribo','ribo'))

col_names <- c("gene_id", "count")

  

count_list <- lapply(count_files, function(f) {
  read_tsv(f, 
           col_names = col_names, 
           col_types = cols(
             gene_id = col_character(),
             count = col_double()))
})

# Create count matrix and sample information for edgeR
genes <- count_list[[1]] %>% pull(gene_id)
genes
count_matrix <- do.call(cbind, lapply(count_list, function(df) df %>% pull(count)))
rownames(count_matrix) <- genes # row names as gene names
colnames(count_matrix) <- sample_info$samples

# View(count_matrix)

# Create DGEList object
y <- DGEList(counts = count_matrix, group = sample_info$condition)
cpm_y <- cpm(y)

colnames(cpm_y)

plot(cpm_y[,'RNA_normal_rep1'], cpm_y[,'RP_normal_rep1'], 
     log="xy", 
     xlim=c(1,10000), ylim=c(1,10000),
     xlab="Transcriptome (CPM)", ylab="Translation (CPM)",
     main = "Normal")
cor(cpm_y[,'RNA_normal_rep1'], cpm_y[,'RP_normal_rep1'])


plot(cpm_y[,'RNA_stress_rep1'], cpm_y[,'RP_stress_rep1'], 
     log="xy", 
     xlim=c(1,10000), ylim=c(1,10000),
     xlab="Transcriptome (CPM)", ylab="Translation (CPM)",
     main = "Stress")
cor(cpm_y[,'RNA_stress_rep1'], cpm_y[,'RP_stress_rep1'])


plot(cpm_y[,'RNA_normal_rep1'], cpm_y[,'RNA_normal_rep2'], 
     log="xy", 
     xlim=c(1,10000), ylim=c(1,10000),
     xlab="Replicate 1 (CPM)", ylab="Replicate 2 (CPM)",
     main="Transcriptome")
cor(cpm_y[,'RNA_normal_rep1'], cpm_y[,'RNA_normal_rep2'])

plot(cpm_y[,'RP_normal_rep1'], cpm_y[,'RP_normal_rep2'], 
     log="xy", 
     xlim=c(1,10000), ylim=c(1,10000),
     xlab="Replicate 1 (CPM)", ylab="Replicate 2 (CPM)",
     main="Translatome")
cor(cpm_y[,'RP_normal_rep1'], cpm_y[,'RP_normal_rep2'])






######################################################
# 
######################################################

library(biomaRt)

# Connect to the Ensembl yeast database
yeast <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")


# Define the list of genes of interest
prlgs <- c("RPL1A","RPL1B","RPL2A","RPL2B","RPL12A","RPL12B","RPL18A","RPL18B","RPL19A","RPL19B",
           "RPL20A","RPL20B","RPL23A","RPL23B","RPL35A","RPL35B","RPL40A","RPL40B","RPL41A","RPL41B",
           "RPL42A","RPL42B","RPL43A","RPL43B","RPS4A","RPS4B","RPS6A","RPS6B","RPS8A","RPS8B","RPS11A",
           "RPS11B","RPS16A","RPS16B","RPS18A","RPS18B","RPS23A","RPS23B","RPS24A","RPS24B","RPS30A",
           "RPS30B")


gene_lengths <- getBM(attributes = c('external_gene_name', 'start_position', 'end_position', 'chromosome_name'),
                      filters = 'external_gene_name', 
                      values = prlgs, 
                      mart = yeast)

gene_lengths$gene_length <- abs(gene_lengths$end_position - gene_lengths$start_position) + 1




# Filter count matrix for specified genes
filtered_counts <- count_matrix[rownames(count_matrix) %in% prlgs,]

# Select only wild type samples
wt_samples <- sample_info$condition == "wt"
filtered_counts_wt <- filtered_counts[, wt_samples]

# Create a DGEList object
dge <- DGEList(counts = filtered_counts_wt)

# Calculate CPM
cpm_values <- cpm(dge)

# Calculate lengths in kilobases for TPM calculation
gene_lengths_kb <- gene_lengths$gene_length / 1000

# Calculate TPM
tpm_values <- sweep(cpm_values, 1, gene_lengths_kb, '/') * 1e6 / rowSums(sweep(cpm_values, 1, gene_lengths_kb, '/'))

write.csv(tpm_values, 'tpm_values_wt_rna_ribo.csv')




# Publication Data

pub_counts <- read.delim('/Users/scampione/Downloads/Blevins_Tavella_etal_tableofcounts (3).txt', sep='\t')
rownames(pub_counts) <- trimws(pub_counts$ID)
pub_genes <- trimws(pub_counts$ID)


# Filter count matrix for specified genes
pub_filtered_counts <- pub_counts[rownames(pub_counts) %in% prlgs,]

pub_filtered_counts_wt <- pub_filtered_counts[,c('RFcontrolRep1', 'RFcontrolRep2', 'RNAcontrolRep1', 'RNAcontrolRep2')]

# output <- pub_filtered_counts[,c('ID', 'RFcontrolRep1', 'RFcontrolRep2', 'RNAcontrolRep1', 'RNAcontrolRep2')]
# colnames(output)
# colnames(output) <- c('gene_id', 'RP_normal_rep1', 'RP_normal_rep2', "RNA_normal_rep1", "RNA_normal_rep2")
# write_csv(output, 'pub_ribo_rna_counts.csv')

# Create a DGEList object
dge <- DGEList(counts = pub_filtered_counts_wt)

# Calculate CPM
cpm_values <- cpm(dge)

# Calculate lengths in kilobases for TPM calculation
gene_lengths_kb <- gene_lengths$gene_length / 1000

# Calculate TPM
tpm_values <- sweep(cpm_values, 1, gene_lengths_kb, '/') * 1e6 / rowSums(sweep(cpm_values, 1, gene_lengths_kb, '/'))

write.csv(tpm_values, 'pub_tpm_values_wt_rna_ribo.csv')


################################
################################
count_matrix
count_matrix_same_as_pub <- count_matrix[rownames(count_matrix) %in% pub_genes,]
count_matrix_to_normalize <- count_matrix_same_as_pub[,c("RNA_normal_rep1", "RNA_normal_rep2", "RP_normal_rep1", "RP_normal_rep2")]

count_df <- data.frame(count_matrix_to_normalize)
count_df$gene_id <- rownames(count_df)


gene_lengths <- getBM(attributes = c('external_gene_name', 'start_position', 'end_position', 'chromosome_name'),
                      filters = 'external_gene_name', 
                      values = count_df$gene_id, 
                      mart = yeast)

gene_lengths$gene_length <- abs(gene_lengths$end_position - gene_lengths$start_position) + 1

length(count_df$gene_id) <- length(gene_lengths$gene_length)

colSums(count_df[,c('RNA_normal_rep1', 'RNA_normal_rep2')])
