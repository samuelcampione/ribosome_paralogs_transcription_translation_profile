library(tidyverse)
library(edgeR)
library(biomaRt)
library(knitr)



##########################################
# Preparing Yeast Ribosome Paralog Gene Names
##########################################
prlgs <- c("RPL1A","RPL1B","RPL2A","RPL2B","RPL12A","RPL12B","RPL18A","RPL18B","RPL19A","RPL19B",
           "RPL20A","RPL20B","RPL23A","RPL23B","RPL35A","RPL35B","RPL40A","RPL40B","RPL41A","RPL41B",
           "RPL42A","RPL42B","RPL43A","RPL43B","RPS4A","RPS4B","RPS6A","RPS6B","RPS8A","RPS8B","RPS11A",
           "RPS11B","RPS16A","RPS16B","RPS18A","RPS18B","RPS23A","RPS23B","RPS24A","RPS24B","RPS30A",
           "RPS30B")

ensembl_yeast <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
# listAttributes(ensembl_yeast)

gene_info <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                   filters = 'external_gene_name',
                   values = prlgs,
                   mart = ensembl_yeast)

standard_prlgs <- gene_info$ensembl_gene_id

id_to_name <- setNames(gene_info$external_gene_name, gene_info$ensembl_gene_id)







##########################################
# Basic Format to get TMM Normalized Values
##########################################
y <- DGEList(count_matrix)
y <- calcNormFactors(y)
cpm(y, log = FALSE)





# GeTMM Normalization


##########################################
# Blevins et al. 2019: GeTMM
##########################################
pub_col_names <- c('ID',
                   'RP_normal_rep1', 'RP_normal_rep2',
                   'RP_stress_rep1', 'RP_stress_rep2', 
                   'RNA_normal_rep1', 'RNA_normal_rep2',
                   'RNA_stress_rep1', 'RNA_stress_rep2') 

pub_matrix <- read_tsv('/Users/scampione/Downloads/Blevins_Tavella_etal_tableofcounts (3).txt', 
                       trim_ws = TRUE,
                       col_names = pub_col_names, 
                       col_types = cols(
                         ID = col_character(),
                         RP_normal_rep1 = col_double(), RP_normal_rep2 = col_double(),
                         RP_stress_rep1 = col_double(), RP_stress_rep2 = col_double(),
                         RNA_normal_rep1 = col_double(), RNA_normal_rep2 = col_double(),
                         RNA_stress_rep1 = col_double(), RNA_stress_rep2 = col_double()),
                       skip=1)


pub_sample_info <- data.frame(
  sra = c('SRR6761669', 'SRR6761670', 'SRR6761667', 'SRR6761668',
          'SRR6761665', 'SRR6761666', 'SRR6761663', 'SRR6761664'),
  samples = c('RP_normal_rep1', 'RP_normal_rep2',
              'RP_stress_rep1', 'RP_stress_rep2', 
              'RNA_normal_rep1', 'RNA_normal_rep2',
              'RNA_stress_rep1', 'RNA_stress_rep2'),
  condition = c('wt', 'wt',
                'stress', 'stress', 
                'wt', 'wt',
                'stress', 'stress'),
  type = c('ribo','ribo','ribo','ribo',
           'rna','rna','rna','rna'))


pub_genes <- pub_matrix$ID

count_matrix <- pub_matrix %>% 
  dplyr::select(-ID) %>% 
  as.data.frame()

rownames(count_matrix) <- pub_genes

gene_info_2 <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'transcript_length'),
                     filters = 'external_gene_name',
                     values = rownames(count_matrix),
                     mart = ensembl_yeast)


df_len <- data.frame(length_kb = gene_info_2$transcript_length / 1000)  # Convert to kb if necessary
rownames(df_len) <- gene_info_2$external_gene_name


# Merge data and handle NA
count_matrix <- merge(count_matrix, df_len, by = "row.names", all.x = TRUE)
count_matrix <- count_matrix[!is.na(count_matrix$length_kb), ]

# Calculate RPK
samples <- colnames(count_matrix)[2:(ncol(count_matrix)-1)]  # Adjust to skip 'Gene' and 'length_kb'

# Calculate RPK for each column of count_matrix corresponding to samples
for (sample in samples) {
  count_matrix[[sample]] <- (count_matrix[[sample]] / count_matrix$length_kb)
}

# TPM Normalization
sub_count_matrix <- count_matrix[c("RNA_stress_rep1", "RNA_stress_rep2", "RNA_normal_rep1", "RNA_normal_rep2", "RP_stress_rep1", "RP_stress_rep2", "RP_normal_rep1", "RP_normal_rep2")]
rownames(sub_count_matrix) <- count_matrix$Row.names

sum_RPK <- colSums(sub_count_matrix)

# Calculate TPM for each sample
for (sample in samples) {
  sub_count_matrix[[sample]] <- (sub_count_matrix[[sample]] / sum_RPK[sample]) * 1e6
}

tpm_matrix <- sub_count_matrix

# Subset for the specific paralog genes
paralogs_TPM_values <- tpm_matrix[rownames(tpm_matrix) %in% prlgs,]

# Sorting paralogs to ensure they match the order in `prlgs` list
sorted_paralogs_TPM <- paralogs_TPM_values[match(prlgs, rownames(paralogs_TPM_values)), ]

sorted_paralogs_TPM <- sorted_paralogs_TPM[c("RNA_normal_rep1", "RNA_normal_rep2", "RP_normal_rep1", "RP_normal_rep2")]

# write.csv(sorted_paralogs_TPM, "Blevins_paralogs_TPM.csv", row.names = TRUE, )


# TMM Normalization
count_matrix$length_kb <- NULL  # remove length_kb column
rownames(count_matrix) <- count_matrix$Row.names
count_matrix$Row.names <- NULL


y <- DGEList(counts = count_matrix)
y <- calcNormFactors(y)
output <- cpm(y, log = FALSE)

paralogs_GeTMM_values <- output[rownames(output) %in% prlgs,]

Blevins_sorted_paralogs_GeTMM <- paralogs_GeTMM_values[match(prlgs, rownames(paralogs_GeTMM_values)), ]
Blevins_sorted_paralogs_GeTMM <- Blevins_sorted_paralogs_GeTMM[,c("RNA_normal_rep1", "RNA_normal_rep2", "RP_normal_rep1", "RP_normal_rep2")]


# Average across replicates
RNA_normal_mean <- rowMeans(Blevins_sorted_paralogs_GeTMM[,c("RNA_normal_rep1", "RNA_normal_rep2")])
Ribo_normal_mean <- rowMeans(Blevins_sorted_paralogs_GeTMM[,c("RP_normal_rep1", "RP_normal_rep2")])
Blevins_avg_paralogs_GeTMM <- data.frame(RNA_normal_mean, Ribo_normal_mean)

View(Blevins_avg_paralogs_GeTMM)

# write.csv(Blevins_avg_paralogs_GeTMM, "Blevins_avg_paralogs_GeTMM.csv", row.names = TRUE, )








##########################################
# Replicated Blevins et al. 2019 (mode: intersection-strict, min quality: 10-—new default) GeTMM
##########################################
setwd('/Users/scampione/Projects/Buck_Institute/Ribo_Profiling/PMID31358845/all_count_files_strict')
count_files <- c("./SRR6761663_gene_counts.txt",
                 "./SRR6761664_gene_counts.txt",
                 "./SRR6761665_gene_counts.txt",
                 "./SRR6761666_gene_counts.txt",
                 "./yes_stranded_Blevins_Tavella_gtf_counts_SRR6761667_full.txt",
                 "./yes_stranded_Blevins_Tavella_gtf_counts_SRR6761668_full.txt",
                 "./yes_stranded_Blevins_Tavella_gtf_counts_SRR6761669_full.txt",
                 "./yes_stranded_Blevins_Tavella_gtf_counts_SRR6761670_full.txt")



sample_info <- data.frame(
  sra = c('SRR6761663', 'SRR6761664', 'SRR6761665', 'SRR6761666',
          'SRR6761667', 'SRR6761668', 'SRR6761669', 'SRR6761670'),
  samples = c('RNA_stress_rep1', 'RNA_stress_rep2', 
              'RNA_normal_rep1', 'RNA_normal_rep2',
              'RP_stress_rep1', 'RP_stress_rep2', 
              'RP_normal_rep1', 'RP_normal_rep2'),
  condition = c('stress', 'stress',  'wt', 'wt',
                'stress', 'stress',  'wt', 'wt'),
  type = c('rna','rna','rna','rna',
           'ribo','ribo','ribo','ribo')
)


col_names <- c("gene_id", "count")

count_list <- lapply(count_files, function(f) {
  read_tsv(f, 
           col_names = col_names, 
           col_types = cols(
             gene_id = col_character(),
             count = col_double()))
})

genes <- count_list[[1]] %>% pull(gene_id)
count_matrix <- do.call(cbind, lapply(count_list, function(df) df %>% pull(count)))
rownames(count_matrix) <- genes 
colnames(count_matrix) <- sample_info$samples

n <- nrow(count_matrix) 
count_matrix <- head(count_matrix, n - 5)


gene_info_2 <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'transcript_length'),
                     filters = 'external_gene_name',
                     values = rownames(count_matrix),
                     mart = ensembl_yeast)

df_len <- data.frame(length_kb = gene_info_2$transcript_length / 1000)  # Convert to kb if necessary
rownames(df_len) <- gene_info_2$external_gene_name

# Merge data and handle NA
count_matrix <- merge(count_matrix, df_len, by = "row.names", all.x = TRUE)
count_matrix <- count_matrix[!is.na(count_matrix$length_kb), ]
colnames(count_matrix)

samples <- colnames(count_matrix)[2:(ncol(count_matrix)-1)]  # Adjust to skip 'Gene' and 'length_kb'

# Calculate RPK for each column of count_data corresponding to samples
for (sample in samples) {
  count_matrix[[sample]] <- (count_matrix[[sample]] / count_matrix$length_kb)
}


# TPM Normalization
sub_count_matrix <- count_matrix[c("RNA_stress_rep1", "RNA_stress_rep2", "RNA_normal_rep1", "RNA_normal_rep2", "RP_stress_rep1", "RP_stress_rep2", "RP_normal_rep1", "RP_normal_rep2")]
rownames(sub_count_matrix) <- count_matrix$Row.names

sum_RPK <- colSums(sub_count_matrix)

# Calculate TPM for each sample
for (sample in samples) {
  sub_count_matrix[[sample]] <- (sub_count_matrix[[sample]] / sum_RPK[sample]) * 1e6
}

tpm_matrix <- sub_count_matrix

# Subset for the specific paralog genes
paralogs_TPM_values <- tpm_matrix[rownames(tpm_matrix) %in% prlgs,]

# Sorting paralogs to ensure they match the order in `prlgs` list
sorted_paralogs_TPM <- paralogs_TPM_values[match(prlgs, rownames(paralogs_TPM_values)), ]

sorted_paralogs_TPM <- sorted_paralogs_TPM[c("RNA_normal_rep1", "RNA_normal_rep2", "RP_normal_rep1", "RP_normal_rep2")]

# write.csv(sorted_paralogs_TPM, "Replication_paralogs_TPM.csv", row.names = TRUE, )



# TMM Normalization
count_matrix$length_kb <- NULL  # remove length_kb column
rownames(count_matrix) <- count_matrix$Row.names
count_matrix$Row.names <- NULL


y <- DGEList(count_matrix)
y <- calcNormFactors(y)
output <- cpm(y, log = FALSE)

paralogs_GeTMM_values <- output[rownames(output) %in% prlgs,]

Replication_sorted_paralogs_GeTMM <- paralogs_GeTMM_values[match(prlgs, rownames(paralogs_GeTMM_values)), ]
Replication_sorted_paralogs_GeTMM <- Replication_sorted_paralogs_GeTMM[,c("RNA_normal_rep1", "RNA_normal_rep2", "RP_normal_rep1", "RP_normal_rep2")]


RNA_normal_mean <- rowMeans(Replication_sorted_paralogs_GeTMM[,c("RNA_normal_rep1", "RNA_normal_rep2")])
Ribo_normal_mean <- rowMeans(Replication_sorted_paralogs_GeTMM[,c("RP_normal_rep1", "RP_normal_rep2")])
Replication_avg_paralogs_GeTMM <- data.frame(RNA_normal_mean, Ribo_normal_mean)
View(Replication_avg_paralogs_GeTMM)

# write.csv(Replication_avg_paralogs_GeTMM, "Replication_avg_paralogs_GeTMM.csv")







##########################################
# Nedialkova et al. 2015: GeTMM
##########################################
setwd('/Users/scampione/Projects/Buck_Institute/Ribosome Paralogs Translaiton Transcription')

rna_count_matrix <- read_tsv("Nedialkova/GSE67387_Scer_mRNA_rawCounts_WT.tsv")
ribo_count_matrix <- read_tsv("Nedialkova/GSE67387_Scer_ribo_rawCounts_WT.tsv")

rna_count_matrix <- rna_count_matrix %>%
  rename_with(~ paste0("RNA_", .), -Gene_ID)

ribo_count_matrix <- ribo_count_matrix %>%
  rename_with(~ paste0("Ribo_", .), -Gene_ID)

merged_data <- inner_join(rna_count_matrix, ribo_count_matrix, by = "Gene_ID")
count_matrix <- data.frame(drop_na(merged_data))
rownames(count_matrix) <- count_matrix$Gene_ID

gene_info_3 <- getBM(attributes = c('ensembl_gene_id', 'transcript_length'),
                     filters = 'ensembl_gene_id',
                     values = count_matrix$Gene_ID,
                     mart = ensembl_yeast)

df_len <- data.frame(length_kb = gene_info_3$transcript_length / 1000)  # Convert to kb if necessary
rownames(df_len) <- gene_info_3$ensembl_gene_id

# Merge data and handle NA
count_matrix <- merge(count_matrix, df_len, by = "row.names", all.x = TRUE)
count_matrix <- count_matrix[!is.na(count_matrix$length_kb), ]
colnames(count_matrix)

samples <- colnames(count_matrix)[3:(ncol(count_matrix)-1)]  # Adjust to skip 'Gene' and 'length_kb'

# Calculate RPK for each column of count_data corresponding to samples
for (sample in samples) {
  count_matrix[[sample]] <- (count_matrix[[sample]] / count_matrix$length_kb)
}

count_matrix$length_kb <- NULL  # remove length_kb column
rownames(count_matrix) <- count_matrix$Gene_ID
count_matrix$Row.names <- NULL
count_matrix$Gene_ID <- NULL


y <- DGEList(counts = count_matrix)
y <- calcNormFactors(y)
output <- cpm(y, log = FALSE)


paralogs_GeTMM_values <- output[rownames(output) %in% standard_prlgs,]

rownames(paralogs_GeTMM_values) <- id_to_name[rownames(paralogs_GeTMM_values)]

Nedialkova_sorted_paralogs_GeTMM <- paralogs_GeTMM_values[match(prlgs, rownames(paralogs_GeTMM_values)), ]

# Average across replicates
RNA_normal_mean <- rowMeans(Nedialkova_sorted_paralogs_GeTMM[,c("RNA_WT_mRNA_YPD_rep1", "RNA_WT_mRNA_YPD_rep2", "RNA_WT_mRNA_YPD_rep3")])
Ribo_normal_mean <- rowMeans(Nedialkova_sorted_paralogs_GeTMM[,c("Ribo_WT_ribo_YPD_rep1", "Ribo_WT_ribo_YPD_rep2", "Ribo_WT_ribo_YPD_rep3")])
Nedialkova_avg_paralogs_GeTMM <- data.frame(RNA_normal_mean, Ribo_normal_mean)

View(Nedialkova_avg_paralogs_GeTMM)

# write.csv(Nedialkova_avg_paralogs_GeTMM, "Nedialkova_avg_paralogs_GeTMM.csv", row.names = TRUE)





##########################################
# Zinshteyn, Gilbert 2013 GeTMM
##########################################
setwd('/Users/scampione/Projects/Buck_Institute/Ribosome Paralogs Translaiton Transcription')

count_matrix <- read_csv("Zinshteyn/GSE45366_WT_counts.csv")

count_matrix <- drop_na(count_matrix)

genes <- count_matrix$gene

count_matrix <- count_matrix %>% 
  dplyr::select(-gene) %>% 
  as.data.frame()

rownames(count_matrix) <- genes

gene_info_3 <- getBM(attributes = c('ensembl_gene_id', 'transcript_length'),
                     filters = 'ensembl_gene_id',
                     values = rownames(count_matrix),
                     mart = ensembl_yeast)

df_len <- data.frame(length_kb = gene_info_3$transcript_length / 1000)  # Convert to kb if necessary
rownames(df_len) <- gene_info_3$ensembl_gene_id

# Merge data and handle NA
count_matrix <- merge(count_matrix, df_len, by = "row.names", all.x = TRUE)
count_matrix <- count_matrix[!is.na(count_matrix$length_kb), ]
colnames(count_matrix)


samples <- colnames(count_matrix)[2:(ncol(count_matrix)-1)]  # Adjust to skip 'Gene' and 'length_kb'

# Calculate RPK for each column of count_data corresponding to samples
for (sample in samples) {
  count_matrix[[sample]] <- (count_matrix[[sample]] / count_matrix$length_kb)
}

count_matrix$length_kb <- NULL  # remove length_kb column
rownames(count_matrix) <- count_matrix$Row.names
count_matrix$Row.names <- NULL


y <- DGEList(count_matrix)
y <- calcNormFactors(y)

cpm(y, log = FALSE)["YDL184C",] #RPL41A
cpm(y, log = FALSE)["YDL133C-A",] #RPL41B
output <- cpm(y, log = FALSE)

paralogs_GeTMM_values <- output[rownames(output) %in% standard_prlgs,]

rownames(paralogs_GeTMM_values) <- id_to_name[rownames(paralogs_GeTMM_values)]

Zinshteyn_sorted_paralogs_GeTMM <- paralogs_GeTMM_values[match(prlgs, rownames(paralogs_GeTMM_values)), ]

# Average across replicates
RNA_normal_mean <- rowMeans(Zinshteyn_sorted_paralogs_GeTMM[,c("WT.T1.p", "WT.T2.p")])
Ribo_normal_mean <- rowMeans(Zinshteyn_sorted_paralogs_GeTMM[,c("WT.FP1.p", "WT.FP2.p")])
Zinshteyn_avg_paralogs_GeTMM <- data.frame(RNA_normal_mean, Ribo_normal_mean)

View(Zinshteyn_avg_paralogs_GeTMM)

# write.csv(Zinshteyn_avg_paralogs_GeTMM, "Zinshteyn_avg_paralogs_GeTMM.csv", row.names = TRUE, )






##########################################
# Heyer et al. 2016 GeTMM
##########################################
setwd('/Users/scampione/Projects/Buck_Institute/Ribosome Paralogs Translaiton Transcription')

count_matrix <- read_tsv("Heyer_Moore/GSE76117_GeneCounts_WT.tsv")
count_matrix <- drop_na(count_matrix)
genes <- count_matrix$Gene

count_matrix <- count_matrix %>% 
  dplyr::select(-Gene) %>% 
  as.data.frame()

rownames(count_matrix) <- genes

gene_info_3 <- getBM(attributes = c('ensembl_gene_id', 'transcript_length'),
                     filters = 'ensembl_gene_id',
                     values = rownames(count_matrix),
                     mart = ensembl_yeast)

df_len <- data.frame(length_kb = gene_info_3$transcript_length / 1000)  # Convert to kb if necessary
rownames(df_len) <- gene_info_3$ensembl_gene_id

# Merge data and handle NA
count_matrix <- merge(count_matrix, df_len, by = "row.names", all.x = TRUE)
count_matrix <- count_matrix[!is.na(count_matrix$length_kb), ]
colnames(count_matrix)


samples <- colnames(count_matrix)[2:(ncol(count_matrix)-1)]  # Adjust to skip 'Gene' and 'length_kb'

# Calculate RPK for each column of count_data corresponding to samples
for (sample in samples) {
  count_matrix[[sample]] <- (count_matrix[[sample]] / count_matrix$length_kb)
}

count_matrix$length_kb <- NULL  # remove length_kb column
rownames(count_matrix) <- count_matrix$Row.names
count_matrix$Row.names <- NULL


y <- DGEList(count_matrix)
y <- calcNormFactors(y)

cpm(y, log = FALSE)["YDL184C",] #RPL41A
cpm(y, log = FALSE)["YDL133C-A",] #RPL41B
output <- cpm(y, log = FALSE)

paralogs_GeTMM_values <- output[rownames(output) %in% standard_prlgs,]
rownames(paralogs_GeTMM_values) <- id_to_name[rownames(paralogs_GeTMM_values)]
Heyer_sorted_paralogs_GeTMM <- paralogs_GeTMM_values[match(prlgs, rownames(paralogs_GeTMM_values)), ]

# Average across replicates
RNA_normal_mean <- rowMeans(Heyer_sorted_paralogs_GeTMM[,c("RNASeq1Counts", "RNASeq2Counts")])
Ribo_normal_mean <- rowMeans(Heyer_sorted_paralogs_GeTMM[,c("global1Counts", "global2Counts")])
Heyer_avg_paralogs_GeTMM <- data.frame(RNA_normal_mean, Ribo_normal_mean)

View(Heyer_avg_paralogs_GeTMM)

# write.csv(Heyer_avg_paralogs_GeTMM, "Heyer_avg_paralogs_GeTMM.csv", row.names = TRUE, )







##########################################
# Chou et al. 2017 GeTMM
##########################################
setwd('/Users/scampione/Projects/Buck_Institute/Ribosome Paralogs Translaiton Transcription')

rna_count_matrix <- read_tsv("Chou/GSE100626_Chou_RNA_WT.tsv")
ribo_count_matrix <- read_tsv("Chou/GSE100626_Chou_Ribo_WT.tsv")

rna_count_matrix <- rna_count_matrix %>%
  rename_with(~ paste0("RNA_", .), -Transcript)

ribo_count_matrix <- ribo_count_matrix %>%
  rename_with(~ paste0("Ribo_", .), -Transcript)

merged_data <- inner_join(rna_count_matrix, ribo_count_matrix, by = "Transcript")
count_matrix <- drop_na(merged_data)

genes <- merged_data$Transcript

count_data <- merged_data %>% 
  dplyr::select(-Transcript) %>% 
  as.data.frame()

rownames(count_data) <- genes


gene_info_3 <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'transcript_length'),
                     filters = 'ensembl_gene_id',
                     values = rownames(count_data),
                     mart = ensembl_yeast)

df_len <- data.frame(length_kb = gene_info_3$transcript_length / 1000)  # Convert to kb if necessary
rownames(df_len) <- gene_info_3$ensembl_gene_id



# Merge data and handle NA
count_data <- merge(count_data, df_len, by = "row.names", all.x = TRUE)
count_data <- count_data[!is.na(count_data$length_kb), ]
colnames(count_data)


samples <- colnames(count_data)[2:(ncol(count_data)-1)]  # Adjust to skip 'Gene' and 'length_kb'

# Calculate RPK for each column of count_data corresponding to samples
for (sample in samples) {
  count_data[[sample]] <- (count_data[[sample]] / count_data$length_kb)
}

count_data$length_kb <- NULL  # remove length_kb column
rownames(count_data) <- count_data$Row.names
count_data$Row.names <- NULL


y <- DGEList(count_data)
y <- calcNormFactors(y)
output <- cpm(y, log = FALSE)

paralogs_GeTMM_values <- output[rownames(output) %in% standard_prlgs,]

rownames(paralogs_GeTMM_values) <- id_to_name[rownames(paralogs_GeTMM_values)]

Chou_sorted_paralogs_GeTMM <- paralogs_GeTMM_values[match(prlgs, rownames(paralogs_GeTMM_values)), ]
colnames(Chou_sorted_paralogs_GeTMM)

# Average across replicates
RNA_normal_mean <- rowMeans(Chou_sorted_paralogs_GeTMM[,c("RNA_WT_1", "RNA_WT_2", "RNA_WT_3", "RNA_WT_4", "RNA_WT_5", "RNA_WT_6", "RNA_WT_7", "RNA_WT_8", "RNA_WT_9", "RNA_WT_10")])
Ribo_normal_mean <- rowMeans(Chou_sorted_paralogs_GeTMM[,c("Ribo_WT_1", "Ribo_WT_2", "Ribo_WT_3", "Ribo_WT_4", "Ribo_WT_5", "Ribo_WT_6", "Ribo_WT_7", "Ribo_WT_8", "Ribo_WT_9", "Ribo_WT_10", "Ribo_WT_11", "Ribo_WT_12", "Ribo_WT_13", "Ribo_WT_14")])
Chou_avg_paralogs_GeTMM <- data.frame(RNA_normal_mean, Ribo_normal_mean)

View(Chou_avg_paralogs_GeTMM)

# write.csv(Chou_avg_paralogs_GeTMM, "Chou_avg_paralogs_GeTMM.csv", row.names = TRUE)





