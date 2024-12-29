# Set working directory
setwd("C:/Users/KIIT/Desktop/bio-project/")

# Decompress the tar.gz file containing the dataset 
untar('brca_tcga_pan_can_atlas_2018.tar.gz')

# Define path for reading the files
folder_path = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
patient_data_path = paste(folder_path,"data_clinical_patient.txt", sep = "/")

# Read using read delim
patient_data = read.delim(patient_data_path)

# Skip first 5 rows as it includes descriptions 
patient_data = patient_data[5:dim(patient_data)[1],]

rna_data_path <- paste(folder_path,"data_mrna_seq_v2_rsem.txt", sep = "/")
rnaseq_data <- read.delim(rna_data_path)

cna_data_path <-paste(folder_path, "data_cna.txt", sep = "/")
cna_data <- read.delim(cna_data_path)

# String manipulation to standardize the patient IDs in CNA data to match other datasets
clean_colnames <- gsub("\\.$", "", gsub("\\.", "-", sub("(\\d{2})$", "", colnames(cna_data)[3:length(colnames(cna_data))])))
clean_colnames<- as.data.frame(clean_colnames)
clean_colnames$id <- substr(clean_colnames$clean_colnames, 1, nchar(clean_colnames$clean_colnames) - 1)
clean_colnames$rna_seq_ids <- colnames(cna_data)[3:length(colnames(cna_data))]

# Match patient IDs between clinical and CNA data
pat_id<- as.data.frame(patient_data$X.Patient.Identifier)
colnames(pat_id)<- 'id'

identifier_patient = merge(clean_colnames, pat_id, by = "id", all = TRUE)
identifier_patient = na.omit(identifier_patient)

# Prepare the RNA-Seq data matrix for analysis
assay = round(as.matrix(rnaseq_data[,-c(1,2)])) 
rownames(assay) = rnaseq_data[,1]

# Find the row where the first column has "ERBB2"
ERBB2_row = which(cna_data[, 1] == "ERBB2")

# Extract the ERBB2 counts (patient IDs start from the 3rd column)
ERBB2_data = matrix(cna_data[ERBB2_row, 3:ncol(cna_data)], 
                    ncol = 1, 
                    dimnames = list(colnames(cna_data)[3:ncol(cna_data)],"ERBB2_Count"))
matching_id <- colnames(assay)

# Subset ERBB2 data to include only matching IDs with RNA-Seq data
ERBB2_filtered <- ERBB2_data[rownames(ERBB2_data) %in% colnames(assay), , drop = FALSE]

# Initialize metadata for HER2 status based on ERBB2 counts
metadata = matrix(0, nrow = nrow(ERBB2_filtered), ncol = 1)
colnames(metadata) = "HER2_Status"

# Determine HER2 status based on ERBB2 counts
metadata[, 1] = ifelse(as.numeric(ERBB2_filtered[, "ERBB2_Count"])>0,1,0)

assay <- assay[, colnames(assay) %in% identifier_patient$rna_seq_ids]

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2
if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)

assay[is.na(assay)] = 0  # Impute with zeros the NA
assay[assay<0] = 0

# Remove genes with too many missing values.
# Filter genes with insufficient expression
smallestGroupSize = 3
keep = rowSums(assay >= 10) >= smallestGroupSize
assay = assay[keep,]

# Create a DESeqDataSet object and perform differential expression analysis
dds =  DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ HER2_Status)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# Extract the results and sort by adjusted p-value
res = results(dds)
res[order(res$padj)[1:10],]

# Perform variance stabilization for PCA and visualization
vsd = vst(dds)

par(mfrow = c(1, 2))
plotPCA(vsd, intgroup=c("HER2_Status"))


if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("enrichplot", quietly = TRUE))
  install.packages("enrichplot")




library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# Gene ontology enrichment analysis for over- and under-expressed genes
res_sig = res[res$padj<0.05,]

# Separate into over and under expressed using log2foldchange
DE_over = rownames(res_sig[res_sig$log2FoldChange>0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,])

go_results_over = enrichGO(
  gene          = DE_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Print and plot results
print(head(go_results_over))
dotplot(go_results_over, showCategory=10) + ggtitle("Gene Ontology Enrichment over Expressed")



go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Print and plot results
#print(head(go_results_under))

dotplot(go_results_under, showCategory=10) + ggtitle("Gene Ontology Enrichment Under Expressed")

# KEGG and Reactome pathway enrichment

if (!requireNamespace("pathview", quietly = TRUE))
  BiocManager::install("pathview")

if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
BiocManager::install("ReactomePA", force = TRUE)

library(ReactomePA)
library(pathview)

gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)


gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_under <- bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(kegg_results_over))

dotplot(kegg_results_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")

print(head(kegg_results_under))

dotplot(kegg_results_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")


gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_under <- bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(kegg_results_over))

dotplot(kegg_results_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")

print(head(kegg_results_under))

dotplot(kegg_results_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")



reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)


print(head(reactome_results_over))

dotplot(reactome_results_over, showCategory=10) + ggtitle("Reactome Pathway Enrichment Over Expressed")

print(head(reactome_results_under))

dotplot(reactome_results_under, showCategory=10) + ggtitle("Reactome Pathway Enrichment Under Expressed")


go_results_under_pw = pairwise_termsim(go_results_under)
treeplot(go_results_under_pw)+ ggtitle("GO Enrichment Under Expressed")

kegg_results_under_pw = pairwise_termsim(kegg_results_under)
treeplot(kegg_results_under_pw)+ ggtitle("KEGG Enrichment Under Expressed")

# Generate heatmap for top differentially expressed genes
# Subset the dataset on differentially expressed gene for this assignment. 
top_DE = order(res$padj)

vsd_DE = assay(vsd)[top_DE[1:20],]
# To simplify the visualization get top most different genes. 


# install packages for nicer heatmap than R's base one. 

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")


library(pheatmap)


annotation_colors = list(HER2_Status = c(Her2_amplified = "#1f78b4", not_amplified = "#33a02c"))



annotation_col = data.frame(HER2_Status = as.matrix(metadata[,1]))
rownames(annotation_col) = colnames(vsd)


pheatmap(
  vsd_DE,
  cluster_rows = TRUE,      
  cluster_cols = TRUE,  
  scale = 'row',
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col)

library(xlsx)
DEG_df <- data.frame(res[res$padj<0.05,])
write.xlsx(DEG_df, "DEG.xlsx", rowNames = TRUE)

#Survival Model (using vst values of DE genes)

library(survival)
library(glmnet)

# Map clinical and RNA-Seq identifiers
mapping_table <- identifier_patient  
id_mapping <- mapping_table[, c("id", "rna_seq_ids")]  
colnames(id_mapping) <- c("clinical_id", "rna_id") 
clinical_ids <- patient_data$X.Patient.Identifier
rna_ids <- colnames(vsd_DE)
id_mapping <- id_mapping[id_mapping$clinical_id %in% clinical_ids & id_mapping$rna_id %in% rna_ids, ]  
# Filter mappings to include only valid IDs

# Filter clinical and RNA-Seq data
clinical_data_filtered <- patient_data[patient_data$X.Patient.Identifier %in% id_mapping$clinical_id, ]  
rna_seq_filtered <- vsd_DE[, colnames(vsd_DE) %in% id_mapping$rna_id] 
rna_seq_filtered <- rna_seq_filtered[, match(id_mapping$rna_id, colnames(rna_seq_filtered))]  
# Reorder columns to match RNA IDs in mapping
all(colnames(rna_seq_filtered) == id_mapping$rna_id)  # Verify alignment of IDs

# Prepare survival data
time <- as.numeric(clinical_data_filtered$Overall.Survival..Months.)  # Extract survival times
status <- ifelse(clinical_data_filtered$Overall.Survival.Status == "1:DECEASED", 1, 0)  
# Encode survival status
valid_indices <- time > 0  # Identify valid indices with positive survival times
time <- time[valid_indices] 
status <- status[valid_indices]
rna_seq_filtered <- rna_seq_filtered[, valid_indices]

# Transpose and scale RNA-Seq data
x <- t(rna_seq_filtered)
x <- x[rowSums(is.na(x)) == 0, ]
x <- x[, colSums(is.na(x)) == 0]
rownames(x) <- id_mapping$clinical_id[match(rownames(x), id_mapping$rna_id)]

# Filter for common identifiers
common_ids <- intersect(rownames(x), clinical_data_filtered$X.Patient.Identifier[valid_indices])  
x <- x[common_ids, ]
x <- scale(x)

# Subset survival data to common IDs
time <- time[clinical_data_filtered$X.Patient.Identifier[valid_indices] %in% common_ids]
status <- status[clinical_data_filtered$X.Patient.Identifier[valid_indices] %in% common_ids]
y <- Surv(time, status)  # Create survival object

# Regularized Cox regression
set.seed(123)  # Set seed for reproducibility
fit <- glmnet(as.matrix(x), y, family = "cox", alpha = 0.5)  # Fit elastic net Cox model
summary(fit)  # Print model details
coef_matrix <- coef(fit, s = min(fit$lambda))  # Extract coefficients at minimum lambda
risk_scores <- as.numeric(x %*% coef_matrix)  # Compute risk scores
risk_groups <- ifelse(risk_scores > median(risk_scores), "High Risk", "Low Risk")  # Categorize risk groups

# Kaplan-Meier survival analysis
km_surv <- Surv(time, status)
library(survminer)
km_fit <- survfit(km_surv ~ risk_groups)
ggsurvplot(km_fit, data = data.frame(time, status, risk_groups),
           pval = TRUE, conf.int = TRUE,  
           risk.table = TRUE, legend.title = "Risk Group",
           main = "Kaplan-Meier Survival Curves by Risk Group",
           xlab = "Time (Months)", ylab = "Survival Probability")  
