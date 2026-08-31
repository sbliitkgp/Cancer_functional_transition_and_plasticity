#This script will run only sample wise. Hence plese change the folder path for each sample run

path <- "/data/subhasis_backup/subhasis/pipeline_testing/data/SMC01-T" #Please give your full path to samples folder where the metadata and gene expression data is present
setwd(path)
library(dplyr)
library(stringr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)

files <- list.files()
tpm_data_file <- files[str_detect(files, ".csv$") & str_detect(files, 'TPM_mod')]
#required_string <- substring(tpm_data_file,1,7)
required_string <- basename(path)
celltype_data_file <- files[str_detect(files, ".tsv$") & str_detect(files, 'required')]

#reading TPM file
data <- read.csv(tpm_data_file,row.names=1)

str <- unlist(strsplit(tpm_data_file, split='_TPM_mod.csv', fixed=TRUE))[1]
cancer <- CreateSeuratObject(counts = data, project = str, min.cells = 0,min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT.")

VlnPlot(cancer, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste(required_string,'feature plot.pdf'),dpi=800,width=300,height=190,units='mm')

plot1 <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(paste(required_string,'feature vs feature plot.pdf'),dpi=800,width=300,height=190,units='mm')

counts.df<- as.data.frame(cancer@assays[["RNA"]]$counts)

#tpm_data <- round(data,1)

# taking only those TPM values where gene name matches to the filtered gene names in raw count dataset
row1 <- rownames(data)
row2 <- rownames(counts.df)
tpm_subset <- data[row1 %in% row2,]

# centering of gene expression values
data1 <- tpm_subset[] - rowMeans(tpm_subset[])[row(tpm_subset[])]

transposed_data <- as.data.frame(t(data1))
colnames(transposed_data) <- NULL
rownames(transposed_data)<- NULL

write.table(transposed_data, file=paste(path,paste(required_string,"_gene_expression_parti_input.tsv", sep=""),sep="/"), sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)

genes <- as.data.frame(rownames(data1))
write.table(genes, file="genes.list", sep=" ",row.names=FALSE,col.names=FALSE, quote=FALSE)