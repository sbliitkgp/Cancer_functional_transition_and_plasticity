library(gplots)
library(glue)
library(stringr)
path <- getwd()
setwd(path)

# Run the Python script and capture its output
output <- system2("python", args = c("../../pipeline_testing/codes/code5.3_distance_cutoff_automate_main_5archetypes.py"), stdout = TRUE)

archetype1_cut_off <- as.numeric(output[1])
archetype2_cut_off <- as.numeric(output[2])
archetype3_cut_off <- as.numeric(output[3])
archetype4_cut_off <- as.numeric(output[4])
archetype5_cut_off <- as.numeric(output[5])
# archetype6_cut_off <- as.numeric(output[6])

archetype1_cut_off
archetype2_cut_off
archetype3_cut_off
archetype4_cut_off
archetype5_cut_off
# archetype6_cut_off

sample_name <- "SMC21-T"

distance_ordering_path <- glue("../../pipeline_testing/distance_and_ordering/{sample_name}")
gene_exp_path <- glue("../../pipeline_testing/data/{sample_name}/")

distance_ordering_files <- list.files(path = distance_ordering_path)
distance_file <- distance_ordering_files[str_detect(distance_ordering_files, ".csv$") & str_detect(distance_ordering_files, "cell_distances")]
ordering_file <- distance_ordering_files[str_detect(distance_ordering_files, ".csv$") & str_detect(distance_ordering_files, "cell_ordering")]

gene_exp_files <- list.files(path = gene_exp_path)
TPM_gene_exp_file <- gene_exp_files[str_detect(gene_exp_files, ".csv$") & str_detect(gene_exp_files, "TPM_mod")]
celltype_data_file <- gene_exp_files[str_detect(gene_exp_files, ".tsv$") & str_detect(gene_exp_files, "modified")]

distance <- read.csv(paste(distance_ordering_path,distance_file,sep="/"),row.names=1)
ordering <- read.csv(paste(distance_ordering_path,ordering_file,sep="/"),row.names=1)
celltype_data <- read.table(paste(gene_exp_path,celltype_data_file,sep="/"),sep="\t",header=TRUE)

#read the main tpm data file and keep only deg genes for further heatmap generation
#Detection of tpm_data file in the folder
gene_exp_data <- read.csv(paste(gene_exp_path,TPM_gene_exp_file,sep="/"),row.names=1, check.names = FALSE) 

#heatmap generation for upregulated genes
#gene_exp_up <- gene_exp_data[rownames(gene_TPM_gene_exp_fileexp_data) %in% deg_superset_up,]

#show the dimension of gene expression data
#print(paste("the up_deg raw matrix dimnesion is", dim(gene_exp_up)))

#separate out specialists based on the distance
#separate out cells based on distance from each archetype
dist_a1 <- distance[1,]
dist_a1 <- dist_a1 [,dist_a1<=archetype1_cut_off]

dist_a2 <- distance[2,]
dist_a2 <- dist_a2[,dist_a2<=archetype2_cut_off]

dist_a3 <- distance[3,]
dist_a3 <- dist_a3[,dist_a3<=archetype3_cut_off]

dist_a4 <- distance[4,]
dist_a4 <- dist_a4[,dist_a4<=archetype4_cut_off]

dist_a5 <- distance[5,]
dist_a5 <- dist_a5[,dist_a5<=archetype5_cut_off]

# dist_a6 <- distance[6,]
# dist_a6 <- dist_a6[,dist_a6<=archetype6_cut_off]


order_specialist_a1 <- ordering[1,1:ncol(dist_a1)]
order_specialist_a2 <- ordering[2,1:ncol(dist_a2)]
order_specialist_a3 <- ordering[3,1:ncol(dist_a3)]
order_specialist_a4 <- ordering[4,1:ncol(dist_a4)]
order_specialist_a5 <- ordering[5,1:ncol(dist_a5)]
# order_specialist_a6 <- ordering[6,1:ncol(dist_a6)]

#extract the specialists from the celltype data
specialist_a1_cells <- celltype_data[match(order_specialist_a1,rownames(celltype_data)),]
specialist_a2_cells <- celltype_data[match(order_specialist_a2,rownames(celltype_data)),]
specialist_a3_cells <- celltype_data[match(order_specialist_a3,rownames(celltype_data)),]
specialist_a4_cells <- celltype_data[match(order_specialist_a4,rownames(celltype_data)),]
specialist_a5_cells <- celltype_data[match(order_specialist_a5,rownames(celltype_data)),]
# specialist_a6_cells <- celltype_data[match(order_specialist_a6,rownames(celltype_data)),]

#extract specialists gene expression data
specialist_a1_gene_exp <- gene_exp_data[,match(specialist_a1_cells$index,colnames(gene_exp_data))]
specialist_a2_gene_exp <- gene_exp_data[,match(specialist_a2_cells$index,colnames(gene_exp_data))]
specialist_a3_gene_exp <- gene_exp_data[,match(specialist_a3_cells$index,colnames(gene_exp_data))]
specialist_a4_gene_exp <- gene_exp_data[,match(specialist_a4_cells$index,colnames(gene_exp_data))]
specialist_a5_gene_exp <- gene_exp_data[,match(specialist_a5_cells$index,colnames(gene_exp_data))]
# specialist_a6_gene_exp <- gene_exp_data[,match(specialist_a6_cells$index,colnames(gene_exp_data))]

#Discard the cells which are common among two or more than two archetypes
all_specialists<-list(specialist_a1_gene_exp,specialist_a2_gene_exp,specialist_a3_gene_exp,specialist_a4_gene_exp,specialist_a5_gene_exp)

total_common_col <- character()

for (i in 1:length(all_specialists)) {
  current_df <- all_specialists[[i]]
  
  for (j in 1:length(all_specialists)) {
    if (j > i) {
        compare_df <- all_specialists[[j]]
        common_cols_iter <- intersect(colnames(current_df), colnames(compare_df))
        total_common_col <- c(total_common_col,common_cols_iter)
    }
  }
}

# Compare column names of one dataframe with the rest
for (i in 1:length(all_specialists)) {
  current_df <- all_specialists[[i]]
  
  for (j in 1:length(all_specialists)) {
    if (j > i) {
      compare_df <- all_specialists[[j]]
      common_cols <- intersect(colnames(current_df), colnames(compare_df))
      
      if (length(common_cols) > 0) {
        print(paste("Common columns between dataframe", i, "and dataframe", j))
        print(length(common_cols))
        
        # Remove common columns from current dataframe
        current_df <- current_df[, !colnames(current_df) %in% common_cols]
        
        # Remove common columns from compare dataframe
        compare_df <- compare_df[, !colnames(compare_df) %in% common_cols]
        
        all_specialists[[i]] <- current_df
        all_specialists[[j]] <- compare_df
      }
    }
  }
}

all_specialists_gene_exp <- cbind(all_specialists[[1]],all_specialists[[2]],all_specialists[[3]],all_specialists[[4]],all_specialists[[5]])
generalists_gene_exp <- gene_exp_data[,!(colnames(gene_exp_data) %in% colnames(all_specialists_gene_exp))]
generalists_gene_exp <- generalists_gene_exp[,!(colnames(generalists_gene_exp) %in% (total_common_col))]  

generalist_cell <- colnames(generalists_gene_exp)

dir.create(glue("../../pipeline_testing/specialist_generalist/{sample_name}"), recursive = TRUE)
write.table(generalist_cell, glue("../../pipeline_testing/specialist_generalist/{sample_name}/generalists_cell_population.tsv"), sep="\t", quote=FALSE)

for (j in 1:length(all_specialists)){
    data <- all_specialists[[j]][,!(colnames(all_specialists[[j]]) %in% colnames(total_common_col))]
    data <- data[,!(colnames(data) %in% colnames(generalists_gene_exp))]
    write.table(colnames(data), glue("../../pipeline_testing/specialist_generalist/{sample_name}/Archetype_{j}_specialists.tsv"),sep="\t", quote=FALSE)
}  

#formatted the gene expression table where all the specialists will be position first and then specialists
formatted_gene_exp_data <- cbind(all_specialists_gene_exp,generalists_gene_exp)

#save the cells information that is going for heatmap generation 
#& save those cells too which are common among archetypes pairwise
#write.csv(formatted_gene_exp_data,"gene_exp_data_formatted_LUNG_N06_heatmap.csv")
write.table(total_common_col,glue("../../pipeline_testing/specialist_generalist/{sample_name}/common_cells_between_archetypes_pairwise.txt"),row.names = FALSE, col.names = FALSE)

#formatted the gene expression table where all the specialists will be position first and then specialists
specialist_generalist_data <- cbind(all_specialists_gene_exp, generalists_gene_exp)

#show the dimension of specialist-generalist data
print(paste("specialist_generalist whole dataframe dimension",dim(specialist_generalist_data)))

#delete rows with rowSums==0
specialist_generalist_data <- specialist_generalist_data[rowSums(specialist_generalist_data)>=1,]
#calculate varience across rows
specialist_generalist_data$varience <- apply(specialist_generalist_data[,-1], 1, var)
min(specialist_generalist_data$varience)
max(specialist_generalist_data$varience)
heatmap_data <- specialist_generalist_data[specialist_generalist_data$varience>0.1,]
heatmap_data <- heatmap_data[,1:ncol(heatmap_data)-1]
print(paste("heatmap input data dimension",dim(heatmap_data)))
matrix_heatmap <- as.matrix(heatmap_data)

#png(filename="heatmap_specialists.png",width=1800,height=1800,units="px")
#my_palette <- colorRampPalette(c("red", "blue", "green"))(n = 100)
#breaks <- c(-50, -20, 20, 50)
dir.create("../../pipeline_testing/heatmap_parti/", recursive = TRUE)

pdf(file=glue("../../pipeline_testing/heatmap_parti/{sample_name}_heatmap_specialist_&_generalist.pdf"),width=40,height=50)
heatmap.2(matrix_heatmap,trace="none",Colv=FALSE,Rowv=TRUE,dendrogram='none',cexRow=1,cexCol=0.5,keysize=0.8,col=bluered(256), margins=c(25,15),ColSideColors=c(rep("royalblue",ncol(all_specialists[[1]])), rep("cyan",ncol(all_specialists[[2]])),rep("olivedrab",ncol(all_specialists[[3]])),rep("seagreen",ncol(all_specialists[[4]])),rep("purple",ncol(all_specialists[[5]])), rep("rosybrown",ncol(generalists_gene_exp))),scale='row',breaks=seq(-2.5,2.5,length.out=257))
legend("topright",legend=c("Archetype 1","Archetype 2","Archetype 3","Archetype 4", "Archetype 5", "Generalists"),fill=c("royalblue", "cyan", "olivedrab", "seagreen", "purple", "rosybrown"),border=FALSE, bty="n", y.intersp = 1.5, cex=3.5)
dev.off()
