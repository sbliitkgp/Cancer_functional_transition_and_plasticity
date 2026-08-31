#This script will run only sample wise. Hence plese change the folder path for each sample run

path_env <- "/data/subhasis_backup/subhasis/pipeline_testing/parti_result/SMC02-T/" # Please give the folder path of the sample where it's ParTI output is present
setwd(path_env)
library(R.matlab)
library(flexclust)
matlab_file <- list.files(path_env,pattern="*.mat")
folder_name <- basename(path_env) 
data <- readMat(matlab_file)
pc <- as.data.frame(data$pc)

#this is for all epithelial cell parti result
#pc <- read.csv("pc_data.csv")

distances_parti <- as.data.frame(data$distances)
ordering_parti <- as.data.frame(data$ordering)
arc_parti <- as.data.frame(data$arc)
dist_calc <- pc[,1:nrow(ordering_parti)-1]
dist_calc <- as.matrix(dist_calc)

if(nrow(ordering_parti)==3){
	colnames(dist_calc) <- c("pc1","pc2")
}else if(nrow(ordering_parti)==4){
	colnames(dist_calc) <- c("pc1","pc2","pc3")
}else if(nrow(ordering_parti)==5){
	colnames(dist_calc) <- c("pc1","pc2","pc3","pc4")
}else if(nrow(ordering_parti)==6){
	colnames(dist_calc) <- c("pc1","pc2","pc3","pc4","pc5")
}else {
	colnames(dist_calc) <- c("pc1","pc2","pc3","pc4","pc5","pc6")
}

#eucleadian distance measurement from all archetypes
x <- seq(1,nrow(ordering_parti),1)
cell_ordering_from_archetype <- list()
distance_from_archetype <- list()
for (i in 1:length(x)){
    vect <- as.data.frame(arc_parti[i,])
    archetype_no <- do.call("rbind", replicate(nrow(pc), vect, simplify = FALSE))
    archetype_no <- as.matrix(archetype_no)
    distances_from_arch <- dist2(dist_calc, archetype_no, method = "euclidean",p=2)
    distances_from_arch <- as.data.frame(distances_from_arch[,i])
    distances_from_arch <- distances_from_arch[order(distances_from_arch$`distances_from_arch[, i]`,decreasing = FALSE), , drop = FALSE]
    colnames(distances_from_arch) <- "dist_from_archetype"
    distance_from_archetype[[i]] <- array(data = c(unlist(distances_from_arch)))
    cell_ordering_from_archetype[[i]] <- rownames(distances_from_arch)
}

#add all the ordering of cells and the distance to a single dataframe
ordering_cells <- data.frame(t(sapply(cell_ordering_from_archetype, c)))
cell_distance_from_vertex <- data.frame(t(sapply(distance_from_archetype, c)))

export_path <- "/data/subhasis_backup/subhasis/pipeline_testing/distance_and_ordering"  #put your destination folder directory path here
dir.create(paste(export_path,folder_name, sep="/"), recursive=T)
path <- paste(export_path,folder_name,sep="/")
file_path_ordering <- file.path(path, "cell_ordering.csv")
write.csv(ordering_cells,file=file_path_ordering)
file_path_distance <- file.path(path, "cell_distances.csv")
write.csv(cell_distance_from_vertex,file=file_path_distance)