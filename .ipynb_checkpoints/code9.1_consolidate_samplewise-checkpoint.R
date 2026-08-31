library(glue)
library(stringr)
library(dplyr)

sample <- 'SMC01-T'

main_dir <- glue("../../pipeline_testing/specialist_generalist_with_cluster/{sample}")

target_files <- list.files(path = main_dir, pattern = "*.tsv$", full.names = TRUE)

sg_data_list <- lapply(target_files, read.table)

filenames <- tools::file_path_sans_ext(basename(target_files))
arc_names <- ifelse(grepl("generalists_cell", filenames, ignore.case = TRUE), "Generalist", gsub("_specialists.*", "", filenames))

sg_data_list <- Map(function(df, n) {
  df$archetype <- n
  return(df)
}, sg_data_list, arc_names)

sg_data_merged <- do.call(rbind, sg_data_list)

rownames(sg_data_merged) <- NULL

cell_removed_before_mellon <- read.table(glue("../../pipeline_testing/specialist_generalist/{sample}/common_cells_between_archetypes_pairwise.txt"), sep="\t")
names(cell_removed_before_mellon) <- 'x'
cell_removed_before_mellon$archetype <- NA
sg_data_merged <- rbind(sg_data_merged, cell_removed_before_mellon)

sample_main_meta_data <- read.table(glue("../../pipeline_testing/data/{sample}/{sample}_required_cell_modified.tsv"), sep="\t", header=T)

sg_data_merged <- merge(sg_data_merged, sample_main_meta_data, by.x = "x", by.y = "index", all = TRUE)
sg_data_merged <- sg_data_merged[,1:2]

cons_full_df <- read.csv("../../pipeline_testing/pareto_front_consolidated_data.csv", row.names=1)
rownames(cons_full_df) <- cons_full_df$index
sample_df <-cons_full_df[cons_full_df$Sample == sample,]
sample_df <- sample_df[,c("PC1", "PC2", "PC3", "Cell_type", "archetype")]

sample_df <- sample_df %>% 
             rename(Arch_after_mellon = archetype)

sample_df$Arch_after_mellon <- gsub("Transitory_cells", NA, sample_df$Arch_after_mellon)
sample_df$Arch_after_mellon <- gsub("generalists", "Generalist", sample_df$Arch_after_mellon)

sample_df_mod <- merge(sg_data_merged, sample_df, by.x = "x", by.y = "row.names", all=TRUE)
sample_df_mod <- sample_df_mod %>%
                  rename(Arch_before_mellon = archetype) %>%
                  rename(index = x)

sample_df_mod  <- sample_df_mod [, c("index", "PC1", "PC2",	"PC3", "Cell_type",	"Arch_before_mellon", "Arch_after_mellon")]

transitory_data <- read.csv(glue("../../pipeline_testing/transitory_cells/{sample}/{sample}_archetype_contribution_score.csv"), row.names=1)

collated_final <- merge(sample_df_mod, transitory_data, by.x = "index", by.y = "row.names", all=TRUE)
# collated_final <- collated_final[, -c('Archetype 1', 'Archetype 2', 'Archetype 3', 'Archetype 4', 'Archetype 5', 'cont_score', 'cont_mod')]
collated_final <- subset(collated_final, select = -c(Archetype_1, Archetype_2, Archetype_3, Archetype_4, Archetype_5, cont, cont_score))
collated_final <- collated_final %>% 
                     rename(Arch_after_transitory = cont_mod)

collated_final$Arch_after_transitory[collated_final$Arch_after_transitory == ""] <- NA
collated_final$Arch_after_transitory <- ifelse(is.na(collated_final$Arch_after_transitory), collated_final$Arch_after_mellon, collated_final$Arch_after_transitory)

collated_final <- collated_final %>%
                      mutate(
                        across(
                          last_col(), 
                          ~ ifelse(
                              str_count(., "Archetype?") == 2,         # Condition: if 'Archetype' appears exactly twice
                              str_replace_all(., "Archetype?", "Transitory"), # True: Replace all occurrences with 'Transitory'
                              .                                       # False: Keep the original text
                          )
                        )
                      )

cent_exp_data <- read.table(glue("../../pipeline_testing/data/{sample}/SMC01-N_gene_expression_parti_input.tsv"))
cell_data <- read.csv(glue("../../pipeline_testing/data/{sample}/SMC01-N_required_cell_modified.tsv"), sep="\t",header=T,row.names=1)
rownames(cent_exp_data) <-  rownames(cell_data)

genes <- read.csv(glue("../../pipeline_testing/data/{sample}/genes.list"), header=FALSE)
colnames(cent_exp_data) <- genes$V1

csc_markers <- c('ABCG2', 'ALCAM', 'ALDH1A1', 'BMI1', 'CD200', 'CD24', 'CD44', 'CEACAM6', 'DCLK1', 'DPP4', 'FLOT2', 'ITGB1', 'LGR5', 'NANOG', 'POU5F1', 'PROM1', 'MSI1', 'SOX9')
epi_markers <- c('MUC2', 'TFF1', 'EPCAM', 'KRT19', 'DCLK1', 'ATP1B3', 'CA1', 'CEACAM1', 'GUCA2A')
myo_markers <- c('CD68', 'CD14', 'FCN1', 'CD11B')
stromal_arkers <- c('TNC', 'VIM', 'COL1A1', 'TAGLN', 'SPARC', 'ACTA2', 'VWF', 'PECAM')
t_markers <- c('CD3D', 'CD3E', 'NKG7', 'CCR6', 'KLRB1', 'CD8A')

all_markers <- unique(c(csc_markers, epi_markers, myo_markers, stromal_arkers, t_markers))

cent_exp_data <- cent_exp_data[,colnames(cent_exp_data) %in% all_markers]
collated_final <- merge(collated_final, cent_exp_data, by.x = "index", by.y = "row.names", all=TRUE)

write.csv(collated_final, glue("../../pipeline_testing/consolidated_data_samplewise/{sample}_consolidated_data.csv"))