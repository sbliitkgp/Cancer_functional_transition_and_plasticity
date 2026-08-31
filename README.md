# Cancer_functional_transition_and_plasticity
Functional transition and plasticity analysis of cell populations in the tumor microenvironment with cancer progression

## Part 1 - Preprocessing
Aim: To check batch effect in the dataset across patient samples and make the data ready for the input to ParTI software.

      A. code1_data_preprocess_and_batch_effect_check.ipynb
      
         Check the batch effect using UMAP based clustering and split the dataframe samplewise
         
      B. code2_preprocessing.R

         It performs centering on the single cell gene expression data along with the metadata and genes list ready as input for ParTI

## Part 2 - Specialist_Generalist_and_Transitory_identification
Aim: To develop a pipeline for the identification of specialist, generalists and transitory cell populations in the tumor 
     microenvironment of each patient sample individually.

      A. code1_distance_measure.R

         Measurement of euclidean distance of each cells with respect to different archetypes

      B. code2_Mixture model fitting modified.ipynb

         Build a mixture model using multiple normal distributions to fit the multimodal distance distribution for each archetype

      C. code3.1_distance_cutoff_automate_main_3archetypes.py

         Extract the statistically validated distance cut-off for each archetype from polytope with 3 archetypes only. 

      D. code3.1_specialists_vs_generalists_extraction_3archetypes.R

         Extract the specialist and generalist cell population based on the distance cut-off for samples with 3 archetypes based 
         
         polytope only. Also check the gene expression pattern of all genes across different specialists and generalist using 
         
         hierarchical clustering based heatmap plot

      E. code3.2_distance_cutoff_automate_main_4archetypes.py

         Extract the statistically validated distance cut-off for each archetype from polytope with 4 archetypes only

      F. code3.2_specialists_vs_generalists_extraction_4archetypes.R

         Extract the specialist and generalist cell population based on the distance cut-off for samples with 4 archetypes based 
         
         polytope only. Also check the gene expression pattern of all genes across different specialists and generalist using 
         
         hierarchical clustering based heatmap plot

      G. code3.3_distance_cutoff_automate_main_5archetypes.py

         Extract the statistically validated distance cut-off for each archetype from polytope with 5 archetypes only

      H. code3.3_specialists_vs_generalists_extraction_5archetypes.R

         Extract the specialist and generalist cell population based on the distance cut-off for samples with 5 archetypes based 
         
         polytope only. Also check the gene expression pattern of all genes across different specialists and generalist using 
         
         hierarchical clustering based heatmap plot

      I. code3.4_distance_cutoff_automate_main_6archetypes.py

         Extract the statistically validated distance cut-off for each archetype from polytope with 6 archetypes only

      J. code3.4_specialists_vs_generalists_extraction_6archetypes.R

         Extract the specialist and generalist cell population based on the distance cut-off for samples with 6 archetypes based 
         
         polytope only. Also check the gene expression pattern of all genes across different specialists and generalist using 
         
         hierarchical clustering based heatmap plot

      K. code4_Umap_clustering_on_gene_expression_data.ipynb

         Cell reclustering and mapping to archetype using UMAP integrated with HDBSCAN based clustering.

      L. code5_seurat_graph_based_clustering.ipynb

         Reclustering and mapping of cells to archetype using seurat SNN (Shared nearest neighbour) graph based clustering which UMAP 
         
         with HDBSCAN failed to cluster.

      M. code6_mellon_density_cut_off_estimation.ipynb

         Calculate physical density of all cells in a sample and build a GMM (Gaussian mixture model) to the density distribution and 
         
         identify the statistically validated density cut-off to robustly capture stable cell states from specialists and generalist

      N. code7_transitory_cell_analysis.ipynb

         Identify the rare transitory cells between archetypes in each sample. Perform archetype to archetype pseudotime trajectory 
         
         analysis using Palantir. 

      O. code8_all_specialist_generalist_consolidated.ipynb

         Consolidation of all data including ParTI output, cell's identity before and after Mellon as specialist, generalist or 
         
         transitory, celltype information along with marker's centered expression.

      P. code9_heatmap_after_mellon.ipynb

         Heatmap plots of gene expression pattern across specialists and generalist in each sample after final identification from   
         
         Mellon.

## Part 3 - Functional_characterization_and_plasticity
Aim: Uncover the functional programs associated with each archetype in each sample and identify the convergence of functional 
     programs with cancer progression
     
     A. code1_Archetype_purity_check_cell_type_wise_acros_samples.ipynb
     
        Identify diversification of functional metaprograms in individual celltype with cancer progression
        
     B. code2_pareto_front_plots.ipynb

        2D and 3D scatter plot of cells inside the pareto front colored according to cell type and mellon density. 

        2D scatter plot of cells inside pareto front colored according to functional programs

     C. code3_functional_characterization_archetypes.ipynb

        Identify the significantly enriched functional programs at each archetype in individual samples

        Identify significantly enriched functional programs in individual cell types with cancer progression
        
     D. code4_archetype_alignment_analysis.ipynb

        Calculate the pairwise similarity of functional programs between all archetypes across all samples

        Perform hierarchical clustering to identify highly similar archetypes in terms of functional programs and their celltype composition

     E. code5_markers_exp_plot_celltype.ipynb

        Stack barplot of expression of marker genes in different cell types with cancer progression
         
         
         
         
   
