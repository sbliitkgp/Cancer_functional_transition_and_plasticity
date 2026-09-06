### Mapping of active TFs in archetypes of all samples 
### Active TF-target map 

library(stringr)
library(Dict)
library(AUCell)
library(GSEABase)

ident_aucell_tfs <- function(my_collat_file,my_grn_connect,my_expr_data,my_cell_ids,my_name,my_samp) {

                #tf_ident = ident_aucell_tfs(collat_file,grn_connect,expr_data,cell_ids,name,arch)

                filt_collat_data = my_collat_file[which(rownames(my_collat_file) %in% my_cell_ids),]
                cell_types = unique(filt_collat_data$cell_type)

		for (ctype in cell_types) { 
			ctt = filt_collat_data[which(filt_collat_data$cell_type == ctype),]

			if(dim(ctt)[1] > 20) { 
				tmp_cell_ids = rownames(ctt) 
				grn_collect = my_grn_connect[[ctype]]
                
				if(length(grn_collect) != 1) { 
		                	
					gene_set = vector()
					gene_set = rownames(my_expr_data) 
	
        	        		tflist = vector()
                			for(target in gene_set) {
                        			search_id = paste0('-',target)
                        			mapped_conn = grn_collect[str_detect(grn_collect,search_id)]
                        			tf = word(mapped_conn,1,sep='-')
                        			tflist = union(tflist,tf)
                        			rm(tf)
                        			rm(search_id)
                        			rm(mapped_conn)
                			}
                			rm(gene_set)

					module = list()
                			for(tf in tflist) {
                        			search_id = paste0(tf,'-')
                        			mapped_conn = grn_collect[str_detect(grn_collect,search_id)]
                        			tg = word(mapped_conn,2,sep='-')
                        			module[[tf]] = tg
                        			rm(search_id)
                        			rm(tg)
                			}

                			rm(tf)


                			filt_expr_data = my_expr_data[,which(colnames(my_expr_data) %in% tmp_cell_ids)]
                			exprMatrix <- as(as.matrix(filt_expr_data), "dgCMatrix")

                			geneSets = vector()

                			for(tf in tflist) {
                         			geneset = GeneSet(unique(module[[tf]]),setName = tf)
                         			geneSets = append(geneSets,geneset)
                         			rm(geneset)
                			}

                			geneSets = GeneSetCollection(geneSets)

                			cells_AUC <- AUCell_run(exprMatrix, geneSets)
                			auc_matrix = getAUC(cells_AUC)

					if(ctype == 'T/NK') { ctype = 'T-NK' } 
                			wrname = paste0('../results/plasticity_aucell_TFs/',my_name,'/',my_name,'_',my_samp,'_',ctype,'_aucMatrix.tsv')
                			write.table(round(auc_matrix,3),wrname,quote=F)
                			rm(wrname)

                			rm(module)
                			rm(cells_AUC)
                			rm(geneSets)

                			## rowMedians(auc_matrix)
                			aa = as.data.frame(rowMedians(auc_matrix))
                			rownames(aa) = rownames(auc_matrix)
                			sorted_aa = aa[order(aa[,1],decreasing=T),,drop=FALSE]
                			selected_TFs =  cbind(rownames(sorted_aa)[1:10],round(sorted_aa[1:10,],3))

                			wrname = paste0('../results/plasticity_aucell_TFs/',my_name,'/',my_name,'_',my_samp,'_',ctype,'_topTFs.tsv')
                			write.table(selected_TFs,wrname,quote=F,col.names=F,row.names=F)
                			rm(wrname)

                			rm(auc_matrix)
                			rm(aa)
                			rm(sorted_aa)
                			rm(tf)
                			rm(tflist)
					rm(selected_TFs) 

                			rm(exprMatrix)
                			rm(filt_expr_data)
				}
				rm(tmp_cell_ids) 
                		rm(grn_collect)

		      }
		      rm(ctt)
		}
		rm(ctype) 
                rm(filt_collat_data)
                rm(cell_types)

                #return(selected_TFs)
}


metamap = list()

metadata_file = read.table("../datasets/Metadata.txt",header=T)
metadata = metadata_file$Sample

cn = 1
for(samp in metadata_file$Sample) {
        metamap[[samp]] = metadata_file$Stage_Sample[cn]
        cn = cn+1
}

print(names(metamap))

collat_path = '../results/collated_data/'

for(name in names(metamap)) {

        print(name)
        collat_file = read.csv(paste0(collat_path,'Collated_',name,'.csv'),row.names=1)
        dir.create(file.path('../results/plasticity_aucell_TFs/', name),showWarnings = F)

        ## Read gene expression for the sample
        gene_expr_path = '/home/subhasis/functional_heterogeneity_pipeline/preprocessing/'
        exprfolder = paste0(name,'_all_genes/')
        exprfile = paste0(name,'_tpm.csv')

        expr_data = read.csv(paste0(gene_expr_path,exprfolder,exprfile),row.names=1)


        ## Read active TF-target map

        cc = metamap[[name]]
        cc = paste0(word(cc,1,sep='_'),'.',word(cc,2,3,sep='_'))
        netpath = '../datasets/samples_filtered_correlations/'
        if((name != 'LUNG_T08') & (name != 'LUNG_T09')) {
            epi_grn = read.csv(paste0(netpath,'Epi_',cc,'_filt.csv'))
        }
        myl_grn = read.csv(paste0(netpath,'Myl_',cc,'_filt.csv'))
        tnk_grn = read.csv(paste0(netpath,'T-NK_',cc,'_filt.csv'))

        grn_connect = list()
        if((name != 'LUNG_T08') & (name != 'LUNG_T09')) {
                grn_connect[['Epithelial']] = paste(epi_grn$Source,epi_grn$Target,sep='-')
        } else {
                grn_connect[['Epithelial']] = ''
        }

	grn_connect[['Myeloid']] = paste(myl_grn$Source,myl_grn$Target,sep='-')
        grn_connect[['T/NK']] = paste(tnk_grn$Source,tnk_grn$Target,sep='-')

        if((name != 'LUNG_T08') & (name != 'LUNG_T09')) {
                rm(epi_grn)
        }
        rm(myl_grn)
        rm(tnk_grn)
        rm(netpath)
        rm(cc)

	archlist = unique(collat_file$Arch_after_transitory)
	select_archlist = archlist[!is.na(archlist) & (str_detect(archlist,'^Archetype'))]

        for(arch in select_archlist) {

                cell_ids = rownames(collat_file[which(collat_file$Arch_after_transitory == arch),])
                print(paste0('    ',arch))
                ident_aucell_tfs(collat_file,grn_connect,expr_data,cell_ids,name,arch)
                rm(cell_ids)

        }
        rm(grn_connect)
}

rm(metamap)
rm(metadata)
rm(metadata_file)
rm(collat_path)
