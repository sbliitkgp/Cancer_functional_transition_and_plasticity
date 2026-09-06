## Mapping the GRN behind archetype transitions 
## Based on active GRN (Ray et al., 2026) 
## For TransMyl and SpecialistMyl within a sample 


library(stringr) 
library(Dict) 
library(AUCell) 
library(GSEABase)
library(data.table) 

## Functions pass: grn_connect[['Myeloid']],target_gene_list,transMyl_data,name,'TransMyl'

ident_aucell_tfs <- function(my_grn_connect,my_targets,my_expr_data,my_name,my_samp) { 

		grn_collect = my_grn_connect

                if(length(grn_collect) == 1) {
                        rm(grn_collect)
                        return('NA')
                }

		my_expr_data = t(my_expr_data)
		mode(my_expr_data) <- 'numeric'
		gene_set = my_targets
 
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
			rm(tg) 
			rm(mapped_conn)
			rm(search_id) 
		}
		
		rm(tf)


		## filt_expr_data = my_expr_data[,which(colnames(my_expr_data) %in% my_cell_ids)]
		## name, word(dir,6,7,sep='/') 

		exprMatrix <- as(as.matrix(my_expr_data), "dgCMatrix")

		geneSets = vector()

		for(tf in tflist) {
			 geneset = GeneSet(unique(module[[tf]]),setName = tf)	
                	 geneSets = append(geneSets,geneset) 
			 rm(geneset) 
		}	

		geneSets = GeneSetCollection(geneSets) 

		cells_AUC <- AUCell_run(exprMatrix, geneSets)
		auc_matrix = getAUC(cells_AUC)

		wrname = paste0('../results/aucell_TFs_Myl/',my_name,'/',my_name,'_',my_samp,'_aucMatrix.tsv')
		write.table(round(auc_matrix,3),wrname,quote=F) 
                rm(wrname) 

		rm(module) 
		rm(cells_AUC)
		rm(geneSets) 


		## rowMeans(auc_matrix)
		#aa = rowMeans(auc_matrix)
                #sorted_aa = as.data.frame(rev(sort(aa)))
		#selected_TFs =  rownames(sorted_aa)[1:5]

		## rowMedians(auc_matrix)
		aa = as.data.frame(rowMedians(auc_matrix))
  		rownames(aa) = rownames(auc_matrix)
		sorted_aa = aa[order(aa[,1],decreasing=T),,drop=FALSE]
		selected_TFs =  cbind(rownames(sorted_aa)[1:10],round(sorted_aa[1:10,],3))

		wrname = paste0('../results/aucell_TFs/',my_name,'/',my_name,'_',my_samp,'_topTFs.tsv')
		write.table(selected_TFs,wrname,quote=F,col.names=F,row.names=F) 
		rm(wrname) 

		rm(auc_matrix)
		rm(aa)
		rm(sorted_aa) 
		rm(tf) 
		rm(tflist) 
		rm(grn_collect) 
		rm(exprMatrix)

		return(selected_TFs) 

		##rm(selected_TFs) 
}


samples = c('BRONCHO_58','EBUS_06','LUNG_T06','LUNG_T08','LUNG_T09','LUNG_T19','LUNG_T20','LUNG_T25','LUNG_T28','LUNG_T30','LUNG_N18','LUNG_N31')  ## Samples with TransMyl


metamap = list()

metadata_file = read.table("../datasets/Metadata.txt",header=T)

cn = 1
for(samp in metadata_file$Sample) {
        metamap[[samp]] = metadata_file$Stage_Sample[cn]
        cn = cn+1
}

print(names(metamap))


datafile = '../MLmodel/MLdata_Myl_LUAD/Merged_myeloid_data.csv'
data = fread(datafile,data.table=F)
rownames(data) = data$V1
data = data[-1]
data[1:10,1:10]


for(name in samples) {

	print(name)
	dir.create(file.path('../results/aucell_TFs_Myl/', name),showWarnings = F)


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

        myl_grn_targets = unique(myl_grn$Target)  ## Only Myl targets 

	rm(myl_grn)
	rm(tnk_grn)
	rm(netpath)
	rm(cc)


        select_data = data[str_detect(rownames(data),name),]
	transMyl_data = select_data[which(select_data$label == 'TransMyl'),]
	transMyl_data = transMyl_data[-ncol(transMyl_data)]
	dim(transMyl_data)
	specMyl_data = select_data[which(select_data$label == 'Myl_DomMyl'),]
	specMyl_data = specMyl_data[-ncol(specMyl_data)]
	dim(specMyl_data)

	### AUcell data format to be created


	print(paste0('    ','TransMyl'))
	tf_ident1 = ident_aucell_tfs(grn_connect[['Myeloid']],myl_grn_targets,transMyl_data,name,'TransMyl')

	print(paste0('    ','SpecialistMyl'))
	tf_ident2 = ident_aucell_tfs(grn_connect[['Myeloid']],myl_grn_targets,specMyl_data,name,'SpecialistMyl')

        only_transMyl <- tf_ident1[!(tf_ident1[,1] %in% tf_ident2[,1])]
        only_specMyl <- tf_ident2[!(tf_ident2[,1] %in% tf_ident1[,1])]

        if(length(only_transMyl) == 0) {  only_transMyl = cbind('NA','NA') }
        if(length(only_specMyl) == 0) {  only_specMyl = cbind('NA','NA') }

        wrfile = paste0('../results/aucell_TFs_Myl/',name,'/','UNIQUE_',name,'_TransMyl_SpecMyl.tsv')
        sink(wrfile)
        cat(c('Only_TransMyl: ',t(only_transMyl),'\n'))
        cat(c('Only_SpecialistMyl: ',t(only_specMyl),'\n'))
        sink()

        rm(wrfile)
        rm(only_transMyl)
	rm(only_specMyl)

	newdf = data.frame(cbind(tf_ident1,tf_ident2))
	colnames(newdf) = c('TF_TransMyl','Score_TransMyl','TF_SpecialistMyl','Score_SpecialistMyl')
	wrfile = paste0('../results/aucell_TFs_Myl/',name,'/','Compiled_',name,'_TransMyl_SpecMyl.tsv')
	write.table(newdf,file=wrfile,row.names=F,quote=F)

	rm(wrfile)
	rm(newdf)
        rm(tf_ident1)
	rm(tf_ident2)

        rm(grn_connect)

}

rm(metamap)
rm(metadata_file)


