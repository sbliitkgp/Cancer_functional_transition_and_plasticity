suppressPackageStartupMessages({
    library(ggplot2)
    library(GGally)
    library(GSEABase)
    library(limma)
    library(reshape2)
    library(data.table)
    library(knitr)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(stringr)
    library(NMF)
    library(rsvd)
    library(RColorBrewer)
    library(MAST)
    library(janitor) 
    library(dplyr) 
    library(showtext)
})

options(mc.cores = 4)
knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE,cache = FALSE,fig.width=8,fig.height=6)
showtext_auto()

freq_expressed <- 0.2
FCTHRESHOLD <- 0 ##log2(1.5)

print('Reading data...') 

data = fread('../../MLmodel2/MLdata_TNK/Merged_DomTNK_data_AllSamples_LUAD.csv',data.table=F)
rownames(data) = data$V1 
col_label = data$label
data = data[-1]

data = data[-ncol(data)]
data = remove_constant(data,na.rm=T)
rownames(data) = paste0(rownames(data),'_',col_label)
data[1:10,1:10]

metadata = as.data.frame(cbind('sample' = rownames(data),'label'= col_label))

scaRaw <- FromMatrix(t(data), metadata)


## Exploratory analysis ========

print('PCA analysis...') 

set.seed(123)
plotPCA <- function(sca_obj){
    projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
    colnames(projection)=c("PC1","PC2","PC3","PC4")
    pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
    print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'PC4'),
            mapping=aes(color=label), upper=list(continuous='blank')))
    invisible(pca)
}

pdf(file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_PCA.pdf',useDingbats=F)
pr = plotPCA(scaRaw)
print(pr)
dev.off()

rm(pr) 


sca = scaRaw

scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
pr = ggplot(flat, aes(x=value))+geom_density() +facet_wrap(~primerid, scale='free_y')

pdf(file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_20genes_distr.pdf',useDingbats=F)
print(pr)
dev.off()

rm(pr)


expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes,]


## ==== Differential expression Hurdle model ========


print('DE analysis...') 

cond<-factor(colData(sca)$label)
cond<-relevel(cond,"NoTransMyl")
colData(sca)$label<-cond
zlmCond <- zlm(~label, sca)


## hypothesis testing

summaryCond <- summary(zlmCond, doLRT='labelWithTransMyl')
#print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
print(summaryCond, n=4, by='D') ## Discrete Z score
print(summaryCond, n=4, by='C') ## Continuous Z score


## Differentially expressed genes

summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='labelWithTransMyl' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='labelWithTransMyl' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.1 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)

write.table(fcHurdleSig,file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_DEresults.txt',quote=F)



## Visualization top 20

print('Visualization top genes...') 

min = 20
tot_genes = dim(fcHurdleSig)[1]
if(tot_genes<min) { min = tot_genes }

entrez_to_plot <- fcHurdleSig[1:min,primerid]
symbols_to_plot <- fcHurdleSig[1:min,primerid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=label, y=value, color=label)) + geom_jitter()+facet_wrap(~primerid, scale='free_y')+ggtitle("DE Genes in TNK Cells with TransMyl")
pr = ggbase+geom_violin()

pdf(file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_topGenes.pdf',useDingbats=F)
print(pr)
dev.off()

rm(pr)



## Mean vs proportion expression

print('Mean vs proportion expression...') 

MM <- model.matrix(~label,unique(colData(sca)[,c("label"),drop=FALSE]))
rownames(MM) <- str_extract(rownames(MM), 'WithTransMyl|NoTransMyl')
predicted <- predict(zlmCond,modelmatrix=MM)

## Avert your eyes...
predicted[, primerid:=as.character(primerid)]
predicted_sig <- merge(mcols(sca), predicted[primerid%in%entrez_to_plot], by='primerid')
predicted_sig <- as.data.table(predicted_sig)


## plot with inverse logit transformed x-axis

pr = ggplot(predicted_sig)+aes(x=invlogit(etaD),y=muC,xse=seD,yse=seC,col=sample)+
    facet_wrap(~primerid,scales="free_y")+theme_linedraw()+
    geom_point(size=0.5)+scale_x_continuous("Proportion expression")+
    scale_y_continuous("Estimated Mean")+
    stat_ell(aes(x=etaD,y=muC),level=0.95, invert='x')

pdf(file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_Prportion_vs_Mean.pdf',useDingbats=F)
print(pr)
dev.off()

rm(pr)


mat_to_plot <- assay(sca[entrez_to_plot,])
rownames(mat_to_plot) <- symbols_to_plot

pdf(file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_heatmap_DEgenes.pdf',useDingbats=F) 
pr = aheatmap(mat_to_plot,annCol=colData(sca)[,"label"],main="DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PuOr",n=10))(20)))
print(pr) 
dev.off() 


## GSEA

print('GSEA analysis...') 

# bootstrap, resampling cells
# R should be set to >50 if you were doing this for real.
boots <- bootVcov1(zlmCond, R = 50)


min_gene_in_module <- 5
module_file = '../datasets/gmt/c5.go.bp.v2026.1.Hs.symbols.gmt'
gene_set <- getGmt(module_file)
gene_ids <- geneIds(gene_set)
##gene_ids <- gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"B cell"]
sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$primerid)
# Only keep modules with at least min_gene_in_module
sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]

gsea <- gseaAfterBoot(zlmCond, boots, sets_indices, CoefficientHypothesis("labelWithTransMyl")) 
z_stat_comb <- summary(gsea, testType='normal')

sigModules <- z_stat_comb[combined_adj<.01][1:10]
gseaTable <- melt(sigModules[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')

write.table(gseaTable, file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_GSEA_enriched_fns.txt',quote=F)


pdf(file = '../results_DiffExprAnalysis/MAST/LUAD_TNKdiff_top10_enriched_modules.pdf',useDingbats=F) 
pr = ggplot(gseaTable, aes(y=set, x=variable, fill=value))+geom_raster() + scale_fill_distiller(palette="PuOr")
print(pr) 
dev.off() 

rm(pr) 


## ================================

rm(data) 
rm(list=ls()) 



