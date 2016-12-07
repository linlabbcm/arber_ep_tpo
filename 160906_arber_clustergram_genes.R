#Load tables
control_v_ep <- read.delim('/grail/projects/arber2/cufflinks/arber_rna_cuffnorm/output/arber_rna_ARBER_CONTROL_REP_vs_ARBER_EP_REP_exprs_matrix.txt', header=TRUE, sep='\t')
control_v_tpo <- read.delim('/grail/projects/arber2/cufflinks/arber_rna_cuffnorm/output/arber_rna_ARBER_CONTROL_REP_vs_ARBER_TPO_REP_exprs_matrix.txt', header=TRUE, sep='\t')
ep_v_tpo <- read.delim('/grail/projects/arber2/cufflinks/arber_rna_cuffnorm/output/arber_rna_ARBER_EP_REP_vs_ARBER_TPO_REP_exprs_matrix.txt', header=TRUE, sep='\t')

#Adjust p-values
pvalues_1 = as.numeric(control_v_ep[,4])
pvalues_2 = as.numeric(control_v_tpo[,4])
pvalues_3 = as.numeric(ep_v_tpo[,4])
padjusted1= p.adjust(pvalues_1,method='BH') #this is fdr correction
padjusted2= p.adjust(pvalues_2,method='BH')
padjusted3= p.adjust(pvalues_3,method='BH')

#insert adjusted p-value column
control_v_ep[,'P-Adjusted'] <- padjusted1
control_v_tpo[,'P-Adjusted'] <- padjusted2
ep_v_tpo[,'P-Adjusted'] <- padjusted3

#filter each table based on adjusted p-value >0.1 and fold change >= 1
sig_genes1 = which(control_v_ep$P_VALUE <=  0.01)
fold_genes1 = which(abs(control_v_ep$LOG2_FOLD_CHANGE) >= 0.5849625)

sig_fold_genes1 = intersect(sig_genes1, fold_genes1)
genes1 = rownames(control_v_ep[sig_fold_genes1,])


sig_genes2 = which(control_v_tpo$P_VALUE <=  0.01)
fold_genes2 = which(abs(control_v_tpo$LOG2_FOLD_CHANGE) >= 0.5849625)

sig_fold_genes2 = intersect(sig_genes2, fold_genes2)
genes2 = rownames(control_v_tpo[sig_fold_genes2,])


sig_genes3 = which(ep_v_tpo$P_VALUE <=  0.01)
fold_genes3 = which(abs(ep_v_tpo$LOG2_FOLD_CHANGE) >= 0.5849625)

sig_fold_genes3 = intersect(sig_genes3, fold_genes3)
genes3 = rownames(ep_v_tpo[sig_fold_genes3,])

#combine lists of found genes and pare down to unique list
genesList=c()

genesList1 <-union(genes1, genes2)
genesList2 <-union(genesList1, genes3)
#genesList <-c(genes1, genes2, genes3)

genesList = unique(genesList2)
#genesList = unique(genesList)


fileConn<-file("/grail/projects/arber2/clustergram_genes_1.5_fold_change.txt")
writeLines(genesList, fileConn)
close(fileConn)

#look at overlaps of control vs ep and control vs tpo
overlaps = intersect(genes1,genes2)

