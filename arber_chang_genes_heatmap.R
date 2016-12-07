expr_table = read.delim('/grail/projects/arber2/cufflinks/arber_rna_cuffnorm/output/arber_rna_all_fpkm_exprs_raw.txt', header=TRUE, sep = '\t' )
enriched_chang_genes = read.delim('/grail/projects/arber2/tables/enriched_chang.txt', header = TRUE, sep = '\t')

enriched_rows = c()
for(gene in enriched_chang_genes[,1]){
  enriched_row = which(rownames(expr_table)==gene)
  enriched_rows = c(enriched_rows,enriched_row)
}
enriched_rows = unique(enriched_rows)
enrichedGenesList = rownames(expr_table[enriched_rows,])

########################################################################

expMatrix = as.matrix(expr_table)
enrichedMatrix = matrix(nrow=length(enrichedGenesList),ncol=ncol(expMatrix))

rownames(enrichedMatrix) = enrichedGenesList
colnames(enrichedMatrix) = colnames(expMatrix)


for(i in 1:length(enrichedGenesList)){
  
  geneName = enrichedGenesList[i]
  expRow = which(rownames(expMatrix)==geneName)
  enrichedMatrix[i,] = expMatrix[expRow,]
  
}

#now we want to make a log2 row median normalized expression matrix

medianVector = apply(enrichedMatrix,1,median) 

medianMatrix = log2(enrichedMatrix/medianVector)


###########################################################################


#===================================================================
#======================CLUSTERING EXPRESSION========================
#===================================================================

expCorDist = as.dist(1-cor(t(medianMatrix)))
expHClust = hclust(expCorDist)

expOrder = expHClust$order
#===================================================================
#=========================MAKING HEATMAPS===========================
#===================================================================


#Set the color spectrum
colorSpectrum <- colorRampPalette(c("blue","white","red"))(100)
#colorSpectrum <- colorRampPalette(c("white","red"))(100)

#setting a color data range
#minValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.025,names=FALSE)
#maxValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.975,names=FALSE)

minValue= -2
maxValue = 2

color_cuts <- seq(minValue,maxValue,length=100)
#color_cuts <- seq(-1,1,length=100)
color_cuts <- c(min(medianMatrix), color_cuts,max(medianMatrix)) # this catches the full dynamic range

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum) #this is something stupid in R that you always have to do

#Making png

#clusterPNGFile = paste(outputFolder,genome,'_',analysisName,'_2d_cluster.png',sep='')
png(filename = '/grail/projects/arber2/figures/enriched_chang_genes_heatmap_median_norm.png',width = 800,height =800)
layout(matrix(data=c(1,1,1,1,1,2,2),ncol= 7))


image(1:ncol(medianMatrix),1:nrow(medianMatrix),t(medianMatrix[expOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')


image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Log2 fold vs median")
dev.off()

write.table(rownames(medianMatrix[expOrder,]), file = "/grail/projects/arber2/tables/enriched_chang_genes_heatmap_order.txt", append = FALSE, quote = FALSE, sep = "\t", eol = '\n', na = "NA", dec = '.', row.names = FALSE, col.names = FALSE,
          qmethod = c("escape", "double"), fileEncoding = "")

########################################################################