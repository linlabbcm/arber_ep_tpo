#setwd('~/Dropbox/arber/')

#=====================================================================
#===============================DATA INPUT============================
#=====================================================================

#load expression table
expTable = read.delim('/grail/projects/arber2/cufflinks/arber_rna_cuffnorm/output/arber_rna_all_fpkm_exprs_raw.txt', header = TRUE, sep ='\t')

#load list of differential genes determines 
diffGenes = read.delim('/grail/projects/arber2/tables/clustergram_genes_1.5_fold_change.txt')

diffGenesList = as.character(diffGenes[,1])

diffGenesList = sort(diffGenesList)
expMatrix = as.matrix(expTable)

#create matrix of differential gene expression
diffMatrix = matrix(nrow=length(diffGenesList),ncol=ncol(expMatrix))

rownames(diffMatrix) = diffGenesList
colnames(diffMatrix) = colnames(expMatrix)


for(i in 1:length(diffGenesList)){
  
  geneName = diffGenesList[i]
  expRow = which(rownames(expMatrix)==geneName)
  diffMatrix[i,] = expMatrix[expRow,]
  
}


#now we want to make a log2 row median normalized expression matrix

medianVector = apply(diffMatrix,1,median) 

medianMatrix = log2(diffMatrix/medianVector)


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
png(filename = '/grail/projects/arber2/Expression_Heatmap_clusters4_5.png',width = 800,height =800)
layout(matrix(data=c(1,1,1,1,1,2,2),ncol= 7))


image(1:ncol(medianMatrix),1:nrow(medianMatrix),t(medianMatrix[expOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')


image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Log2 fold vs median")
dev.off()


#===================================================================
#===========================DENDROGRAMS=============================
#===================================================================

#plot(expHClust)
clusterCut <- cutree(expHClust, 10)
#plot(expHClust, h = -.5)


clusterList <-c(clusterCut)


#===================================================================
#===========================ZOOM IN 4 & 5===========================
#===================================================================

#add list of clusters 1 thru 10 to diffMatrix
diffMatrix <- cbind(diffMatrix,clusterList)

#isolate clusters of interest
cluster2genes_rows <- which(diffMatrix[,9] == 2)
cluster4genes_rows <- which(diffMatrix[,9] == 4)
cluster1genes_rows <- which(diffMatrix[,9] == 1)
cluster5genes_rows <- which(diffMatrix[,9] == 5)


#isolate gene names
cluster2genes <- rownames(diffMatrix[cluster2genes_rows,])
cluster4genes <- rownames(diffMatrix[cluster4genes_rows,])
cluster1genes <- rownames(diffMatrix[cluster1genes_rows,])
cluster5genes <- rownames(diffMatrix[cluster5genes_rows,])

all_genes = c(cluster1genes_rows, cluster2genes_rows,cluster5genes_rows,cluster4genes_rows)

clusterTable = diffMatrix[all_genes,]

#creates a table of just selected clusters and expression
for(i in 1:length(clusterTable[,1])){
  if(clusterTable[i,9]==1){
    clusterTable[i,9]='C'
  }
  if(clusterTable[i,9]==2){
    clusterTable[i,9]='A'
  }
  if(clusterTable[i,9]==4){
    clusterTable[i,9]='B'
  }
  if(clusterTable[i,9]==5){
    clusterTable[i,9]='D'
  }
}

write.table(clusterTable,  file = "/grail/projects/arber2/cluster_A-D_members_and_expr.txt", append = FALSE, quote = FALSE, sep = "\t", eol = '\n', na = "NA", dec = '.', row.names = TRUE, col.names = TRUE,
            qmethod = c("escape", "double"), fileEncoding = "")

#create list output for cluster1 genes
fileConn<-file("/grail/projects/arber2/cluster_1_genes_list.txt")
writeLines(cluster1genes, fileConn)
close(fileConn)

# create list output for cluster2 genes
fileConn<-file("/grail/projects/arber2/cluster_2_genes_list.txt")
writeLines(cluster2genes, fileConn)
close(fileConn)

# create list output for cluster4 genes
fileConn<-file("/grail/projects/arber2/cluster_4_genes_list.txt")
writeLines(cluster4genes, fileConn)
close(fileConn)

#create list output for cluster7 genes
fileConn<-file("/grail/projects/arber2/cluster_5_genes_list.txt")
writeLines(cluster5genes, fileConn)
close(fileConn)
