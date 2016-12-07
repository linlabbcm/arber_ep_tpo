
#=====================================================================
#===============================DATA INPUT============================
#=====================================================================


expTable = read.delim('/grail/projects/arber2/cufflinks/arber_rna_cuffnorm/output/arber_rna_all_fpkm_exprs_raw.txt', header = TRUE, sep ='\t')

diffGenes = read.delim('/grail/projects/arber2/clustergram_genes_1.5_fold_change.txt')

diffGenesList = as.character(diffGenes[,1])

diffGenesList = sort(diffGenesList)

expMatrix = as.matrix(expTable)

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


#zero hour mean normalization

meanVector = apply(diffMatrix[,1:2],1,mean)

meanMatrix = log2(diffMatrix/meanVector)

#meanMatrix <- meanMatrix[,-1]

#===================================================================
#======================CLUSTERING EXPRESSION========================
#===================================================================

expCorDist = as.dist(1-cor(t(meanMatrix)))
expHClust = hclust(expCorDist)

expOrder = expHClust$order
#===================================================================
#=========================MAKING EXEMPLARS==========================
#===================================================================



enhancerMatrix = diffMatrix

#medianMatrix = diffMatrix

meanVector = apply(enhancerMatrix[,1:2],1,mean) 

meanMatrix = log2(enhancerMatrix/meanVector)

sampleDist = as.dist(1-cor(meanMatrix))
sampleHC = hclust(sampleDist)
sampleOrder = sampleHC$order

#distances by SEs
seDist = as.dist(1-cor(t(meanMatrix)))
seHC = hclust(seDist)
#plot(seHC)

seOrder = seHC$order

seClusterOrder = cutree(seHC,k=10)

exemplarMatrix = matrix(nrow=max(seClusterOrder),ncol = ncol(meanMatrix[,]))

#Plot cluster exemplar plots
pdf('/grail/projects/arber2/mean_normalized_exemplar_plots.pdf')
for(i in 1:max(seClusterOrder)){
  
  cluster = matrix(meanMatrix[which(seClusterOrder==i),],ncol=ncol(meanMatrix[,]))
  plot(1:ncol(meanMatrix),rep(0,ncol(meanMatrix[,])),ylim =quantile(cluster,c(0.025,0.975)),cex=0,
      xlab='Samples',xaxt='n',ylab='Control log2 row normalized mean signal',main = paste('Cluster',i))
  axis(1,1:ncol(meanMatrix),colnames(meanMatrix[,]),las=2)
  for(j in 1:nrow(cluster)){
    lines(1:ncol(meanMatrix),cluster[j,],col = rgb(0.0,0.5,1,0.1),lwd=1)
    
  }
  lines(1:ncol(meanMatrix),apply(cluster,2,mean),col='blue',lwd=4)
  legend(1,quantile(cluster,.90),paste("n =",nrow(cluster)))
  if(nrow(cluster) == 1){
    exemplarMatrix[i,] = cluster[1,]
  }else{
    exemplarMatrix[i,] = apply(cluster[,],2,mean)
  }
  print(nrow(cluster))
}

colnames(meanMatrix) = colnames(meanMatrix)
rownames(exemplarMatrix) = 1:max(seClusterOrder)
dev.off()