#Load tables from NES output
#These tables have a column that was added with python, where a 1 signifies the gene set
#is a cell cycle, growth, or proliferation gene set and may be highlighted red
nes_1 = read.delim('/grail/projects/arber2/tables/neg_nes_red_column.txt')

nes_2 = read.delim('/grail/projects/arber2/tables/neg_pos_red_column.txt')



fdr_vector = c(nes_1$FDR.q.val,nes_2$FDR.q.val)
nes_vector = c(nes_1$NES,nes_2$NES)
red_vector = c(nes_1$RED,nes_2$RED)
red_points = intersect(which(red_vector == 1),which(fdr_vector<0.1))

#plot nes vs fdr
pdf(file='/grail/projects/arber2/figures/drug_tpo_nes.pdf',width = 6,height =6)
plot(nes_vector,fdr_vector,pch=16,col='grey',xlab='Normalized Enrichment Score',ylab='FDR Q value')
points(nes_vector[red_points],fdr_vector[red_points],col='red',pch=16,cex=1.5)
abline(h=0.1,lty=2)
dev.off()



