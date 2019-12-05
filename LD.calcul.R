---
  title: "LD.calcul"
---
  
## R
PATH.data = paste0(getwd(),"/data")
PATH.script = paste0(getwd(),"/scripts")
source(paste0(PATH.script,'/cor.test.matrix.r'))
source(paste0(PATH.script,'/cone.proj.invsqrt.r'))

geneA =read.csv(paste0(PATH.data,"/example.gene.1.txt"), header=T, row.names=1, sep=',')
geneB =read.csv(paste0(PATH.data,"/example.gene.2.txt"), header=T, row.names=1, sep=',')

kinship =read.csv(paste0(PATH.data,"/kinship.matrix.txt"), header=T, row.names=1, sep=',')
kinship.inv = cone.proj.sdp.invsqrt(kinship)



## Calcul PCA

if (nrow(geneA) != 0) {
  # calcul PCA
  geneA.transpo <- t(geneA)
  #PC1 and PC1 % kinship.inv
  acp.A <- prcomp(geneA.transpo ,scale=FALSE)
  PC1.gene.A <- as.data.frame(acp.A$x[,1])
  PC1.gene.A.cor <- as.data.frame(kinship.inv%*%(acp.A$x[,1]))
}

if (nrow(geneB) != 0) {
  # calcul PCA
  geneB.transpo <- t(geneB)
  #PC1 and PC1 % kinship.inv
  acp.B <- prcomp(geneB.transpo ,scale=FALSE)
  PC1.gene.B <- as.data.frame(acp.B$x[,1])
  PC1.gene.B.cor <- as.data.frame(kinship.inv%*%(acp.B$x[,1]))
}



## Calcul LD

LD_corPC1 = cor(PC1.gene.A,PC1.gene.B)
LD_corPC1.pval = cor.test.p2(PC1.gene.A,PC1.gene.B)
LD_corPC1v = cor(PC1.gene.A.cor,PC1.gene.B.cor)
LD_corPC1v.pval = cor.test.p2(PC1.gene.A.cor,PC1.gene.B.cor)

result = c(LD_corPC1,LD_corPC1.pval,LD_corPC1v,LD_corPC1v.pval)
result
