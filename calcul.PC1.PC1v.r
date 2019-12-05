calcul.PC1.PC1v<-function(gene,kinship.inv){
  if (nrow(gene) != 0) {
  # calcul PCA
  gene.transpo <- t(gene)
  #PC1 and PC1 % kinship.inv
  acp.A <- prcomp(gene.transpo ,scale=FALSE)
  PC1.gene <- as.data.frame(acp.A$x[,1])
  PC1.gene.cor <- as.data.frame(kinship.inv%*%(acp.A$x[,1]))
}
}