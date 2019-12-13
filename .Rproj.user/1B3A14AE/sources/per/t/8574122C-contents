cone.proj.sdp.invsqrt<-function(matrix,cst=0.000001){
  #projection of the matrix on definite positive matrix cone and V-1/2
  mat.decomp<-eigen(matrix,symmetric=TRUE)
  valpp.mat<-mat.decomp$values
  valpp.mat[which(valpp.mat<cst)]<-cst  # transform too small or negative value
  valpp.mat.inv<-1/sqrt(valpp.mat)
  res<-mat.decomp$vectors %*% (diag(valpp.mat.inv)) %*% t(mat.decomp$vectors)
  colnames(res)<-colnames(matrix)
  rownames(res)<-rownames(matrix)
  res
}
