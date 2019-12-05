cor.test.p2 <- function(x,z){
  ## x and z must have same length
  FUN1 <- function(x, y) cor.test(x, y,alternative="greater")[["p.value"]]
  FUN2 <- function(x, y) cor.test(x, y,alternative="less")[["p.value"]]
  y <- outer(
    colnames(x),
    colnames(z),
    Vectorize(function(i,j) ifelse(cor(x[,i], z[,j])>=0, FUN1(x[,i], z[,j]), FUN2(x[,i], z[,j])))
  )
  names(y) <- list(colnames(z))
  y
}
