LD.corPC1.corPC1v <- function(PC1.A,PC1.B,PC1.Av,PC1.Bv) {

LD_corPC1 = as.numeric(0)
LD_corPC1.pval = as.numeric(0)
LD_corPC1v = as.numeric(0)
LD_corPC1.vpval = as.numeric(0)


## Calcul LD

LD_corPC1 = cor(PC1.A,PC1.B)
LD_corPC1.pval = cor.test.p2(PC1.A,PC1.B)
LD_corPC1v = cor(PC1.Av,PC1.Bv)
LD_corPC1v.pval = cor.test.p2(PC1.Av,PC1.Bv)

res=t(rbind(LD_corPC1,
	LD_corPC1.pval,
	LD_corPC1v,
	LD_corPC1v.pval))
colnames(res)= c("LD.cor","LD.cor.pval","LD.corv","LD.corv.pval")
res
}
