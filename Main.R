###
  title: "LD.calcul.with.corPC1"
###
  
## R
##for R version R-3.6.0
PATH.data = paste0(getwd(),"/data")
PATH = paste0(getwd(),"")
source(paste0(PATH,'/cor.test.matrix.r'))
source(paste0(PATH,'/cone.proj.invsqrt.r'))
source(paste0(PATH,'/calcul.PC1.PC1v.r'))
source(paste0(PATH,'/LD.corPC1.corPC1v.r'))
source(paste0(PATH,'/GWESS.manhattan.plot.r'))
library(dplyr)
####################################EXAMPLE 1####################################
##Example 1 illustrates how to calculate LD with corPC1 and corPC1v statistics between a bait gene and two other genes
## SNP in line and id in columns 
geneA =read.csv(paste0(PATH.data,"/example.gene.1.txt"), header=T, row.names=1, sep=',')
geneB =read.csv(paste0(PATH.data,"/example.gene.2.txt"), header=T, row.names=1, sep=',')
geneC =read.csv(paste0(PATH.data,"/example.gene.3.txt"), header=T, row.names=1, sep=',')

kinship =read.csv(paste0(PATH.data,"/kinship.matrix.txt"), header=T, row.names=1, sep=',')
kinship.inv = cone.proj.sdp.invsqrt(kinship)

############CALCULATE corPC1, corPC1v
##gene: SNPs in lines and id in columns
PC1.A = calcul.PC1(geneA)
colnames(PC1.A)=("PC1.A")
PC1.Av = calcul.PC1v(geneA,kinship.inv)
colnames(PC1.Av)=("PC1.Av")
PC1.B = calcul.PC1(geneB)
colnames(PC1.B)=("PC1.B")
PC1.Bv = calcul.PC1v(geneB,kinship.inv)
colnames(PC1.Bv)=("PC1.Bv")
PC1.C = calcul.PC1(geneC)
colnames(PC1.C)=("PC1.C")
PC1.Cv = calcul.PC1v(geneC,kinship.inv)
colnames(PC1.Cv)=("PC1.Cv")

##create PC1 data
data.frame.PC1 = cbind(PC1.B,PC1.C)
data.frame.PC1v = cbind(PC1.Bv,PC1.Cv)
## Calculate LD
result = LD.corPC1.corPC1v (PC1.A,data.frame.PC1,PC1.Av,data.frame.PC1v)
result

####################################EXAMPLE 2####################################
##Example 2 illustrates how to calculate LD with corPC1 and corPC1v statistics between a bait gene and all genes in the M. truncatula genome and represents the results as a manhattan plot.
PC1.data = get(load(paste0(PATH.data,"/PC1_data.RData")))
PC1v.data = get(load(paste0(PATH.data,"/PC1v_data.RData")))
list.genes = read.csv(paste0(PATH.data,"/list.genes.txt"), header=T, sep=',')
##choose one bait gene
bait.MtSUNN = as.data.frame(PC1.data[which(rownames(PC1.data) == "Medtr4g070970"),])
bait.MtSUNNv = as.data.frame(PC1v.data[which(rownames(PC1v.data) == "Medtr4g070970"),])

##calculate p-value
#take few minutes depending on data
result = LD.corPC1.corPC1v (t(bait.MtSUNN),t(PC1.data),t(bait.MtSUNNv),t(PC1v.data))
result = as.data.frame(result)

##create file for manhattan plot
## columns: position.start,chr,pvalue
#corPC1
result.corPC1 = list.genes
result.corPC1$pvalue = result$LD.cor.pval[match(result.corPC1$ID_Mt4,rownames(result))]
result.corPC1 = result.corPC1[,c("position_start","chromosome_file","pvalue")]

position_bait_genes = list.genes$position_start[which(list.genes$ID_Mt4 == "Medtr4g070970")]
NB.CHR = 8
position_bait_genes = list.genes$position_start[which(list.genes$ID_Mt4 == "Medtr4g070970")]
chr_bait_gene = 4
statistic = "corPC1"
GWESS.manhattan.plot(result.corPC1,NB.CHR,position_bait_genes,chr_bait_gene,statistic)

#corPC1v
result.corPC1v = list.genes
result.corPC1v$pvalue = result$LD.corv.pval[match(result.corPC1v$ID_Mt4,rownames(result))]
result.corPC1v = result.corPC1v[,c("position_start","chromosome_file","pvalue")]

position_bait_genes = list.genes$position_start[which(list.genes$ID_Mt4 == "Medtr4g070970")]
NB.CHR = 8
position_bait_genes = list.genes$position_start[which(list.genes$ID_Mt4 == "Medtr4g070970")]
chr_bait_gene = 4
statistic = "corPC1v"
GWESS.manhattan.plot(result.corPC1v,NB.CHR,position_bait_genes,chr_bait_gene,statistic)
