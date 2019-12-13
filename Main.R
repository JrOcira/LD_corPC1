---
  title: "LD.calcul"
---
  
## R
PATH.data = paste0(getwd(),"/data")
PATH = paste0(getwd(),"/")
source(paste0(PATH,'/cor.test.matrix.r'))
source(paste0(PATH,'/cone.proj.invsqrt.r'))
source(paste0(PATH,'/calcul.PC1.PC1v.r'))
source(paste0(PATH,'/LD.corPC1.corPC1v.r'))
source(paste0(PATH,'/GWESS.manhattant.plot.r'))

####################################EXAMPLE 1####################################
## SNP in line and id in columns 
geneA =read.csv(paste0(PATH.data,"/example.gene.1.txt"), header=T, row.names=1, sep=',')
geneB =read.csv(paste0(PATH.data,"/example.gene.2.txt"), header=T, row.names=1, sep=',')
geneC =read.csv(paste0(PATH.data,"/example.gene.3.txt"), header=T, row.names=1, sep=',')

kinship =read.csv(paste0(PATH.data,"/kinship.matrix.txt"), header=T, row.names=1, sep=',')
kinship.inv = cone.proj.sdp.invsqrt(kinship)

############CALCUL corPC1, corPC1v
##gene: SNPs in line and id in columns
PC1.A = calcul.PC1(geneA)
PC1.Av = calcul.PC1v(geneA,kinship.inv)
PC1.B = calcul.PC1(geneB)
PC1.Bv = calcul.PC1v(geneB,kinship.inv)
PC1.C = calcul.PC1(geneC)
PC1.Cv = calcul.PC1v(geneC,kinship.inv)

##creat PC1 data
data.frame.PC1 = cbind(PC1.B,PC1.C)
data.frame.PC1v = cbind(PC1.Bv,PC1.Cv)
## Calcul LD
result = LD.corPC1.corPC1v (PC1.A,data.frame.PC1,PC1.Av,data.frame.PC1v)
result

####################################EXAMPLE 2####################################
PC1.data = get(load(paste0(PATH,"/PC1_data.RData")))
PC1v.data = get(load(paste0(PATH,"/PC1v_data.RData")))
list.genes = read.csv(paste0(PATH,"/list.genes.txt"), header=T, sep=',')
##choose one bait gene
bait.MtSUNN = as.data.frame(PC1.data[which(rownames(PC1.data) == "Medtr4g070970"),])
bait.MtSUNNv = as.data.frame(PC1v.data[which(rownames(PC1v.data) == "Medtr4g070970"),])

##calcul p-value
#take few minutes depending on data
result = LD.corPC1.corPC1v (t(bait.MtSUNN),t(PC1.data),t(bait.MtSUNNv),t(PC1v.data))
result = as.data.frame(result)

##creat file manhattan plot
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
GWESS.manhattant.plot(result.corPC1,NB.CHR,position_bait_genes,chr_bait_gene,statistic)

#corPC1v
result.corPC1v = list.genes
result.corPC1v$pvalue = result$LD.corv.pval[match(result.corPC1v$ID_Mt4,rownames(result))]
result.corPC1v = result.corPC1v[,c("position_start","chromosome_file","pvalue")]

position_bait_genes = list.genes$position_start[which(list.genes$ID_Mt4 == "Medtr4g070970")]
NB.CHR = 8
position_bait_genes = list.genes$position_start[which(list.genes$ID_Mt4 == "Medtr4g070970")]
chr_bait_gene = 4
statistic = "corPC1v"
GWESS.manhattant.plot(result.corPC1v,NB.CHR,position_bait_genes,chr_bait_gene,statistic)
