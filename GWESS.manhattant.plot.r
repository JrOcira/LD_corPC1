GWESS.manhattant.plot <- function(bait,nb.chr,position_bait_genes,chr_bait_gene,statistic){
library(dplyr)
library(magrittr)
library(ggplot2)
library("RColorBrewer")
colnames(bait) = c("position_start","CHR","P_corPC1v")
bait$P_corPC1v = -log10(bait$P_corPC1v)
bait$position_start = as.numeric(as.character(bait$position_start))
bait$CHR = as.numeric(as.character(bait$CHR))

#replace Inf values by values max
bait$P_corPC1v[bait$P_corPC1v == Inf] <- NA
df2 <- bait
colourCount = nb.chr
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

cbPalette <- getPalette(colourCount)
don2 <- df2 %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(as.numeric(as.character(position_start)))) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df2, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(as.numeric(as.character(CHR)), as.numeric(as.character(position_start))) %>% #arrange() = rérodonne data frame en fonction des groupes
  mutate( BPcum=as.numeric(as.character(position_start))+tot) #mutate() = add new variable (here new column)

don2$CHR = as.numeric(as.character(don2$CHR))
axisdf2 = don2 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

colourCount = nb.chr
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cbPalette <- getPalette(colourCount)
q <- ggplot(don2, aes(x=BPcum, y=P_corPC1v)) +
	# Show all points
	geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
	scale_color_manual(values=cbPalette, nb.chr ) +
	# custom X axis:
	scale_x_continuous( label = as.numeric(as.character(axisdf2$CHR)), breaks= axisdf2$center ) +
	scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
	ylim(low=0, high=max(don2$P_corPC1v,na.rm=T)) +
	# Custom the theme:
	theme_bw() +
	theme( 
	  legend.position="none",
	  panel.border = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank()
	)  +xlab("chromosomes")+ylab("-log10(p-value)") + 
	geom_vline(xintercept = don2$BPcum[which(don2$CHR == chr_bait_gene & don2$position_start == position_bait_genes)], linetype="dotted", color = "black", size=0.5)


png(paste0(getwd(),"/manhattant.plot.",statistic,".png"), width = 1000, height = 600)
print(q) 
dev.off()



}
