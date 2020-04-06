# foldingnets Figure 2

setwd('~/foldingnets')

library(ggplot2)
library(ggExtra)
library(matrixStats)
library(ggpubr)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# structural alignment


data <- read.table("data/processed/structaln.txt", header=T)
data <- as.data.frame(data)
data$pcaln = (data$Nfit / data$Aling) 


postscript("figures/Figure2/B_structaln.ps", width=4, height=4, horizontal=T, paper="special", onefile=F)

p <- ggplot(data, aes(x=SeqID, y=RMS)) + 
  geom_point(aes(color=pcaln)) + 
  geom_point(aes(size=pcaln, color=pcaln), show.legend = FALSE) + 
  labs(x="% sequence ID", y="RMSD") +
  theme_classic() + 
  scale_color_gradient(low="#222222", high="#29ABE2") + 
  theme( 
    axis.text=element_text(size=15), 
    axis.title=element_text(size=16), 
    legend.position = c(0.8, 0.5)
    )


ggExtra::ggMarginal(
  p,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = '#222222',
  fill = '#555555' ##ffa500' #'#29ABE2'
)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



postscript("figures/Figure2/C_alignment.ps", width=6, height=3, paper="special", horizontal=T, onefile=F)

par(mar=c(4,6,2,1))

aln <- read.table("data/processed/c3718_aln.txt", sep='\t') 
names <- aln[,1]
aln <- aln[,-1]
nr <- nrow(aln)
nc <- ncol(aln)

plot(-1, -1, xlim=c(0, nc), ylim=c(0, nr), xlab="", ylab="", main="", axes=F)
axis(1, line=0, cex.axis=0.8)
axis(2, at=c(1:16)+0.5, labels=names, las=2, tick=F, cex.axis=0.8)

polar <- "#5caaf1"
hydro <- "#f4b042"

for (i in 1:nc){
  for (j in 1:nr){
    if (j <= 11){offset <- 0}
    else {offset <- 0.5}
    
    current_aa = aln[j,i]
    if (current_aa == '-'){ colo <- "white"}		
    if (current_aa == 'I'){ colo <- hydro}
    if (current_aa == 'V'){ colo <- hydro}
    if (current_aa == 'L'){ colo <- hydro}
    if (current_aa == 'F'){ colo <- hydro}
    if (current_aa == 'C'){ colo <- hydro}
    if (current_aa == 'M'){ colo <- hydro}
    if (current_aa == 'A'){ colo <- hydro}
    if (current_aa == 'G'){ colo <- hydro}
    if (current_aa == 'D'){ colo <- "black"}
    if (current_aa == 'E'){ colo <- "black"}
    if (current_aa == 'H'){ colo <- "black"}
    if (current_aa == 'K'){ colo <- "black"}
    if (current_aa == 'R'){ colo <- "black"}
    if (current_aa == 'T'){ colo <- polar}
    if (current_aa == 'S'){ colo <- polar}
    if (current_aa == 'W'){ colo <- polar}
    if (current_aa == 'Y'){ colo <- polar}
    if (current_aa == 'P'){ colo <- polar}
    if (current_aa == 'Q'){ colo <- polar}
    if (current_aa == 'N'){ colo <- polar}
    
    polygon(x=c(i, i+1, i+1, i), y=c(j, j, j+1, j+1), col=colo, border=F)
  }
}

abline(h=6, lwd=3, col="white")
legend(440, 5, legend=c("hydrophobic", "polar", "charged"), pch=15, col=c("#f4b042", "#5caaf1", "black"), bty='n', xpd=T, cex=0.8)
mtext("Alignment position", side=1, line=2, cex=0.8)

dev.off()



# ASA profile, use ggsave and cairo engine to save transparency to postscript! 

par(mar=c(6,6,2,1))

asa <- read.table("data/processed/c3718.asa", sep='\t')
asa <- as.matrix(asa)

asa.mean <- colMeans(asa, na.rm=T)
asa.std <- colSds(asa, na.rm=T)
asa.std[is.na(asa.std)] <- 0

asa2 <- data.frame(idx=c(1:length(asa.mean)), avrg=asa.mean, stnd=asa.std)

p <- ggplot(asa2, aes(x=idx, y=avrg)) + geom_ribbon(aes(ymin=avrg-stnd, ymax=avrg+stnd), alpha=0.3) + geom_line() + theme_classic() + theme(axis.text.y = element_text(size=12) ) + ylab("ASA") + scale_color_manual(values=c("#666666")) + scale_y_continuous(breaks=c(0,100, 200))

ggsave(filename = "figures/Figure2/asa_alignment.ps",
       plot = print(p),
       width=6, height=1,
       device = cairo_ps)













#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D+E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# properties of struct align
hyd_buried_aligned <- as.matrix( read.table("data/processed/fig2_hyd_buried_alinged.txt", sep='\t') )
hyd_buried_notaln <- as.matrix( read.table("data/processed/fig2_hyd_buried_notalinged.txt", sep='\t') )
hyd_exposed_aligned <- as.matrix( read.table("data/processed/fig2_hyd_exposed_alinged.txt", sep='\t') )
hyd_exposed_notaln <- as.matrix( read.table("data/processed/fig2_hyd_exposed_notalinged.txt", sep='\t') )

agg_buried_aligned <- as.matrix( read.table("data/processed/fig2_agg_buried_alinged.txt", sep='\t') )
agg_buried_notaln <- as.matrix( read.table("data/processed/fig2_agg_buried_notalinged.txt", sep='\t') )
agg_exposed_aligned <- as.matrix( read.table("data/processed/fig2_agg_exposed_alinged.txt", sep='\t') )
agg_exposed_notaln <- as.matrix( read.table("data/processed/fig2_agg_exposed_notalinged.txt", sep='\t') )





postscript("figures/Figure2/D_core_hydrophobicity.ps", height=3, width=2.5, paper="special", horizontal=T, onefile=F)

p <- hyd_buried_aligned[,1] == 1
n <- hyd_buried_aligned[,1] == 0

hyd_df <- data.frame(hyd=c(as.numeric(hyd_buried_aligned[p,-1]), as.numeric(hyd_buried_aligned[n,-1])), substr=c(rep(1, length(as.numeric(hyd_buried_aligned[p,-1]))), rep(0, length(as.numeric(hyd_buried_aligned[n,-1])))) )
hyd_df <- hyd_df[!is.na(hyd_df[,1]),]
hyd_df$substr <- as.factor(hyd_df$substr)

ggplot(hyd_df, aes(x=substr, y=hyd, fill=substr)) + 
  geom_violin(trim=F) + 
  labs(x="", y="Hydrophobicity") +
  scale_fill_manual(values=c("#FF0000", "#29ABE2")) + 
  theme_classic() + 
  theme(
    text = element_text(size=16), 
    axis.text = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.position="none"
  ) + 
  scale_x_discrete(breaks=c("0","1"), labels=c("NI", "I") )

dev.off()




postscript("figures/Figure2/E_core_agg.ps", height=3, width=2.5, paper="special", horizontal=T, onefile=F)

p <- agg_buried_aligned[,1] == 1
n <- agg_buried_aligned[,1] == 0

agg_df <- data.frame(agg=c(as.numeric(agg_buried_aligned[p,-1]), as.numeric(agg_buried_aligned[n,-1])), substr=c(rep(1, length(as.numeric(agg_buried_aligned[p,-1]))), rep(0, length(as.numeric(agg_buried_aligned[n,-1])))) )
agg_df <- agg_df[!is.na(agg_df[,1]),]
agg_df$substr <- as.factor(agg_df$substr)

ggplot(agg_df, aes(x=substr, y=agg, fill=substr)) + 
  geom_violin(trim=F) + 
  labs(x="", y="Agg. propensity") +
  scale_fill_manual(values=c("#FF0000", "#29ABE2")) + 
  theme_classic() + 
  theme(
    text = element_text(size=16), 
    axis.text = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.position="none"
    ) + 
  scale_x_discrete(breaks=c("0","1"), labels=c("NI", "I") )

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# contact map

postscript("figures/Figure2/F_contactmap.ps", width=5, height=4, paper="special", horizontal=T, onefile=F)

cm <- read.table("data/processed/consensus_cm_df.txt", sep='\t', header=T)
cm$n <- as.factor(cm$n)

ggplot(cm, aes(x=i, y=j)) + 
  geom_point(aes(color=n), shape=20) + 
  theme_classic() + 
  scale_colour_viridis_d(option = "plasma") + 
  xlab("Alignment position") + 
  ylab("Alignment position") + 
  theme(
    text=element_text(size=20),
    axis.text=element_text(size=20),  
    legend.position = 'right' )

dev.off()


postscript("figures/Figure2/F_contactmap_hist.ps", width=2.5, height=4, paper="special", horizontal=T, onefile=F)

cm2 <- read.table("data/processed/consensus_cm_df.txt", sep='\t', header=T)
ggplot(cm2, aes(x=n)) + 
  geom_histogram(binwidth=1, color="white", fill="darkblue") + 
  theme_classic() + 
  xlab("Contact") + 
  ylab("Count") +
  theme(
    text=element_text(size=20),
    axis.text=element_text(size=20)  
    )

dev.off()


