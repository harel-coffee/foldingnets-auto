# foldingnets Figure 1

setwd("~/foldingnets/")


library(ggplot2)
library(reshape2)
library("ape")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg(file = "figures/Figure1/B_histogram.svg", width=5, height=4)

sf <- as.data.frame(read.table("data/fprocessed/ssb_sfam.txt", header=T))
sf <- sf[sort(sf$total, decreasing=T, index.return=T)$ix,]
sf$sf <- factor(sf$sf, levels=sf$sf)
sf <- sf[,-2]
sf_m <- melt(sf)

ggplot(sf_m, aes(x=sf, y=value, fill=variable)) + 
  geom_col(position="dodge2") + 
  theme_classic() + 
  scale_fill_manual(values=c("#29ABE2", "red")) + 
  labs(x="SCOP superfamily", y="Number of proteins") +
  theme(
    axis.text = element_text(size=24),
    axis.text.x = element_text(size=24, angle = 90, hjust = 1), 
    axis.ticks.x=element_blank(),
    axis.line.x = element_blank(), 
    text = element_text(size=24), 
    legend.position = c(0.8, 0.8), 
    legend.title = element_blank(),
    legend.text = element_text(size=24)
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# dendrogram

t <- as.data.frame( read.table("data/processed/Q_DF.txt", header=T) )
names <- t[,1]
names <- do.call(paste, c(t[,c(1,2,3)], sep="_"))

data <- as.matrix( t[, 4:ncol(t)] )
row.names(data) <- names
colnames(data) <- names

d <- as.dist( 1 - data)
hc <- hclust(d)


postscript("figures/Figure1/C_tree.ps", width=6, height=20, onefile=F, paper="special", horizontal=F)

plot(as.phylo(hc), cex = 1.2, label.offset = 0.1)

dev.off()




# secondary structure diagrams in order of structural phylogeny

ss <- read.table("data/processed/dssp_c371.txt", header=T, sep='\t')
ss.names <- ss[,1]
ss.data <- as.matrix( ss[,-1] )
ss.order <- hc$order

plot_ss <- function(data){
  
  xmax <- ncol(data)	
  n <- nrow(data)
  
  par(mar=c(4, 6, 2, 1))
  plot(-10, -10, xlim=c(0, xmax), ylim=c(0,n+1), xlab="", ylab="", main="", axes=F )
  
  axis(2, at=c(1:n), labels=ss.names[ss.order], las=2, tick=F, line=0)
  axis(1, at=c(0, 100, 200), cex.axis=1.2)
  
  for (j in 1:n){	
    current_data <- as.matrix(data[j,])
    current_data <- current_data[!is.na(current_data)]
    
    data.c <- current_data
    data.c[data.c==1] <- "#66AA29"
    data.c[data.c==2] <- "#222222"
    data.c[data.c==0] <- "white"
    
    lines(x=c(1,length(current_data)), y=c(j, j))
    
    for (i in 1:length(current_data)){		
      if (current_data[i] > 0 ) {
        polygon(x=c(i, i+1, i+1, i), y=c(j-0.2, j-0.2, j+0.2, j+0.2), col=data.c[i], border=data.c[i])
      }	
    }	
  }	
}

postscript("figures/Figure1/C_c371_ss.ps", width=15, height=22, onefile=F, paper="special", horizontal=F)

plot_ss(ss.data[ss.order,])

legend(250, 6, legend="beta sheet", pch=15, col="#222222", bty='n', pt.cex=1.4, xpd=T)
legend(250, 5, legend="alpha helix", pch=15, col="#66AA29", bty='n', pt.cex=1.4, xpd=T)

dev.off()


