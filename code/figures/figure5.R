# foldingnets Figure 5


library(ggplot2)
library(reshape2)
library(pROC)
library(RColorBrewer)
library(viridis)
library(scales)
library(cowplot)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pax <- as.data.frame(read.table("data/processed/yeast_abundance.txt", header=T))
pax$substrate <- as.factor(pax$substrate)



svg(file = "figures/Figure5/A_abundance.svg", height = 4, width = 2.5)

ggplot(pax) + 
  geom_boxplot(aes(x=substrate, y=log(abundance), fill=substrate), outlier.shape = NA) + 
  theme_classic() + 
  scale_fill_manual( values = c("red", "#29ABE2") ) + 
  scale_x_discrete(breaks=c("0","1"), labels=c("NI", "I") ) +
  labs(x="", y="Protein abundance (log)") +
  theme(
    axis.text.y = element_text(size=20),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size=20), 
    legend.position = 'none'
  ) 
  
dev.off()



#wilcox.test(pax$abundance[pax$substrate==1], pax$abundance[nopt$substrate==0])
#Wilcoxon rank sum test with continuity correction
#data:  pax$abundance[pax$substrate == 1] and pax$abundance[nopt$substrate == 0]
#W = 1975, p-value = 0.0002181
#alternative hypothesis: true location shift is not equal to 0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nopt <- as.data.frame(read.table("data/processed/nopt_clusters.txt",header=T))


svg(file = "figures/Figure5/B_pctnopt.svg", height = 4, width = 2.5)

ggplot(nopt) + 
  geom_boxplot(aes(x=as.factor(substrate), y=nopt*100, fill=as.factor(substrate)), outlier.shape = NA) + 
  theme_classic() + 
  scale_fill_manual( values = c("red", "#29ABE2") ) + 
  scale_x_discrete(breaks=c("0","1"), labels=c("NI", "I") ) +
  labs(x="", y="% nonoptimal codons") +
  theme(
    axis.text.y = element_text(size=20),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size=20), 
    legend.position = 'none'
  ) 

dev.off()


#wilcox.test(nopt$nopt[nopt$substrate==1], nopt$nopt[nopt$substrate==0])
#Wilcoxon rank sum test with continuity correction
#data:  nopt$nopt[s] and nopt$nopt[-s]
#W = 2784.5, p-value = 0.3674
#alternative hypothesis: true location shift is not equal to 0



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rand_nopt  <- as.matrix( read.table("data/processed/rarecodon_rand_nopt_10000.txt") )
rand_Cnopt <- as.matrix( read.table("data/processed/rarecodon_rand_Cnopt_10000.txt") )


result.nopt <- c(rep(NA, 1000))
result.Cnopt <- c(rep(NA, 1000))
for (i in 1:1000){
  result.nopt[i] <- roc(nopt$substrate, rand_nopt[,i], direction=">")$auc
  result.Cnopt[i] <- roc(nopt$substrate, rand_Cnopt[,i], direction=">")$auc
}

roc.nopt <- roc(nopt$substrate, nopt$nopt, direction=">")
roc.Cnopt <- roc(nopt$substrate, nopt$Cnopt, direction=">")

p.nopt <- sum(roc.nopt$auc < result.nopt)/1000
p.Cnopt <- sum(roc.Cnopt$auc < result.Cnopt)/1000


result.nopt <- data.frame("AUC"=result.nopt)
result.Cnopt <- data.frame("AUC"=result.Cnopt)



svg(file = "figures/Figure5/C_rand_nopt.svg", height = 4, width = 3.5)

plot.nopt <- ggplot(result.nopt) + 
                geom_histogram(aes(x=AUC), binwidth=0.01, fill="grey20", color='white', alpha=0.5) + 
                scale_x_continuous( lim=c(0.35, 0.75) ) +
                scale_y_continuous( lim=c(0, 100) ) + 
                labs(x="AUC (% nopt. codons)", y="Count") +
                theme_classic() + 
                geom_segment(inherit.aes=F, aes(x=roc.nopt$auc, xend=roc.nopt$auc, y=25, yend=0), size=1.5, color="red", arrow=arrow(length = unit(0.02, "npc"), ends="last", type="closed") ) +
                #geom_hline(aes(yintercept = mean(Median)), size=0.5, linetype="twodash" ) + 
                theme(
                  text = element_text(size=16), 
                  legend.position = 'none'
                )
  
plot.Cnopt <- ggplot(result.Cnopt) + 
                  geom_histogram(aes(x=AUC), binwidth=0.01, fill="grey20", color='white', alpha=0.5) + 
                  scale_x_continuous( lim=c(0.35, 0.75) ) +
                  scale_y_continuous( lim=c(0, 100) ) + 
                  labs(x="AUC (% nopt. codon clusters)", y="Count") +
                  theme_classic() + 
                  geom_segment(inherit.aes=F, aes(x=roc.nopt$auc, xend=roc.nopt$auc, y=25, yend=0), size=1.5, color="red", arrow=arrow(length = unit(0.02, "npc"), ends="last", type="closed") ) +
                  #geom_hline(aes(yintercept = mean(Median)), size=0.5, linetype="twodash" ) + 
                  theme(
                    text = element_text(size=16), 
                    legend.position = 'none'
                  )

plot_grid(plot.nopt, plot.Cnopt, labels ="", ncol = 1, align = 'v')

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid <- as.data.frame(read.table("data/processed/rarecodon_grid_CC.txt", header=T))

result.grid <- c(rep(NA, ncol(grid)))
result.i <- c(rep(NA, ncol(grid)))
result.j <- c(rep(NA, ncol(grid)))
for (i in 1:ncol(grid)){
  result.grid[i] <- roc(nopt$substrate, grid[,i], direction=">")$auc
  tmp <- strsplit(colnames(grid)[i], '_')
  result.i[i] <- as.numeric(tmp[[1]][2])
  result.j[i] <- as.numeric(tmp[[1]][3])
}
result.grid2 <- data.frame("theta_cont"=result.i, "theta_dist"=result.j, "AUC"=result.grid)



svg(file = "figures/Figure5/D_heatmap.svg", height = 4, width = 5)

ggplot(result.grid2, aes(x = theta_cont, y = theta_dist)) + 
  geom_tile(aes(fill=AUC),color="white", size=0.1) + 
  labs(x="Minimum contacts", y="Distance") + 
  scale_y_continuous(breaks=c(5:14), labels=c(5:14)) + 
  scale_x_continuous(breaks=c(0:4), labels=c(1:5)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20)
  ) +
  scale_fill_viridis_c(option = "plasma")

dev.off()








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


svg(file = "figures/Figure5/E_hdclusters.svg", height=4, width=2.5)

ggplot(nopt) + 
  geom_boxplot(aes(x=factor(substrate), y=as.numeric(CC), fill=factor(substrate))) + 
  theme_classic() + 
  scale_fill_manual( values = c("red", "#29ABE2") ) + 
  scale_x_discrete(breaks=c("0","1"), labels=c("NI", "I") ) +
  scale_y_continuous(breaks=c(0,3,6,9), labels=c(0,3,6,9)) + 
  labs(x="", y="# hd clusters") +
  theme(
    axis.text.y = element_text(size=20),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size=20), 
    legend.position = 'none'
  ) 

dev.off()


#wilcox.test(nopt$CC[nopt$substrate==1], nopt$CC[nopt$substrate==0])
#W = 1097, p-value = 0.04593


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rand_CC <- as.matrix( read.table("data/processed/rarecodon_rand_CC_10000.txt") )
result.CC <- c(rep(NA, 1000))
for (i in 1:1000){
  result.CC[i] <- roc(nopt$substrate, rand_CC[,i], direction=">")$auc
}
roc.CC <- roc(nopt$substrate, nopt$CC, direction=">")
p.CC <- sum(roc.CC$auc < result.CC)/1000
result.CC <- data.frame("AUC"=result.CC)



plot.CC <- ggplot(result.CC) + 
              geom_histogram(aes(x=AUC), binwidth=0.01, fill="grey20", color='white', alpha=0.5) + 
              scale_x_continuous( lim=c(0.35, 0.75) ) +
              scale_y_continuous( lim=c(0, 100) ) + 
              labs(x="AUC (hd clusters)", y="Count") +
              theme_classic() + 
              geom_segment(inherit.aes=F, aes(x=roc.CC$auc, xend=roc.CC$auc, y=25, yend=0), size=1.5, color="red", arrow=arrow(length = unit(0.02, "npc"), ends="last", type="closed") ) +
              theme(
                text = element_text(size=16), 
                legend.position = 'none'
              )


reverse_rocdf <- function(roc){
  data <- data.frame(spec=roc$specificities, sens=roc$sensitivities)
  newspec <- c()
  newsens <- c()
  unique_spec <- unique(data$spec)
  for (i in 1:length(unique_spec)){
    sel <- data$spec == unique_spec[i]
    newspec <- c(newspec, rep(1-unique_spec[i], sum(sel)))
    newsens <- c(newsens, rev(data$sens[sel]) )
  }
  newdata <- data.frame(specificities=newspec, sensitivities=newsens)
}
roc_CC_rev <- reverse_rocdf(roc.CC)



plot.roc <- ggplot(roc_CC_rev, aes(x=specificities, y=sensitivities)) + 
              geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
              geom_line(color="red", size=2, show.legend = FALSE) + 
              theme_classic() + 
              scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(1, 0.75, 0.5, 0.25, 0)) + 
              scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1)) + 
              theme(
                axis.text = element_text(size=16),
                text = element_text(size=16)
              ) + 
              labs(x="Specificity", y="Sensitivity")


svg(file = "figures/Figure5/F_CC_ROC.svg", height = 4, width = 3.5)

plot_grid(plot.CC, plot.roc, labels ="", ncol = 1, align = 'v')

dev.off()






