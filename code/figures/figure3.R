# Analysis and figure code for Figure 3 panels (July 2019)
# Updated March 2020
# contact density

library(ggplot2)
library(pROC)
library(MASS)
library(rpart)
library(reshape2)

setwd("~/foldingnets")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

svg(file = "figures/Figure3/A_corecontacts.svg", height = 4, width = 3.5)

cd <- as.data.frame( read.table("data/processed/c3718_structaln_contacts.txt", sep='\t', header=T) )
cd$substrate <- as.factor(cd$substrate)

cd <- data.frame(cd$ORF, cd$substrate, cd$mASA20, cd$mASA30, cd$mASA50)
colnames(cd) <- c("ORF", "substrate", "mASA<20", "mASA<30", "mASA<50")
cd_m <- melt(cd)

ggplot(cd_m, aes(x=variable, y=value, fill=substrate)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#FF0000", "#29ABE2")) + 
  labs(x="", y="Contact density") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size=18), 
    axis.title = element_text(size=20),
    axis.text.x = element_text(size=18, angle=35, hjust = 1), 
    legend.position="right",
    legend.text = element_text(size=16), 
    legend.title = element_blank()
    )
  
dev.off()


#wilcox.test( cd$mASA20[cd$substrate==0], cd$mASA20[cd$substrate==1])
#W = 47.5, p-value = 0.02683

#wilcox.test( cd$mASA30[cd$substrate==0], cd$mASA30[cd$substrate==1])
#W = 43, p-value = 0.08902

#wilcox.test( cd$mASA50[cd$substrate==0], cd$mASA50[cd$substrate==1])
#W = 41.5, p-value = 0.1259



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

deg <- as.matrix( read.table("data/processed/c3718_structaln_buried_degree.txt", sep='\t', header=T) )

n <- deg[,1] == 0
p <- deg[,1] == 1
nsubstr <- as.numeric(deg[n,-1])
substr  <- as.numeric(deg[p,-1])

mat <- matrix(NA, ncol=3,nrow=2)
	mat[1,] <- c( sum(nsubstr == 0), sum(nsubstr == 1), sum(nsubstr > 1) )
	mat[2,] <- c( sum(substr == 0),  sum(substr == 1),  sum(substr > 1) )

#fisher.test(mat[,c(1,2)])
#p-value = 0.02197

#fisher.test(mat[,c(1,3)])
#p-value = 0.001952

mat <- as.data.frame( mat / rowSums(mat) )
colnames(mat) <- c("n=0", "n=1", "n>1")
mat$class <- factor(c("NI", "I"), levels=c("NI", "I"))
mat2 <- melt(mat)



svg(file = "figures/Figure3/B_core_degree.svg", height = 4, width = 2.5)

ggplot(mat2, aes(x=variable, y=value, fill=class) ) + 
  geom_bar(stat="identity", position=position_dodge() ) +
  theme_classic() + 
  labs(x="Number of contacts", y="Frequency") + 
  scale_fill_manual(values=c("red", "#29ABE2")) +
  theme(
    axis.text = element_text(size=18), 
    axis.title = element_text(size=20), 
    legend.position = c(0.9, 0.6), 
    legend.text = element_text(size=16),
    legend.title = element_blank()
    ) 

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


svg(file = "figures/Figure3/C_contact_distances.svg", height = 4, width = 2)

cd <- as.data.frame( read.table("data/processed/c3718_structaln_contacts.txt", sep='\t', header=T) )
cd$substrate <- as.factor(cd$substrate)
cd_dists <- data.frame(ORF=cd$ORF, substrate=cd$substrate, dists=cd$dists)

ggplot(cd_dists, aes(x=factor(substrate), y=dists, fill=factor(substrate) )) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#FF0000", "#29ABE2")) + 
  scale_x_discrete(breaks=c("0","1"), labels=c("NI", "I") ) +
  labs(x="", y="Average distance (aa)") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size=18), 
    axis.title = element_text(size=20),
    axis.text.x = element_text(size=18), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = c(0.27, 0.8),
    legend.text = element_text(size=16), 
    legend.title = element_blank()
  )

dev.off()

#wilcox.test(cd$dists[cd$substrate==0], cd$dists[cd$substrate==1])
#W = 4, p-value = 0.005495


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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


c3718 <- as.data.frame( read.table("data/processed/c3718_structaln_contacts.txt", sep='\t', header=T) )
roc.c3718.cd <- roc(c3718$substrate, c3718$mASA20, direction=">")
roc.c3718.cd.rev <- reverse_rocdf(roc.c3718.cd)
roc.c3718.cd.rev$class <- rep("c.37.1.8", length(roc.c3718.cd.rev$specificities))
roc.c3718.cd.rev$weight <- rep(0.9, length(roc.c3718.cd.rev$specificities))
roc.c3718.cd.rev$feature <- rep("contact dens.", length(roc.c3718.cd.rev$specificities))

roc.c3718.dist <- roc(c3718$substrate, c3718$dists, direction="<")
roc.c3718.dist.rev <- reverse_rocdf(roc.c3718.dist)
roc.c3718.dist.rev$class <- rep("c.37.1.8", length(roc.c3718.dist.rev$specificities))
roc.c3718.dist.rev$weight <- rep(0.9, length(roc.c3718.dist.rev$specificities))
roc.c3718.dist.rev$feature <- rep("distance", length(roc.c3718.dist.rev$specificities))

all <- as.data.frame( read.table("data/processed/all_contact_density_normalized.txt", header=T) )
roc.all.cd <- roc(all$substrate, all$normCD, direction=">")
roc.all.cd.rev <- reverse_rocdf(roc.all.cd)
roc.all.cd.rev$class <- rep("all", length(roc.all.cd.rev$specificities)) 
roc.all.cd.rev$weight <- rep(0.8, length(roc.all.cd.rev$specificities)) 
roc.all.cd.rev$feature <- rep("contact dens.", length(roc.all.cd.rev$specificities))

roc.all.dist <- roc(all$substrate, all$buried_ldens, direction="<")
roc.all.dist.rev <- reverse_rocdf(roc.all.dist)
roc.all.dist.rev$class <- rep("all", length(roc.all.dist.rev$specificities))
roc.all.dist.rev$weight <- rep(0.9, length(roc.all.dist.rev$specificities))
roc.all.dist.rev$feature <- rep("distance", length(roc.all.dist.rev$specificities))

roc_df <- rbind(roc.c3718.cd.rev, roc.c3718.dist.rev, roc.all.cd.rev, roc.all.dist.rev)



svg(file = "figures/Figure3/D_roc.svg", height = 4, width = 4)

ggplot(roc_df, aes(x=specificities, y=sensitivities, color=class)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
  geom_line(size=1.5, aes(linetype=feature)) + 
  theme_classic() + 
  scale_color_manual(values = c(  "#555555", "orange") ) + 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(1, 0.75, 0.5, 0.25, 0)) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1)) + 
  theme(
    axis.text = element_text(size=18),
    text = element_text(size=20), 
    legend.position=c(0.8, 0.3), 
    legend.text = element_text(size=16),
    legend.title = element_blank()
  ) + 
  labs(x="Specificity", y="Sensitivity")

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sfstats <- as.data.frame( read.table("data/processed/superfamily_contact_density.txt", header=T))
sfclassify <- data.frame(sf=sfstats$sf, correct=rowSums(cbind(sfstats$TP, sfstats$TN)), false=rowSums(cbind(sfstats$FP, sfstats$FN)) )
sfclassify_m <- melt(sfclassify)

svg(file = "figures/Figure3/E_classify.svg", width=3, height=4)

ggplot(sfclassify_m, aes(x=sf, y=value, fill=variable)) + 
  geom_col(position="dodge2") + 
  theme_classic() + 
  scale_fill_manual(values=c("#222222", "#BBBBBB")) + 
  labs(x="SCOP SF", y="Num. classified") +
  coord_flip() + 
  theme(
    axis.text = element_text(size=18),
    axis.text.x = element_text(size=16, hjust = 1), 
    axis.ticks.x=element_blank(),
    axis.line.x = element_blank(), 
    text = element_text(size=20), 
    legend.position = c(0.8, 0.8), 
    legend.title = element_blank(),
    legend.text = element_text(size=16)
  )

dev.off()

