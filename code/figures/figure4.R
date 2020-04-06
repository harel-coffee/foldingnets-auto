# foldingnets Figure 4

library(ggplot2)
library(reshape2)
library(pROC)
library(cowplot)

setwd("~/foldingnets")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

peaks.pred1 <- as.data.frame(read.table("data/processed/result_peaks_hyd05_R1/predict.txt", header=T))
pp1 <- melt( peaks.pred2[,c(1,7,8,9)] , id="N")
pp1$R <- rep("R1", nrow(pp1))

peaks.pred2 <- as.data.frame(read.table("data/processed/result_peaks_hyd05_R2/predict.txt", header=T))
pp2 <- melt( peaks.pred2[,c(1,7,8,9)] , id="N")
pp2$R <- rep("R2", nrow(pp2))

peaks.pred3 <- as.data.frame(read.table("data/processed/result_peaks_hyd05_R3/predict.txt", header=T))
pp3 <- melt( peaks.pred3[,c(1,7,8,9)] , id="N")
pp3$R <- rep("R3", nrow(pp3))

peaks.pred <- rbind(pp1, pp2, pp3)


svg(file = "figures/Figure4/B_peaks_auc.svg", height=4, width=5.5)

ggplot(peaks.pred) + 
  geom_line(aes(x=N, y=value, color=variable, linetype=factor(R)) ) + 
  scale_x_reverse() + 
  labs(x="Number of clusters", y="Discriminative power (AUC)") +
  scale_colour_viridis_d(option = "viridis") +
  scale_y_continuous(limits=c(0.5, 0.8)) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=20),
    text = element_text(size=20), 
    legend.title = element_blank()
    ) 

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary.1 <- as.data.frame(read.table('data/processed/result_peaks_hyd05_R1/clustering_summary.txt', header=T))
summary.2 <- as.data.frame(read.table('data/processed/result_peaks_hyd05_R2/clustering_summary.txt', header=T))
summary.3 <- as.data.frame(read.table('data/processed/result_peaks_hyd05_R3/clustering_summary.txt', header=T))

summary.size.1 <- as.data.frame(summary.1[,c(1,2,3)])
summary.size.1$R <- rep("R1", nrow(summary.1))
summary.size.2 <- as.data.frame(summary.2[,c(1,2,3)])
summary.size.2$R <- rep("R2", nrow(summary.2))
summary.size.3 <- as.data.frame(summary.3[,c(1,2,3)])
summary.size.3$R <- rep("R3", nrow(summary.3))
summary.size <- rbind(summary.size.1, summary.size.2, summary.size.3)

summary.ent.1 <- as.data.frame(summary.1[,c(1,6,7)])
summary.ent.1$R <- rep("R1", nrow(summary.1))
summary.ent.2 <- as.data.frame(summary.2[,c(1,6,7)])
summary.ent.2$R <- rep("R2", nrow(summary.2))
summary.ent.3 <- as.data.frame(summary.3[,c(1,6,7)])
summary.ent.3$R <- rep("R3", nrow(summary.3))
summary.ent <- rbind(summary.ent.1, summary.ent.2, summary.ent.3)

summary.hyd.1 <- as.data.frame(summary.1[,c(1,10,11)])
summary.hyd.1$R <- rep("R1", nrow(summary.1))
summary.hyd.2 <- as.data.frame(summary.2[,c(1,10,11)])
summary.hyd.2$R <- rep("R2", nrow(summary.2))
summary.hyd.3 <- as.data.frame(summary.3[,c(1,10,11)])
summary.hyd.3$R <- rep("R3", nrow(summary.3))
summary.hyd <- rbind(summary.hyd.1, summary.hyd.2, summary.hyd.3)



p.size <- ggplot(summary.size, aes(x=Iter)) + 
    geom_line(aes(y=size_mean, color=factor(R)) ) + 
    geom_ribbon(aes(ymin=size_mean-size_sd, ymax=size_mean+size_sd, fill=factor(R)), alpha=0.3) + 
    geom_line(data=summary.1, inherit.aes=F, aes(x=Iter, y=size_best3), color="#222222", linetype=1, size=0.8) + 
    geom_line(data=summary.2, inherit.aes=F, aes(x=Iter, y=size_best3), color="#222222", linetype=2, size=0.8) + 
    geom_line(data=summary.3, inherit.aes=F, aes(x=Iter, y=size_best3), color="#222222", linetype=3, size=0.8) + 
    geom_line(data=summary.1, inherit.aes=F, aes(x=Iter, y=size_best3_hyd), color="orange", linetype=1, size=0.8) + 
    geom_line(data=summary.2, inherit.aes=F, aes(x=Iter, y=size_best3_hyd), color="orange", linetype=2, size=0.8) + 
    geom_line(data=summary.3, inherit.aes=F, aes(x=Iter, y=size_best3_hyd), color="orange", linetype=3, size=0.8) + 
    scale_color_manual(values=c("grey20", "grey50", "grey70")) + 
    scale_fill_manual(values=c("grey20", "grey50", "grey70")) + 
    theme_classic() + 
    labs(x="", y="Size")+
    scale_y_log10() + 
    theme(
      text = element_text(size=16),
      axis.text = element_text(size=12), 
      legend.position = 'none'
    )

p.ent <- ggplot(summary.ent, aes(x=Iter)) + 
    geom_line(aes(y=ent_mean, color=factor(R)) ) + 
    geom_ribbon(aes(ymin=ent_mean-ent_sd, ymax=ent_mean+ent_sd, fill=factor(R)), alpha=0.3) + 
    geom_line(data=summary.1, inherit.aes=F, aes(x=Iter, y=ent_best3), color="#222222", linetype=1, size=0.8) + 
    geom_line(data=summary.2, inherit.aes=F, aes(x=Iter, y=ent_best3), color="#222222", linetype=2, size=0.8) + 
    geom_line(data=summary.3, inherit.aes=F, aes(x=Iter, y=ent_best3), color="#222222", linetype=3, size=0.8) + 
    geom_line(data=summary.1, inherit.aes=F, aes(x=Iter, y=ent_best3_hyd), color="orange", linetype=1, size=0.8) + 
    geom_line(data=summary.2, inherit.aes=F, aes(x=Iter, y=ent_best3_hyd), color="orange", linetype=2, size=0.8) + 
    geom_line(data=summary.3, inherit.aes=F, aes(x=Iter, y=ent_best3_hyd), color="orange", linetype=3, size=0.8) + 
    scale_color_manual(values=c("grey20", "grey50", "grey70")) + 
    scale_fill_manual(values=c("grey20", "grey50", "grey70")) + 
    theme_classic() +
    labs(x="", y="Entropy") +
    theme(
      text = element_text(size=16),
      axis.text = element_text(size=12), 
      legend.position = 'none'
    )

p.hyd <- ggplot(summary.hyd, aes(x=Iter)) + 
  geom_line(aes(y=hyd_mean, color=factor(R)) ) + 
  geom_ribbon(aes(ymin=hyd_mean-hyd_sd, ymax=hyd_mean+hyd_sd, fill=factor(R)), alpha=0.3) + 
  geom_line(data=summary.1, inherit.aes=F, aes(x=Iter, y=hyd_best3), color="#222222", linetype=1, size=0.8) + 
  geom_line(data=summary.2, inherit.aes=F, aes(x=Iter, y=hyd_best3), color="#222222", linetype=2, size=0.8) + 
  geom_line(data=summary.3, inherit.aes=F, aes(x=Iter, y=hyd_best3), color="#222222", linetype=3, size=0.8) + 
  geom_line(data=summary.1, inherit.aes=F, aes(x=Iter, y=hyd_best3_hyd), color="orange", linetype=1, size=0.8) + 
  geom_line(data=summary.2, inherit.aes=F, aes(x=Iter, y=hyd_best3_hyd), color="orange", linetype=2, size=0.8) + 
  geom_line(data=summary.3, inherit.aes=F, aes(x=Iter, y=hyd_best3_hyd), color="orange", linetype=3, size=0.8) + 
  scale_color_manual(values=c("grey20", "grey50", "grey70")) + 
  scale_fill_manual(values=c("grey20", "grey50", "grey70")) + 
  theme_classic() +
  labs(x="Clustering iterating", y="Hydro") +
  theme(
    text = element_text(size=16),
    axis.text = element_text(size=12), 
    legend.position = 'none'
  )

plot.joint <- plot_grid(p.size, p.ent, p.hyd, labels ="", ncol = 1, align = 'v')


ggsave(filename = "figures/Figure4/C_summary.ps",
       plot = print(plot.joint),
       width=4, height=4,
       device = cairo_ps)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct.pred <- as.data.frame(read.table("data/processed/result_structures_hyd05/predict.txt", header=T))
struct.pred.hyd <- as.data.frame(read.table("data/processed/result_structures_hyd05/predict_hyd.txt", header=T))
tmp.all <- melt( struct.pred[,c(1,7,8,9)] , id="N")
tmp.all$feature = rep("all", nrow(tmp.all))
tmp.hyd <- melt( struct.pred.hyd[,c(1,7,8,9)] , id="N")
tmp.hyd$feature = rep("hyd", nrow(tmp.hyd))
structures.pred <- rbind(tmp.all, tmp.hyd)


svg(file = "figures/Figure4/D_structures_auc.svg", height=4, width=5.5)

ggplot(structures.pred) + 
  geom_line(aes(x=N, y=value, color=variable, linetype=factor(feature)), size=0.3 ) + 
  scale_x_reverse() + 
  labs(x="Number of clusters", y="Discriminative power (AUC)") +
  scale_colour_viridis_d(option = "viridis") +
  scale_y_continuous(limits=c(0.5, 0.9)) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=20),
    text = element_text(size=20), 
    legend.title = element_blank()
  ) 

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

roc1 <- read.table("data/processed/result_structures_hyd05/roc1.txt")
roc1 <- data.frame(substrate=roc1[,1], ftc=roc1[,2])
roc3 <- read.table("data/processed/result_structures_hyd05/roc3.txt")
roc3 <- data.frame(substrate=roc3[,1], ftc=rowSums(roc3[,2:4]))
roc5 <- read.table("data/processed/result_structures_hyd05/roc5.txt")
roc5 <- data.frame(substrate=roc5[,1], ftc=rowSums(roc5[,2:6]))

roc1h <- read.table("data/processed/result_structures_hyd05/roc1_hyd.txt")
roc1h <- data.frame(substrate=roc1h[,1], ftc=roc1h[,2])
roc3h <- read.table("data/processed/result_structures_hyd05/roc3_hyd.txt")
roc3h <- data.frame(substrate=roc3h[,1], ftc=rowSums(roc3h[,2:4]))
roc5h <- read.table("data/processed/result_structures_hyd05/roc5_hyd.txt")
roc5h <- data.frame(substrate=roc5h[,1], ftc=rowSums(roc5h[,2:6]))

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


roc.1 <- roc(roc1$substrate, roc1$ftc, direction="<")
roc.1.rev <- reverse_rocdf(roc.1)
roc.1.rev$class <- rep("AUC1", length(roc.1.rev$specificities))
roc.1.rev$weight <- rep(0.9, length(roc.1.rev$specificities))
roc.1.rev$feature <- rep("all", length(roc.1.rev$specificities))

roc.3 <- roc(roc3$substrate, roc3$ftc, direction="<")
roc.3.rev <- reverse_rocdf(roc.3)
roc.3.rev$class <- rep("AUC3", length(roc.3.rev$specificities))
roc.3.rev$weight <- rep(0.9, length(roc.3.rev$specificities))
roc.3.rev$feature <- rep("all", length(roc.3.rev$specificities))

roc.5 <- roc(roc5$substrate, roc5$ftc, direction="<")
roc.5.rev <- reverse_rocdf(roc.5)
roc.5.rev$class <- rep("AUC5", length(roc.5.rev$specificities))
roc.5.rev$weight <- rep(0.9, length(roc.5.rev$specificities))
roc.5.rev$feature <- rep("all", length(roc.5.rev$specificities))

roc.1h <- roc(roc1h$substrate, roc1h$ftc, direction="<")
roc.1h.rev <- reverse_rocdf(roc.1h)
roc.1h.rev$class <- rep("AUC1", length(roc.1h.rev$specificities))
roc.1h.rev$weight <- rep(0.9, length(roc.1h.rev$specificities))
roc.1h.rev$feature <- rep("hyd", length(roc.1h.rev$specificities))

roc.3h <- roc(roc3h$substrate, roc3h$ftc, direction="<")
roc.3h.rev <- reverse_rocdf(roc.3h)
roc.3h.rev$class <- rep("AUC3", length(roc.3h.rev$specificities))
roc.3h.rev$weight <- rep(0.9, length(roc.3h.rev$specificities))
roc.3h.rev$feature <- rep("hyd", length(roc.3h.rev$specificities))

roc.5h <- roc(roc5h$substrate, roc5h$ftc, direction="<")
roc.5h.rev <- reverse_rocdf(roc.5h)
roc.5h.rev$class <- rep("AUC5", length(roc.5h.rev$specificities))
roc.5h.rev$weight <- rep(0.9, length(roc.5h.rev$specificities))
roc.5h.rev$feature <- rep("hyd", length(roc.5h.rev$specificities))

roc_df <- rbind(roc.1.rev, roc.3.rev, roc.5.rev, roc.1h.rev, roc.3h.rev, roc.5h.rev)



svg(file = "figures/Figure4/F_roc.svg", height = 4, width = 4)

ggplot(roc_df, aes(x=specificities, y=sensitivities, color=class)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
  geom_line(aes(linetype=feature), size=1.2) + 
  theme_classic() + 
  scale_colour_viridis_d(option = "viridis") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(1, 0.75, 0.5, 0.25, 0)) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1)) + 
  theme(
    axis.text = element_text(size=18),
    text = element_text(size=20), 
    legend.position=c(0.85, 0.3), 
    legend.text = element_text(size=16),
    legend.title = element_blank()
  ) + 
  labs(x="Specificity", y="Sensitivity")

dev.off()



