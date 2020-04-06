# Figure 6 | SP | 4/2020


library(ggplot2)
library(reshape2)
library(rpart)

setwd("~/foldingnets/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


pred <- as.data.frame(read.table("data/processed/pred_tradeoffs.txt", header=T))
pred$feature <- factor(pred$feature, rev(pred$feature))
pred_m <- melt(pred)



svg(file = "figures/Figure6/A_predict.svg", width=5, height=3)

ggplot(pred_m, aes(x=feature, y=value*100, fill=variable)) + 
  geom_col(position="dodge2") + 
  theme_classic() + 
  scale_fill_manual(values=c("#222222", "orange")) + 
  labs(x="", y="% accuracy") +
  coord_flip() + 
  theme(
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20), 
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank(), 
    text = element_text(size=20), 
    legend.position = c(0.88, 0.44), 
    legend.title = element_blank(),
    legend.text = element_text(size=16)
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


data.all <- as.data.frame(read.table("data/processed/data_tradeoffs.txt", header=T))
data.hyd <- as.data.frame(read.table("data/processed/data_tradeoffs_hyd.txt", header=T))
tree.all <- rpart(substrate ~ cd + counts + nopt, data=data.all, method='class')
tree.all.nocounts <- rpart(substrate ~ cd + nopt, data=data.all, method='class')
tree.hyd <- rpart(substrate ~ cd + counts + nopt, data=data.hyd, method='class')

svg(file = "figures/Figure6/B_tree.svg", width=5, height=1.5)
par(mar=c(0,1,0,1))

par(mfrow=c(1,2))

plot(tree.all)
text(tree.all, pretty = 0)

plot(tree.all.nocounts)
text(tree.all.nocounts, pretty = 0)

dev.off()





