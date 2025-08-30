#test some script
setwd("D:/research/mens_cycle/demo")
load("./mc2u_test.rda")
mc2su=mc2su_test
library(SCpubr)
library(paletteer)
library(ggplot2)
Idents(mc2su)<-mc2su$disease;mc2su<-subset(mc2su,idents=c("other"),invert=T)
col_cluster0=paletteer_d("awtools::spalette",6)
sc=mc2su@meta.data;sc$x=mc2su@reductions$metacell@cell.embeddings[,1];sc$y=mc2su@reductions$metacell@cell.embeddings[,2]
p1 = ggplot(sc)+
  geom_point(aes(x = x,y = y,fill = cluster0),shape = 21,color = "grey20",size = 3,stroke = 0.05)+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank())+
  scale_fill_manual(values = col_cluster0)
ggsave("Cluster0_2d1.png",p1,width = 6,height = 5,dpi = 300) # save the image


  Idents(mc2su)<-mc2su$disease;mc2su<-subset(mc2su,idents=c("CE","healthy"))
  pdf("./Number_of_UMI_count_per_cell.pdf",height=5,width=5)
  h = hist(mc2su$nCount_RNA,ylab="Cell frequency",xlab="UMI count" ,
      breaks=seq(0, max(mc2su$nCount_RNA) + 100, by=100),col="#69AB32",xlim=c(0,20000),main="Number of UMI count per cell")
  abline(v=median(mc2su$nCount_RNA), col='red', lty=2)
  text(max(mc2su$nCount_RNA), max(h$count)/1, "mean=1004", pos=2, col='red')
  dev.off()
  pdf("./Number_of_detected_genes_per_cell.pdf",height=5,width=5)
  h = hist(mc2su$nFeature_RNA,ylab="Cell frequency",xlab="Number of detected genes" ,
       breaks=seq(0, max(mc2su$nFeature_RNA) + 50, by=50),col="#F0E356",xlim=c(0,6000),main="Number of detected genes per cell")
  abline(v=median(mc2su$nFeature_RNA), col='red', lty=2)
  text(max(mc2su$nFeature_RNA), max(h$count)/1, "mean=703", pos=2, col='red')
  dev.off()
