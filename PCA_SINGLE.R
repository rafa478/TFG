rm(list=ls(all=TRUE)) 

a1<- read.table("matriz_G3.csv", header = TRUE, sep=",", row.names = 1, stringsAsFactors=FALSE)
a2<- read.table("matriz_G4.csv", header = TRUE, sep=",", row.names = 1, stringsAsFactors=FALSE)

colnames(a1)<-gsub("á","a", colnames(a1),fixed = T)
colnames(a1)<-gsub("ó","o", colnames(a1),fixed = T)
colnames(a2)<-gsub("á","a", colnames(a2),fixed = T)
colnames(a2)<-gsub("ó","o", colnames(a2),fixed = T)

library(FactoMineR)
library(ggplot2)
library(ggrepel)

res.pca_G3 <- PCA(a1, quali.sup=1)

PC1 <- res.pca_G3$ind$coord[,1]-res.pca_G3$ind$coord["WT",1]
PC2 <- res.pca_G3$ind$coord[,2]-res.pca_G3$ind$coord["WT",2]
labs <- rownames(res.pca_G3$ind$coord)
PCs <- data.frame(cbind(PC1,PC2))
rownames(PCs) <- labs

vPC1 <- res.pca_G3$var$coord[,1]
vPC2 <- res.pca_G3$var$coord[,2]
vlabs <- rownames(res.pca_G3$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- colnames(PCs)

PCcat<-PCs
PCcat$TYPE <- "SINGLE"
PCcat$TYPE[rownames(PCcat)=="WT"] <- "WT"

cols <- c("DOUBLE"="grey80","SINGLE"="dark red","WT"="orange", "repel"=NA)
vPCs2<-vPCs*5

S<-subset(PCcat, TYPE=="SINGLE")
W<-subset(PCcat, TYPE=="WT")
PCcat2<-PCcat[!(tail(sort(PCcat$PC1),2)[1]==PCcat$PC1 | tail(sort(PCcat$PC1),2)[2]==PCcat$PC1),]
p_G3<-ggplot()+ theme_bw(base_size = 20) +
  geom_point(aes(x=PCcat2[PCcat2$TYPE=="DOUBLE",]$PC1,y=PCcat2[PCcat2$TYPE=="DOUBLE",]$PC2, colour="DOUBLE"),size=4)+
  geom_point(aes(x=PCcat2[PCcat2$TYPE=="SINGLE",]$PC1,y=PCcat2[PCcat2$TYPE=="SINGLE",]$PC2, colour="SINGLE"),size=4)+
  geom_text_repel(data=vPCs2, aes(x=vPCs2$PC1,y=vPCs2$PC2,label=rownames(vPCs2)), size=5, colour="black",segment.color=NA)+
  geom_text_repel(aes(x=PCcat2[(PCcat2$TYPE=="SINGLE" & (PCcat2$PC1>1 | PCcat2$PC2>1 | PCcat2$PC1 <= -1 | PCcat2$PC2<= -1)),]$PC1,y=PCcat2[(PCcat2$TYPE=="SINGLE" & (PCcat2$PC1>1 | PCcat2$PC2>1 | PCcat2$PC1 <= -1 | PCcat2$PC2<= -1)),]$PC2,label=rownames(PCcat2[(PCcat2$TYPE=="SINGLE" & (PCcat2$PC1>1 | PCcat2$PC2>1 | PCcat2$PC1 <= -1 | PCcat2$PC2<= -1)),])), size=4.7, colour="dark red",segment.color=NA) +
  geom_segment(data=vPCs2, aes(x = 0, y = 0, xend = vPCs2$PC1, yend = vPCs2$PC2), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30") +
  geom_point(aes(x=PCcat2[PCcat2$TYPE=="WT",]$PC1,y=PCcat2[PCcat2$TYPE=="WT",]$PC2, colour="WT"),size=4) +
  scale_colour_manual(name="Mutants",values=cols) +
  scale_y_continuous(limits=c(-5, 5)) +
  ylab("PC2_WTcentered") + xlab("PC1_WTcentered") + labs(title="PCA on G3 mutants") +
  theme(plot.title = element_text(hjust = 0.5))
p_G3

ggsave("PCA_G3.png", dpi=300, width = 9, units = "in")

library(FactoMineR)
library(ggplot2)
library(ggrepel)

res.pca_G4 <- PCA(a2, quali.sup=1)

PC1 <- res.pca_G4$ind$coord[,1]-res.pca_G4$ind$coord["WT",1]
PC2 <- res.pca_G4$ind$coord[,2]-res.pca_G4$ind$coord["WT",2]
labs <- rownames(res.pca_G4$ind$coord)
PCs <- data.frame(cbind(PC1,PC2))
rownames(PCs) <- labs

vPC1 <- res.pca_G4$var$coord[,1]
vPC2 <- res.pca_G4$var$coord[,2]
vlabs <- rownames(res.pca_G4$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- colnames(PCs)

PCcat<-PCs
PCcat$TYPE <- "SINGLE"
PCcat$TYPE[rownames(PCcat)=="WT"] <- "WT"

cols <- c("DOUBLE"="grey80","SINGLE"="dark red","WT"="orange", "repel"=NA)
vPCs2<-vPCs*5

S<-subset(PCcat, TYPE=="SINGLE")
W<-subset(PCcat, TYPE=="WT")
PCcat2<-PCcat[!(tail(sort(PCcat$PC1),2)[1]==PCcat$PC1 | tail(sort(PCcat$PC1),2)[2]==PCcat$PC1),]
p_G4<-ggplot()+ theme_bw(base_size = 20) +
  geom_point(aes(x=PCcat2[PCcat2$TYPE=="DOUBLE",]$PC1,y=PCcat2[PCcat2$TYPE=="DOUBLE",]$PC2, colour="DOUBLE"),size=4)+
  geom_point(aes(x=PCcat2[PCcat2$TYPE=="SINGLE",]$PC1,y=PCcat2[PCcat2$TYPE=="SINGLE",]$PC2, colour="SINGLE"),size=4)+
  geom_text_repel(data=vPCs2, aes(x=vPCs2$PC1,y=vPCs2$PC2,label=rownames(vPCs2)), size=5, colour="black",segment.color=NA)+
  geom_text_repel(aes(x=PCcat2[(PCcat2$TYPE=="SINGLE" & (PCcat2$PC1>1 | PCcat2$PC2>1 | PCcat2$PC1 <= -1 | PCcat2$PC2<= -1)),]$PC1,y=PCcat2[(PCcat2$TYPE=="SINGLE" & (PCcat2$PC1>1 | PCcat2$PC2>1 | PCcat2$PC1 <= -1 | PCcat2$PC2<= -1)),]$PC2,label=rownames(PCcat2[(PCcat2$TYPE=="SINGLE" & (PCcat2$PC1>1 | PCcat2$PC2>1 | PCcat2$PC1 <= -1 | PCcat2$PC2<= -1)),])), size=4.7, colour="dark red",segment.color=NA) +
  geom_segment(data=vPCs2, aes(x = 0, y = 0, xend = vPCs2$PC1, yend = vPCs2$PC2), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30") +
  geom_point(aes(x=PCcat2[PCcat2$TYPE=="WT",]$PC1,y=PCcat2[PCcat2$TYPE=="WT",]$PC2, colour="WT"),size=4) +
  scale_colour_manual(name="Mutants",values=cols) +
  scale_y_continuous(limits=c(-5, 5)) +
  ylab("PC2_WTcentered") + xlab("PC1_WTcentered") + labs(title="PCA on G4 mutants") +
  theme(plot.title = element_text(hjust = 0.5))
p_G4

ggsave("PCA_G4.png", dpi=300, width = 9, units = "in")

library(tidyverse)

res.pca_G3 <- PCA(a1, quali.sup=1)
G3_PC1 <- res.pca_G3$ind$coord[,1]-res.pca_G3$ind$coord["WT",1]
G3_PC2 <- res.pca_G3$ind$coord[,2]-res.pca_G3$ind$coord["WT",2]
G3_labs <- rownames(res.pca_G3$ind$coord)
G3_PCs <- data.frame(cbind(G3_PC1,G3_PC2))
G3_PCs$Grupo <- "G3"
colnames(G3_PCs)<-gsub("G3_", "", colnames(G3_PCs),fixed = T)
G3_PCs$mutante <- rownames(G3_PCs)

res.pca_G4 <- PCA(a2, quali.sup=1)
G4_PC1 <- res.pca_G4$ind$coord[,1]-res.pca_G4$ind$coord["WT",1]
G4_PC2 <- res.pca_G4$ind$coord[,2]-res.pca_G4$ind$coord["WT",2]
G4_labs <- rownames(res.pca_G4$ind$coord)
G4_PCs <- data.frame(cbind(G4_PC1,G4_PC2))
G4_PCs$Grupo <- "G4"
colnames(G4_PCs)<-gsub("G4_", "", colnames(G4_PCs),fixed = T)
G4_PCs$mutante <- rownames(G4_PCs)

PCs_ambos <- rbind(G3_PCs,G4_PCs)
PCs_ambos2 <- PCs_ambos[(PCs_ambos$PC1>1 | PCs_ambos$PC2>1 | PCs_ambos$PC1 <= -1 | PCs_ambos$PC2<= -1),] %>% filter(Grupo == "G3")

ggplot()+ theme_bw(base_size = 20) +
  geom_point(data=PCs_ambos, aes(x=PC1, y=PC2, colour = Grupo),size=4, alpha=.6)+
  geom_line(data=PCs_ambos, aes(x=PC1, y=PC2, group = mutante)) +
  geom_text_repel(data=PCs_ambos2,aes(x=PC1,y=PC2,
    label=mutante),
    size=4.7,
    colour="dark red",
    segment.color=NA)+
  scale_y_continuous(limits=c(-5, 5)) +
  ylab("PC2_WTcentered") + xlab("PC1_WTcentered") + labs(title="PCA difference between G3 and G4 mutants") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("PCA_diff.png", dpi=300, width = 9, units = "in")

mutantes_comunes <- intersect(G3_PCs$mutante, G4_PCs$mutante)
mutantes_comunes <- mutantes_comunes[mutantes_comunes != "WT"]

distancias <- sapply(mutantes_comunes, function(m) {
  x1 <- as.numeric(G3_PCs[G3_PCs$mutante == m, c("PC1", "PC2")])
  x2 <- as.numeric(G4_PCs[G4_PCs$mutante == m, c("PC1", "PC2")])
  sqrt(sum((x1 - x2)^2))
})

distancias_df <- data.frame(
  Mutante = mutantes_comunes,
  Distancia = distancias
)
distancias_df <- distancias_df[order(-distancias_df$Distancia),]

ggplot(distancias_df, aes(x = reorder(Mutante, Distancia), y = Distancia)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +
  coord_flip() +
  labs(
    title = "Distance between G3 and G4 KO",
    x = "",
    y = "G3-G4 distance"
  ) +
  geom_text(aes(label = round(Distancia, 2)), hjust = -0.2, size = 5) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

ggsave("DistanciasG3-G4.png", dpi=300, width = 9, units = "in")
