rm(list=ls(all=TRUE)) 

a1 <- read.table("matriz_G3_Doble.csv", header = TRUE, sep=",", row.names = 1, stringsAsFactors=FALSE)
a2 <- read.table("matriz_G4_Doble.csv", header = TRUE, sep=",", row.names = 1, stringsAsFactors=FALSE)

colnames(a1) <- gsub("á","a", colnames(a1), fixed = T)
colnames(a1) <- gsub("ó","o", colnames(a1), fixed = T)
colnames(a2) <- gsub("á","a", colnames(a2), fixed = T)
colnames(a2) <- gsub("ó","o", colnames(a2), fixed = T)

rownames(a1) <- trimws(rownames(a1))
rownames(a2) <- trimws(rownames(a2))

library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(tidyverse)

# PCA G3 

res.pca_G3 <- PCA(a1, quali.sup=1)
PC1 <- res.pca_G3$ind$coord[,1] - res.pca_G3$ind$coord["WT",1]
PC2 <- res.pca_G3$ind$coord[,2] - res.pca_G3$ind$coord["WT",2]
labs <- rownames(res.pca_G3$ind$coord)
PCs <- data.frame(PC1=PC1, PC2=PC2)
rownames(PCs) <- labs

vPC1 <- res.pca_G3$var$coord[,1]
vPC2 <- res.pca_G3$var$coord[,2]
vlabs <- rownames(res.pca_G3$var$coord)
vPCs <- data.frame(PC1=vPC1, PC2=vPC2)
rownames(vPCs) <- vlabs

PCcat <- PCs
PCcat$TYPE <- "DOUBLE"
PCcat$TYPE[rownames(PCcat)=="WT"] <- "WT"
cols <- c("DOUBLE"="royalblue", "WT"="orange", "repel"=NA)
vPCs2 <- vPCs * 5

PCcat2 <- PCcat
PCcat2$mutante <- rownames(PCcat2)
PCcat2$label <- ifelse(
  (PCcat2$TYPE == "DOUBLE") & 
    (PCcat2$PC1 > 1 | PCcat2$PC2 > 1 | PCcat2$PC1 <= -1 | PCcat2$PC2 <= -1), 
  PCcat2$mutante, "")

p_G3 <- ggplot() + theme_bw(base_size = 20) +
  geom_point(data=PCcat2[PCcat2$TYPE=="DOUBLE",], aes(x=PC1, y=PC2, colour="DOUBLE"), size=4) +
  geom_text_repel(data=vPCs2, aes(x=PC1, y=PC2, label=rownames(vPCs2)), size=5, colour="black", segment.color=NA) +
  geom_text_repel(data=PCcat2[PCcat2$label != "",], 
                  aes(x=PC1, y=PC2, label=label), size=4.7, colour="royalblue", segment.color=NA) +
  geom_segment(data=vPCs2, aes(x=0, y=0, xend=PC1, yend=PC2), 
               arrow = arrow(length = unit(1/2, "picas")), color = "grey30") +
  geom_point(data=PCcat2[PCcat2$TYPE=="WT",], aes(x=PC1, y=PC2, colour="WT"), size=4) +
  scale_colour_manual(name="Mutants", values=cols) +
  scale_y_continuous(limits=c(-5, 5)) +
  ylab("PC2_WTcentered") + xlab("PC1_WTcentered") + 
  labs(title="PCA on G3 double mutants") +
  theme(plot.title = element_text(hjust = 0.5))
p_G3
ggsave("PCA_G3_Doble.png", dpi=300, width = 9, units = "in")

# PCA G4 (invertir el eje Y para poder compararla con G3)

res.pca_G4 <- PCA(a2, quali.sup=1)
PC1 <- res.pca_G4$ind$coord[,1] - res.pca_G4$ind$coord["WT",1]
PC2 <- -(res.pca_G4$ind$coord[,2] - res.pca_G4$ind$coord["WT",2])
labs <- rownames(res.pca_G4$ind$coord)
PCs <- data.frame(PC1=PC1, PC2=PC2)
rownames(PCs) <- labs

vPC1 <- res.pca_G4$var$coord[,1] 
vPC2 <- -(res.pca_G4$var$coord[,2]) 
vlabs <- rownames(res.pca_G4$var$coord)
vPCs <- data.frame(PC1=vPC1, PC2=vPC2)
rownames(vPCs) <- vlabs

PCcat <- PCs
PCcat$TYPE <- "DOUBLE"
PCcat$TYPE[rownames(PCcat)=="WT"] <- "WT"
cols <- c("DOUBLE"="royalblue", "WT"="orange", "repel"=NA)
vPCs2 <- vPCs * 5

PCcat2 <- PCcat
PCcat2$mutante <- rownames(PCcat2)
PCcat2$label <- ifelse(
  (PCcat2$TYPE == "DOUBLE") & 
    (PCcat2$PC1 > 1 | PCcat2$PC2 > 1 | PCcat2$PC1 <= -1 | PCcat2$PC2 <= -1), 
  PCcat2$mutante, "")

p_G4 <- ggplot() + theme_bw(base_size = 20) +
  geom_point(data=PCcat2[PCcat2$TYPE=="DOUBLE",], aes(x=PC1, y=PC2, colour="DOUBLE"), size=4) +
  geom_text_repel(data=vPCs2, aes(x=PC1, y=PC2, label=rownames(vPCs2)), size=5, colour="black", segment.color=NA) +
  geom_text_repel(data=PCcat2[PCcat2$label != "",], 
                  aes(x=PC1, y=PC2, label=label), size=4.7, colour="royalblue", segment.color=NA) +
  geom_segment(data=vPCs2, aes(x=0, y=0, xend=PC1, yend=PC2), 
               arrow = arrow(length = unit(1/2, "picas")), color = "grey30")+
  geom_point(data=PCcat2[PCcat2$TYPE=="WT",], aes(x=PC1, y=PC2, colour="WT"), size=4) +
  scale_colour_manual(name="Mutants", values=cols) +
  scale_y_continuous(limits=c(-5, 5)) +
  ylab("PC2_WTcentered") + xlab("PC1_WTcentered") + 
  labs(title="PCA on G4 double mutants") +
  theme(plot.title = element_text(hjust = 0.5))
p_G4
ggsave("PCA_G4_Doble.png", dpi=300, width = 9, units = "in")

# PCA comparación G3 y G4

res.pca_G3 <- PCA(a1, quali.sup=1)
G3_PC1 <- res.pca_G3$ind$coord[,1] - res.pca_G3$ind$coord["WT",1]
G3_PC2 <- res.pca_G3$ind$coord[,2] - res.pca_G3$ind$coord["WT",2]
G3_labs <- rownames(res.pca_G3$ind$coord)
G3_PCs <- data.frame(PC1=G3_PC1, PC2=G3_PC2)
G3_PCs$Grupo <- "G3"
G3_PCs$mutante <- rownames(G3_PCs)

res.pca_G4 <- PCA(a2, quali.sup=1)
G4_PC1 <- res.pca_G4$ind$coord[,1] - res.pca_G4$ind$coord["WT",1]
G4_PC2 <- -(res.pca_G4$ind$coord[,2] - res.pca_G4$ind$coord["WT",2])  # Multiplicamos por -1 para invertir el eje Y
G4_labs <- rownames(res.pca_G4$ind$coord)
G4_PCs <- data.frame(PC1=G4_PC1, PC2=G4_PC2)
G4_PCs$Grupo <- "G4"
G4_PCs$mutante <- rownames(G4_PCs)

PCs_ambos <- rbind(G3_PCs, G4_PCs)

mutantes_comunes <- intersect(G3_PCs$mutante, G4_PCs$mutante)
mutantes_comunes <- setdiff(mutantes_comunes, "WT") # Elimina WT

PCs_ambos_filtrado <- rbind(
  G3_PCs[G3_PCs$mutante %in% mutantes_comunes, ],
  G4_PCs[G4_PCs$mutante %in% mutantes_comunes, ]
)

PCs_ambos2 <- PCs_ambos_filtrado[
  (PCs_ambos_filtrado$PC1 > 1 | PCs_ambos_filtrado$PC2 > 1 |
     PCs_ambos_filtrado$PC1 <= -1 | PCs_ambos_filtrado$PC2 <= -1) &
    PCs_ambos_filtrado$Grupo == "G3", ]

ggplot() + theme_bw(base_size = 20) +
  geom_point(data=PCs_ambos_filtrado, aes(x=PC1, y=PC2, colour = Grupo), size=4, alpha=.6) +
  geom_line(data=PCs_ambos_filtrado, aes(x=PC1, y=PC2, group = mutante), color="black") +
  geom_text_repel(data=PCs_ambos2, aes(x=PC1, y=PC2, label=mutante), 
                  size=4.7, colour="royalblue", segment.color=NA) +
  scale_colour_manual(values=c("G3"="pink", "G4"="cyan3")) +
  scale_y_continuous(limits=c(-5, 5)) +
  ylab("PC2_WTcentered") + xlab("PC1_WTcentered") +
  labs(title="PCA difference between G3 and G4 double mutants") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("PCA_diff_Doble.png", dpi=300, width = 9, units = "in")

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

ggsave("DistanciasG3-G4_Dobles.png", dpi=300, width = 9, units = "in")
