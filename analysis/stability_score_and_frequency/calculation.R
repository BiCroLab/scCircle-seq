library(ggplot2)
library(tidyverse)
library(colorspace)  
library(ggExtra)
HeLa <- readRDS('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_stability/HeLa/HeLa_frequency_jaccard_index.rds')
Colo320DM <- readRDS('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_stability/Colo320DM/Colo320DM_frequency_jaccard_index.rds')
PC3 <- readRDS('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_stability/PC3/PC3_frequency_jaccard_index.rds')
K562 <- readRDS('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_stability/K562/K562_frequency_jaccard_index.rds')
H293T <- readRDS('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_stability/293T/293T_frequency_jaccard_index.rds')


HeLa$jaccard_index <- HeLa$jaccard_index/0.8835216986
#HeLa, frequency>0.208
HeLa[which(HeLa$jaccard_index > 1),]$jaccard_index <- 1
HeLa$frequency <- HeLa$frequency/0.91666667
HeLa[which(HeLa$frequency > 1),]$frequency <- 1
HeLa$stability_score <- HeLa$frequency*HeLa$jaccard_index
HeLa <- HeLa[,c(1,2,3,6,7,8)]
HeLa$length <- HeLa$V3 - HeLa$V2
HeLa$length <- log(HeLa$length+1)
HeLa$type <- 'Low frequency unstable circle(LU)'
HeLa[which(HeLa$stability_score > 0.3 & HeLa$frequency <= 0.65),]$type <- 'Low frequency stable circle(LS)'
HeLa[which(HeLa$stability_score > 0.3 & HeLa$frequency > 0.65),]$type <- 'High frequency stable circle(HS)'
HeLa[which(HeLa$stability_score <= 0.3 & HeLa$frequency > 0.65),]$type <- 'High frequency unstable circle(HU)'
HeLa[which(HeLa$V1 == 'chrM'),]$type <- 'mtDNA'

# HeLa[order(HeLa$stability_score),]
# HeLa <- HeLa[which(!HeLa$V1 == 'chrM'),]

Colo320DM$jaccard_index <- Colo320DM$jaccard_index/9.410328e-01
#Colo320DM, frequency>0.102; jaccard_index > 
Colo320DM[which(Colo320DM$jaccard_index > 1),]$jaccard_index <- 1
Colo320DM$frequency <- Colo320DM$frequency/0.71428571
Colo320DM[which(Colo320DM$frequency > 1),]$frequency <- 1
Colo320DM$stability_score <- Colo320DM$frequency*Colo320DM$jaccard_index
Colo320DM <- Colo320DM[,c(1,2,3,6,7,8)]
Colo320DM$length <- Colo320DM$V3 - Colo320DM$V2
Colo320DM$length <- log(Colo320DM$length+1)
Colo320DM$type <- 'Low frequency unstable circle(LU)'
Colo320DM[which(Colo320DM$stability_score > 0.3 & Colo320DM$frequency <= 0.65),]$type <- 'Low frequency stable circle(LS)'
Colo320DM[which(Colo320DM$stability_score > 0.3 & Colo320DM$frequency > 0.65),]$type <- 'High frequency stable circle(HS)'
Colo320DM[which(Colo320DM$stability_score <= 0.3 & Colo320DM$frequency > 0.65),]$type <- 'High frequency unstable circle(HU)'
Colo320DM[which(Colo320DM$V1 == 'chrM'),]$type <- 'mtDNA'

# Colo320DM[order(Colo320DM$stability_score),]
# Colo320DM <- Colo320DM[which(!Colo320DM$V1 == 'chrM'),]

PC3$jaccard_index <- PC3$jaccard_index/0.8288986119
#PC3, frequency>0.217
PC3[which(PC3$jaccard_index > 1),]$jaccard_index <- 1
PC3$frequency <- PC3$frequency/0.86956522
PC3$stability_score <- PC3$frequency*PC3$jaccard_index
PC3 <- PC3[,c(1,2,3,6,7,8)]
PC3$length <- PC3$V3 - PC3$V2
PC3$length <- log(PC3$length+1)
PC3$type <- 'Low frequency unstable circle(LU)'
PC3[which(PC3$stability_score > 0.3 & PC3$frequency <= 0.65),]$type <- 'Low frequency stable circle(LS)'
PC3[which(PC3$stability_score > 0.3 & PC3$frequency > 0.65),]$type <- 'High frequency stable circle(HS)'
PC3[which(PC3$stability_score <= 0.3 & PC3$frequency > 0.65),]$type <- 'High frequency unstable circle(HU)'
PC3[which(PC3$V1 == 'chrM'),]$type <- 'mtDNA'
# PC3[order(PC3$stability_score),]
# PC3 <- PC3[which(!PC3$V1 == 'chrM'),]

K562$jaccard_index <- K562$jaccard_index/8.986955e-01
#K562, frequency>0.172
K562[which(K562$jaccard_index > 1),]$jaccard_index <- 1
K562$frequency <- K562$frequency/0.89655172
K562[which(K562$frequency > 1),]$frequency <- 1
K562$stability_score <- K562$frequency*K562$jaccard_index
K562 <- K562[,c(1,2,3,6,7,8)]
K562$length <- K562$V3 - K562$V2
K562$length <- log(K562$length+1)
K562$type <- 'Low frequency unstable circle(LU)'
K562[which(K562$stability_score > 0.3 & K562$frequency <= 0.65),]$type <- 'Low frequency stable circle(LS)'
K562[which(K562$stability_score > 0.3 & K562$frequency > 0.65),]$type <- 'High frequency stable circle(HS)'
K562[which(K562$stability_score <= 0.3 & K562$frequency > 0.65),]$type <- 'High frequency unstable circle(HU)'
K562[which(K562$V1 == 'chrM'),]$type <- 'mtDNA'
# K562[order(K562$stability_score),]
# K562 <- K562[which(!K562$V1 == 'chrM'),]

H293T$jaccard_index <- H293T$jaccard_index/9.650736e-01
#293T, frequency>0.192
H293T[which(H293T$jaccard_index > 1),]$jaccard_index <- 1
H293T$frequency <- H293T$frequency/1
H293T[which(H293T$frequency > 1),]$frequency <- 1
H293T$stability_score <- H293T$frequency*H293T$jaccard_index
H293T <- H293T[,c(1,2,3,6,7,8)]
H293T$length <- H293T$V3 - H293T$V2
H293T$length <- log(H293T$length+1)
H293T$type <- 'Low frequency unstable circle(LU)'
H293T[which(H293T$stability_score > 0.3 & H293T$frequency <= 0.65),]$type <- 'Low frequency stable circle(LS)'
H293T[which(H293T$stability_score > 0.3 & H293T$frequency > 0.65),]$type <- 'High frequency stable circle(HS)'
H293T[which(H293T$stability_score <= 0.3 & H293T$frequency > 0.65),]$type <- 'High frequency unstable circle(HU)'
H293T[which(H293T$V1 == 'chrM'),]$type <- 'mtDNA'
# H293T[order(H293T$stability_score),]
# H293T <- H293T[which(!H293T$V1 == 'chrM'),]


p_Colo320DM <- ggplot(data = Colo320DM,aes(x=frequency,y=stability_score,color = type))+
geom_point() + theme_classic()+ scale_colour_manual(values=c("#FC4E07","#56B4E9","#999999","#E69F00","#8642CC")) +
xlab("Frequency") + ylab("Stability score") + theme(legend.position = "none")  
ggsave('Colo320DM_Frequency-Stability.svg')

p_PC3 <- ggplot(data = PC3,aes(x=frequency,y=stability_score,color = type))+
geom_point() + theme_classic()+ scale_colour_manual(values=c("#FC4E07","#56B4E9","#999999","#E69F00","#8642CC")) +
xlab("Frequency") + ylab("Stability score") + theme(legend.position = "none")
ggsave('PC3_Frequency-Stability.svg')

p_HeLa <- ggplot(data = HeLa,aes(x=frequency,y=stability_score,color = type))+
geom_point() + theme_classic()+ scale_colour_manual(values=c("#FC4E07","#56B4E9","#999999","#E69F00","#8642CC")) +
xlab("Frequency") + ylab("Stability score") + theme(legend.position = "none") 
ggsave('HeLa_Frequency-Stability.svg')

p_K562 <- ggplot(data = K562,aes(x=frequency,y=stability_score,color = type))+
geom_point() + theme_classic()+ scale_colour_manual(values=c("#FC4E07","#56B4E9","#E69F00","#8642CC")) +
xlab("Frequency") + ylab("Stability score") + theme(legend.position = "none")
ggsave('K562_Frequency-Stability.svg')

p_H293T <- ggplot(data = H293T,aes(x=frequency,y=stability_score,color = type))+
geom_point() + theme_classic()+ scale_colour_manual(values=c("#56B4E9","#999999","#E69F00","#8642CC")) +
xlab("Frequency") + ylab("Stability score") + theme(legend.position = "none")
ggsave('293T_Frequency-Stability.svg')

ggplot(data = Colo320DM,aes(x=frequency,y=stability_score,color = type))+
geom_point() + theme_classic()+ scale_colour_manual(values=c("#FC4E07","#56B4E9","#999999","#E69F00","#8642CC")) +
xlab("Frequency") + ylab("Stability score")
ggsave('Label.svg')

pdf('Colo320DM_Frequency-Stability.pdf')
ggMarginal(p_Colo320DM, type = "density", groupColour = TRUE, groupFill = TRUE)
dev.off()

pdf('PC3_Frequency-Stability.pdf')
ggMarginal(p_PC3, type = "density", groupColour = TRUE, groupFill = TRUE)
dev.off()

pdf('K562_Frequency-Stability.pdf')
ggMarginal(p_K562, type = "density", groupColour = TRUE, groupFill = TRUE)
dev.off()

pdf('HeLa_Frequency-Stability.pdf')
ggMarginal(p_HeLa, type = "density", groupColour = TRUE, groupFill = TRUE)
dev.off()

pdf('293T_Frequency-Stability.pdf')
ggMarginal(p_H293T, type = "density", groupColour = TRUE, groupFill = TRUE)
dev.off()

write.table(Colo320DM,'Colo320DM_calculated.bed',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(PC3,'PC3_calculated.bed',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(HeLa,'HeLa_calculated.bed',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(K562,'K562_calculated.bed',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(H293T,'293T_calculated.bed',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
