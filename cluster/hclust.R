#Adapted from https://uc-r.github.io/hc_clustering

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(data.table)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

# Define settings 
io <- list()
#metadata of the samples - added the germlayer level
io$metadata <-  "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/sample_metadata_bulk_pseudobulk3.txt"
#output directory
io$outdir <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/refilter/differential/Ectoderm_refilter/"

enhancers<-("id_met")
enhancers <- c(enhancers, read.csv(paste0(io$outdir,"/E7.5_Ect_Enh_edit.refiltered_E11_E16_P0_8w_2.bed"),stringsAsFactors=FALSE,header=FALSE,sep='\t')$V4)
enhancers <- c(enhancers, read.csv(paste0(io$outdir,"/E7.5_Mesend_Enh_edit.refiltered_E11_E16_P0_8w_2.bed"),stringsAsFactors=FALSE,header=FALSE,sep='\t')$V4)
#enhancers <- c(enhancers, read.csv(paste0(io$outdir,"/E7.5_Mes_Enh_edit.refiltered_E11_E16_P0.bed"),stringsAsFactors=FALSE,header=FALSE,sep='\t')$V4)
stats <- read.csv(paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Ect_E11_E16_P0_individual_enhancers.txt"),stringsAsFactors=FALSE,sep='\t') # %>% .[,mean:=mean*100]
stats2 <- read.csv(paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Mesend_E11_E16_P0_individual_enhancers.txt"),stringsAsFactors=FALSE,sep='\t') # %>% .[,mean:=mean*100]
#stats3 <- read.csv(paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Mes_E11_E16_P0_individual_enhancers.txt"),stringsAsFactors=FALSE,sep='\t') # %>% .[,mean:=mean*100]
#stats <- stats %>% merge(stats2, by="id_met") %>% merge(stats3, by="id_met") %>% unique()
stats <- stats %>% merge(stats2, by="id_met") %>% unique()


# remove samples that are not relavent for clustering 
exclude_samples<-c("Hon_2013_thymus_GRCm38_bismark_bt2.deduplicated.bismark","E5.5_Epiblast_pseudobulk","E4.5_Primitive_endoderm_pseudobulk","E4.5_Epiblast_pseudobulk","E6.5_Anterior_Primitive_Streak_pseudobulk","E3.5_ICM","GSE101478_Craniofacial_E10.5_cpg_meth","GSE82724_Craniofacial_E14.5_cpg_meth","GSE82812_Craniofacial_E13.5_cpg_meth","GSE82828_Craniofacial_E12.5_cpg_meth","GSE82907_Craniofacial_E15.5_cpg_meth","GSE83041_Craniofacial_E11.5_cpg_meth")
stats <- stats[!(stats$id_met %in% exclude_samples),]
#stats = subset(stats, select = -c(anno,anno.x,anno.y) )
#stats = subset(stats, select = -c(anno.x,anno.y) )
stats<- stats[c(enhancers)]

#read in metadata 
metadata <- read.csv(io$metadata,stringsAsFactors=FALSE,sep='\t') %>% .[,c("id_met","stage","tissue","germ_layer")] 
#assign colours to tissues
blood <- c("Blood_progenitors_1","Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3","Haematoendothelial_progenitors","Macrophage","macrophage","Treg","blood")
brain <- c("Rostral_neurectoderm","Forebrain_Midbrain_Hindbrain","Neural_crest","Spinal_cord","Hindbrain","Neural_tube","cerebellum","olfactory_bulb","Neural","Photoreceptors","CNU","glia","OLF","HPF","cerebrum","Midbrain","Isocortex","Forebrain","cortex","NPC")
skin <- setdiff(unique(metadata[metadata$germ_layer=="Ectoderm",]$tissue),brain)
endoderm <- unique(metadata[metadata$germ_layer=="Endoderm",]$tissue)
mesoderm <- setdiff(unique(metadata[metadata$germ_layer=="Mesoderm",]$tissue),blood)
metadata$Col<-NA
metadata$Col[metadata$tissue %in%blood] <- "Red"
metadata$Col[metadata$tissue %in% skin] <- "#5373B9"
metadata$Col[metadata$tissue %in% brain] <- "#619CFF"
metadata$Col[metadata$tissue %in% endoderm] <- "#2FBA37"
metadata$Col[metadata$tissue %in% mesoderm] <- "#F8766D"

metadata$seperate[metadata$tissue %in%blood] <- "Blood"
metadata$seperate[metadata$tissue %in% skin] <- "Skin"
metadata$seperate[metadata$tissue %in% brain] <- "Brain"
metadata$seperate[metadata$tissue %in% endoderm] <- "Endoderm"
metadata$seperate[metadata$tissue %in% mesoderm] <- "Mesoderm"

metadata <- metadata[!(metadata$id_met %in% exclude_samples),]
stats_new<- stats
stats<- stats[c(1:9,11:29,33,34,38,49,50,60,67:163,170:182,201:207),]

stats <- stats %>% merge(metadata, by="id_met")
stats$name2<- make.unique(paste(stats$stage,stats$tissue,sep="_"))
#names<-stats[,1]
#colors<-stats[,47]
stats <- stats[complete.cases(stats),]
#rename the rows by the sample name

rownames(stats)<-stats[,dim(stats)[2]]

#get only a numeric matrix
methylation<- stats[,2:(dim(stats)[2]-6)]
#methylation <- scale(methylation) # you dont neeed to scale in this instance 


# get the average methylation across all site 
methylation2<-rowSums(methylation)/nrow(methylation)

distSamples <- dist(methylation, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc <- hclust(distSamples, method = "complete" )

# generate a dendrogram
dend <- as.dendrogram(hc)

# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) <- stats$Col[order.dendrogram(dend)]
# Plotting the new dendrogram
par(mar=c(17,1,1,1))
plot(dend)
#library(dunn.test)
#dunn.test(methylation)

pdf(paste0(io$outdir,"/dendrogram_ect_mesend_E11_E16_P0_8w_2.pdf"), width=15, height=10)
par(mar=c(17,1,1,1))
plot(dend)
dev.off()



new_distances<- cor(t(methylation), method="pearson")
new_distances<-as.matrix(new_distances, labels=TRUE)

pdf(paste0(io$outdir,"/dendrogram_heatmap_ect_mesend_E11_E16_P0_8w_2.pdf"), width=17, height=10)

ha = HeatmapAnnotation(germ_layer = stats$seperate, col = list(germ_layer = c("Blood" = "Red", "Endoderm" = "#2FBA37", "Brain" = "#619CFF", "Skin"="#5373B9", "Mesoderm"="#F8766D")), show_legend = FALSE )
heatmap_plot= Heatmap(new_distances, cluster_columns = hc,cluster_rows=hc,column_dend_height = unit(8, "cm"),bottom_annotation=ha, show_heatmap_legend = FALSE,
       colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),show_row_dend=FALSE,show_row_names = FALSE) 

lgd = packLegend( Legend(title = "germ layer", at = c("Blood","Endoderm","Brain","Skin","Mesoderm"), legend_gp = gpar(fill = c("Red","#2FBA37","#619CFF","#5373B9","#F8766D"))),
  Legend(col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), title = "Pearson's correlation") )
draw(heatmap_plot, padding = unit(c(2, 2, 2, 20), "mm"))
draw(lgd, x = unit(0.95, "npc"), y = unit(0.8, "npc"))

dev.off()



stats_new <- stats_new %>% merge(metadata, by="id_met")
stats_new$name2<- make.unique(paste(stats_new$stage,stats_new$tissue,sep="_"))
#names<-stats[,1]
#colors<-stats[,47]
#stats_new <- stats_new[complete.cases(stats_new),]
#rename the rows by the sample name
rownames(stats_new)<-stats_new[,dim(stats)[2]]
#get only a numeric matrix
methylation_new<- stats_new[,2:(dim(stats_new)[2]-6)]

# calculate the euclidean distance of the new samples to the old reference samples 

new_distances<- -dist(methylation_new, method="euclidean")
new_distances<-as.matrix(new_distances, labels=TRUE)

new_samples<-c("E13.5_Fibrobast","E13.5_Fibrobast.1","2dR_1_TransDiffFibrobast", "2dR_1_TransDiffFibrobast.1",
"5dR_1_TransDiffFibrobast", "5dR_1_TransDiffFibrobast.1", "22dR_1_TransDiffFibrobast", "22dR_1_TransDiffFibrobast.1",
"5dR_3_TransDiffFibrobast", "5dR_3_TransDiffFibrobast.1", "22dR_3_TransDiffFibrobast", "22dR_3_TransDiffFibrobast.1",
"0dD_NPC", "0dD_NPC.1", "7dD_NPC", "7dD_NPC.1", "14dD_NPC", "21dD_NPC")

#stats_new<- stats_new[c(176:193),]
matrix_to_plot<-new_distances[c(new_samples),]
#matrix_to_plot<-matrix_to_plot[, stats$name2[order.dendrogram(dend)]]
matrix_to_plot<-matrix_to_plot[, stats$name2]

# mycolors <- c("Red","#2FBA37","#619CFF","#5373B9","#F8766D")
# names(mycolors) <- c("Blood","Endoderm","Brain","Skin","Mesoderm")
# mycolors <- list(mycolors = mycolors)


pdf(paste0(io$outdir,"/dendrogram_with_reprogramming_heatmap_euclidean.pdf"), width=25, height=10)

row_split=c(rep("1", 2),rep("2", 6),rep("3", 4),rep("4", 6))
ha = HeatmapAnnotation(germ_layer = stats$seperate, col = list(germ_layer = c("Blood" = "Red", "Endoderm" = "#2FBA37", "Brain" = "#619CFF", "Skin"="#5373B9", "Mesoderm"="#F8766D")), show_legend = FALSE )
heatmap_plot= Heatmap(matrix_to_plot, cluster_columns = hc,cluster_rows=FALSE,column_dend_height = unit(8, "cm"),row_split = row_split,bottom_annotation=ha, show_heatmap_legend = FALSE)

lgd = packLegend( Legend(title = "germ layer", at = c("Blood","Endoderm","Brain","Skin","Mesoderm"), legend_gp = gpar(fill = c("Red","#2FBA37","#619CFF","#5373B9","#F8766D"))),
  Legend(col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), title = "Similarity", at = c(0, 0.5, 1), labels = c("Low", "Median", "High")) )
draw(heatmap_plot, padding = unit(c(2, 2, 2, 20), "mm"))
draw(lgd, x = unit(0.9, "npc"), y = unit(0.8, "npc"))


dev.off()



# calculate the correlation of the new samples to the old reference samples 
new_distances<- cor(t(methylation_new), method="pearson")
new_distances<-as.matrix(new_distances, labels=TRUE)

new_samples<-c("E13.5_Fibrobast","E13.5_Fibrobast.1","2dR_1_TransDiffFibrobast", "2dR_1_TransDiffFibrobast.1",
"5dR_1_TransDiffFibrobast", "5dR_1_TransDiffFibrobast.1", "22dR_1_TransDiffFibrobast", "22dR_1_TransDiffFibrobast.1",
"5dR_3_TransDiffFibrobast", "5dR_3_TransDiffFibrobast.1", "22dR_3_TransDiffFibrobast", "22dR_3_TransDiffFibrobast.1",
"0dD_NPC", "0dD_NPC.1")

new_samples<-c("E13.5_Fibrobast","E13.5_Fibrobast.1","0dD_NPC", "0dD_NPC.1")

#stats_new<- stats_new[c(176:193),]
matrix_to_plot<-new_distances[c(new_samples),]
#matrix_to_plot<-matrix_to_plot[, stats$name2[order.dendrogram(dend)]]
matrix_to_plot<-matrix_to_plot[, stats$name2]

# mycolors <- c("Red","#2FBA37","#619CFF","#5373B9","#F8766D")
# names(mycolors) <- c("Blood","Endoderm","Brain","Skin","Mesoderm")
# mycolors <- list(mycolors = mycolors)


pdf(paste0(io$outdir,"/dendrogram_with_reprogramming_heatmap_pearsoncor.pdf"), width=25, height=10)


row_split=c(rep("1", 2),rep("2", 6),rep("3", 4),rep("4", 2))
ha = HeatmapAnnotation(germ_layer = stats$seperate, col = list(germ_layer = c("Blood" = "Red", "Endoderm" = "#2FBA37", "Brain" = "#619CFF", "Skin"="#5373B9", "Mesoderm"="#F8766D")), show_legend = FALSE )
heatmap_plot= Heatmap(matrix_to_plot, cluster_columns = hc,cluster_rows=FALSE,column_dend_height = unit(8, "cm"),row_split = row_split,bottom_annotation=ha, show_heatmap_legend = FALSE,
       colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

lgd = packLegend( Legend(title = "germ layer", at = c("Blood","Endoderm","Brain","Skin","Mesoderm"), legend_gp = gpar(fill = c("Red","#2FBA37","#619CFF","#5373B9","#F8766D"))),
  Legend(col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), title = "Pearson's correlation") )
draw(heatmap_plot, padding = unit(c(2, 2, 2, 20), "mm"))
draw(lgd, x = unit(0.95, "npc"), y = unit(0.8, "npc"))


dev.off()



# calculate the correlation of the new samples to the old reference samples 
new_distances<- cor(t(methylation_new), method="pearson",use="pairwise.complete.obs")
new_distances<-as.matrix(new_distances, labels=TRUE)

#new_samples<-c("E13.5_fibroblast", "day3_reprogramming","day6_reprogramming","day9_reprogramming","day12_reprogramming","0_iPSC","E5.5_Epi","E5.5_Epi.1","E5.5_Epi.2","E6.5_Epi","E6.5_Epi.1","E6.5_Epi.2","E5.5_Epi.3","E5.5_Epi.4","E6.5_Epi.3","E6.5_Epiblast","E6.5_Primitive_Streak","0_EpiSC")
new_samples<-c("E13.5_fibroblast", "day3_reprogramming","day6_reprogramming","day9_reprogramming","day12_reprogramming","0_iPSC","E6.5_Epiblast","E6.5_Primitive_Streak")

#stats_new<- stats_new[c(176:193),]
matrix_to_plot<-new_distances[c(new_samples),]
#matrix_to_plot<-matrix_to_plot[, stats$name2[order.dendrogram(dend)]]
matrix_to_plot<-matrix_to_plot[, stats$name2]

# mycolors <- c("Red","#2FBA37","#619CFF","#5373B9","#F8766D")
# names(mycolors) <- c("Blood","Endoderm","Brain","Skin","Mesoderm")
# mycolors <- list(mycolors = mycolors)


pdf(paste0(io$outdir,"/dendrogram_with_iPSC_reprogramming_heatmap_pearsoncor.pdf"), width=25, height=10)

ha = HeatmapAnnotation(germ_layer = stats$seperate, col = list(germ_layer = c("Blood" = "Red", "Endoderm" = "#2FBA37", "Brain" = "#619CFF", "Skin"="#5373B9", "Mesoderm"="#F8766D")), show_legend = FALSE )
heatmap_plot= Heatmap(matrix_to_plot, cluster_columns = hc,cluster_rows=FALSE,column_dend_height = unit(8, "cm"),bottom_annotation=ha, show_heatmap_legend = FALSE,
       colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

lgd = packLegend( Legend(title = "germ layer", at = c("Blood","Endoderm","Brain","Skin","Mesoderm"), legend_gp = gpar(fill = c("Red","#2FBA37","#619CFF","#5373B9","#F8766D"))),
  Legend(col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), title = "Pearson's correlation") )
draw(heatmap_plot, padding = unit(c(2, 2, 2, 20), "mm"))
draw(lgd, x = unit(0.95, "npc"), y = unit(0.8, "npc"))


dev.off()

## see if it is all the sites or just some that change 


reprogramming_samples<-c("E13.5_Fibrobast","2dR_1_TransDiffFibrobast", "2dR_1_TransDiffFibrobast.1",
"5dR_1_TransDiffFibrobast", "5dR_1_TransDiffFibrobast.1", "22dR_1_TransDiffFibrobast", "22dR_1_TransDiffFibrobast.1",
"5dR_3_TransDiffFibrobast", "5dR_3_TransDiffFibrobast.1", "22dR_3_TransDiffFibrobast", "22dR_3_TransDiffFibrobast.1","0dD_NPC")



sample1="E15.5_Limb"
sample2="E13_cortex"

pdf(paste0(io$outdir,"/correaltion_",sample1,"_",sample2,".pdf"), width=15, height=10)

par(mar=c(4,4,4,4))
par(mfrow = c(4, 6))

for (sample in reprogramming_samples){
  plot(c(methylation_new[sample1,]),c(methylation_new[sample,]), xlab="E15.5 Limb", ylab=sample, pch=19, col=c(rep("blue",29),rep("red",52)),ylim=c(0,100),xlim=c(0,100))
  #Ect_cor=paste("Ect=",signif(cor(unlist(methylation_new[sample1,1:29]), unlist(methylation_new[sample,1:29]), method = c("pearson")),2),sep="")
  #Mesend_cor=paste("Mesend=",signif(cor(unlist(methylation_new[sample1,30:81]), unlist(methylation_new[sample,30:81]), method = c("pearson")),2),sep="")
  #mtext(paste(Ect_cor,Mesend_cor,sep=","))
  mtext(paste("corr=",signif(cor(unlist(methylation_new[sample1,]), unlist(methylation_new[sample,]), method = c("pearson")),2),sep=""))
}
for (sample in reprogramming_samples){
  plot(c(methylation_new[sample2,]),c(methylation_new[sample,]), xlab=sample2, ylab=sample, pch=19,col=c(rep("blue",29),rep("red",52)),ylim=c(0,100),xlim=c(0,100))
  #Ect_cor=paste("Ect=",signif(cor(unlist(methylation_new[sample2,1:29]), unlist(methylation_new[sample,1:29]), method = c("pearson")),2),sep="")
  #Mesend_cor=paste("Mesend=",signif(cor(unlist(methylation_new["",30:81]), unlist(methylation_new[sample,30:81]), method = c("pearson")),2),sep="")
  #mtext(paste(Ect_cor,Mesend_cor,sep=","))
  mtext(paste("corr=",signif(cor(unlist(methylation_new[sample2,]), unlist(methylation_new[sample,]), method = c("pearson")),2),sep=""))
}
dev.off()


## plot a pca

       
train <-  methylation                        # training set
test  <- methylation_new                       # test set

pca <- prcomp(train, scale.=T, retx=T)        # pca on training set
summary(pca)$importance[3,1:4]                # 84% of variation in first 2 PC
# PC1     PC2     PC3     PC4 
# 0.60268 0.83581 0.89643 0.92139             
scores <- data.frame(pca$x)[1:2]              # so use first two PC
distSamples <- dist(scores, method = "euclidean") # hclust cluster analysis
# Hierarchical clustering using Complete Linkage
km <- hclust(distSamples, method = "complete" )
km <- cutree(km, k = 2)

pc.test  <- predict(pca,test)[,1:2]           # transform the test set

library(ggplot2)
# plot first two PC with cluster id
gg.train  <- data.frame(cluster=factor(km), scores)
centroids <- aggregate(cbind(PC1,PC2)~cluster,data=gg.train,mean)
gg.train  <- merge(gg.train,centroids,by="cluster",suffixes=c("",".centroid"))
# label cluster 1 as Mesendoderm and cluster 2 as ectoderm 
gg.train$germ_layer<-"NA"
gg.train$germ_layer[gg.train$cluster=="1"] <- "Mesendoderm"
gg.train$germ_layer[gg.train$cluster=="2"] <- "Ectoderm"


gg.test   <- data.frame(pc.test[,1:2])
gg.test   <- gg.test[c(new_samples),]
#specify which samples of the new samples are which colours
gg.test_col <-  c(rep("black",2),rep("orchid4",10),rep("grey",6))
#specify which samples are in the same group to draw lines between
gg.test$paired=c(1,1,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4)

pdf(paste0(io$outdir,"/pca_sample_clustering.pdf"), width=7, height=5)

# generate cluster plot...
cluster.plot <- ggplot(gg.train, aes(x=PC1, y=PC2, color=germ_layer)) +
  geom_point(size=2,alpha=0.7) +
  geom_point(data=gg.test,color=gg.test_col,size=3)+
  geom_line(data=gg.test, aes(group = paired,color="trajectory")) + 
  scale_color_manual(values = c("Mesendoderm"="firebrick","Ectoderm"="skyblue","trajectory"="black")) +
  scale_fill_manual( values = c("Mesendoderm"="firebrick","Ectoderm"="skyblue","trajectory"="black")) +
  #geom_text(data=gg.train,label=rownames(gg.train),color="purple")+
  theme_bw() +
  coord_fixed()

plot(cluster.plot)

dev.off()






