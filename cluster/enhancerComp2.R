#after merging outputs from enhancerComp1 eg:
#awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $2 } END { for(i=1;i<=FNR;i++) print a[i] }' $(ls -1v E7.5_End_*) | tr " " "\t" | paste <(cut -f 1 E7.5_End_1_enhancer_meth_for_clustering.txt) -  > all_E7.5_End_enhancers.txt
#cat <( head -n1 Ect/all_E7.5_Ect_enhancers.txt ) <(tail -n +2 Ect/all_E7.5_Ect_enhancers.txt) <(tail -n +2 End/all_E7.5_End_enhancers.txt) <(tail -n +2 Mes/all_E7.5_Mes_enhancers.txt) > all_E7.5_lineage_specific_enhancers.txt


library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(data.table)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library("GenomicRanges")

args<-c("1","Ect")

# Define settings 
io <- list()
io$data_dir<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/"
#metadata of the samples - added the germlayer level
io$metadata <-  "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/sample_metadata_bulk_pseudobulk.txt"
#input the cpg methlation calls
#annotation directory whicha has the regions you would like to look at for
io$annos_dir  <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/"
#output directory
io$outdir <- paste("/bi/group/reik/Jannat/nmt/analysis/methacc/meth/meth_comparison/",args[2],"/",sep="")
io$all_enhancers <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/genomic_contexts/mm10_encode_cRES_mainonly.txt"


opts <- list()
# Define which cells to use (use NULL for all cells)
opts$stage_lineage <- c(
  "Mesoderm",
  "Ectoderm",
  "Endoderm",
  "Epiblast",
  "ICM"
)
#opts$stage_lineage <- NULL

# Define genomic contexts (use NULL for no genomic context filtering)
opts$annos <- c(
  "E7.5_Ect_Enh_edit",
  "E7.5_End_Enh_edit",
  "E7.5_Mes_Enh_edit"
)

# Define which cells to use ie the plate and position of the cells to use 
tmp <- fread(io$metadata) %>% .[!is.na(id_met)] 
if (is.null(opts$stage_lineage)) {
  opts$cells <- tmp[,id_met]
} else {
  opts$cells <- tmp %>%
    .[,stage_lineage:=paste(germ_layer)] %>%
    .[pass_metQC=="TRUE",id_met]
}

# Load sample metadata
metadata <- fread(io$metadata,stringsAsFactors=F) %>% 
  .[,c("id_met","stage","germ_layer","tissue")] %>%
  .[id_met%in%opts$cells] %>%
  .[,stage_lineage:=as.factor(paste(germ_layer,sep="_"))]

pseudobulk_notE8.5 <-  tmp %>% .[stage=="E7.5" | stage=="E6.5" | stage=="E5.5" | stage=="E4.5" | stage=="E3.5" | id_met=="101weeks_Mesoderm_pseudobulk"  | id_met=="77weeks_Mesoderm_pseudobulk"  | id_met=="36weeks_Mesoderm_pseudobulk",id_met]
pseudobulk_E8.5 <-  tmp %>% .[stage=="E8.5" | id_met=="Neural_Ectoderm_8weeks_pseudobulk" ,id_met]
cortes_morales <-  tmp %>% .[tissue =="Photoreceptors" | tissue =="Treg" | tissue =="Pancreas",id_met]
lister<-grep("Lister_2013",opts$cells, value=T)
luo<-grep("Luo",opts$cells, value=T)
Charlton<-grep("Charlton",opts$cells, value=T)
knaupp<-grep("Knaupp",opts$cells, value=T)
zhang<-grep("Zhang",opts$cells, value=T)
cortes_morales_lister<-c(cortes_morales, lister,luo,Charlton,knaupp,zhang)


merged_meth <- fread(paste0("/bi/group/reik/Jannat/nmt/analysis/methacc/meth/meth_comparison/all_E7.5_lineage_specific_enhancers.txt"),stringsAsFactors=FALSE, fill=TRUE,check.names=FALSE, header=TRUE) # %>% .[,mean:=mean*100]
merged_meth<-merged_meth[,c(1:6,8:44,46:49,51:61,63:70,72:82,83:103,140:195)]
#merged_meth<-merged_meth[,c(1:4,6,8,10,18,29,59,86,87,162,163,187:195,196:213)]


merged_meth_num <-  data.frame(apply(merged_meth[,2:ncol(merged_meth)], 2, function(x) as.numeric(as.character(x))))

correlation<-cor(merged_meth_num,use="pairwise.complete.obs",method="pearson")
correlation<-as.matrix(correlation, labels=TRUE)

names<- colnames(merged_meth_num) %>% gsub("MethylC\\.seq", "MethylC-seq", .) %>% gsub("MethylC\\.Seq", "MethylC-Seq", .) %>% gsub("X101weeks_Mesoderm_pseudobulk", "101weeks_Mesoderm_pseudobulk", .)  %>% gsub("X36weeks_Mesoderm_pseudobulk", "36weeks_Mesoderm_pseudobulk", .)  %>% gsub("X77weeks_Mesoderm_pseudobulk", "77weeks_Mesoderm_pseudobulk", .) %>% data.frame() %>% setnames("id_met") 

names <- names %>% merge(metadata, by="id_met", all.x=TRUE )
names$id<-make.unique(paste(names$stage,names$tissue,sep="_"))
colnames(correlation)<- names$id
rownames(correlation)<- names$id

#order_vec=c(names[names$stage_lineage=="Ectoderm",]$id, names[names$stage_lineage=="Endoderm",]$id,names[names$stage_lineage=="Mesoderm",]$id, names[names$stage_lineage=="ICM",]$id, names[names$stage_lineage=="Epiblast",]$id)
#germlayer_vec=c(names[names$stage_lineage=="Ectoderm",]$stage_lineage, names[names$stage_lineage=="Endoderm",]$stage_lineage,names[names$stage_lineage=="Mesoderm",]$stage_lineage, names[names$stage_lineage=="ICM",]$stage_lineage, names[names$stage_lineage=="Epiblast",]$stage_lineage)


#correlation<-correlation[match(order, row.names(correlation)),1:nrow(correlation),drop=FALSE]
#correlation<-correlation[1:nrow(correlation),match(order, colnames(correlation)),drop=FALSE]

pdf(paste0("/bi/group/reik/Jannat/nmt/analysis/methacc/meth/meth_comparison/dendrogram_heatmap_all_lineage_specific_enhancers.pdf"), width=19, height=10)

ha = HeatmapAnnotation(germ_layer = names$germ_layer, col = list(germ_layer = c("Endoderm" = "#2FBA37", "Ectoderm" = "#619CFF", "Mesoderm"="#F8766D")), show_legend = FALSE )
heatmap_plot= Heatmap(correlation,column_dend_height = unit(8, "cm"), show_heatmap_legend = FALSE,bottom_annotation=ha,
       colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),show_row_dend=FALSE,show_row_names = FALSE)


lgd = packLegend( Legend(title = "germ layer", at = c("Endoderm","Ectoderm","Mesoderm"), legend_gp = gpar(fill = c("#2FBA37","#619CFF","#F8766D"))),
  Legend(col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), title = "Pearson's correlation") )
draw(heatmap_plot, padding = unit(c(2, 2, 2, 20), "mm"))
draw(lgd, x = unit(0.95, "npc"), y = unit(0.8, "npc"))

dev.off()



