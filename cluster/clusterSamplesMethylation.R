
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
io$data_dir<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/"
#metadata of the samples - added the germlayer level
io$metadata <-  "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/sample_metadata_bulk_pseudobulk.txt"
#input the cpg methlation calls
io$data <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/"
#annotation directory whicha has the regions you would like to look at for
io$annos_dir  <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/refilter/differential/Ectoderm_refilter/"
#output directory
io$outdir <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/refilter/differential/Ectoderm_refilter/"
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



# Define which cells to use ie the plate and position of the cells to use 
tmp <- fread(io$metadata) %>% .[!is.na(id_met)] 
if (is.null(opts$stage_lineage)) {
  opts$cells <- tmp[,id_met]
} else {
  opts$cells <- tmp %>%
    .[,stage_lineage:=paste(germ_layer)] %>%
    .[pass_metQC=="TRUE" & stage_lineage%in%opts$stage_lineage,id_met]
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
cortes_morales_lister<-c(cortes_morales, lister)


sample1="GSE82692_Forebrain_E16.5_cpg_meth"

data1=fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data_dir,sample1), sep="\t", header=TRUE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
      setnames(c("chr","pos1","pos2","rate","met_cpgs","nonmet_cpgs"))
data1$id=paste(data1$chr,data1$pos2,sep="_")

merged_meth<-data.frame(data1$id) %>%  setnames("id")

for (cell in opts$cells) {
    if (file.exists(sprintf("%s/%s.cov.gz",io$data,cell))) {
        # print(sprintf("Loading %s methylation...",cell))
        if (cell  %in%  pseudobulk_notE8.5 ) {
        # Load sample methylation data
        data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=TRUE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
        setnames(c("chr","pos","met_cpgs","nonmet_cpgs","rate"))
        data <- data.table(chr=paste("chr",data_tmp$chr,sep=""),pos=data_tmp$pos,rate=data_tmp$rate)
        } else if (cell  %in%  pseudobulk_E8.5) {
        # Load sample methylation data
        data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=TRUE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
        setnames(c("chr","pos","met_cpgs","nonmet_cpgs","rate"))
        data <- data.table(chr=data_tmp$chr,pos=data_tmp$pos,rate=data_tmp$rate)
        } else if (cell  %in%  cortes_morales_lister ) {
        # Load sample methylation data
        data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=FALSE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
        setnames(c("chr","pos","pos","rate","met_cpgs","nonmet_cpgs"))
        data <- data.table(chr=paste("chr",data_tmp$chr,sep=""),pos=data_tmp$pos,rate=data_tmp$rate)
        } else {
        # Load sample methylation data
        data <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", verbose=F, showProgress=F, stringsAsFactors=T) %>%
        setnames(c("chr","pos","pos2","rate","met_reads","nonmet_reads"))
        } 
        id<-paste(data$chr,data$pos2,sep="_")
        data2_input_to_matrix= data.frame(id, data$rate) %>% setnames(c("id",cell))
        merged_meth<- merged_meth %>%  merge(data2_input_to_matrix, by="id",all.x=TRUE)
        print(dim(merged_meth))
    } else {
        print(sprintf("Sample %s not found for methylation",cell))
    }
}
fwrite(merged_meth, file=paste0(io$outdir,"/all_CpG_meth_for_clustering.txt"), sep="\t", row.names=F, col.names=T)


correlation<-cor(merged_meth[2:ncol(merged_meth)])
correlation<-as.matrix(correlation, labels=TRUE)

names<-data.frame(colnames(merged_meth)[2:length(colnames(merged_meth))]) %>% setnames("id_met")
names <- names %>% merge(metadata, by="id_met" )
names$id<-make.unique(paste(names$stage,names$tissue,sep="_"))
colnames(correlation)<- names$id
rownames(correlation)<- names$id


heatmap_plot= Heatmap(correlation,cluster_rows=FALSE,cluster_columns=FALSE,column_dend_height = unit(8, "cm"), show_heatmap_legend = FALSE,
       colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
heatmap_plot


ha = HeatmapAnnotation(germ_layer = stats$seperate, col = list(germ_layer = c("Blood" = "Red", "Endoderm" = "#2FBA37", "Brain" = "#619CFF", "Skin"="#5373B9", "Mesoderm"="#F8766D")), show_legend = FALSE )
heatmap_plot= Heatmap(matrix_to_plot, cluster_columns = hc,cluster_rows=FALSE,column_dend_height = unit(8, "cm"),bottom_annotation=ha, show_heatmap_legend = FALSE,
       colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

lgd = packLegend( Legend(title = "germ layer", at = c("Blood","Endoderm","Brain","Skin","Mesoderm"), legend_gp = gpar(fill = c("Red","#2FBA37","#619CFF","#5373B9","#F8766D"))),
  Legend(col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), title = "Pearson's correlation") )
draw(heatmap_plot, padding = unit(c(2, 2, 2, 20), "mm"))
draw(lgd, x = unit(0.95, "npc"), y = unit(0.8, "npc"))



