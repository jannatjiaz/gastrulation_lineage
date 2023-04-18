
args <- commandArgs(trailingOnly = TRUE)

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



# Load bed coordinates of enhancer regions
if (!is.null(opts$annos)) {
  anno_dt <- lapply(opts$annos, function(anno) fread(sprintf("%s/%s.bed",io$annos_dir,anno), stringsAsFactors=T)) %>%  
    rbindlist %>% setnames(c("chr","start","end","id","anno","call","strand"))
}
# add a column for which embryonic stage this enhancer is specific to 
anno_dt$stage = paste(as.character(lapply(strsplit(as.character(anno_dt$anno), split="_"), "[", c(3))), as.character(lapply(strsplit(as.character(anno_dt$anno), split="_"), "[", c(4))), sep="_")
anno_dt$id<-paste(anno_dt$chr,anno_dt$start,anno_dt$end,sep="_")


cells_arg<-opts$cells[as.numeric(args[1])]
print(cells_arg)

merged_meth <- data.table(expand.grid(anno_dt[anno_dt$stage==paste("E7.5_",args[2],sep=""),id])) %>% setnames(c("id")) %>% .[,c(cells_arg):=as.numeric(NA)]

for (cell in cells_arg) {
  if (file.exists(sprintf("%s/%s.cov.gz",io$data,cell))) {
    print(cell)
    # print(sprintf("Loading %s methylation...",cell))
     if (cell  %in%  cortes_morales_lister ) {
    # Load sample methylation data
    data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=FALSE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
    setnames(c("chr","pos","pos","rate","met_cpgs","nonmet_cpgs"))
    data <- data.table(chr=paste("chr",data_tmp$chr,sep=""),pos=data_tmp$pos,rate=data_tmp$rate)
    } else if (cell  %in%  pseudobulk_E8.5) {
    # Load sample methylation data
    data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=TRUE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
    setnames(c("chr","pos","met_cpgs","nonmet_cpgs","rate"))
    data <- data.table(chr=data_tmp$chr,pos=data_tmp$pos,rate=data_tmp$rate)
    } else if (cell  %in%  pseudobulk_notE8.5 ) {
    # Load sample methylation data
    data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=TRUE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
    setnames(c("chr","pos","met_cpgs","nonmet_cpgs","rate"))
    data <- data.table(chr=paste("chr",data_tmp$chr,sep=""),pos=data_tmp$pos,rate=data_tmp$rate)
    } else {
    # Load sample methylation data
    data <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", verbose=F, showProgress=F, stringsAsFactors=T) %>%
    setnames(c("chr","pos","pos2","rate","met_reads","nonmet_reads"))
    } 
    data_gr = GRanges(seqnames=data$chr, ranges = IRanges(start = data$pos,end = data$pos+1), rate=data$rate)
    for (layer in merged_meth$id){
        germ_layer = makeGRangesFromDataFrame( anno_dt[anno_dt$id==layer,c("chr","start","end")] )
        data2 = as.data.table(subsetByOverlaps(data_gr,germ_layer))
        merged_meth[id==layer, paste(cell):=mean(data2$rate)]
    }
  } else {
      print(sprintf("Sample %s not found for methylation",cell))
  }
}

fwrite(merged_meth, file=paste0(io$outdir,"/E7.5_",args[2],"_",args[1],"_enhancer_meth_for_clustering.txt"), sep="\t", row.names=F, col.names=T)
