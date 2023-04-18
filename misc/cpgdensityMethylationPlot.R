#"Gastrulation: general methylation statistics of cells"
# Adapted from Ricards scripts 
# Currently taking the entire enhancer peak

library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library("GenomicRanges")
library(gridExtra)
library(grid)
library(dplyr)

theme_pub <- function() {
  theme(
    axis.title.y = element_text(colour="black", size=rel(1.3)),
    axis.title.x = element_text(colour="black", size=rel(1.3)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.3)),
    axis.ticks = element_line(colour="black"),
  )
}

# Define settings 

io <- list()
#metadata of the samples - added the germlayer level
io$metadata <-  "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/sample_metadata_bulk_pseudobulk.txt"
#input the cpg methlation calls
io$data <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/"
#annotation directory whicha has the regions you would like to look at for
io$annos_dir  <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/refilter/"
#output directory
io$outdir <- "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/refilter/"
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

# Define genomic contexts (use NULL for no genomic context filtering)
opts$annos <- c( "E7.5_Ect_Enh_edit.refiltered2_withcounts")

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

pseudobulk_notE8.5 <-  tmp %>% .[stage=="E7.5" | stage=="E6.5" | stage=="E5.5" | stage=="E4.5" | stage=="E3.5",id_met]
pseudobulk_E8.5 <-  tmp %>% .[stage=="E8.5" | id_met=="Neural_Ectoderm_8weeks_pseudobulk" ,id_met]
cortes_morales <-  tmp %>% .[tissue =="Photoreceptors" | tissue =="Treg",id_met]

# Load bed coordinates of enhancer regions
if (!is.null(opts$annos)) {
  anno_dt <- lapply(opts$annos, function(anno) fread(sprintf("%s/%s.bed",io$annos_dir,anno), stringsAsFactors=T)) %>%  
    rbindlist %>% setnames(c("chr","start","end","id","anno","call","strand","pos","counts"))
}
anno_dt$stage = paste(as.character(lapply(strsplit(as.character(anno_dt$anno), split="_"), "[", c(3))), as.character(lapply(strsplit(as.character(anno_dt$anno), split="_"), "[", c(4))), sep="_")


# Load methylation data and filter genomic contexts and calculate stats

# general plots for all methylation per sample


# stage specific stats
stats <- data.table(expand.grid(opts$cells,c(names(opts$annos),"all"))) %>% setnames(c("id_met","anno")) %>% .[,as.vector(anno_dt[stage=="E7.5_Ect",id]):=as.numeric(NA)]
germ_layers = c("Ect")

#stats <- data.table(expand.grid(opts$cells,c(names(opts$annos),"all"))) %>% setnames(c("id_met","anno")) %>% .[,c("coverage_Ect","mean_Ect","base_met","enh_base_met"):=as.numeric(NA)]
#germ_layers = c("Ect")

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
    } else if (cell  %in%  cortes_morales ) {
      # Load sample methylation data
      data_tmp <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", header=FALSE, verbose=F, showProgress=F, stringsAsFactors=T) %>%
      setnames(c("chr","pos","pos2","rate","met_cpgs","nonmet_cpgs"))
      data <- data.table(chr=paste("chr",data_tmp$chr,sep=""),pos=data_tmp$pos,rate=data_tmp$rate)
    } else {
     # Load sample methylation data
      data <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", verbose=F, showProgress=F, stringsAsFactors=T) %>%
      setnames(c("chr","pos","pos2","rate","met_reads","nonmet_reads"))
    } 
    
    data_gr = GRanges(seqnames=data$chr, ranges = IRanges(start = data$pos,end = data$pos+1), rate=data$rate)

    # Compute genome-wide statistics
    for (layer in as.vector(anno_dt[stage=="E7.5_Ect",id])){
      germ_layer = makeGRangesFromDataFrame(anno_dt[id==paste(layer,sep=""),c("chr","start","end")], keep.extra.columns=T)
      data2 = as.data.table(subsetByOverlaps(data_gr,germ_layer))
      stats[id_met==cell & anno=="all", c(layer):=list((mean(data2$rate)))]
    }
  } 
  else {
    print(sprintf("Sample %s not found for methylation",cell))
  }
}
fwrite(stats, file=paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Ect_all_individual_enhancers.txt"), sep="\t", row.names=F, col.names=T)


# Load pre-computed statistics
stats <- fread(paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Ect_all_individual_enhancers.txt"),stringsAsFactors=FALSE) # %>% .[,mean:=mean*100]


for (sample in as.vector(stats$id)){
    # Get the methylation and counts for each enhancer in sample of interest
    methylation=c()
    counts=c()
    sample_id=sample
    for (enhancer in anno_dt$id){
    met<-stats[stats$id_met==sample_id]%>%pull(enhancer)
    methylation <- c(methylation,met)
    enhancer_size<- anno_dt[anno_dt$id==enhancer]$end-anno_dt[anno_dt$id==enhancer]$start
    counts<- c(counts,anno_dt[anno_dt$id==enhancer]$counts/enhancer_size)
    }

    pdf(paste0(io$outdir,"/CpG_density_plot/",sample,".pdf"), width=6, height=6)
    plot(counts,methylation,pch = 19,xlab = "CpG density", ylab = "Methylation (%)",frame = FALSE,ylim=c(0,100),main=sample_id)
    dev.off()

}
