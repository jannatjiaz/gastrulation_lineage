#"Gastrulation: general methylation statistics of cells"
# Cell types have been annotated by ricard and stevens transcription PCA method 
# Currently taking the entire enhancer peak could be more useful to take the center of the peak and 500bp around it.

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

# Define genomic contexts (use NULL for no genomic context filtering)
opts$annos <- c(
 "E7.5_Ect_Enh_edit.refiltered_E11_E16_P0",
 "E7.5_End_Enh_edit.refiltered_E11_E16_P0",
 "E7.5_Mes_Enh_edit.refiltered_E11_E16_P0",
 "E7.5_Mesend_Enh_edit.refiltered_E11_E16_P0"
)

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
luo<-grep("Luo",opts$cells, value=T)
Charlton<-grep("Charlton",opts$cells, value=T)
knaupp<-grep("Knaupp",opts$cells, value=T)
zhang<-grep("Zhang",opts$cells, value=T)
cortes_morales_lister<-c(cortes_morales, lister,luo,Charlton,knaupp,zhang)

# Load bed coordinates of enhancer regions
if (!is.null(opts$annos[1:4])) {
  anno_dt <- lapply(opts$annos[1:4], function(anno) fread(sprintf("%s/%s.bed",io$annos_dir,anno), stringsAsFactors=T)) %>%  
    rbindlist %>% setnames(c("chr","start","end","id","anno","call","strand"))
}
# add a column for which embryonic stage this enhancer is specific to 
anno_dt$stage = paste(as.character(lapply(strsplit(as.character(anno_dt$anno), split="_"), "[", c(3))), as.character(lapply(strsplit(as.character(anno_dt$anno), split="_"), "[", c(4))), sep="_")



# Load bed coordinates of enhancer regions for all  enhancers
if (!is.null(io$all_enhancers)) {
  all_enhancers  <- fread(io$all_enhancers, stringsAsFactors=T)  %>% setnames(c("chr","start","end","name","score","strand","thichStart","thickEnd","reserved","ccre","encodeLabel","zscore","ucscLabel","accessionLabel","description")) %>% .[,c("chr","start","end","ucscLabel")]
}
all_enhancers <- all_enhancers[all_enhancers$ucscLabel=="enhD" | all_enhancers$ucscLabel=="enhP"]


# Load methylation data and filter genomic contexts and calculate stats

# general plots for all methylation per sample


# stage specific stats
stats <- data.table(expand.grid(opts$cells,c(names(opts$annos),"all"))) %>% setnames(c("id_met","anno")) %>% .[,as.vector(anno_dt[stage=="E7.5_Mesend",id]):=as.numeric(NA)]


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
      setnames(c("chr","pos","pos2","rate","met_cpgs","nonmet_cpgs"))
      data <- data.table(chr=paste("chr",data_tmp$chr,sep=""),pos=data_tmp$pos,rate=data_tmp$rate)
    } else {
     # Load sample methylation data
      data <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",io$data,cell), sep="\t", verbose=F, showProgress=F, stringsAsFactors=T) %>%
      setnames(c("chr","pos","pos2","rate","met_reads","nonmet_reads"))
    } 
    
    data_gr = GRanges(seqnames=data$chr, ranges = IRanges(start = data$pos,end = data$pos+1), rate=data$rate)

    # Compute genome-wide statistics
    for (layer in as.vector(anno_dt[stage=="E7.5_Mesend",id])){
      germ_layer = makeGRangesFromDataFrame(anno_dt[id==paste(layer,sep=""),c("chr","start","end")], keep.extra.columns=T)
      data2 = as.data.table(subsetByOverlaps(data_gr,germ_layer))
      stats[id_met==cell & anno=="all", c(layer):=list((mean(data2$rate)))]
    }
  } 
  else {
    print(sprintf("Sample %s not found for methylation",cell))
  }
}
fwrite(stats, file=paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Mesend_E11_E16_P0_individual_enhancers.txt"), sep="\t", row.names=F, col.names=T)

stats <- fread(paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Mesend_E11_E16_P0_individual_enhancers.txt"),stringsAsFactors=FALSE) # %>% .[,mean:=mean*100]
#opts$cells <- opts$cells[186:195]
#stats2 <- data.table(expand.grid(opts$cells,c(names(opts$annos),"all"))) %>% setnames(c("id_met","anno")) %>% .[,as.vector(anno_dt[stage=="E7.5_Mesend",id]):=as.numeric(NA)]
#stats <-  rbind(stats, stats2)


# Load pre-computed statistics
stats <- fread(paste0(io$outdir,"/sample_stats_bygermlayer_pseudobulk_Mesend_E11_E16_P0_individual_enhancers.txt"),stringsAsFactors=FALSE) # %>% .[,mean:=mean*100]
stats <- stats[complete.cases(stats)]
stats <- stats %>% merge(metadata, by="id_met")

blood <- c("Blood_progenitors_1","Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3","Haematoendothelial_progenitors","Macrophage","macrophage","Treg")
brain <- c("Rostral_neurectoderm","Forebrain_Midbrain_Hindbrain","Neural_crest","Spinal_cord","Hindbrain","Midbrain","Forebrain","Neural_tube","cerebellum","cortex","olfactory_bulb","Neural","glia","Photoreceptors")
skin <- c("Surface_ectoderm","skin" )
endoderm <- unique(stats[stats$stage_lineage=="Endoderm"]$tissue)
mesoderm <- setdiff(unique(stats[stats$stage_lineage=="Mesoderm"]$tissue),blood)
stats$seperate <- stats$tissue
stats$seperate[stats$tissue %in%blood] <- "Blood"
stats$seperate[stats$tissue %in% skin] <- "Skin"
stats$seperate[stats$tissue %in% brain] <- "Brain"
stats$seperate[stats$tissue %in% endoderm] <- "Endoderm"
stats$seperate[stats$tissue %in% mesoderm] <- "Mesoderm"
stats$seperate[stats$tissue %in% c("ICM")] <- "ICM"
stats$seperate[stats$tissue %in% c("Epiblast","Primitive_endoderm","Anterior_Primitive_Streak","Primitive_Streak")] <- "Epiblast"


# Boxplots with mean methylation rate per split by sublieage
#order.lineages <- c("Ectoderm","Mesoderm","Endoderm")
order.stage <- c("E3.5","E4.5","E5.5","E6.5","E7.5","E8.5","E10.5","E11.5","E12.5","E13.5", "E14.5","E15.5","E16.5","P0","4weeks","8weeks","12weeks","48weeks","72weeks","88weeks","96weeks")
stats_filt <- stats %>% copy %>%
  .[, tissue:=stringr::str_replace_all( tissue,"_"," ")] %>%
  .[,tissue:=factor(tissue)]

stats_filt=stats_filt[stats_filt[,tissue!="Craniofacial"]]
stats_filt=stats_filt[stats_filt[,tissue!="thymus"]]
stats_filt=stats_filt[stats_filt[,stage!="E3.5"]]
stats_filt=stats_filt[stats_filt[,stage!="E4.5"]]
stats_filt=stats_filt[stats_filt[,stage!="E5.5"]]
stats_filt=stats_filt[stats_filt[,stage!="E6.5"]]
stats_filt=stats_filt[stats_filt[,tissue!="Erythroid"]]


stats_filt$stage <- factor(stats_filt$stage, levels=c(order.stage))
#aggregate(mean_Ect~stage+germ_layer, stats_filt, median)

stats_filt$germ_layer <- paste(stats_filt$germ_layer,"derived",sep=" ")

stats_filt$stage2 <- 0
for (row in 1:nrow(stats_filt)){
  age <- as.character(stats_filt$stage[row])
  if (grepl("weeks", age, fixed = TRUE)){
    stats_filt$stage2[row] = as.numeric(gsub('weeks', '', age))*7
  } else if (grepl("E", age, fixed = TRUE)){
    stats_filt$stage2[row] = as.numeric(gsub('E', '', age))-21
  }
}

#ect_average <- stats_filt %>% group_by(germ_layer,stage2)%>%summarise(mean_val=mean(mean_Ect))
enhancers<- as.vector(anno_dt[stage=="E7.5_Mesend",id])
order.ticks <- c("88weeks","8weeks","12weeks","48weeks","72weeks","96weeks","E7.5","E8.5","E10.5","E12.5","E15.5","P0","E11.5","E16.5","E14.5","E13.5","4weeks")

breaks <-c(0.0,28.0,56.0,84.0,336.0,504.0,616.0,672.0)
order.ticks <- c("P0","4weeks","8weeks","12weeks","48weeks","72weeks","88weeks","96weeks")

for (enh in 1:length(enhancers)){
  data<-stats_filt
  names(data)[names(data) == enhancers[enh]] <- "to.plot"

  ect <- ggplot( data, aes(x=stage2, y=to.plot,color=seperate,fill=seperate)) + 
    stat_summary(geom = "line", fun = mean,size=1) +
    geom_point(size=0.5) +
    #stat_smooth(method=loess,se=FALSE) + 
    #geom_errorbar(data=data,aes(ymin=mean_Ect+enh_base_met-sd(mean_Ect+enh_base_met), ymax=mean_Ect+enh_base_met+sd), width=.2) + 
    stat_summary(fun.data=mean_sdl, geom="errorbar", fun.args=list(mult=1)) +
    #geom_ribbon(data=data, aes(y = mean_Ect+enh_base_met, ymin = mean_Ect+enh_base_met-sd(mean_Ect+enh_base_met), ymax = mean_Ect+enh_base_met+sd(mean_Ect+enh_base_met), fill = seperate), alpha = .2) +
    coord_cartesian(ylim=c(0,100),xlim=c(-14,700)) +
    ylab("Methylation (%)") +
    xlab("Days post fertilization") +
    scale_x_continuous(breaks = breaks, labels=order.ticks) +
    ggtitle("E7.5 Mesendoderm Enhancers") + 
    scale_color_manual(values = c("Skin"="#5373B9","Brain"="#619CFF","Endoderm"="#2FBA37","Blood"="red","Mesoderm"="#F8766D")) +
    scale_fill_manual( values = c("Skin"="#5373B9","Brain"="#619CFF","Endoderm"="#2FBA37","Blood"="red","Mesoderm"="#F8766D")) +
    theme_bw() +
    theme_pub() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(colour="black",size=rel(1.2), angle=45, hjust=1, vjust=1),
      strip.background = element_blank(),
      strip.text = element_text(color="white", size=rel(0.1))
    )

  pdf(paste0(io$outdir,"/enhancer_qc_plot_E11_E16_P0/mesendoderm_trajectory_pseudobulk_aged_enhancer_ect_",enhancers[enh],".pdf"), width=20, height=6)
  grid.arrange(ect, vp=viewport(width=0.9, height=0.9))
  dev.off()
}

