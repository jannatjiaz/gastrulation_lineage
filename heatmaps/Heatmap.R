# Heatmaps to show cooperativity of methylation in enhancers
#while read SAMPLE; do 
#    while read GERMLAYER; do 
#        ssub -o logs/${SAMPLE}_${GERMLAYER}.log -e logs/${SAMPLE}_${GERMLAYER}.err  -c 1 Rscript enhancer_plot.R ${GERMLAYER} ${SAMPLE}
#    done < germ_layers.txt
#done < samples.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(purrr)
library(ggplot2)
library(pheatmap)
library(viridis)


io <- list()
#metadata of the samples - added the germlayer level
#io$sample <- "GSE82725_Neural_tube_E15.5_cpg_meth"
#io$germ_layer <- "Ectoderm"
io$sample <- args[2]
io$germ_layer <- args[1]
if (io$germ_layer=="Ectoderm"){
    io$enhancer_positions<-"../../methylation_states_at_enhancers/E7.5_Ect_Enh_edit.filt.bed"
    io$data<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/ectoderm_enhancer/"
    io$files<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/ectoderm_enhancer/files.txt"
} else if (io$germ_layer=="Endoderm"){
    io$enhancer_positions<-"../../methylation_states_at_enhancers/E7.5_End_Enh_edit.filt.bed"
    io$data<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/endoderm_enhancer/"
    io$files<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/endoderm_enhancer/files.txt"
} else if (io$germ_layer=="Mesoderm"){
    io$enhancer_positions<-"../../methylation_states_at_enhancers/E7.5_Mes_Enh_edit.filt.bed"
    io$data<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/mesoderm_enhancers/"
    io$files<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/mesoderm_enhancers/files.txt"
}
io$outdir<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/enhancer_plots/"


files <- c(fread(io$files,header=F,stringsAsFactors=F)$V1)
enhancer_positions <- fread(io$enhancer_positions,header=F,stringsAsFactors=F) %>% setnames(c("chr","start","end","id","anno","call","strand"))
enhancer_positions$location_id <- NA
for (x in 1:nrow(enhancer_positions)){
    enhancer_positions$location_id[x]<-paste(strsplit(enhancer_positions$id[x],("_"))[[1]][7],".bed",sep="")
}

#just the methylated bases
num_sites=c()
for (enh in 1:length(files)){
    enhancer<-files[enh]
    meth_calls<-fread(paste(io$data,io$sample,"/Mes_Overlap_",io$sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F) %>% setnames(c("chr","start","end","meth"))
    num_sites<-c(num_sites,nrow(meth_calls))
}

columns<-paste("p",1:max(num_sites),sep="")
methylation <- data.table(expand.grid(files)) %>% setnames(c("id_met")) %>% .[,c(columns):=as.numeric(NA)]
for (enh in 1:length(files)){
    enhancer<-files[enh]
    meth_calls<-fread(paste(io$data,io$sample,"/Mes_Overlap_",io$sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F) %>% setnames(c("chr","start","end","meth"))
    meth_calls$matrix_pos<-0
    matrix_start_pos<-(max(num_sites)-nrow(meth_calls))/2
    for (x in 1:nrow(meth_calls)){
        methylation[enh,x+matrix_start_pos+1]=meth_calls$meth[x]
    }
}

methylation = methylation[order(rowSums(is.na(methylation))), ]

methylation_num = methylation[,2:ncol(methylation)]


pdf(paste(io$outdir,io$sample,"_",io$germ_layer,"_heatmap.pdf",sep=""), width=17, height=8)
pheatmap(methylation_num,cluster_rows=F,cluster_cols=F,annotation_names_row=F,annotation_names_col=F,show_rownames=F,show_colnames=F,na.rm = TRUE,na_col = "white",main=paste(io$sample,"_",io$germ_layer, sep=""))
dev.off()


# all the bases across the regions

#just the methylated bases
#size=c()
#for (enh in 1:nrow(enhancer_positions)){
#    size_enh = enhancer_positions$end[enh]-enhancer_positions$start[enh]
#    size<-c(size,size_enh)
#}

#columns<-paste("p",1:max(size)+1,sep="")
#methylation <- data.table(expand.grid(files)) %>% setnames(c("id_met")) %>% .[,c(columns):=as.numeric(NA)]

#for (enh in 1:length(files)){
#    enhancer<-files[enh]
#    meth_calls<-fread(paste(io$data,io$sample,"/Mes_Overlap_",io$sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F) %>% setnames(c("chr","start","end","meth"))
#    start<-enhancer_positions[enhancer_positions$location_id==enhancer]$start
#    end<-enhancer_positions[enhancer_positions$location_id==enhancer]$end    
#    meth_calls$matrix_pos<-0
#    if (end-start+1==1000){
#        meth_calls$matrix_pos<-meth_calls$end-start
#    } else {
#        size_adj<-(3000-(end-start+1))/2
#        meth_calls$matrix_pos<-meth_calls$end-start+size_adj
#    }
#    for (x in 1:nrow(meth_calls)){
#        methylation[enh,meth_calls$matrix_pos[x]+1]=meth_calls$meth[x]
#    }
#}

#methylation_num = methylation[,2:ncol(methylation)]
#pheatmap(methylation_num,cluster_rows=F,cluster_cols=F,annotation_names_row=F,annotation_names_col=F,show_rownames=F,show_colnames=F,na.rm = TRUE)
