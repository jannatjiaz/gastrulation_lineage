# Heatmaps to show cooperativity of methylation in enhancers
#while read STAGE; do 
#    ssub -o logs/${STAGE}.log -e logs/${STAGE}.err  -c 1 Rscript enhancer_plot_multiple_plots2.R ${STAGE}
#done < stages.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(purrr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(grid)
library(gridExtra)
library(stringr)

io <- list()
io$germ_layer <- c("Ectoderm","Endoderm","Mesoderm")
io$outdir<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/enhancer_plots/"
io$metadata <-  "/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/sample_metadata_all_data_postqc.txt"
io$stage <- args[1]


opts<-list()
# Define which cells to use ie the plate and position of the cells to use 
tmp <- fread(io$metadata) %>% .[!is.na(id_met)] 
opts$cells_tmp <- tmp %>% .[pass_metQC=="TRUE" & stage%in%io$stage & id_met!="GSE82481_Hindbrain_P0_cpg_meth",id_met]


# Load sample metadata
metadata <- fread(io$metadata,stringsAsFactors=F) %>% 
  .[,c("id_met","stage","germ_layer","tissue")] %>%
  .[id_met%in%opts$cells_tmp] %>%
  .[,stage_lineage:=as.factor(paste(germ_layer,sep="_"))]

opts$cells<-str_sub(sort(paste(metadata$stage_lineage,metadata$id_met,sep="_")),10,-1)

plot_list=list()
n=1
for (samp in 1:length(opts$cells)){
    sample = opts$cells[samp]
    for (germ_lay in 1:length(io$germ_layer)){
        germ=io$germ_layer[germ_lay]
        if (germ=="Ectoderm"){
            io$enhancer_positions<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/methylation_states_at_enhancers/refilter/E7.5_Ect_Enh_edit.refiltered.bed"
            io$data<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/ectoderm_enhancer/"
            io$files<-"ect_files_2.tmp"
        } else if (germ=="Endoderm"){
            io$enhancer_positions<-"../../methylation_states_at_enhancers/E7.5_End_Enh_edit.filt.bed"
            io$data<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/endoderm_enhancer/"
            io$files<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/endoderm_enhancer/files.txt"
        } else if (germ=="Mesoderm"){
            io$enhancer_positions<-"../../methylation_states_at_enhancers/E7.5_Mes_Enh_edit.filt.bed"
            io$data<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/mesoderm_enhancers/"
            io$files<-"/bi/group/reik/Jannat/nmt/analysis/methacc/meth/all_data/mesoderm_enhancers/files.txt"
        }
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
            if (file.size(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""))!=0 && nrow(fread(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F))>10 ){            
                nrows = nrow(fread(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F))
                new_start_pos <- nrows-round(nrows*0.7)
                new_end_pos <- round(nrows*0.7)
                meth_calls<-fread(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F)[new_start_pos:new_end_pos,] %>% setnames(c("chr","start","end","meth"))
                num_sites<-c(num_sites,nrow(meth_calls))
            }
        }
        columns<-paste("p",1:max(num_sites),sep="")
        methylation <- data.table(expand.grid(files)) %>% setnames(c("id_met")) %>% .[,c(columns):=as.numeric(NA)]
        for (enh in 1:length(files)){
            enhancer<-files[enh]
            if (file.size(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""))!=0 && nrow(fread(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F))>10){
                nrows = nrow(fread(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F))
                new_start_pos <- nrows-round(nrows*0.7)
                new_end_pos <- round(nrows*0.7)                
                meth_calls<-fread(paste(io$data,sample,"/Mes_Overlap_",sample,"_",enhancer,sep=""),header=F,stringsAsFactors=F)[new_start_pos:new_end_pos,] %>% setnames(c("chr","start","end","meth"))
                meth_calls$matrix_pos<-0
                matrix_start_pos<-(max(num_sites)-nrow(meth_calls))/2
                for (x in 1:nrow(meth_calls)){
                    methylation[enh,x+matrix_start_pos+1]=meth_calls$meth[x]
                }
            }
        }
        methylation = methylation[order(rowSums(is.na(methylation))), ]
        methylation_num = methylation[,2:ncol(methylation)]
        plot_list[[n]] = pheatmap(as.matrix(methylation_num),cluster_rows=F,cluster_cols=F,annotation_names_row=F,annotation_names_col=F,show_rownames=F,show_colnames=F,na_col = "white",main=paste(sample,"_",germ, sep=""))[[4]] 
        n=n+1
    }
}

#pheatmap(methylation_num,cluster_rows=F,cluster_cols=F,annotation_names_row=F,annotation_names_col=F,show_rownames=F,show_colnames=F,na.rm = TRUE,na_col = "white",main=paste(sample,"_",germ, sep=""))
#g<-do.call(grid.arrange,c(plot_list, ncol=3))

#pdf(paste(io$outdir,"/enhancer_heatmaps_",io$stage,".pdf",sep=""), width=17, height=(10*length(io$germ_layer)))
pdf(paste(io$outdir,"/enhancer_heatmaps_",io$stage,".pdf",sep=""), width=17, height=(3*length(opts$cells)))
g<-do.call(grid.arrange,c(plot_list, ncol=3))
dev.off()
