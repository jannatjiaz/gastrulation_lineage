## Add germ layer info to metadata

import pandas as pd
import numpy as np

#load sv calls - since this is huge keep it as a pandas df
metadata = pd.read_csv("/bi/group/reik/Jannat/nmt/analysis/methacc/meth/sample_metadata_postqc.txt",sep="\t",skiprows=1,header=None)
metadata.columns = ["cell","sample","id_met","id_acc","plate","id_rna","method","embryo","stage","nCount_RNA","nFeature_RNA","mit_percent_RNA","rib_percent_RNA","celltype","celltype2","celltype3","celltype.score","closest.cell","nCG","met_rate","nGC","acc_rate","pass_rnaQC","pass_metQC","pass_accQC","germ_layer"]
metadata['germ_layer'] = "NA"


celltype_dict = {"Allantois": "Extra-embryonic", "Blood_progenitors_1": "Mesoderm", "Blood_progenitors_2": "Mesoderm", "Cardiomyocytes": "Mesoderm",
"Endothelium": "Mesoderm", "Erythroid1": "Mesoderm", "Erythroid2": "Mesoderm", "Erythroid3": "Mesoderm", "ExE_ectoderm": "Extra-embryonic",
"ExE_endoderm": "Extra-embryonic", "ExE_mesoderm": "Extra-embryonic", "Forebrain_Midbrain_Hindbrain": "Ectoderm", "Gut": "Endoderm", 
"Haematoendothelial_progenitors": "Mesoderm", "Intermediate_mesoderm": "Mesoderm", "Mesenchyme": "Mesoderm", "Nascent_mesoderm": "Mesoderm",
"Neural_crest": "Ectoderm", "NMP": "NA", "Paraxial_mesoderm": "Mesoderm", "Parietal_endoderm": "Extra-embryonic", "Pharyngeal_mesoderm": "Mesoderm",
"Rostral_neurectoderm": "Ectoderm", "Somitic_mesoderm": "Mesoderm", "Spinal_cord": "Ectoderm", "Surface_ectoderm": "Ectoderm", 
"Visceral_endoderm": "Extra-embryonic", "NA":"NA","Anterior_Primitive_Streak":"Epiblast","Caudal_epiblast":"Epiblast","Caudal_neurectoderm": "Epiblast",
"Def._endoderm": "Endoderm","Epiblast":"Epiblast","ICM":"Epiblast","Mixed_mesoderm":"Mesoderm","Notochord":"Mesoderm","PGC":"Yolk_Sac",
"Primitive_endoderm":"Epiblast","Primitive_Streak":"Epiblast","none":"NA"}

 "Blood_progenitors_1","Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3","Haematoendothelial_progenitors","Heart","Macrophage"

for x in range(len(metadata)):
    #if metadata["stage"][x]=="E8.5":
    if pd.isna(metadata["celltype"][x])==False:
        germ_layer = celltype_dict[metadata["celltype"][x]]
        metadata["germ_layer"][x]=germ_layer


order=["cell","sample","id_met","id_acc","plate","id_rna","method","embryo","stage","nCount_RNA","nFeature_RNA","mit_percent_RNA","rib_percent_RNA","celltype","celltype2","celltype3","celltype.score","closest.cell","nCG","met_rate","nGC","acc_rate","pass_rnaQC","pass_metQC","pass_accQC","germ_layer"]
metadata.to_csv("{}".format('/bi/group/reik/Jannat/nmt/analysis/methacc/meth/sample_metadata_postqc_withgermlayers_all.txt'), index=False, sep='\t', na_rep="none", header=order)
