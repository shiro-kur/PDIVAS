import argparse
import os
import sys
from cyvcf2 import VCF, Writer
import pickle
from statistics import mean
import re
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np

def scoring_to_vcf(var,vep_info_index,spliceai_info_index,clf,features,filt):
	#Filt1: Return annotaitons without annotations
	if (var.INFO.get("CSQ") is None
		or var.INFO.get("SpliceAI") is None
	   ) :
		if filt == "off" :
			var.INFO["PDIVAS"] = "None|wo_annots"
			return var
		else :
			return None

   
	# arrange VEP annotations and filtering out annotations outside of Deep intron
	vep_annot = {dict(zip(vep_info_index, x.strip().split("|")))["Gene"]:
				dict(zip(vep_info_index, x.strip().split("|")))
				for x in var.INFO["CSQ"].strip().split(",")
				}
	
	vep_annot = ({k: v for k, v in vep_annot.items()
				  if ('intron_variant' in v["Consequence"]) and 
				  abs(int("".join(re.findall('.+c\.-?\*?[0-9].*([+,-][0-9]+)', v["HGVSc"])) or "0"))>=50})

	# arrange SpliceAI annotations and filtering out annotatins with "." predictions
	spliceai_annot = {dict(zip(spliceai_info_index, x.strip().split("|")))["GENE_ID"]:
					dict(zip(spliceai_info_index, x.strip().split("|")))
					for x in var.INFO["SpliceAI"].strip().split(",")
					}

	spliceai_annot = ({k: v for k, v in spliceai_annot.items() 
						if v["DS_AG"]!="."}) 

	#Filt2: Return annotaitons without candidates
	if (var.CHROM=="Y") | (var.CHROM=="chrY") | (vep_annot == {}) | (spliceai_annot == {}) :
		if filt == "off" :
			var.INFO["PDIVAS"] = "None|out_of_scope"
			return var 
		else :
			return None

	#Get ConSplice and MES_max
	for gene, annot in vep_annot.items():
		#Prepare ConSplice
		vep_annot[gene]["ConSplice"] = (max(annot["ConSplice"].split("&")) if "&" in annot["ConSplice"] 
								  else -1 if annot["ConSplice"]==""
								  else annot["ConSplice"] )
		#Prepare MES
		annot["MES-SWA_acceptor_alt"] = float(annot["MES-SWA_acceptor_alt"] or "0")
		annot["MES-SWA_donor_alt"] = float(annot["MES-SWA_donor_alt"] or "0")
		annot["MES-SWA_acceptor_diff"] = float(annot["MES-SWA_acceptor_diff"] or "0")
		annot["MES-SWA_donor_diff"] = float(annot["MES-SWA_donor_diff"] or "0")
		vep_annot[gene]["MES_SA"] = annot["MES-SWA_acceptor_alt"] if (annot["MES-SWA_acceptor_alt"] > 0 and annot["MES-SWA_acceptor_diff"]) < 0 else 0
		vep_annot[gene]["MES_SD"] = annot["MES-SWA_donor_alt"] if (annot["MES-SWA_donor_alt"] > 0 and annot["MES-SWA_donor_diff"]) < 0 else 0
		vep_annot[gene]["MES_max"] = max([annot["MES_SA"],annot["MES_SD"]])

		del (vep_annot[gene]["MES-SWA_acceptor_diff"], vep_annot[gene]["MES-SWA_donor_diff"], 
			 vep_annot[gene]["MES-SWA_acceptor_alt"], vep_annot[gene]["MES-SWA_donor_alt"])

	#Get SpliceAI_del_gain_max,mean,raw_gain_mean
	for gene, annot in spliceai_annot.items():
		spliceai_annot[gene]["DS_G_MAX"] = max(float(annot["DS_AG"]),float(annot["DS_DG"])) if annot["DS_AG"] != "." else -1
		spliceai_annot[gene]["DS_G_MEAN"] = mean([float(annot["DS_AG"]),float(annot["DS_DG"])]) if annot["DS_AG"] != "." else -1
		spliceai_annot[gene]["RS_G_MEAN"] = mean([float(annot["RS_AG"]),float(annot["RS_DG"])]) if annot["DS_AG"] != "." else -1

	#Gene match between SpliceAI and ConSplice
	combined_score_list = []
	for gene_id in spliceai_annot :
		
		if (gene_id in vep_annot) & (spliceai_annot[gene_id]["DS_G_MAX"]!=-1) :
			combined_score_list.append(
				{"GENE_ID" : gene_id,
				 "AI_del_gain_max" : spliceai_annot[gene_id]["DS_G_MAX"],
				 "AI_del_gain_mean" : spliceai_annot[gene_id]["DS_G_MEAN"],
				 "AI_raw_gain_mean" : spliceai_annot[gene_id]["RS_G_MEAN"],
				 "ConSplice" : vep_annot[gene_id]["ConSplice"],
				 "MES_max" : vep_annot[gene_id]["MES_max"]	
				})

	#Return annotations with matched genes
	if len(combined_score_list)>0 :
		new_anno = []
		var_df = pd.DataFrame(combined_score_list)
		var_df["PDIVAS"] = clf.predict_proba(var_df[features])[:,1]
		new_anno = ["{}|{}".format(d["GENE_ID"], round(d["PDIVAS"],3) \
					if (var.CHROM != "Y") and (var.CHROM != "Mt") and (d["ConSplice"]!=-1) else "None") \
					for d in var_df.to_dict("records")]
		var.INFO["PDIVAS"] = ",".join(new_anno) 
		return var

	#Filt3: Return annotations without matched genes
	else :
		if filt == "off" :
			var.INFO["PDIVAS"] = "None|no_gene_match"
			return var 
		else :
			return None

