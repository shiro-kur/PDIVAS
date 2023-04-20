from cyvcf2 import VCF, Writer
import csv
import sys

def vep_editor(var,vep_info_index):
    vep_annot = {dict(zip(vep_info_index, x.strip().split("|")))["Gene"]:
                        x.strip().split("|")
                        for x in var.INFO["CSQ"].strip().split(",")}
    
    base_anno = [var.CHROM,var.start+1,var.ID,var.REF,var.ALT[0],var.QUAL,(var.FILTER or "PASS")]
    oup_list = []
    
    for gene_id in vep_annot:    
        #output formatting
        oup_list += [base_anno+vep_annot[gene_id]+\
                     ["." for i in range(len(spliceai_info_index))]+["." for i in range(len(PD_info_index))]]
    
    return oup_list

def vep_ai_matcher(var,vep_info_index,spliceai_info_index,PD_info_index):
    vep_annot = {dict(zip(vep_info_index, x.strip().split("|")))["Gene"]:
                        x.strip().split("|")
                        for x in var.INFO["CSQ"].strip().split(",")}

    spliceai_annot = {dict(zip(spliceai_info_index, x.strip().split("|")))["GENE_ID"]:
                        x.strip().split("|")
                        for x in var.INFO["SpliceAI"].strip().split(",")}

    PD_annot = {dict(zip(PD_info_index, x.strip().split("|")))["PD_ID"]:
                        x.strip().split("|")
                        for x in var.INFO["PDIVAS"].strip().split(",")}

    #Gene match & output formatting
    base_anno = [var.CHROM,var.start+1,var.ID,var.REF,var.ALT[0],var.QUAL,(var.FILTER or "PASS")]
    oup_list = []

    for gene_id in vep_annot:
        if (gene_id in spliceai_annot) :
            #output formatting
            if (gene_id in PD_annot) :
                oup_list += [base_anno+vep_annot[gene_id]+\
                             spliceai_annot[gene_id]+PD_annot[gene_id]]
            else :
                oup_list += [base_anno+vep_annot[gene_id]+\
                             spliceai_annot[gene_id]+list(PD_annot.values())[0]]

        else :
            oup_list += [base_anno+vep_annot[gene_id]+\
                         ["." for i in range(len(spliceai_info_index))]+["." for i in range(len(PD_info_index))]]     

    return oup_list