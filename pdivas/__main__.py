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
import csv

from .__init__ import __cur_path__, __version__
from pdivas.scoring_to_vcf import scoring_to_vcf
from pdivas.vcf2tsv import vep_editor
from pdivas.vcf2tsv import vep_ai_matcher

def get_options():
    parser = argparse.ArgumentParser(
        prog="pdivas",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
                " ========================================================================================\n"
                "|| PDIVAS: Pathogenicity predictor for Deep-Intronic Variants causing Aberrant Splicing ||\n"
                " ========================================================================================\n"),
                        )
    
    # Calculator of PDIVAS scores
    subparsers = parser.add_subparsers(dest="subcommand")
    
    pred_parser = subparsers.add_parser("predict")
    pred_parser.add_argument(
        "-I",
        required=True,
        metavar="input.vcf/vcf.gz",
        type=str,
        help="The path to the vcf(.gz) file to add PDIVAS annotation")
        
        
    pred_parser.add_argument(
        "-O",
        required=True,
        metavar="output.vcf/vcf.gz",
        type=str,
        help="The path to output vcf(.gz) file name and pass")
        
        
    pred_parser.add_argument(
        "-F",
        metavar="filtering:off/on",
        required=False,
        default="off",
        type=str,
        choices=["off","on"],
        help="Output all variants (-F off; default) or only deep-intronic variants with PDIVAS scores (-F on)")
        
    pred_parser.set_defaults(handler=predict)
    
    
    
    tsv_parser = subparsers.add_parser("vcf2tsv")
    tsv_parser.add_argument(
        "-I",
        required=True,
        metavar="input.vcf/vcf.gz",
        type=str,
        help="The path to the vcf(.gz) file with PDIVAS annotation")
        
        
    tsv_parser.add_argument(
        "-O",
        required=True,
        metavar="output.tsv",
        type=str,
        help="The path to output tsv file name and pass")
        
    tsv_parser.set_defaults(handler=vcf2tsv)
    args = parser.parse_args()
    
    return args,parser
   

def predict(args):
        
    try:
        vcf = VCF(args.I)
    
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
        
    vep_info_index = vcf["CSQ"]["Description"].strip().split("Format:")[1].strip().replace('"',"").split("|")
    spliceai_info_index = vcf["SpliceAI"]["Description"].strip().split("Format:")[1].strip().replace('"',"").split("|")
    paths = os.path.join(os.path.dirname(str(__cur_path__)),'model/PDIVAS.sav')
    clf = pickle.load(open(paths, 'rb'))
    features = ['AI_del_gain_mean','AI_del_gain_max','AI_raw_gain_mean',"ConSplice","MES_max"]
    
    vcf.add_info_to_header(
        {"ID" :"PDIVAS",
        "Description" : "Predictor of Deep-Intronic Variant causing Aberrant Splicing (PDIVAS). Format: GENE_ID|PDIVAS",
        "Type":"String",
        "Number":"."})

    try:
        oup_file = args.O
        scored_vcf = (
            Writer(oup_file, vcf, "w") if oup_file.endswith("vcf")
            else Writer(oup_file, vcf, "wz") if oup_file.endswith("vcf.gz")
            else Writer(oup_file, vcf, "w")
                    )

    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    itera = 0
    print("Start calculating PDIVAS")
    for var in vcf:
        itera += 1
        var_oup = scoring_to_vcf(var,vep_info_index,spliceai_info_index,clf,features,args.F)

        if var_oup != None :
            scored_vcf.write_record(var_oup)
        
        else :
            print("scoring error at row",itera,"(",var.CHROM+":"+str(int(var.start)+1)+":"+var.REF+":"+var.ALT[0],")")
            continue
            
    print("PDIVAS successfully executed")
    vcf.close()
    scored_vcf.close()



def vcf2tsv(args):
    
    try:
        vcf = VCF(args.I)
    
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    print("Start converting VCF file to TSV file (1 gene annotation per 1 line)")
    
    #Set input/output files
    inp_file = args.I
    vcf = VCF(inp_file)

    fir_var = next(vcf, None)
    fmt_tags = fir_var.FORMAT if fir_var and fir_var.FORMAT else []

    for var in vcf:
        current_tags = var.FORMAT if var.FORMAT else []
        for tag in current_tags:
            if tag not in fmt_tags:
                fmt_tags.append(tag)
    vcf.close()

    vep_info_index = vcf["CSQ"]["Description"].strip().split("Format:")[1].strip().replace('"',"").split("|")
    spliceai_info_index = vcf["SpliceAI"]["Description"].strip().split("Format:")[1].strip().replace('"',"").split("|")
    PD_info_index = ["PD_ID","PDIVAS"]
    SM_index = vcf.samples if vcf.samples else []

    oup_col = ["chr","coordinate","ID","ref","alt","qual","filt"]+vep_info_index+spliceai_info_index+PD_info_index+SM_index
    hgvs_idx = [i for i, x in enumerate(oup_col) if x == 'Gene'][0]

    oup_file = args.O
    csvfile = open(oup_file, "w")
    csv_writer = csv.writer(csvfile, delimiter='\t')
    csv_writer.writerow(["chr","coordinate","ID","ref","alt","qual","filt"]+\
                        vep_info_index+spliceai_info_index+PD_info_index+SM_index)

    itera = 0
    filt_vep_error = 0
    filt_ai_error = 0
    uniq_annots = set([])
    uniq_vars = set([])

    vcf = VCF(inp_file)
    for var in vcf:
        itera += 1
        
        if (var.INFO.get("CSQ") is None) :
            print(var.CHROM, var.start, var.end, var.ID)
            filt_vep_error += 1
            continue

        if (var.INFO.get("SpliceAI") is None) :
            filt_ai_error += 1
            elem = vep_editor(var,vep_info_index,spliceai_info_index,PD_info_index)
        
        else :
            elem = vep_ai_matcher(var,vep_info_index,spliceai_info_index,PD_info_index)
            
        if len(elem) == 0 :
            print("Output error")
            break

        if len(elem) != 0:
            var_ID = ":".join([str(k) for k in [var.CHROM,var.start+1,var.REF,var.ALT[0]]])
            uniq_vars.add(var_ID)

        #Sample columns
        SMs = []
        if len(SM_index)>=1:
            for idx, sample in enumerate(SM_index):
                SM = []
                if "GT" in fmt_tags:
                    GT =str(var.gt_types[idx]) if var.gt_types[idx]!=2 else "NA"
                    GT = "2" if GT =="3" else GT
                    SM.append(GT)
                if "AD" in fmt_tags:
                    AD =str(var.gt_alt_depths[idx])+","+str(var.gt_ref_depths[idx]) if var.gt_alt_depths[idx] != -1 else "NA,NA" 
                    SM.append(AD)
                for tag in [i for i in fmt_tags if i not in ["GT","AD"]] :
                    sm = var.format(tag)[idx] if var.format(tag) is not None else "NA"
                    if isinstance(sm, np.ndarray) :
                        sm = ",".join([str(i) if i > -200000000 else "NA" for i in sm])
                    elif isinstance(sm, np.str_) :
                        sm = sm if sm != "." else "NA"
                    SM.append(sm)
                SMs.append(":".join(SM))

        for gen in range(len(elem)) :
            annot_ID = ":".join([str(k) for k in [var.CHROM,var.start+1,var.REF,var.ALT[0],elem[gen][hgvs_idx]]])
            if annot_ID in uniq_annots :
                print("duplication detected @",annot_ID)
                continue
            uniq_annots.add(annot_ID)
        
            oup_list = elem[gen]+SMs
            csv_writer.writerow(oup_list)


    vcf.close()
    csvfile.close()

    print("# of input variants : ",itera)
    print("# of VEP-lacked variants :",filt_vep_error,";",round(100*filt_vep_error/itera,2),"(%)")
    print("# of AI-lacked variants :",filt_ai_error,";",round(100*filt_ai_error/itera,2),"(%)")
    print("# of output unique variants :",len(uniq_vars),";",round(100*len(uniq_vars)/itera,2),"(%)")
    print("# of output unique annotations : ",len(uniq_annots))

def main(args=None):
    args,parser = get_options()
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        # show help for unknown subcommands
        parser.print_help()

if __name__ == '__main__':
    main()
