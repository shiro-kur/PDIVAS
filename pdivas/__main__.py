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

from .__init__ import __cur_path__, __version__
from pdivas.scoring_to_vcf import scoring_to_vcf

#https://qiita.com/kzkadc/items/e4fc7bc9c003de1eb6d0
#http://www.yamamo10.jp/yamamoto/comp/Python/library/argparse/sample/prog/index.html
#https://qiita.com/tanabe13f/items/6c09f8f71eb2efb1ac75

def get_options():
	parser = argparse.ArgumentParser(
		prog="PDIVAS",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=(
		"\n\t=====================================================\n"
		"\t|| PDIVAS, Pathogenic predictor of Deep-Intronic Variants causing Aberrant Splicing ||\n"
		"\t=====================================================\n\n"),
						)
						
	parser.add_argument(
		"-I",
		required=True,
		metavar="input",
		type=str,
		help="The path to the vcf(.gz) file to add PDIVAS annotation")
		
		
	parser.add_argument(
		"-O",
		required=True,
		metavar="output",
		type=str,
		help="The path to output file name and pass")
		
		
	parser.add_argument(
		"-F",
		metavar="filtering",
		required=False,
		default="off",
		type=str,
		choices=["off","on"],
		help="Output all variants (-F off; default) or only deep-intronic variants with PDIVAS scores (-F on)")

	args = parser.parse_args()
	return args
    
#https://www.javadrive.jp/python/userfunc/index5.html
def main(args=None):
	args = get_options()
	if None in [args.I, args.O]:
		logging.error('Usage: pdivas [-h] [-I [input]] [-O [output]]  [-F [True/False]]')
		exit()
		
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
		"Description" : "Predictor of Deep-Intronic Variant causing Aberrant Splicing (PDIVAS). Format: GENE_ID|PDIVAS_score",
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
	for var in vcf:
		itera += 1
		var_oup = scoring_to_vcf(var,vep_info_index,spliceai_info_index,clf,features,args.F)

		if var_oup != None :
			scored_vcf.write_record(var_oup)
		
		else :
			print("scoring error at row",itera,"(",var.CHROM+":"+str(int(var.start)+1)+":"+var.REF+":"+var.ALT[0],")")
			continue
			
	print("Succeeded")
	vcf.close()
	scored_vcf.close()

if __name__ == '__main__':
    main()
