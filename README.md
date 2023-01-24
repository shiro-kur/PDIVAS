# PDIVAS : Pathogenic Predictor of Deep-Intronic Variants causing Aberrant Splicing
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![PDIVAS image](/PDIVAS.png)

## Sumary
- PDIVAS is a pathogenic predictor of deep-intronic variants causing aberrant splicing.
- The deep-intronic variants can cause pathogenic pseudoexons or extending exons which disturb the normal gene expression and can be the causal of patiens with Mendelian diseases. 
- PDIVAS efficiently prioritizes the causal candidates from enourmous deep-intronic variants detected by whole-genome sequencing. 
- The scope of PDIVAS prediction is variants in protein-coding genes on autosomes and X chromosome. 
- This command line interface is compatible with variant file on VCF format. 
 
PDIVAS is modeled on random forest algorism using features from  
1) **Splicing predictors** of [SpliceAI](https://github.com/Illumina/SpliceAI) ([Jaganathan et al., Cell 2019](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub)) and [MaxEntScan](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) ([Yao and Berge, j. Comput. Biol. 2004](https://www.liebertpub.com/doi/10.1089/1066527041410418?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed))  
(*)The output module of SpliceAI wass customed for PDIVAS features (see the Option2, for the details).
          
 2) **Human splicing constraint score** of [ConSplice](https://github.com/mikecormier/ConSplice) ([Cormier et al., BMC Bioinfomatics 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05041-x)).

## Reference & contact

## Option1. Prediction with the PDIVAS-precomputed files (SNV+ short indels(1~4nt))
For the quick implementation of PDIVAS, please use the score-precomputed file [file here]().
Possible rare SNVs and short indels(1~4nt) in genes (n=4,512) of Mendelian diseases were comprehensively annotated in the file.
To annotate your VCF file, please run the command below,for example.

**0. Installation**
```sh
conda install -c bioconda tabix bcftools
(or conda install pdivas (including xsamtools))
```

**1. Peform PDIVAS prediction**
```sh
bgzip -c examples/input.vcf > example/input.vcf.gz
tabix examples/input.vcf.gz
#The site file enables quick annotation
bcftools query -f'%CHROM\t%POS\n' examples/input.vcf.gz > examples/input_sites.txt
bcftools annotate -c 'INFO/PDIVAS' -a PDIVAS_snv_precomputed_GRCh38.vcf.gz -R examples/input_sites.txt examples/input.vcf.gz | bgzip -c > examples/output_precomp.vcf.gz
#Compare the output_precomp.vcf.gz with output_precomp_expect.vcf.gz to validate the succcessful annotation.
```

## Option2. Peform annotation of individual features and calculation of PDIVAS scores 
For more complehensive annotation than pre-computed files, run PDIVAS by following the description below.

**0-1. Installation**
```sh
conda install PDIVAS_env.yml
conda activate PDIVAS
pip install PDIVAS
```
The successful installation was verified on anaconda(version=.....)

**0-2. Setting custom usages**
For output-customed SpliceAI
```sh
git clone git@github.com:shiro-kur/PDIVAS.git
cd PDIVAS/Customed_SpliceAI
cp ./__main__for_customed_SpliceAI.py path_to_your_installed_path/__main__.py
cp ./utils_for_customed-SpliceAI.py path_to_your_installed_path/__utils__.py
cp -r ./annotations_for_customed_SpliceAI path_to_your_installed_path/annotations

# check the successful custom by comparing the output file between ~~.vcf
```
For VEP custom usage,
- Donwload VEP cache file (version>=107, should correspond to VEP tool version).  
Follow the instruction of "Manually downloading caches" part below.  
(https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
- To implement MaxEntScan plugin, follow the instruction below.  
(https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#maxentscan)
- Download ConSplice score file from [here]().  
The file was editted from the original score file of ([Cormier et al., BMC Bioinfomatics 2022](https://home.chpc.utah.edu/~u1138933/ConSplice/best_splicing_constraint_model/)).

**1. Preprocessing VCF format (resolve the mullti-allelic site to biallelic sites)**
```sh
bcftools norm -m - multi.vcf > bi.vcf
```

**2. Add gene annotations, MaxEntScan scores and ConSplice scores with VEP.**
```sh
vep \
--cache --offline --cache_version 107 --assembly GRCh38 --hgvs --pick_allele_gene \
--fasta ../../reference/hg38.fa.gz --vcf --force \
--custom ../../reference_output/ConSplice.50bp_region.inverse_proportion_refor.bed.gz,ConSplice,bed,overlap,0 \
--plugin MaxEntScan,../../reference/MaxEntScan/fordownload,SWA,NCSS \
--fields "Consequence,SYMBOL,Gene,INTRON,HGVSc,STRAND,ConSplice,MES-SWA_acceptor_diff,MES-SWA_acceptor_alt,MES-SWA_donor_diff,MES-SWA_donor_alt" \
--compress_output bgzip
-i examples/input.vcf.gz -o examples/input_vep.vcf.gz
```

**3. Add output-customed SpliceAI scores**
```sh
spliceai -I examples/input_vep.vcf.gz -O examples/input_vep_AI.vcf -R hg38.fa -A grch38 -D 300 -M 1
```

**4. Perform the detection of deep-intronic variants and PDIVAS prediction**
```sh
pdivas -I input_vep_AI.vcf -O input_vep_AI_PD.vcf.gz -F off
```
## Usage of PDIVAS command line
Required parameters:
 - ```-I```: Input VCF(.vcf/.vcf.gz) with variants of interest.
 - ```-O```: Output VCF(.vcf/.vcf.gz) with PDIVAS predictions `GENE_ID|PDIVAS_score` Variants in multiple genes have separate predictions for each gene.
Optional parameters:
 - ```-F```: filtering function (off/on) : Output all variants (-F off; default) or only deep-intronic variants with PDIVAS scores (-F on)")
 
 Details of PDIVAS INFO field:

|    ID    | Description |
| -------- | ----------- |
|  GENE_ID  | Ensembl gene ID based on GENCODE V41(GRCh38) or V19(GRCh37) |
|  PDIVAS  | \<Predicted result\> <br> **Pattern 1 : 0.000-1.000 float value**  (The higher, the more deleterious) <br> \<Exceptions\> <br> - Output with '-F off'. Filtered with '-F on'. <br> **Pattern 2 : 'wo_annots'**, variants out of VEP or SpliceAI annotations : <br>**Pattern 3 : 'out_of_scope'**, variants without PDIVAS annotation scope<br>       (chrY, non-coding gene or non-deep-intronic variants)　<br>**Pattern 4 :'no_gene_match'**, variants without matched gene annotation between VEP and SpliceAI|

## Interpretation of PDIVAS score

