# PDIVAS : Pathogenic Predictor of Deep-Intronic Variants causing Aberrant Splicing
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![PDIVAS image](/PDIVAS.png)

## Sumary
PDIVAS is a pathogenic predictor of deep-intronic variants causing aberrant splicing. The deep-intronic variants can cause pathogenic pseudoexons or extending exons which disturb the normal gene expression and can be the causal of patiens with Mendelian diseases. PDIVAS efficiently prioritizes the causal candidates from enourmous deep-intronic variants detected by whole-genome sequencing. This command line interface is compatible with variant file on VCF format.
 
PDIVAS is modeled on random forest algorism using features from 
 1. Splicing predictors of [SpliceAI](https://github.com/Illumina/SpliceAI) ([Jaganathan et al., Cell 2019](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub)) and [MaxEntScan](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) ([Yao and Berge, j. Comput. Biol. 2004](https://www.liebertpub.com/doi/10.1089/1066527041410418?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed))
 2. Human splicing constraint score of [ConSplice](https://github.com/mikecormier/ConSplice) ([Cormier et al., BMC Bioinfomatics 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05041-x)).

## Availability
To perform PDIVAS prediction on your environment, there are two options.
### 1.Prediction with the PDIVAS-precomputed files (SNV+ short indels(1~4nt))
File link to the PDIVAS-precomputed files
To annotate your VCF file with the precomputed file, please run the command below,for example.
**1-0. Installation**
```sh
conda install -c bioconda tabix bcftools
(or conda install pdivas (including xsamtools))
```

**1-1. Peform PDIVAS prediction**
```sh
bgzip input.vcf.gz
tabix input.vcf.gz
bcftools annotate -c 'INFO/PDIVAS' -a PDIVAS_snv_precomputed_GRCh38.vcf.gz input.vcf.gz | bgzip -c > output.vcf.gz
```

### 2.Peform annotation of individual features and calculation of PDIVAS scores 
**2-0-1. Installation**
```sh
conda install PDIVAS_env.yml
conda activate PDIVAS
pip install PDIVAS
```
The successful installation was verified on anaconda(version=.....)

**2-0-2. Setting custom usages**
For output-customed SpliceAI
```sh
cp ./
cp ./
# check the successful custom by comparing the output file between ~~.vcf

```
For VEP custom usage,
- Donwload VEP cache file (version>=107, should correspond to VEP tool version).
Follow the instruction of "Manually downloading caches" part below.
https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html
- Download ConSplice.50bp_region.inverse_proportion_refor.bed.gz from ....
- To implement MaxEntScan plugin, follow the instruction below.
https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#maxentscan
- 

**2-1.Preprocessing VCF format (resolve the mullti-allelic site to single sites)**
```sh
bcftools norm 
```

**2-2.Add gene annotations, MaxEntScan scores and ConSplice scores with VEP.**
```sh
vep \
--cache --offline --cache_version 107 --assembly GRCh38 --hgvs --pick_allele_gene \
--fasta ../../reference/hg38.fa.gz --vcf --force \
--custom ../../reference_output/ConSplice.50bp_region.inverse_proportion_refor.bed.gz,ConSplice,bed,overlap,0 \
--plugin MaxEntScan,../../reference/MaxEntScan/fordownload,SWA,NCSS \
--fields "Consequence,SYMBOL,Gene,INTRON,HGVSc,STRAND,ConSplice,MES-SWA_acceptor_diff,MES-SWA_acceptor_alt,MES-SWA_donor_diff,MES-SWA_donor_alt" \
--compress_output bgzip
-i ../../input/ex_inp.vcf.gz -o ../../data_output/ex_inp_vep.vcf.gz
```

**2-3. Add output-customed SpliceAI scores**
```sh
spliceai -I ../../data_output/ex_inp_vep.vcf.gz -O ../../data_output/ex_inp_vep_AI.vcf -R ../../reference/hg38.fa -A grch38 -D 300 -M 1
bgzip ../../data_output/ex_inp_vep_AI.vcf
```

**2-4. Perform the detection of deep-intronic varaints and PDIVAS prediction**
```sh
pdivas -I ../data_output/ex_inp_vep_AI.vcf.gz -O ../data_output/ex_inp_vep_AI_PD.vcf.gz
```
### Usage of PDIVAS command line
Required parameters:
 - ```-I```: Input VCF with variants of interest.
 - ```-O```: Output VCF with PDIVAS predictions `GENE_ID|PDIVAS_score` Variants in multiple genes have separate predictions for each gene.
Optional parameters:
 - ```-F```: filtering function (off/on) : Output all variants (-F off; default) or only deep-intronic variants with PDIVAS scores (-F on)")
 
 Details of SpliceAI INFO field:

|    ID    | Description |
| -------- | ----------- |
|  GENE_ID  | Ensembl gene ID based on GENCODE V...(GRCh38) V...(GRCh37) When  |
|  PDIVAS_score  | \<Predicted result\> <br> **Pattern 1 : 0.000-1.000 float value**  (The higher, the more deleterious) <br> \<Exceptions\> <br> - Output with '-F off'. Filtered with '-F on'. <br> **Pattern 2 : 'wo_annots'**, variants out of VEP or SpliceAI annotations : <br>**Pattern 3 : 'out_of_scope'**, variants without PDIVAS annotation scope (chrY, non-coding gene or non-deep-intronic variants)　<br>**Pattern 4 :'no_gene_match'**, variants without matched gene annotation between VEP and SpliceAI|


