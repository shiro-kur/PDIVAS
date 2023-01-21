# PDIVAS : Pathogenic Predictor of Deep-Intronic Variants causing Aberrant Splicing
Predict the deleterious effect of variant on RNA splicing from VCF (variant call format) file.

![PDIVAS image](/PDIVAS.png)
## Availability
To perform PDIVAS prediction on your environment, there are two options.
### 1.Refer to the PDIVAS-precomputed files (SNV+ short indels(1~4nt))
File link to the PDIVAS-precomputed files
To annotate your VCF file with the precomputed file, please run the command below,for example.

### 2.Peform individual feature annotation and PDIVAS prediction for your VCF 
**2-0. Installation**
```sh
conda install PDIVAS_env.yml
conda activate PDIVAS
pip install PDIVAS
```
The successful installation was verified on anaconda(version=.....)

output-customed SpliceAI

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

**2-3. Annotate output-customed SpliceAI scores**
```sh
spliceai -I ../../data_output/ex_inp_vep.vcf.gz -O ../../data_output/ex_inp_vep_AI.vcf -R ../../reference/hg38.fa -A grch38 -D 300 -M 1
bgzip ../../data_output/ex_inp_vep_AI.vcf
```

**2-4. Perform the detection of deep-intronic varaint and PDIVAS prediction**
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
|  PDIVAS_score  | \<Predicted result\> <br> **Pattern 1 : 0.000-1.000 float value**  (The higher, the more deleterious) <br> \<Exceptions\> <br> - Output if '-F off'. Filtered if '-F on' <br> **Pattern 2 : 'wo_annots'**, variants without VEP or SpliceAI annotations : <br>**Pattern 3 : 'out_of_scope'**, variants without PDIVAS annotation scope (chrY, non-coding gene or non-deep-intronic variants)　<br>**Pattern 4 :'no_gene_match'**, variants without matched gene annotation between VEP and SpliceAI|


