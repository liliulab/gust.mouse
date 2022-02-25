## GUST.mouse (Genes Under Selection in Tumors, rebuilt for mouse)
This is the GUST model rebuilt to make predictions in mouse.
GUST predicts oncogenes (OGs), tumor suppressor gens (TSGs) and passenger genes (PGs) using somatic SNVs detected in a collection of tumors. It can classify single genes from targetted sequencing or multiple genes from whole-exome sequencing. 

## Install. 
````
# Be patient. It may take some time.
 devtools::install_github('liliulab/gust.mouse')
````

## Usage
```` 
library(gust.mouse)

# to make predictions. 
# A sample input can be downloaded from the "examples" folder. 
 gust.mouse(input.file.name='./examples/ACC.mm.txt', output.folder='./examples/', output.prefix='ACC.mm');

# Compare outputs from the above command with files in the "examples" folder.
 
# to plot distribution
 m=plot.this.one('Ctnnb1', './examples/', 'ACC.mm') 
````

The input file needs to be in the VCF format. Required fields include Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, and Entrez_Gene_Id. We recommend using SnpEff to annotate somatic variants, which will generate a VCF file containing the required fields.
Additional usages can be found via using "?gust.mouse" after installing the package.

## Reference
Please cite our publication in Bioinformatics.

Chandrashekar P, Ahmadinejad N, Sekulic A, Wang J, Kumar S, Maley C, Liu L. (2019) Somatic selection distinguishes oncogenes and tumor suppressor genes. Bioinformatics. 36(6):1712-1717 
Open Access Article is available at https://academic.oup.com/bioinformatics/article/36/6/1712/5625621

## Contributors
Algorithm of GUST was developed by Li Liu. Database was developed by Pramod Chandrashekar. Please contact liliu at asu.edu for any questions or suggestions.