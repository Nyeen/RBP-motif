# RBP-motif
RNA binding proteins (RBPs) have unique features for RNA binding regions. Often, they form protein clusters which are biologically significant for cellular mechanism of proteins and our primary objective to identify. Here, our target domain is 28 distinct RBPs across 5 biological categories in human species. We have extended the co-occurrance concept of DNA transcription factor binding sites to learn specificity of RBP pairs in clusters using CNN, achieving a validation accuracy of 92.61%. Finally, the trained CNN prepares a scored seed-corpus for seed and extend method to search co-regulatory RBPs on RNA sequence alignment.

## Datasets
For data analysis, feature engineering, training CNN and clustering process, we have used the following public datasets.
### POSTAR3 [1]
It is an updated CLIP-seq database consisting of RBP binding sites information for 7 species, including human, mouse, zebra fish, fly, worm, Arabidopsis, and yeast. This is the source of "human.txt" file.
### UCSC genome browser [2]
It is a very well-known bio-informatics resource or database for genomics research, containing various annotations along with whole sequence data regularly updated. This is the source of "hg38.fa" file.
### GENCODE [3]
It is also an up-to-date human and mouse gene annotations database. This is the source of "gencode.v38lift37.annotation.gtf" file.
Collect data (.bigWig) from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSunyRipSeq/ and convert them to .bed files.
To be updated...
