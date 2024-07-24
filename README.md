# RBP-motif
RNA binding proteins (RBPs) have unique features for RNA binding regions. Often, they form protein clusters which are biologically significant for cellular mechanism of proteins and our primary objective to identify. Here, our target domain is 28 distinct RBPs across 5 biological categories in human species. We have extended the co-occurrance concept of DNA transcription factor binding sites to learn specificity of RBP pairs in clusters using CNN, achieving a validation accuracy of 92.61%. Finally, the trained CNN prepares a scored seedcorpus for seed and extend method to search co-regulatory RBPs on RNA sequence alignment.

## Training Data
Collect data (.bigWig) from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSunyRipSeq/ and convert them to .bed files.
To be updated...
