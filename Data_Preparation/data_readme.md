The scripts in this directory are used to prepare data for CNN training and the clustering process. The following open source data are being imported in different stages of the process.
* gencode.v38lift37.annotation.gtf
* hg38.fa
* Individual motif files of RBPs from POSTAR3/human.txt


The steps are as follows, run the scripts in the following order.
1. First script to run, 'extract_gene.py' uses 'gencode.v38lift37.annotation.gtf' file to extract gene start and end position to identify RBP location.
2. Secondly, run 'rbp_pairs_genes.py' file.
3. Next, 'get_positive.py' & 'get_negative.py' to accumulate 2 class data as per definition; however it is only the positions of the binding sites.
4. Finally, run 'get_sequence.py' to acquire nucleotide sequences which is the defined input data form.
