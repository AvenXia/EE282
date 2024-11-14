mamba activate ee282
cd HW3
#file integrety
md5sum dmel-all-r6.48.gtf.gz
#Total number of features of each type, sorted from the most common to the least common
bioawk -c gff ' { print $feature } ' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -k1,1nr
#Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)
bioawk -c gff ' $feature == "gene" { print $seqname } ' dmel-all-r6.48.gtf.gz | sort | uniq -c

