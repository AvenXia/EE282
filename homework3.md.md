---
title: Homework3
Name: Qingyuan Xia
---

# Part 1: Summarize Genome Assembly
## File Integrity
The file integrity is checked by function md5sum. We can check the digit we get from the code with txt file on the web. (Check the code from hw3_genome_summary.sh)
``` sh
 md5sum dmel-all-chromosome-r6.48.fasta.gz
```
It turns out that the digit is correct
## Calculate Summaries of the Genome
The genome can be calculated and summarized by faSize. (Check the code from hw3_genome_summary.sh)
``` sh
faSize dmel-all-chromosome-r6.48.fasta.gz
```
The results are: 
Total number of nucleotides: 143726002
Total number of Ns: 1152978
Total number of sequences: 1870

# Part2: Summarize an Annotation File
## File Integrity
The file integrity is checked by function md5sum. We can check the digit we get from the code with txt file on the web. (Check the code from hw3_annotation_summary.sh)
``` sh
 md5sum dmel-all-r6.48.gtf.gz
```
It turns out that the digit is correct
## Compile a Report Summarizing the Annotation
Summarizing the annotation can be achieved by bioawk.
### Total number of features of each type
Here $feature is printed. (Check the code from hw3_annotation_summary.sh)
``` sh
bioawk -c gff ' { print $feature } ' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -k1,1nr
```
The results are:
  190050 exon
 163242 CDS
  46802 5UTR
  33738 3UTR
  30885 start_codon
  30825 stop_codon
  30799 mRNA
  17896 gene
   3053 ncRNA
    485 miRNA
    365 pseudogene
    312 tRNA
    300 snoRNA
    262 pre_miRNA
    115 rRNA
     32 snRNA
Here another  | sort -k1,1nr is added to re-order the features from high number to low number.
```
### Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)
Here let's take "gene" as feature and print seq name, and include uniq and sort to find the numbers of genes from each chromosome type. (Check the code from hw3_annotation_summary.sh)
```sh
bioawk -c gff ' $feature == "gene" { print $seqname } ' dmel-all-r6.48.gtf.gz | sort | uniq -c
```
The results are:
3515 2L
   3653 2R
   3489 3L
   4227 3R
    114 4
    38 mitochondrion_genome
    21 rDNA
    2 Unmapped_Scaffold_8_D1580_D1567
   2708 X
    113 Y


 



