---
title: Homework4
Name: Qingyuan Xia
---

# Part 1: Summarize partitions of a genome assembly
## Calculate the following for all sequences ≤ 100kb and all sequences > 100kb: 1.Total number of nucleotides 2. Total number of Ns 3. Total number of sequences
The file of fly all chromosome from homework3 is used again here for analysis. The file integrity is checked by function md5sum. 
The first step is to filter the seq based on the seq length. Here faFilter is utilized to set max or minsize of the length so that we can seperate the whole data into two: small.fa and large.fa
``` sh
zcat dmel-all-chromosome-r6.48.fasta.gz | faFilter -maxSize=100000 stdin small.fa
zcat dmel-all-chromosome-r6.48.fasta.gz | faFilter -minSize=100001 stdin large.fa
faSize small.fa
faSize large.fa
```
By doing this, we can get the answer:
for ≤ 100kb:
Total number of nucleotides: 6178042
Total number of Ns: 662593
Total number of sequences: 1873

for > 100kb:
Total number of nucleotides: 137547960
Total number of Ns: 490385
Total number of sequences: 7

## Plots of the following for for all sequences ≤ 100kb and all sequences > 100kb:
To filter the seq with specific length, the data is summarized based on bioawk, where we can filter length($seq) and calculage both length and gc proportion by length($seq) and gc()
``` sh
bioawk -c fastx 'length($seq) <=100000 { print $name "\t" length($seq) "\t" gc($seq)}' dmel-all-chromosome-r6.48.fasta.gz > small.txt
bioawk -c fastx 'length($seq) >100000 { print $name "\t" length($seq) "\t" gc($seq)}' dmel-all-chromosome-r6.48.fasta.gz > large.txt
```
After getting the summary table: small and large txt, we can use R to do hist to draw the plot:
``` sh
Rscript -e 'small <- read.table("small.txt");
small_length <- small[,2];
small_gc<- small[,3];
png(file = "small_length_hist.png")
hist(small_length, breaks=100, main="Length Dist (≤100kb)", xlab="Length (bp)")
dev.off();
png(file = "small_gc_hist.png")
hist(small_gc, breaks=100, main="gc% (≤100kb)", xlab="Length (bp)")
dev.off()'

Rscript -e 'large <- read.table("large.txt");
large_length <- large[,2];
large_gc<- large[,3];
png(file = "large_length_hist.png")
hist(large_length, breaks=7, main="Length Dist (>100kb)", xlab="Length (bp)")
dev.off();
png(file = "large_gc_hist.png")
hist(large_gc, breaks=7, main="gc% (>100kb)", xlab="Length (bp)")
dev.off()'
```
The figure is showed in github. It turns out that the seq length is close to around 30-40% for seq <= 100kb.
To do the CDF plot, our strategy is to first order seq length into decreasing order. Then use the function cumsum from R to do the cumulation, and finally plot the line.
``` sh
Rscript -e 'small <- read.table("small.txt");
small_lengths_descending <- sort(small$V2, decreasing = TRUE);
small_cum_lengths <- cumsum(small_lengths_descending);
png("small_cum_length.png")
plot(small_cum_lengths, 
     type="l",               # 'l' for a line plot
     xlab="Sequence Rank (largest to smallest)", 
     ylab="Cumulative Size (bp)", 
     main="Cumulative Sequence Size Distribution (<=100kb)",
     las=1)                  # 'las=1' to make axis labels horizontal
dev.off()'
Rscript -e 'large <- read.table("large.txt");
large_lengths_descending <- sort(large$V2, decreasing = TRUE);
large_cum_lengths <- cumsum(large_lengths_descending);
png("large_cum_length.png")
plot(large_cum_lengths, 
     type="l",               # 'l' for a line plot
     xlab="Sequence Rank (largest to smallest)", 
     ylab="Cumulative Size (bp)", 
     main="Cumulative Sequence Size Distribution (>100kb)",
     las=1)                  # 'las=1' to make axis labels horizontal
dev.off()'
```
For CDF, since the length is orded in decreasing order, the line will first increase and then stay the same to a specific value as maximum. The figures are shown in the github.

# Part2: Summarize an Annotation File
## Assemble a genome using Pacbio HiFi reads
To do the assembly, we first copy the data from ee282 file and use hifiasm to do the assembly. SRUN is used to make code running faster.
``` sh
srun -c 16 -A class_ee282 --pty bash -i
hifiasm -o ios1.asm -t 16 -f 0 ISO1_Hifi_AdaptorRem.40X.fasta.gz 2> log
```
A gfa file is obtained after the running
## Assembly assessment
### Calculate the N50
This can be done with faSize+awk+sort:
``` sh
awk '$1 == "S" { print ">" $2; print $3 }' iso1.asm.bp.p_ctg.gfa > iso1.asm.bp.p_ctg.fasta
bioawk -c fastx '{ print length($seq) }'  iso1.asm.bp.p_ctg.fasta \
    | sort -nr \
    | awk '{
        total += $1
        lengths[NR] = $1
    }
    END {
        half = total/2
        csum = 0
        for (i = 1; i <= NR; i++) {
            csum += lengths[i]
            if (csum >= half) {
                print lengths[i]
                exit
            }
        }
    }'
```
The final result of N50 is 21715751, which is a little bit different from the reference from the link.
### Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot.
To do this, since our previous data from flybase is scaffold assembly, contig assembly should be obtained by converting previous data to contig one. This can be achieved by using faSplitByN 
``` sh
faSplitByN dmel-all-chromosome-r6.48.fasta contig.fasta 10
```
Then we can summarize both types of assembly by bioawk
``` sh
bioawk -c fastx '{ print $name "\t" length($seq) "\t" gc($seq)}' dmel-all-chromosome-r6.48.fasta.gz > scaffold.txt
bioawk -c fastx '{ print $name "\t" length($seq) "\t" gc($seq)}' contig.fasta > contig.txt
bioawk -c fastx '{ print $name "\t" length($seq) "\t" gc($seq)}' iso1.asm.bp.p_ctg.fasta > myassembly.txt
```
Using R as mentioned above, we can plot CDF of three assemblies together
``` R
scaffold <- read.table("scaffold.txt")
scaffold_lengths_descending <- sort(scaffold$V2,decreasing = TRUE)
scaffold_cum_lengths <- cumsum(scaffold_lengths_descending)
contig <- read.table("contig.txt")
contig_lengths_descending <- sort(contig$V2,decreasing = TRUE)
contig_cum_lengths <- cumsum(contig_lengths_descending)
my <- read.table("myassembly.txt")
my_lengths_descending <- sort(my$V2,decreasing = TRUE)
my_cum_lengths <- cumsum(my_lengths_descending)

png("compare_cum_length.png")
plot(scaffold_cum_lengths, 
     type="l",               # 'l' for a line plot
     xlab="Sequence Rank (largest to smallest)", 
     ylab="Cumulative Size (bp)", 
     main="Cumulative Sequence Size Distribution",
     las=1,    # 'las=1' to make axis labels horizontal
     log='x',
     col="red")
lines(contig_cum_lengths,col="blue")
lines(my_cum_lengths,col="green")
dev.off()
```
The figure will be shown in the github. Red is scaffold assembly, blue contig assembly while green is our own assembly. Interestingly, Our assembly is quite different from the true ones as our cumulative mamimum is larger and slower to be reached.
### busco score
We use compleasm here to calculate the busco score.
``` sh
compleasm run -a ISO1.p_ctg.fa -o compleasm/ -l diptera_odb10 -L compleasm_my/
compleasm run -a contig.fa -o compleasm/ -l diptera_odb10 -L compleasm_my/
```
For the reference:
S:99.70%, 3275
    D:0.30%, 10
    F:0.00%, 0
    I:0.00%, 0
    M:0.00%, 0
    N:3285

For my assembly:
S:99.63%, 3273
    D:0.24%, 8
    F:0.00%, 0
    I:0.00%, 0
    M:0.12%, 4
    N:3285