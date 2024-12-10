mamba activate ee282
cd HW4
#file Integrety
md5sum dmel-all-chromosome-r6.48.fasta.gz
# check with  md5sum.txt
#Calculate Genome
zcat dmel-all-chromosome-r6.48.fasta.gz | faFilter -maxSize=100000 stdin small.fa
zcat dmel-all-chromosome-r6.48.fasta.gz | faFilter -minSize=100001 stdin large.fa
#Total number of nucleotides & sequences
faSize small.fa
faSize large.fa
#
bioawk -c fastx 'length($seq) <=100000 { print $name "\t" length($seq) "\t" gc($seq)}' dmel-all-chromosome-r6.48.fasta.gz > small.txt
bioawk -c fastx 'length($seq) >100000 { print $name "\t" length($seq) "\t" gc($seq)}' dmel-all-chromosome-r6.48.fasta.gz > large.txt
Rscript -e 'small <- read.table("small.txt");
small_length <- small[,2];
small_gc<- small[,3];
png(file = "small_length_hist.png")
hist(small_length, breaks=100, main="Length Dist (¡Ü100kb)", xlab="Length (bp)")
dev.off();
png(file = "small_gc_hist.png")
hist(small_gc, breaks=100, main="gc% (¡Ü100kb)", xlab="Length (bp)")
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

bioawk -c fastx 'length($seq) <=100000 { print length($seq)}'  dmel-all-chromosome-r6.48.fasta.gz | sort -nr > small_descending.txt
bioawk -c fastx 'length($seq) >100000 { print length($seq)}'  dmel-all-chromosome-r6.48.fasta.gz | sort -nr > large_descending.txt
awk '{sum+=$1; print sum}' small_descending.txt > small_cumulative_sizes.txt
awk '{sum+=$1; print sum}' large_descending.txt > large_cumulative_sizes.txt
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

