mamba activate ee282
cd HW4
#assembly
srun -c 16 -A class_ee282 --pty bash -i
hifiasm -o ios1.asm -t 16 -f 0 ISO1_Hifi_AdaptorRem.40X.fasta.gz 2> log
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
# compare CDF
faSplitByN dmel-all-chromosome-r6.48.fasta contig.fasta 10
bioawk -c fastx '{ print $name "\t" length($seq) "\t" gc($seq)}' dmel-all-chromosome-r6.48.fasta.gz > scaffold.txt
bioawk -c fastx '{ print $name "\t" length($seq) "\t" gc($seq)}' contig.fasta > contig.txt
bioawk -c fastx '{ print $name "\t" length($seq) "\t" gc($seq)}' iso1.asm.bp.p_ctg.fasta > myassembly.txt
# busco score
compleasm run -a contig.fa -o compleasm/ -l diptera_odb10 -L compleasm_my/
compleasm run -a ISO1.p_ctg.fa -o compleasm/ -l diptera_odb10 -L compleasm_my/