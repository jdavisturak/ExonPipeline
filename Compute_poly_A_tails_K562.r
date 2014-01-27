## In human hg19:
cd /home/jeremy/RNAseq/Tilgner/Cytosol
sh ~/Code/filterPolyAtags.sh wgEncodeCshlLongRnaSeqK562CytosolPapFastqRd2Rep1.fastq Cyt1_rd2
sh ~/Code/filterPolyAtags.sh wgEncodeCshlLongRnaSeqK562CytosolPapFastqRd1Rep1.fastq Cyt1_rd1

bowtie2 --end-to-end --fast -p 20 -x  /home/RNAseq/indices/Homo_sapiens/UCSC/hg19/hg19 -U  Cyt1_rd1T.fastq,Cyt1_rd2A.fastq -S aligned1.sam --un unaligned1.fastq
                                                                                       
samtools view -bS -o aligned1.bam aligned1.sam
samtools sort aligned1.bam aligned1.sorted
bamToBed -i aligned1.sorted.bam > aligned1.all.bed
mergeBed -s -i aligned1.all.bed > aligned1.bed

## Now trim the unaligned1.fastq
awk -f ~/Code/trimPolyAtags.awk unaligned1.fastq > unaligned.trimmed.fastq

bowtie2 --end-to-end --fast -p 20 -x  /home/RNAseq/indices/Homo_sapiens/UCSC/hg19/hg19 -U  unaligned.trimmed.fastq -S aligned2.sam 

samtools view -bS -o aligned2.bam aligned2.sam
samtools sort aligned2.bam aligned2.sorted
bamToBed -i aligned2.sorted.bam > aligned2.all.bed
mergeBed -s -i aligned2.all.bed > K562_polyA_sites_rep1.bed

####
# Different data:
