sed '/chrM/d;/random/d;/chrUn/d' $1.sam > $1_filtered.sam

samtools view -Sb $1_filtered.sam > $bam

samtools view -c $1.bam #count reads of bam files
