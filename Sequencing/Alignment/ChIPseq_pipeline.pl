#usr/bin/perl
use strict;
use warnings;

###Declare variation
print("Write the project directory :\n");
my $p_directory = <STDIN>;
chomp($p_directory);
`mkdir $p_directory/trimmed`;
`mkdir $p_directory/trimmed/Fastqc`;

###Adaptor trimming
#Extract each fastq file drictory from $fastq_dir and save to @file_list
my @file_list = ();
push(@file_list, `find $p_directory/rawdata -name '*.fastq.gz'`);
@file_list = sort(@file_list); #sort in suffix
#check the @file_list
for (@file_list){
    chomp($_);
}
for (@file_list){
    print("$_\n");
}
#Run Trimgalore
for (my $i=0; $i<$#file_list+1; $i+=2){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/\.fastq\./, $Oname_1[-1]);
    print("Trimming file named $Oname_2[0]....\n");
#    `skewer -t 4 -z -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $file_list[$i] $file_list[$i+1]`;
    `/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/TrimGalore-0.6.10/trim_galore --paired -j 4 -o $p_directory/trimmed/ --fastqc_arg "-4 --outdir $p_directory/trimmed/Fastqc/" $file_list[$i] $file_list[$i+1]`;
}

###Move trimmed files to /trimmed/
#`mv $p_directory/rawdata/*trimmed* $p_directory/trimmed/`;

###Alignment
`mkdir $p_directory/aligned`;
#Declare fastq file directory list
@file_list = ();
push(@file_list, `find $p_directory/trimmed -name '*.fq.gz'`);
@file_list = sort(@file_list); #sort in suffix
for (@file_list){
    chomp($_);
}
for (@file_list){
    print("$_\n");
}
#Run bowtie2 and samtools sorting
for (my $i=0; $i<$#file_list+1; $i+=2){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/\.fq\./, $Oname_1[-1]);
    print("Making alignment file named $Oname_2[0]....\n");
    `bowtie2 -p 6 --very-sensitive -x /media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/KHM/align_data/Bowtie2indx/hg38/GRCh38_noalt_as -X 1000 -1 $file_list[$i] -2 $file_list[$i+1] 2>>$p_directory/aligned/bt2.$Oname_2[0].log | samtools view -@ 2 -S -bh -q 10 -f 2 -F 1804 - | samtools sort -@ 2 -m 24G -o $p_directory/aligned/$Oname_2[0].bam 2>>$p_directory/aligned/bt2.$Oname_2[0].log`; 
}

###Run picard_MarkDup -> Remove duplicate, chrM, blacklist. And run samtools indexing
#my $specific_name = "*.bam"; #specific_name for fastq file
`mkdir $p_directory/duprm`;
my $blacklist = "/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/KHM/align_data/GRCh38_unified_blacklist.bed.gz";
#Extract each fastq file directory from $fastq_dir and save to @file_list
@file_list = ();
push(@file_list, `find $p_directory/aligned -name '*.bam'`);
@file_list = sort(@file_list); #sort in suffix
#check if the @file_list
for (@file_list){
    chomp($_);
}
for (@file_list){
    print("$_\n");
}
#Run picard MarkDuplicate with removing chrM, blacklist, samtools indexing
for (my $i=0; $i<$#file_list+1; $i+=1){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/\.fastq-/, $Oname_1[-1]);
    print("Making alignment file named $Oname_2[0]....\n");
    `java -jar /media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/picard.jar MarkDuplicates I=$file_list[$i] O=$p_directory/duprm/$Oname_2[0]_duprm.bam M=$p_directory/duprm/$Oname_2[0]_marked_dup_metrics.txt`;
    `samtools view -h -@ 2 $p_directory/duprm/$Oname_2[0]_duprm.bam | egrep -v chrM | samtools view -bh -@ 2 - | bedtools intersect -abam stdin -b $blacklist -v > $p_directory/duprm/$Oname_2[0]_duprm_rmChrM_rmBlack.bam`;
    `samtools index -@ 2 $p_directory/duprm/$Oname_2[0]_duprm_rmChrM_rmBlack.bam`;
    `rm $p_directory/duprm/$Oname_2[0]_duprm.bam`;
}

