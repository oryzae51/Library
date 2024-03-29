#usr/bin/perl
use strict;
use warnings;

#Project directory and declare files
print("Write the project directory :\n");
my $p_directory = <STDIN>;
chomp($p_directory);
my $specific_name = "*.bam"; #specific_name for fastq file
`mkdir $p_directory/duprm`;
my $blacklist = "/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/KHM/align_data/GRCh38_unified_blacklist.bed.gz";

#Extract each fastq file drictory from $fastq_dir and save to @file_list
my @file_list = ();
push(@file_list, `find $p_directory/aligned -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix
#check if the @file_list is made correctly
for (@file_list){
    chomp($_);
}
for (@file_list){
    print("$_\n");
}

for (my $i=0; $i<$#file_list+1; $i+=1){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/\.fastq-/, $Oname_1[-1]);
    print("Making alignment file named $Oname_2[0]....\n");
    `java -jar /media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/picard.jar MarkDuplicates I=$file_list[$i] O=$p_directory/duprm/$Oname_2[0]_duprm.bam M=$p_directory/duprm/$Oname_2[0]_marked_dup_metrics.txt`;
    `samtools view -h -@ 2 $p_directory/duprm/$Oname_2[0]_duprm.bam | egrep -v chrM | samtools view -bh -@ 2 - | bedtools intersect -abam stdin -b $blacklist -v > $p_directory/duprm/$Oname_2[0]_duprm_rmChrM_rmBlack.bam`;
    `samtools index -@ 2 $p_directory/duprm/$Oname_2[0]_duprm_rmChrM_rmBlack.bam`;
    #`rm $p_directory/duprm/$Oname_2[0]_duprm.bam`;
}
