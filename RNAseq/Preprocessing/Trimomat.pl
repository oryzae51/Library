#!/usr/bin/perl
use warnings;
use strict;

print("fastq raw data dir\n");
my $fastq_dir = <STDIN>; #RAW fastq file directory
chomp($fastq_dir);
print("output dir\n");
my $output_dir = <STDIN>; #output file directory for FastQC
chomp($output_dir);
print("specific name for fastq\n");
my $specific_name = <STDIN>; #specific_name for fastq file
chomp($specific_name);
#Extract each fastq file drictory from $fastq_dir and save to @file_list

#Declare fastq file directory list
my @file_list = ();
push(@file_list, chomp(`find $fastq_dir -name '$specific_name'`));
@file_list = sort(@file_list); #sort in suffix

#Iterate each file in @file_list and run fastq

my $i = 0;
for (my $i=0; $i<$#file_list+1; $i+=2) {
    print("File1: $file_list[$i]\nFile2: $file_list[$i+1]\n");
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_1_str = split(/\./, $Oname_1[-1]);
    my @Oname_2 = split(/\//, $file_list[$i+1]);
    my @Oname_2_str = split(/\./, $Oname_2[-1]);

    print("$Oname_1_str[0]\n$Oname_2_str[0]\n");
    print("Trim-O-mating.........\nForward specific name : $Oname_1_str[0] \nReverse specific name : $Oname_2_str[0] \n");
    print("java -jar /media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 -phred33 
        $file_list[$i] 
        $file_list[$i+1] 
        $output_dir/$Oname_1_str[0]_f_paired.fastq.gz
        $output_dir/$Oname_1_str[0]_f_unpaired.fastq.gz 
        $output_dir/$Oname_2_str[0]_r_paired.fastq.gz 
        $output_dir/$Oname_2_str[0]_r_unpaired.fastq.gz 
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n");
    `java -jar /media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 -phred33 $file_list[$i] $file_list[$i+1] $output_dir/$Oname_1_str[0]_f_paired.fastq.gz$output_dir/$Oname_1_str[0]_f_unpaired.fastq.gz $output_dir/$Oname_2_str[0]_r_paired.fastq.gz $output_dir/$Oname_2_str[0]_r_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`;
}