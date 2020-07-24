#!/usr/bin/perl
use warnings;
use strict;

#print("Name of RAW data dir : ");
my $fastq_dir = <STDIN>; #RAW fastq file directory
chomp($fastq_dir);
my $output_dir = <STDIN>; #output file directory for FastQC
chomp($output_dir);
my $specific_name = <STDIN>; #specific_name for fastq file
chomp($specific_name);
#Extract each fastq file drictory from $fastq_dir and save to @file_list

#Declare fastq file directory list
my @file_list = ();
push(@file_list, `find $fastq_dir -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix

#Iterate each file in @file_list and run fastq

my $i = 0;
foreach (@file_list) {
	print("Quality checking of fastq file \n");
	print("/media/bm/ETL4TiB/Tools/FastQC/fastqc -o $output_dir --noextract $_");
	#`/media/bm/ETL4TiB/Tools/FastQC/fastqc -o '$output_dir' --noextract $_`;
}