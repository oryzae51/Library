#!/usr/bin/perl
use warnings;
use strict;

#print("Name of RAW data dir : ");
my $fastq_dir = shift; #RAW fastq file directory
my $output_dir = shift; #output file directory for FastQC
my $specific_name = "*fastq.gz"; #specific_name for fastq file
#Extract each fastq file drictory from $fastq_dir and save to @file_list

#Declare fastq file directory list
my @file_list = ();
push(@file_list, `find $fastq_dir -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix

#Iterate each file in @file_list and run fastq

my $i = 0;
foreach (@file_list) {
	print("Quality checking of fastq file \n");
	`~/data/FastQC/fastqc -o '$output_dir' --noextract $_`;
}