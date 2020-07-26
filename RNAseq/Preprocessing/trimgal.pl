#!/usr/bin/perl
use warnings;
use strict;

#Write fastq directory
print("Write fastq file dir:\n");
my $fastq_dir = <STDIN>;
chomp($fastq_dir);

#find fastq files
my @fastq_list = ();
push(@fastq_list, `find $fastq_dir -name "*.fastq" -print`);

@fastq_list = sort(@fastq_list);

#Write output pathway you want
print("Write output pathway you want:\n");
my $out_dir = <STDIN>;
chomp($out_dir);

#dummping with 
foreach (@fastq_list){
    my @Oname = split(/\n/, $_);
    my @Oname_1 = split(/\//, $Oname[-1]);
    print("Making counting file named $Oname_1[-1]....\n");
    print("$_\n");
    print("printing command:\n/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/TrimGalore-0.6.5/trim_galore -q 28 -j 8 -o $out_dir $_\n\n");
    `/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/TrimGalore-0.6.5/trim_galore -q 28 -j -o $out_dir $_`;
}