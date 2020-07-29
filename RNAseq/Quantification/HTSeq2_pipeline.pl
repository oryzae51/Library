#!/usr/bin/perl
use strict;
use warnings;

print("Write SAM file directory\n");
my $directory = <STDIN>;
chomp($directory);

print("Write gtf directory\n");
my $gtf_dir = <STDIN>;
chomp($gtf_dir);

print("Write outpu directory\n")
my $output_dir = <STDIN>;
chomp($output_dir);

print("Choose number of strandedness\n
	0 : non\n
	1 : stranded\n
	2 : reverse stranded\n
	");
my $strandedness = <STDIN>;
chomp($strandedness);
if ($strandedness == 0){
	$strandedness = "no";
}elsif($strandedness == 1){
	$strandedness = "yes";
}elsif($strandedness == 2){
	$strandedness = "reverse";
}

my $file = `find $directory -name "SRR*"`;
my @file_list = split(/\n/, $file);
@file_list = sort(@file_list);

my $i = 0;
for (@file_list){
    my @file_name = split(/\//, $file_list[$i]);
    print("$file_name[-1]\n");
    print("printing command:\nhtseq-count -s $strandedness $file_list[$i] $gtf_dir > $output_dir/$file_name[-1]\n\n");
    `htseq-count -s $strandedness $file_list[$i] $gtf_dir > $output_dir/$file_name[-1]`;
    #`touch $output_dir/$file_name[-1]`;
    $i++;
}
