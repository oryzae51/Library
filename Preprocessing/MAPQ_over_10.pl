#!/usr/bin/perl
use warnings;
use strict;

my $directory = shift; #SAM file directory
my $specific_name = shift; #Specific name for SAM file
#Extract each SAM file directory from $directory and save to @file_list

#Declare SAM file directory list
my @file_list = ();
push(@file_list, `find $directory -name '$specific_name'`);
@file_list = sort(@file_list);

#Extract file name from @file_list directory
my @filename_list = ();
foreach (@file_list) {
    push(split(/\//, $_)[-1], @filename_list);
};

#Cut of with samtools
my $i = 0;
foreach (@file_list) {
    print("Cutting SAM file $_ ...\n");
	`samtools view -q 10 -h -o $directory$filename_list[$i] $_`
}