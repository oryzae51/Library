#!/usr/bin/perl
use warnings;
use strict;

#Collect fastq file name from sample directory
#my $dir = shift;
#my $format = shift;
my @filename_list = ();
push(`find ~/sequence_data/raw_seq -name '*.fastq.gz'`, @filename_list); 

print(@filename_list);
#my @filename_list = split(/\n/, $filename);
my @line_list = ();
my @half = ();
my @quartile = ();
print("@filename_list\n");
#Count the number of line of fastq file
for (@filename_list) {
    print("Counting and pushing line\n");
    my $count = `wc -l < $_`;
    push(@line_list, $count);
}
print("@line_list");

#Declare the dividing count variable
for (@line_list){
    @half = sprintf("%.0f", $_/2);
    @quartile = sprintf("%.0f", $_/4);
}
print("@half\n");
print("@quartile\n");
