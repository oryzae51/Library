#!/usr/bin/perl
use strict;
use warnings;

my $directory = shift;
my $gtf_dir = shift;
my $output_dir = shift;
my $iattr = shift;

my $file = `find $directory -name "SRR*"`;
my @file_list = split(/\n/, $file);
@file_list = sort(@file_list);

my $i = 0;
for (@file_list){
    my @file_name = split(/\//, $file_list[$i]);
    print("$file_name[-1]\n");
    `htseq-count -s reverse -i $iattr $file_list[$i] $gtf_dir > $output_dir/$file_name[-1]`;
    #`touch $output_dir/$file_name[-1]`;
    $i++;
}
