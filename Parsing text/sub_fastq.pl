#!/usr/bin/perl
use warnings;
use strict;

my $sample = shift;
open(F, "gunzip -c $sample |") or die "gunzip $sample: $!";
my $subset_file1 = shift;
my $subset_file2 = shift;
#`touch ~/code/$subset_file1.txt`;
#`touch ~/code/$subset_file2.txt`;
open(O, ">$subset_file1");
open(P, ">$subset_file2");
my $reads = shift;
my $i = 0;
while(<F>){
    my $line = $_;
    if ($i<$reads*4){
        print O "$_";
    }
    elsif ($i>=$reads*4 && $i<$reads*4*2){
        print P "$_";
    }
    $i++;
    if ($i == $reads*4*2){
        last;
    }
}


close(F);
close(O);
close(P);
