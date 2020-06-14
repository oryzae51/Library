#!/usr/bin/perl
use strict;
use warnings;

my $samfile = shift;
my $FPKMfile = shift;
open(F, "$samfile");
open(O, ">$FPKMfile");

while(<F>){
    my $line = $_;
    my @line_split = split(/\s+|\t/, $line);
    my $FPKM = $line_split[-9];
    print O $FPKM;
    print O "\n";
#    print("$FPKM\n");
}

close(O);
close(F);
