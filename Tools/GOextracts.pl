#!/usr/bin/perl
use strict;
use warnings;

my $GOlist = "/Users/hmkim/data/GOlist.csv";
my $GOextracts="/Users/hmkim/data/GOextracts.csv";
open(O, ">$GOextracts");
open(F, "$GOlist");
my $i = 0;
while(<F>){
	chomp;	
	my @line = split(/GO:/, $_);
	@line = split(/\)/, $line[1]);
	# print("@line\n");
	# print("@line\n");
	#print O "$line[0]\n";
	print O "GO:$line[0]\n";
	# print("$_\n");
	
	# $i++;
	# if($i==10){
	# 	last;
	# }
}
close(F);
close(O);
