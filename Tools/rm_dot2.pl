#!/usr/bin/perl
use strict;
use warnings;

my $file = "/Users/hmkim/data/test.csv";
my $rm_dot="/Users/hmkim/data/rm_dot.csv";
open(O, ">$rm_dot");
open(F, "$file");
my $i = 0;
while(<F>){
	chomp;	
	my @line = split(/"/, $_);
	my @line = split(/\./, $line[2]);
	@line = split(/"/, $line[0]);
	# print("@line\n");
	# print("@line\n");
	#print O "$line[0]\n";
	print O "$line[1]\n";
	# print("$_\n");
	
	# $i++;
	# if($i==10){
	# 	last;
	# }
}
close(F);
close(O);
