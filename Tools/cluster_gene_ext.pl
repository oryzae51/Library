#/usr/bin/perl
use strict;
use warnings;

my $KvsH = "/Users/hmkim/data/M124W_clustering/KvsH_K4.csv";
my $gene_ext_K4="/Users/hmkim/data/M124W_clustering/gene_ext_H4K1.csv";
open(O, ">$gene_ext_K4");
open(F, "$KvsH");
my $i = 0;
while(<F>){
	chomp;	
	my @line = split(/,/, $_);
	# print("@line\n");
	if($line[1]==4 & $line[2]==1){
		print O "$_\n";
		# print("$_\n");
	}
	# $i++;
	# if($i==10){
	# 	last;
	# }
}
close(F);
close(O);