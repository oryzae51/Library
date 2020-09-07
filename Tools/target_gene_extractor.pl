#usr/bin/perl
use strict;
use warnings;

#This code is for extract target gene symbol from data
#Data is downloaded from 
#https://amp.pharm.mssm.edu/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets
#Target gene measurement : Transcription factory DNA-binding by Chip-seq

print("Write the dir of transcription factor target data frim ENCODE database:\n");
my $file = <STDIN>; chomp($file);
print("Transcription factor name for output file name:\n");
my $TFactor = <STDIN>; chomp($TFactor);
my @gene_list = ();
my $output_file = "/Users/hmkim/data/deg_data/target_gene_list/$TFactor.csv";


open(O, ">$output_file");
open(F, "$file");
my $i = 0;
while(<F>){
	chomp;
	my @splitOne = split(/{"gene":{"symbol":"/, $_);
	foreach(@splitOne){
		my @splitTwo = split(/"/, $_);
		push(@gene_list, $splitTwo[0]);
	}
}
shift(@gene_list);
foreach(@gene_list){
	print O "$_\n";
}
close(F);
close(O);
