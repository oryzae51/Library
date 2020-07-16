#/usr/bin/perl
use warnings;
use strict;

print("CDS sequence file dir : \n");
my $fasta_not_chomp = <STDIN>;
print("Front CDS sequence of domain : \n");
my $aa_region_f = <STDIN>;
print("Reverse CDS sequence of domain : \n");
my $aa_region_r = <STDIN>;

open(F, $fasta_not_chomp);
my $fasta = "";
while (<F>) {
    chomp($_);
    $fasta = $fasta . $_;
}
#print("$fasta\n");
my $base_region_f = ($aa_region_f - 1) * 3;
my $base_region_r = ($aa_region_r - 1) * 3;
my $i = 0;
foreach my $fasta_char (split //, $fasta) {
    if ($i >= $base_region_f && $i < $base_region_r){
        print($fasta_char);
    }
    $i++;
}
print("\n");

close(F);
