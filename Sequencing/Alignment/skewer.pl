#usr/bin/perl
use strict;
use warnings;

#Project directory
print("Write the project directory :\n");
my $p_directory = <STDIN>;
chomp($p_directory);
my $specific_name = "*.fastq.gz"; #specific_name for fastq file
`mkdir $p_directory/trimmed`;

#Extract each fastq file drictory from $fastq_dir and save to @file_list
my @file_list = ();
push(@file_list, `find $p_directory/rawdata -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix
#check if the @file_list is made correctly
for (@file_list){
    chomp($_);
}
for (@file_list){
    print("$_\n");
}

for (my $i=0; $i<$#file_list+1; $i+=2){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/\.fastq\./, $Oname_1[-1]);
    print("Making alignment file named $Oname_2[0]....\n");
    `skewer -t 4 -z -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $file_list[$i] $file_list[$i+1]`;
}
