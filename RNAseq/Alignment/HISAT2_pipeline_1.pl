#/usr/bin/perl
use strict;
use warnings;

#Write fastq file dir and index dir
print("Write the directory of raw fastq file: \n");
my $directory = <STDIN>;
chomp($directory); #raw fastq file directory
print("Write the directory of index file and basename: \n");
my $hisat2_indx = <STDIN>;
chomp($hisat2_indx); #Index file directory and basename
print("Write the specific name for fastq file: \n");
my $specific_name = "*fastq.gz"; #specific_name for fastq file
print("Write the directory for output files: \n");
my $hit_file = <STDIN>;	
chomp($hit_file); #file directory of SAM

#Extract each fastq file drictory from $fastq_dir and save to @file_list
#my $file = `find $directory -name "*.fastq"`;
#my @file_list = split(/\n/, $file);
#Declare fastq file directory list
my @file_list = ();
push(@file_list, `find $directory -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix
#shift(@file_list);
#check if the @file_list is made correctly
for (@file_list){
    chomp($_);
}
for (@file_list){
    print($_);
}

for (my $i=0; $i<$#file_list+1; $i+=2){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/_/, $Oname_1[-1]);
    print("Making alignment file named $Oname_2[0]....\n");
    print("$file_list[$i]\n");
    print("$file_list[$i+1]\n");
    print("hisat2 -x $hisat2_indx -p 8 -1 $file_list[$i] -2 $file_list[$i+1] -S $hit_file/$Oname_2[0]\n");
    `hisat2 -x $hisat2_indx -p 8 -1 $file_list[$i] -2 $file_list[$i+1] -S $hit_file/$Oname_2[0]`;
}
