#usr/bin/perl
use strict;
use warnings;

#Project directory
print("Write the project directory :\n");
my $p_directory = <STDIN>;
chomp($p_directory);
#Write fastq file dir and index dir
#print("Write fastq file dir :\n");
#my $directory = ;
#chomp($directory); #raw fastq file directory
#print("Write the basename of files :\n");
my $specific_name = "*.fastq.gz"; #specific_name for fastq file
#chomp($specific_name);
#print("Write the file directory for output :\n");
#my $hit_file = <STDIN>;
#chomp($hit_file); #file directory for save output
`mkdir $p_directory/aligned`;

#Extract each fastq file drictory from $fastq_dir and save to @file_list
#my $file = `find $directory -name "*.fastq"`;
#my @file_list = split(/\n/, $file);
#Declare fastq file directory list
my @file_list = ();
push(@file_list, `find $p_directory/trimmed -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix
#shift(@file_list);
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
    print("$file_list[$i]\n");
    `bowtie2 -p 8 --very-sensitive -x /media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/KHM/align_data/Bowtie2indx/hg38 -X 2000 -1 $file_list[$i] -2 $file_list[$i+1] | samtools view -Sb - | samtools view -b -q 10 -f 2 - | samtools sort - $p_directory/aligned/$Oname_2[0]`;
}
