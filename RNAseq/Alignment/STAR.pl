#usr/bin/perl
use strict;
use warnings;

#Write fastq file dir and index dir
print("Write fastq file dir :\n");
my $directory = <STDIN>;
chomp($directory); #raw fastq file directory
print("Write index file dir: \n");
my $indx = <STDIN>;
chomp($indx); #Index file directory and basename
print("Write the basename of files :\n");
my $specific_name = <STDIN>; #specific_name for fastq file
chomp($specific_name);
print("Write the file directory for SAM output :\n");
my $hit_file = <STDIN>;
chomp($hit_file); #file directory for SAM output

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

	print("printing command:\n\n/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir $indx --readFilesCommand zcat --readFilesIn $file_list[$i] $file_list[$i+1] --outFileNamePrefix $hit_file/$Oname_2[0]/\n");

    `/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir $indx --readFilesCommand zcat --readFilesIn $file_list[$i] $file_list[$i+1] --outFileNamePrefix $hit_file/$Oname_2[0]/`;
}

