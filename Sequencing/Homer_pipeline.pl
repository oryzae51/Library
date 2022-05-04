#usr/bin/perl
use strict;
use warnings;

#Project directory
print("Write the project directory :\n");
my $p_directory = <STDIN>;
chomp($p_directory);
my $specific_name = "*_duprm.bam"; #specific_name for dup removed bam file
`mkdir $p_directory/Tagdir_dir`;

#Extract each bam file drictory from bam file dir and save to @file_list
my @file_list = ();
push(@file_list, `find $p_directory/aligned -name '$specific_name'`);
@file_list = sort(@file_list); #sort in suffix
#shift(@file_list);
#check if the @file_list is made correctly
for (@file_list){
    chomp($_);
}
for (@file_list){
    print("$_\n");
}

#Extracting file name for output and run tools(bowtie2, Homer, etc)
for (my $i=0; $i<$#file_list+1; $i+=1){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/\.fastq-/, $Oname_1[-1]);
    print("Making Tagdir with Homer makeTagDirectory, file named $Oname_2[0]....\n");
    `mkdir $p_directory/Tagdir_dir/$Oname_2[0]`;
    #print("$file_list[$i]\n");
    #print("/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/HOMER/bin/makeTagDirectory $p_directory/Tagdir_dir/$Oname_2[0] $file_list[$i]\n");
    `/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/HOMER/bin/makeTagDirectory $p_directory/Tagdir_dir/$Oname_2[0] $file_list[$i]` ;
    #print("Making bed file with Homer findPeaks, options -style factor, histone\n");
    #`/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/HOMER/bin/findPeaks $p_directory/Tagdir_dir/$Oname_2[0] -style factor -o auto`;
    #`/media/bm/790240e4-2887-451f-ad02-1b19c4b4e120/Tools/HOMER/bin/findPeaks $p_directory/Tagdir_dir/$Oname_2[0] -style histone -o auto`; 
}
