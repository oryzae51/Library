#/usr/bin/perl
use strict;
use warnings;

my $directory = shift;
my $idx_dir = shift;
#my $mode = shift;
my $file = `find $directory -name "*.fastq"`;
my @file_list = split(/\n/, $file);
@file_list = sort(@file_list);
for (@file_list){
    print("$_\n");
}

my $i = 0;
for (@file_list){
    my @Oname = split(/\//, $file_list[$i]);
    @Oname = split(/\.f/, $Oname[-1]);
    print("Making aligning file named $Oname[0]....\n");
    `hisat2 -x $idx_dir -1 $file_list[$i] -2 $file_list[$i+2] -S ~/sequence_data/hisat2_aligned/$Oname[0]`;
    if ($i%2==0){
        $i+=1;
    }
    elsif ($i%2!=0){
        $i+=3;
    }
    if ($i >=16){last;}
}
