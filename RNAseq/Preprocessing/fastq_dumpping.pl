#!/usr/bin/perl
use warnings;
use strict;

#Write SRR directory
print("Write SRR file dir:\n");
my $srr_dir = <STDIN>;
chomp($srr_dir);

#find srr files
my @srr_list = ();
push(@srr_list, `find $srr_dir -name "SRR*" -print`);
@srr_list = sort(@srr_list);

#output directory
#Write output pathway you want
print("Write output pathway you want:\n");
my $out_dir = <STDIN>;
chomp($out_dir);


#dummping with 
foreach (@srr_list){
    my @Oname = split(/\//, $_);
    my @Oname_1 = split(/\n/, $Oname[-1]);
    print("Making counting file named $Oname_1[-1]....\n");
    print("$_\n");
    print("printing command:\n/media/bm/ETL4TiB/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump --gzip $out_dir/$Oname_1[-1].txt\n\n");
    #`/media/bm/ETL4TiB/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump --gzip $out_dir/$Oname_1[-1].txt`;
}