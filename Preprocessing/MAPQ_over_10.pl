#!/usr/bin/perl
use warnings;
use strict;

print("Sam file directory : ");
my $directory = <STDIN>; #SAM file directory
#chomp($directory);
print("Specific name for SAM file : ");
my $specific_name = <STDIN>; #Specific name for SAM file
#chomp($specific_name);
print("Output file directory with last / : ");
my $outdir = <STDIN>;
chomp($outdir);
#Extract each SAM file directory from $directory and save to @file_list

#Declare SAM file directory list
my @file_list = ();
push(@file_list, `find $directory -name '$specific_name'`);
@file_list = sort(@file_list);
for (@file_list){
    chomp($_);
}
shift(@file_list);
#Extract file name from @file_list directory
for (@file_list){
    my @Oname = split(/\//, $_);
    print("\n\n--------\nCutting SAM file named $Oname[-1]...\n--------\n");
    print("samtools view -q 10 -h -b $outdir$Oname[-1] $_");
    `samtools view -q 10 -h -b $outdir$Oname[-1] $_`;

}

# #Cut of with samtools
# my $i = 0;
# foreach (@file_list) {
#     print("Cutting SAM file $_ ...\n");
#     `samtools view -q 10 -h -o $directory$filename_list[$i] $_`;
# }
