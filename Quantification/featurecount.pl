#/usr/bin/perl
use strict;
use warnings;

#Write annotation gtf file dir
print("Write annotation gtf file dir:\n")
my $gtf_dir = <STDIN>;
chomp($gtf_dir);

#Write input BAM/SAM dir
print("Write input BAM/SAM dir\n");
my $sam_dir = <STDIN>;
chomp($sam_dir);

#add sam file dir to list
my @sam_list = ();
push(@sam_list, `find $sam_dir`);
shift(@sam_list);
@sam_list = sort(@sam_list);
for (@sam_list){
    chomp($_);
}

for (my $i=0; $i<$#file_list+1; $i+=2){
    my @Oname_1 = split(/\//, $file_list[$i]);
    my @Oname_2 = split(/_/, $Oname_1[-1]);
    print("Making alignment file named $Oname_2[0]....\n");
    print("$file_list[$i]\n");
    print("$file_list[$i+1]\n");
    print("printing command:\nfeatureCounts -T 4 -p -a $gtf_dir -t exon -g gene id -o $_.txt $sam_dir\n\n");
    #`featureCounts -T 4 -p -a $gtf_dir -t exon -g gene id -o $_.txt $sam_dir`;
}
