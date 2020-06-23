#/usr/bin/perl
use strict;
use warnings;

#Write annotation gtf file dir
print("Write annotation gtf file dir:\n");
my $gtf_dir = <STDIN>;
chomp($gtf_dir);

#Write input align_data dir
print("Write input align_data dir\n");
my $align_data_dir = <STDIN>;
chomp($align_data_dir);

#Write output pathway you want
print("Write output pathway you want:\n");
my $out_dir = <STDIN>;
chomp($out_dir);

# #Write input BAM/SAM dir
# print("Write input BAM/SAM dir\n");
# my $sam_dir = <STDIN>;
# chomp($sam_dir);

#add sam folder dir to list
my @align_data_list = ();
push(@align_data_list, `find $align_data_dir`);
@align_data_list = sort(@align_data_list);
shift(@align_data_list)
for (@align_data_list){
    chomp($_);
}
for (@align_data_list){
    print("$_\n");
}

#iterate @align_data_list
foreach(@align_data_list){
	my @sam_list = ();
	push(@sam_list, `find "$(pwd)`);#find absolute path of file
	shift(@sam_list);
	@sam_list = sort(@sam_list);
	for (@sam_list){
	    chomp($_);
	}
	for (@sam_list){
		print("$_\n");
	}

	foreach (@sam_list){
	    my @Oname_1 = split(/\//, $_);
	    print("Making counting file named $Oname_1[-1]....\n");
	    print("$_\n");
	    print("printing command:\nfeatureCounts -T 4 -p -a $gtf_dir -t exon -g gene id -o $out_dir/$Oname_1[-1].txt $_\n\n");
	    #`featureCounts -T 4 -p -a $gtf_dir -t exon -g gene id -o $_.txt $sam_dir`;
	}
}
