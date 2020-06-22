#/usr/bin/perl
use strict;
use warnings;

#Write annotation gtf file dir
print("Write annotation gtf file dir:\n")
my $gtf_dir = <STDIN>;
chomp($gtf_dir);
my @gtf_list = ();
push(@gtf_list, `find $gtf_dir`);
shift(@gtf_list);
@gtf_list = sort(@gtf_list);
for (@gtf_list){
    chomp($_);
}
for (@gtf_list){
    print("$_\n");
}

#Write input align_data dir
print("Write input align_data dir\n");
my $align_data_dir = <STDIN>;
chomp($align_data_dir);

# #Write input BAM/SAM dir
# print("Write input BAM/SAM dir\n");
# my $sam_dir = <STDIN>;
# chomp($sam_dir);

#add sam folder dir to list
my @align_data_list = ();
push(@align_data_list, `find $align_data_dir`);
shift(@align_data_list);
@align_data_list = sort(@align_data_list);
for (@align_data_list){
    chomp($_);
}
for (@align_data_list){
    print("$_\n");
}

#iterate @align_data_list
foreach(@align_data_list){
	my @sam_list = ();
	push(@sam_list, `find $_`);
	shift(@sam_list);
	@sam_list = sort(@sam_list);
	for (@sam_list){
	    chomp($_);
	}
	for (@sam_list){
		print("$_\n");
	}

	foreach (@sam_list){
	    my @Oname_1 = split(/\//, $sam_list[$i]);
	    print("Making counting file named $Oname_1[-1]....\n");
	    print("$sam_list[$i]\n");
	    print("printing command:\nfeatureCounts -T 4 -p -a $gtf_dir -t exon -g gene id -o $Oname_1[-1].txt $_\n\n");
	    #`featureCounts -T 4 -p -a $gtf_dir -t exon -g gene id -o $_.txt $sam_dir`;
	}
}
