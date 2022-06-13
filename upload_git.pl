#usr/bin/perl
use strict;
use warnings;

print("Commit : ");
my $comm = <STDIN>; chomp($comm);

`git init`;
`git add . `;
`git commit -m $comm`;
`git push -u origin master`;
