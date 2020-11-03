#usr/bin/perl
use strict;
use warnings;

my @union = @intersection = @difference = ();
my %count = ();

my @array1 = (1, 2, 3, 5, 7, 23, 8, 14, 95, 19);
my @array2 = (3, 14, 6, 22, 88, 19, 100);

foreach my $element (@array1, @array2) {
  $count{$element}++
 }

foreach my $element (keys %count) {
   push @union, $element;
   push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
}

foreach my $k ( keys %count ) {
  if ( $count{$k} > 1 ) {
    print "$k exist on both the arrays\n";
  }
}