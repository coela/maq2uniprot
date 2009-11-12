use strict;

open my $FILE, $ARGV[0];

my $i=0;
while(<$FILE>){
  chomp;
  $i++;
  s/\s//g if $i%2;
  print $_,"\n";
}

