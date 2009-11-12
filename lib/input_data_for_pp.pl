use strict;
use Bio::SeqIO;
use Data::Dumper;

#KEGG map
open my $LIST, "../data/KEGG/sce_sgd-sce.list" or die;

my %sgd2kegg;
while(<$LIST>){
  chomp;
  my ($gene, $sgd) = split /\t/;
  if ($sgd =~ /sgd-sce\:(\S+)/){
    $sgd2kegg{$1} = $gene;
  }
}

open my $HOGE, $ARGV[0] or die "no input file";
#PP things
 my $stream = Bio::SeqIO->new(-file => '../data/KEGG/S.cerevisiae.ent', -format => 'KEGG');

print "Organism:sce\n";
my %hash;

while ( my $seq = $stream->next_seq() ) {
print Dumper $seq;
  $hash{ $seq->primary_seq->accession_number } = $seq;
}

print Dumper \%hash;
while(<$HOGE>){
  chomp; 
  my @line = split /\t/;
  my $sgd =  $line[1];
  my $kegg_gene = $sgd2kegg{$sgd};
  $kegg_gene = /sce\:(\S+)/;
  my $kegg_gene_simple = $1;
  my $seq = $hash{$kegg_gene_simple};

  if ($seq->primary_seq->display_id){  
    print $seq->primary_seq->display_id;
    print ",red\n";
  }
}

