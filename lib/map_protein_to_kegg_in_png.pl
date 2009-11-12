use strict;
use Data::Dumper;
use SOAP::Lite;
my $wsdl = 'http://soap.genome.jp/KEGG.wsdl';
my $serv = SOAP::Lite -> service($wsdl);

my $var;

open my $LIST, "../data/KEGG/sce_sgd-sce.list" or die;

my %sgd2kegg;
while(<$LIST>){
  chomp;
  my ($gene, $sgd) = split /\t/;
  if ($sgd =~ /sgd-sce\:(\S+)/){
    $sgd2kegg{$1} = $gene;
  }
}

open my $HOGE, $ARGV[0] or die "no input file\n";
my %check;
while(<$HOGE>){
  chomp; 
  my @line = split /\t/;
  my $sgd =  $line[1];
  my $kegg_gene = $sgd2kegg{$sgd};
  next if $check{$kegg_gene};
  $check{$kegg_gene} = 1;
  for my $path (@{$serv->get_pathways_by_genes([$kegg_gene])}){
    push @{$var->{$path}},$kegg_gene;
  }
}

for my $path (keys %$var){
  print $serv->mark_pathway_by_objects($path, $var->{$path});
  print "\n";
}


# sub routines implicitly used in the above code
sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}
