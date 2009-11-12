use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::Tools::IUPAC;
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::LiveSeq::Mutation;

use Getopt::Long;
my $up_data;
my $sgd_dir;
my $maq_snp_file;

my $result = GetOptions ("up=s" => \$up_data,    # numeric
    "sgd_dir=s"   => \$sgd_dir,      # string
    "maq=s"  => \$maq_snp_file);  # flag

## STEP1 Uniprot things ###


my %SGD2Uniprot;      # SGD <-> uniprot id translation hash table

#create Bio::SeqIO object for uniprot CHAGE THE FILENAME if you want to load a different uniprot file;
my $stream = Bio::SeqIO->new(-file => $up_data, -format=>"swiss");   

my %uniprot_place;    #uniprot_id <->[]object translation hash table

#look at every protein entry in the uniprot file. $seq is the object of a protein entry
while ( my $seq = $stream->next_seq() ) {

#look at every database_links in protein entry
  for my $dblink ($seq->{_annotation}->get_Annotations('dblink')){

#look at SGD annotation
    if($dblink->database eq 'SGD'){

      if (defined $SGD2Uniprot{$dblink->primary_id}){
#this part is for the error handling, I'm asking questions to uniprot 
      }

      $SGD2Uniprot{$dblink->primary_id} = $seq->primary_seq->display_id; # SGD<->uniprot id mapping table
        $uniprot_place{$seq->primary_seq->display_id} = $seq;              # uniprot id <-> uniprot object translation
    }
  }
}


## STEP1 done ##


## STEP2 MAQ snp file parsing and get the protein information to do with the snp ##

open my $RES1, ">Protein_details.txt" or die;
#open file of MAQ snp results CHAGE THE FILENAME if you want to load a different file
open my $FILE, $maq_snp_file or die;

my $old_chr = "";
my $seqio_object;
my $seq_object;
my %gene_list;

# looking at the snp file line by line
while(<$FILE>){
  chomp;

# MAQ details;
  my @line = split /\t/;            # information in a line is stored in @line
    my $chr = $line[0];               # chromosome number
    my $mutation_position = $line[1]; # mutation position in chromosome
    my $ref_nuc = $line[2];           # reference nuculeotide
    my $mut_nuc = $line[3];           # detected snp nucleotide (in IUB symbol)

    my $chr_num;                      # current chromosome num

# this part desides which GenBank file to look at
    if ($chr =~ /ref\|NC_0011(\d+)\|/){
      $chr_num = $1 - 32;
      $chr_num = "0".$chr_num if (length $chr_num == 1);
    }else{
      die;
    }

# for speedup, load the GenBank file only when the snp is in the defferent chr from the privious
  if ($old_chr ne $chr_num){
    $seqio_object = Bio::SeqIO->new(-file => "${sgd_dir}/chr${chr_num}.gbf");   
    $seq_object = $seqio_object->next_seq;
  }

  $old_chr = $chr_num;              # remember the privious chromosome


# look at each features in genbank $feat_object is the [] object
    for my $feat_object ($seq_object->get_SeqFeatures) {        

      next unless $feat_object->primary_tag eq "CDS";   # only look for CDS
        my $location = $feat_object->location;            # get the []Object


# check if the mutation is in the CDS
        if ($location->contains($mutation_position)){

          my $original_AA = $feat_object->seq->translate->seq;      # get the reference AA sequence of the CDS 

# this part is to change the IUB symbles( $mut_nuc ) to ATGC 
            my $ambiseq = Bio::Seq->new(-seq => $mut_nuc, -alphabet => 'dna');
          my $stream  = Bio::Tools::IUPAC->new(-seq => $ambiseq);

# for each nucleotide in IUB symble
          while (my $uniqueseq = $stream->next_seq()) {
# get the mutated nucletide sequence of the genome (not the CDS sequence!)
            Bio::SeqUtils->mutate($feat_object->entire_seq,
                Bio::LiveSeq::Mutation->new(-seq => $uniqueseq->seq,
                  -pos => $mutation_position
                  ));
# check if the mutated CDS AA is different from the reference AA
            if ($original_AA ne $feat_object->seq->translate->seq){

# if so, get the uniprot id from the SGD id (using the hash from STEP1)
              my $uniprot_id;
              my $SGD_id;
              for my $tag ($feat_object->get_all_tags) {        
                if ($tag eq  "db_xref"){        
                  for my $value ($feat_object->get_tag_values($tag)) {        
                    if ($value =~ /SGD:(.+)/){ 
                      $uniprot_id = $SGD2Uniprot{$1} or print "no up id for $1\n";
                      $SGD_id = $1;
                    } 
                  }
                }
              }    
              die "no SGD found" unless $SGD_id;                      #error check
                my $uniprot_feat = $uniprot_place{$uniprot_id} or die;         

# get the gene name, if there is no gene name, get the locus_tag 
              my $gene = ${$feat_object->{_gsf_tag_hash}->{gene}}[0] || ${$feat_object->{_gsf_tag_hash}->{locus_tag}}[0];
              $gene = "" unless $gene;
# error check
              if (defined $gene_list{$SGD_id}){
                print "\n\n",$gene_list{$SGD_id},"\t",$gene,"\n\n" unless $gene_list{$SGD_id} eq $gene; 
              }

# store the protein name for STEP 3
              $gene_list{$SGD_id} = $gene;

              my $j;  # positions of AA
                my @letters = split (//,$feat_object->seq->translate->seq);

# find where AA chage occured by itterating the positions in AA ($j)
              for my $MA (@letters){
                $j++;   
                my $OA = substr ($original_AA,$j-1,1);   # $OA is the original AA $MA is the mutated AA 

                  if( $OA ne $MA){
                    my $upseq = $uniprot_feat->seq."*";
                    my @ft_description = ();

# check if AA in SGD is same as uniprot
                    if ($original_AA eq $upseq){

# look for FT tables (which contains things like domains) in uniprot
                      for my $ft (@{$uniprot_feat->{_as_feat}}){
                        next if $ft->{_primary_tag} eq "CHAIN";                 # CHAIN is not a fun thing to watch at 
                        next unless ($ft->start =~ /^\d+$/ and $ft->end =~ /^\d+$/); 
                          if ($ft->contains($j)){                                 # check if the mutated AA position is in the FT
                            my $des = ${$ft->{_gsf_tag_hash}->{description}}[0] || "";  # get the annotation of FT
# die if scalar @{$ft->{_gsf_tag_hash}->{description}};
                              push @ft_description, $ft->{_primary_tag}.":".$des;   # push the information in to array for the print
                          }
                      }
                    }else{
                      push @ft_description, "couldn't see details, becuase AA in SGD & UniProt is different";            # if the sequence is different between SGD and UniProt
                    }
                    
#print the information you have found so far!
                    print $RES1 $chr_num,"\t", $mutation_position,"\t",$gene,"\t",$SGD_id,"\t",$ref_nuc,"\t",$uniqueseq->seq,"\t",$uniprot_id,"\t",$j,"\t",$OA,"\t",$MA,"\t";
                    print $RES1 join ", ", @ft_description;
                    print $RES1 "\n";
                  } 
              }
            }
#fix the mutated genome to the reference genome
            Bio::SeqUtils->mutate($feat_object->entire_seq,
                Bio::LiveSeq::Mutation->new(-seq => $ref_nuc,
                  -pos => $mutation_position
                  ));
          }
        }
    }
}
## STEP 2 end ##

open my $RES2, ">Protein_List.txt" or die;
## STEP 3 print the information of the prottein ##
for my $sgdid (keys %gene_list){
  my $up_id = $SGD2Uniprot{$sgdid} or die;
  my $up_entry = $uniprot_place{$up_id} or die;

  my @kegg_id = grep {$_->database eq 'KEGG'} @{$up_entry->{_annotation}->{_annotation}->{dblink}};
  my @go_id = grep {$_->database eq 'GO'} @{$up_entry->{_annotation}->{_annotation}->{dblink}};

  print $RES2 $gene_list{$sgdid},"\t";
  print $RES2 $sgdid,"\t";
  print $RES2 $up_id,"\t";
  print $RES2 join "," , map{$_->primary_id}@kegg_id;
  print $RES2 "\t";
  print $RES2 join "," , map{$_->primary_id}@go_id;
  print $RES2 "\t";
  print $RES2 $up_entry->primary_seq->desc;
  print $RES2 "\n";
}
## STEP 3 and script end ##

__END__

end of script

coela

