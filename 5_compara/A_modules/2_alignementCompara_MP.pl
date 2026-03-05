#!/usr/bin/env perl
use strict;
use warnings;
use Bio::AlignIO;
use Bio::EnsEMBL::Registry;
use Time::HiRes qw( time );


## Mandatory variable
my $query_species = $ARGV[0];
my $filename = $ARGV[1];

## State of art
my $start_time = time();


print("--------------------------------------------\n");
print("            Module/Adaptor loading          \n");
## Connexion to ensembl
print("• Ensembl connexion ... ");
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
  );
print("DONE \n");
## Get the Compara Adaptor for gene in the query species
print("• Gene adaptor ... ");
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $query_species, "core", "gene" );
print("DONE \n");
## Get the Compara Adaptor for slice in the query species
print("• Slice adaptor ... ");
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($query_species, "core", "Slice");
print("DONE \n");
## Get the Compara Adaptor for MethodLinkSpeciesSet
print("• MethodLinkSpeciesSet adaptor ... ");
my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
print("DONE \n");
## Use of the Mercator Pecan Algortihm for the amniotes group
print("• Mercator Pecan module ... ");
my $methodLinkSpeciesSet = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name("PECAN", "amniotes");
print("DONE \n");
# Get the Compara Adaptor for GenomicAlignBlocks
print("• GenomicAlignBlock adaptor ... ");
my $genomic_align_block_adaptor =   Bio::EnsEMBL::Registry->get_adaptor("Multi", "compara", "GenomicAlignBlock");
print("DONE \n");
print("--------------------------------------------\n\n");

my $total = `wc -l < ${filename}`;
chomp($total);
print("Alignement of $total lncRNAs from $query_species ... : \n");

## Opening of the file containing all the lncRNA id that we want to map
open(FH, '<', $filename) or die $!;

## Progress indicator
my $count = 1;

## For each of this lncRNA
while(<FH>){
    ## Extraction of hte lncRNA id
    my $id_LNC = $_;
    chomp $id_LNC;

    ## Naming of the file and put it in the directory asscoiated 
    my $destinationFile = "${query_species}_alignement_MP_63amniotes/${id_LNC}_${query_species}_alignementAmniotes.txt";
    open(DES,'>',$destinationFile) or die $!; # Opening in writing of the file that will contain alignement information

    ## Extraction of gene information
    my $gene = $gene_adaptor->fetch_by_stable_id($id_LNC);
    # Define the start and end positions for the alignment + extraction of information :
    my $seqname = $gene->seq_region_name();
    my $start   = $gene->seq_region_start();
    my $end     = $gene->seq_region_end();
    my $strand = $gene->seq_region_strand();
    #print "$query_species\t$id_LNC\t$seqname:$start-$end:$strand\n";
    print DES "$query_species\t$id_LNC\t$seqname:$start-$end:$strand\n";

    #The fetch_all_by_MethodLinkSpeciesSet_Slice() returns a ref to an array of GenomicAlingBlock objects
    # Creation of a slice for the alignement (can't take the gene directly - caution the strans is always "+")
    my $specie_slice = $slice_adaptor->fetch_by_gene_stable_id($id_LNC,0);
    print "[$count/$total] - $id_LNC or $seqname:$start-$end:$strand\n";

    my $all_genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($methodLinkSpeciesSet, $specie_slice);

    my $all_aligns;
    foreach my $this_genomic_align_block(@$all_genomic_align_blocks) {
        my $skip_empty_GenomicAligns = 1;
        my $restricted_gab = $this_genomic_align_block->restrict_between_reference_positions($start, $end);
        # fetch the GenomicAligns and move through
        my @genomic_aligns = @ { $restricted_gab->get_all_GenomicAligns};
        foreach my $genomic_align (@genomic_aligns) {
            # Extraction of the scientific_name (short name no longer supported)
            my $species = $genomic_align->genome_db->get_scientific_name;
            # Assembly information to be sure of the version
            my $genebuild = $genomic_align->genome_db->assembly;
            # Extraction of the slice corresponding to the alignement
            my $slice = $genomic_align->get_Slice;
            # If there is an alignement
            if (defined($slice)) {
                my $seqRegionName = $slice->seq_region_name();
                my $start = $slice->start();
                my $end = $slice->end();
                my $strand = $slice->strand();
                print DES $species, "\t", $genebuild, "\t", $seqRegionName, ":", $start, "-", $end, ":", $strand,  "\n";
            } else { # If no alignemnt is found
                print DES $species, "\t", $genebuild, "\t","undefined", "\n";
            }
        }
        print DES "\n";
    }
    close(DES);
    $count ++;
}

close(FH);
print("\n");

my $end_time = time();
printf("Execution Time: %0.02f s\n", $end_time - $start_time);