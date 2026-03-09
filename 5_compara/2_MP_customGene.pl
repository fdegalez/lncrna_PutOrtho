#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Getopt::Long qw(GetOptions);

############################################################
# Arguments
############################################################

my $speciesName_source = 'NA';
my $input  = 'NA';
my $output = 'NA';

GetOptions(
    'name=s'   => \$speciesName_source,
    'input=s'  => \$input,
    'output=s' => \$output,
) or die "Usage: $0 --name species --input file.tsv --output file.tsv\n";

die "Missing --name\n"  if $speciesName_source eq 'NA';
die "Missing --input\n" if $input eq 'NA';
die "Missing --output\n" if $output eq 'NA';

############################################################
# Registry loading
############################################################

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -verbose => 0
);

############################################################
# Adaptors
############################################################

my $slice_adaptor =
    $registry->get_adaptor($speciesName_source, "core", "Slice");

my $method_link_species_set_adaptor =
    $registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");

my $methodLinkSpeciesSet =
    $method_link_species_set_adaptor
    ->fetch_by_method_link_type_species_set_name("PECAN", "amniotes");

my $genomic_align_block_adaptor =
    $registry->get_adaptor("Multi", "compara", "GenomicAlignBlock");

############################################################
# Input / Output
############################################################

open(my $IN,  "<", $input)  or die "Cannot open $input\n";
open(my $OUT, ">", $output) or die "Cannot write $output\n";

my $count = 0;
my $start_run = time();

############################################################
# Processing
############################################################

while (my $line = <$IN>) {

    chomp $line;

    $count++;

    my ($seq_region, $start, $end, $strand, $seq_id) = split("\t", $line);

    print "[$count] $seq_id\n";

    my $species_slice =
        $slice_adaptor->fetch_by_region(
            'toplevel',
            $seq_region,
            $start,
            $end,
            $strand
        );

    next unless defined $species_slice;

    my @gab_list =
        @{ $genomic_align_block_adaptor
            ->fetch_all_by_MethodLinkSpeciesSet_Slice(
                $methodLinkSpeciesSet,
                $species_slice
            )
        };

    foreach my $gab (@gab_list) {

        my $restricted_gab =
            $gab->restrict_between_reference_positions($start, $end);

        next unless defined $restricted_gab;

        my @genomic_aligns =
            @{ $restricted_gab->get_all_GenomicAligns };

        foreach my $genomic_align (@genomic_aligns) {

            my $genome_db = $genomic_align->genome_db;

            my $species_scientificName =
                $genome_db->get_scientific_name();

            my $species_shortName =
                $genome_db->get_short_name();

            my $species_displayName =
                $genome_db->display_name();

            my $genebuild = $genome_db->genebuild();
            my $assembly  = $genome_db->assembly();

            my $slice = $genomic_align->get_Slice;

            next unless defined $slice;

            my $seqRegionName = $slice->seq_region_name;
            my $slice_start   = $slice->start;
            my $slice_end     = $slice->end;
            my $slice_strand  = $slice->strand;

            my @genes = @{ $slice->get_all_Genes };

            my $genes_inSlice = "";

            foreach my $gene (@genes) {
                $genes_inSlice .= $gene->stable_id . ";";
            }

            chop($genes_inSlice) if $genes_inSlice ne "";

            print $OUT join(
                "\t",
                $speciesName_source,
                $seq_id,
                $start,
                $end,
                $strand,
                $species_scientificName,
                $species_shortName,
                $species_displayName,
                $genebuild,
                $assembly,
                $seqRegionName,
                $slice_start,
                $slice_end,
                $slice_strand,
                $genes_inSlice
            ), "\n";
        }
    }
}

close($IN);
close($OUT);

############################################################
# Runtime
############################################################

my $end_run = time();
my $run_time = $end_run - $start_run;

print "Job took $run_time seconds\n";