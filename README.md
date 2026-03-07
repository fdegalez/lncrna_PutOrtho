# lncRNA Orthology Inference Pipeline

# Overview

This repository provides a modular workflow to infer **putative orthologous relationships between long non‑coding RNAs (lncRNAs)** across multiple species.

Because lncRNAs evolve rapidly and often lack strong sequence conservation, classical orthology inference methods designed for protein‑coding genes are frequently insufficient. This pipeline integrates **three complementary strategies** to detect lncRNA orthologs:

1.  **Synteny-based inference**
2.  **Genomic context conservation (FEELnc classification)**
3.  **Genome alignment conservation (Ensembl Compara, here: Mercator--Pecan alignments)**

The workflow provides a flexible framework allowing users to detect candidate lncRNA orthologs using complementary evidence.
For a more general overview, you could point to [the associated paper](https://www.biorxiv.org/content/10.1101/2024.10.03.616473v1).

------------------------------------------------------------------------

# Repository Structure

    .
    ├── data
    │   └── config.txt
    │
    ├── 1_extractionGenes
    │   └── extract_genes.sh
    │
    ├── 2_extractionOrthologyPCG
    │   └── run_OrthologyExtraction.sh
    │
    ├── 3_synteny
    │   └── run_synteny.sh
    │
    ├── 4_FEELnc
    │   └── run_orthoFEELnc.sh
    │
    ├── 5_compara
    │   └── run_compara_pipeline.sh
    │
    ├── 6_summary
    │   └── TODO
    │
    └── run_all_methods.sh

------------------------------------------------------------------------

# Pipeline Overview

The workflow is composed of **five main modules**:

    GTF annotations
    │
    ├── 1_extractionGenes
    │
    ├── 2_extractionOrthologyPCG
    │
    ├── 3_synteny
    │
    ├── 4_FEELnc
    │
    └── 5_compara

The two first modules are mandatory.
Each module then can be executed independently or through the global pipeline script (run_all_methods.sh).

------------------------------------------------------------------------

# Requirements / Software Dependencies

## R

Minimum version:
    [R](https://cran.r-project.org/) ≥ 4.0

Required R packages:
    - [BiomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) v2.38.0 or more - Interface to BioMart databases (i.e Ensembl). 
    - [stringr](https://cloud.r-project.org/web/packages/stringr/index.html) v1.5.0 or more - Simple, Consistent Wrappers for Common String Operations. 
    - [stringi](https://cran.r-project.org/web/packages/stringi/index.html) v1.8 or more  - Character String Processing Facilities.  
    - [pbapply](https://cran.r-project.org/web/packages/pbapply/index.html) v1.7 or more  - Adding Progress Bar to '*apply' Functions.  
    - [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html) v1.3.1 or more  - Tidy Messy Data.
    - [dplyr](https://cran.r-project.org/web/packages/stringi/index.html) v1.1.4 or more  - A Grammar of Data Manipulation.


------------------------------------------------------------------------

## Perl

Required for the Ensembl API.
    [Perl5+](https://www.perl.org/) ≥ 5.32

Required modules:
    - [Ensembl API](https://www.ensembl.org/info/docs/api/api_installation.html) : tested with e! v104  
    - [Bioperl](http://www.bioperl.org/wiki/Main_Page) : tested with version 1.7.8  
    - [TimeHiRes](https://metacpan.org/pod/Time::HiRes) : tested with version 1.9764      


------------------------------------------------------------------------

## FEELnc

Used for genomic context classification of lncRNAs.

Repository:

- [FEELnc](https://github.com/tderrien/FEELnc) : test with version 0.2.1 - 2022-07-20
    - [FEELnc_classifier.pl](https://github.com/tderrien/FEELnc#3--feelnc_classifierpl) : Classify lncRNAs based on their genomic localization with others transcripts.
    - [FEELnc_tpLevel2gnLevelClassifcation.R](https://github.com/tderrien/FEELnc/blob/master/scripts/FEELnc_tpLevel2gnLevelClassification.R) : Transformation of transcript-level configurations to gene-level models. 


------------------------------------------------------------------------

# Input Data

## Species configuration file

All species analyzed must be listed in:

    data/config.txt

Example format:

    shortName    ensemblName    completeName    pathToGTF
    human         hsapiens       Homo_sapiens        /PATH/TO/Homo_sapiens.GRCh38.109.gtf
    mouse         mmusculus      Mus_musculus       /PATH/TO/Mus_musculus.GRCm39.109.gtf 
    dog         clfamiliaris        Canis_lupus_familiaris         /PATH/TO/Canis_lupus_familiaris.ROS_Cfam_1.0.109.gtf

Columns description:
  - shortName: short identifier
  - ensemblName: Ensembl species name
  - completeName: Full species identifier 
  - pathToGTF: path to the genome annotation (GTF)

N.B : These names follow the Ensembl nomenclature. 

------------------------------------------------------------------------

# Module Descriptions

## 1. Gene Extraction

Directory:

    1_extractionGenes

This step extracts standardized gene information from GTF annotations.

Extracted fields:

-   gene_id
-   gene_name
-   gene_biotype
-   genomic coordinates
-   strand

Output:

    results_gnInfo/species_gnInfo.tsv

Run module:

``` bash
bash 1_extractionGenes/extract_genes.sh data/config.txt
```

------------------------------------------------------------------------

## 2. Protein‑Coding Gene Orthology Extraction

Directory:

    2_extractionOrthologyPCG

Step:
1.  Orthology relationships between protein‑coding genes are retrieved automatically from **Ensembl BioMart**.

The user can indicate the Ensembl version he wants to use, if nothing is indicated, the latest version will be use.

Output:

    results/speciesA-speciesB_homology.tsv

Run module:

``` bash
bash 2_extractionOrthologyPCG/run_OrthologyExtraction.sh [Ensembl Version]
```

------------------------------------------------------------------------

## 3. Synteny‑Based lncRNA Orthology

Directory:

    3_synteny

This module identifies lncRNA orthology based on conserved genomic neighborhoods of protein‑coding genes.

Steps:
1.  Identify flanking PCGs around each lncRNA
2.  Identify orthologous PCG pairs
3.  Detect lncRNAs located between orthologous PCGs

Output:

    results_synteny/speciesA-speciesB_synteny.tsv
    results_syntenyMerged/speciesA-speciesB_synteny.tsv

Run module:

``` bash
bash 3_synteny/run_synteny.sh data/config.txt
```

------------------------------------------------------------------------

## 4. FEELnc Genomic Context Orthology

Directory:

    4_FEELnc

lncRNAs are classified according to their genomic configuration relative to nearby PCGs using FEELnc.
Orthology is inferred when configurations are conserved between species.

Steps:
1.  FEELnc classification
2.  Transcript → gene level conversion
3.  Cross‑species comparison

Outputs:

    results_orthoFeelnc/
    results_orthoFeelnc_agg/
    results_orthoFeelnc_merged/

Run module:

``` bash
bash 4_FEELnc/run_orthoFEELnc.sh data/config.txt
```

------------------------------------------------------------------------

## 5. Genome Alignment Orthology (Ensembl Compara)

Directory:

    5_compara

This module queries **Ensembl Compara genomic alignments** generated with the Mercator--Pecan pipeline.

For each lncRNA:

-   the genomic region is queried
-   aligned regions across species are retrieved
-   conserved regions are recorded

Outputs:

    output/
    output_isMatching/

Run module:

``` bash
bash 5_compara/run_compara_pipeline.sh data/config.txt
```

------------------------------------------------------------------------

## Running the Full Pipeline

To run the entire workflow:

``` bash
bash run_all_methods.sh [options] [config_file]
```

Options:
-   --all        Run the full pipeline (default)
-   --synteny    Run method 1 only
-   --feelnc     Run method 2 only
-   --compara    Run method 3 only

------------------------------------------------------------------------

# Output Overview

The pipeline produces orthology predictions from **three independent methods**.

  Method                 Output directory
  ---------------------- ------------------------------
  Synteny                3_synteny/results_syntenyMerged
  FEELnc configuration   4_FEELnc/results_orthoFeelnc_merged
  Genome alignment       5_compara/results_isMatching

These outputs can be integrated to obtain **high‑confidence candidate lncRNA orthologs**.

------------------------------------------------------------------------

# Typical Workflow

1.  Prepare `config.txt`
2.  Provide genome annotation files (GTF)
3.  Run each methods
4.  Integrate orthology predictions from the three methods.

------------------------------------------------------------------------

# Comments / Questions / Bugs / TODO

- Fabien Degalez
- Sandrine Lagarrigue

------------------------------------------------------------------------

# Authors

Fabien Degalez

------------------------------------------------------------------------

# License

MIT License
