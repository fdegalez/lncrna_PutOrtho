# lncRNA Orthology Inference Pipeline


# Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Pipeline Overview](#pipeline-overview)
- [Requirements / Software Dependencies](#requirements--software-dependencies)
  - [R](#r)
  - [Perl](#perl)
  - [FEELnc](#feelnc)
- [Input Data](#input-data)
- [Module Descriptions](#module-descriptions)
  - [1. Gene Extraction](#1-gene-extraction)
  - [2. Protein-Coding Gene Orthology Extraction](#2-protein-coding-gene-orthology-extraction)
  - [3. Synteny-Based lncRNA Orthology](#3-syntenybased-lncrna-orthology)
  - [4. FEELnc Genomic Context Orthology](#4-feelnc-genomic-context-orthology)
  - [5. Genome Alignment Orthology](#5-genome-alignment-orthology-ensembl-compara)
  - [(6.) Summary Scripts](#6-summary-scripts)
- [Running the Full Pipeline](#running-the-full-pipeline)
- [Summary Scripts](#summary-scripts)
- [Output Overview](#output-overview)
- [Typical Workflow](#typical-workflow)
- [Notes and Limitations](#notes-and-limitations)
- [Authors](#authors)
- [License](#license)


# Overview

This repository provides a modular workflow to infer **putative orthologous relationships between long non‑coding RNAs (lncRNAs)** across multiple species.

Because lncRNAs evolve rapidly and often lack strong sequence conservation, classical orthology inference methods designed for protein‑coding genes are frequently insufficient. This pipeline integrates **three complementary strategies** to detect lncRNA orthologs:

1.  **Synteny-based inference**
2.  **Genomic context conservation (FEELnc classification)**
3.  **Genome alignment conservation (Ensembl Compara, here: Mercator--Pecan alignments)**

The workflow provides a flexible framework allowing users to detect candidate lncRNA orthologs using complementary evidence.
For a more general overview, you could point to [the associated paper](https://www.biorxiv.org/content/10.1101/2024.10.03.616473v1).



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
    │   └── run_summary.sh
    │
    └── run_all_methods.sh



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



## Perl

Required for the Ensembl API.
    [Perl5+](https://www.perl.org/) ≥ 5.32

Required modules:
    - [Ensembl API](https://www.ensembl.org/info/docs/api/api_installation.html) : tested with e! v104  
    - [Bioperl](http://www.bioperl.org/wiki/Main_Page) : tested with version 1.7.8  
    - [TimeHiRes](https://metacpan.org/pod/Time::HiRes) : tested with version 1.9764      



## FEELnc

Used for genomic context classification of lncRNAs.

Repository:

- [FEELnc](https://github.com/tderrien/FEELnc) : test with version 0.2.1 - 2022-07-20
    - [FEELnc_classifier.pl](https://github.com/tderrien/FEELnc#3--feelnc_classifierpl) : Classify lncRNAs based on their genomic localization with others transcripts.
    - [FEELnc_tpLevel2gnLevelClassifcation.R](https://github.com/tderrien/FEELnc/blob/master/scripts/FEELnc_tpLevel2gnLevelClassification.R) : Transformation of transcript-level configurations to gene-level models. 




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
  - `shortName`: short identifier
  - `ensemblName`: Ensembl species name
  - `completeName`: Full species identifier 
  - `pathToGTF`: path to the genome annotation (GTF)

N.B : These names follow the Ensembl nomenclature. 



# Module Descriptions

## 1. Gene Extraction

Directory:

    1_extractionGenes

This step extracts standardized gene information from GTF annotations.

Extracted fields:

-   `gene_id`
-   `gene_name`
-   `gene_biotype`
-   `genomic coordinates`
-   `strand`

Output:

    results_gnInfo/species_gnInfo.tsv

Run module:

``` bash
bash 1_extractionGenes/extract_genes.sh data/config.txt
```
**Note**: The pipeline was designed using **Ensembl GTF annotations**, however, custom GTF files can also be used provided that they follow the same conventions, in particular:
- standard GTF column structure
- presence of the following attributes:
  - `gene_id`
  - `gene_name`
  - `gene_biotype`

If these attributes follow the standard Ensembl naming scheme, the pipeline should function without modification.

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



**Note**: The default behavior is to automatically retrieve protein-coding gene orthology relationships from Ensembl BioMart using Ensembl gene identifiers.
When using **custom annotations**, two alternatives are possible:
1. **Provide a custom orthology file**  
   The user can directly place orthology files in:

       2_extractionOrthologyPCG/results/

   These files must follow the usual naming convention:

       speciesA-speciesB_homology.tsv

   Species names must correspond to those defined in `config.txt`.

2. **Convert custom gene identifiers to Ensembl identifiers (recommended)**

   Alternatively, a preliminary step can be performed to map custom gene identifiers to Ensembl gene IDs and convert the GTFs according to this nomenclature.  
   This allows the pipeline to run normally for all modules

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


## 4. FEELnc Genomic Context Orthology

The pipeline assumes that `FEELnc_classifier.pl` is available in the system PATH.
If this is not the case, the path to the executable must be manually modified in the script: `4_FEELnc/1_extractionFEELnc.sh`

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

Note: this step may take a significant amount of time depending on the number of species and the size of the genome annotations, as FEELnc classification can be computationally intensive.

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

Note: this step can be time-consuming.  
Each lncRNA query is performed independently through the Ensembl Compara API.  
To avoid excessive requests that could lead to temporary API blocking, queries are performed sequentially, which can make this step relatively slow when large numbers of lncRNAs are analyzed.

### Notes on the Compara alignment

At present, the pipeline queries the **Mercator–Pecan multiple genome alignment of 65 amniote species**, which is the alignment set currently used by default in the Ensembl Compara database.
In future updates, the script will be extended to allow users to specify the **Compara alignment version or species set** they wish to query directly from the Ensembl database.
Using a fully custom multiple genome alignment is not currently supported within the pipeline itself. However, this can be achieved externally by the user. In that case, results obtained from the custom alignment can be integrated into the workflow by generating tables following the **same format as those produced in `results_isMatching/`**, which can then be used by the downstream summary scripts.

## (6.) Summary Scripts

Two additional scripts are provided in the `6_summary` directory to  generate summary tables of orthology detection across the three methods and can be launched through the `run_summary.sh` script :


1. **Pairwise summary (`1_summaryPairwise.R`)**  
   Generates a table summarizing orthology detection for each species pair.

2. **Species-centered summary (`2_summaryAll.R`)**  
   Generates a merged table considering one reference species against all other species.

These scripts are **not integrated into the global pipeline (`run_all_methods.sh`)** because their correct execution depends on which modules of the pipeline have been run beforehand. 
The optimal case is when the three methods are used howerver users may therefore want to run summaries after running **only a subset of the methods**.

------------------------------------------------------------------------

## Running the Full Pipeline

To run the entire workflow:

``` bash
bash run_all_methods.sh [options] [config_file]
```

Options:
-   `--all
`: Run the full pipeline (default)
-   `--synteny`: Run method 1 only
-   `--feelnc`: Run method 2 only
-   `--compara`: Run method 3 only



# Output Overview

The pipeline produces orthology predictions from **three independent methods**:

    Method                 Output directory
    Synteny                3_synteny/results_syntenyMerged
    FEELnc configuration   4_FEELnc/results_orthoFeelnc_merged
    Genome alignment       5_compara/results_isMatching

These outputs can be integrated to obtain **high‑confidence candidate lncRNA orthologs**.


# Typical Workflow

1.  Prepare `config.txt`
2.  Provide genome annotation files (GTF)
3.  Run each methods
4.  Integrate orthology predictions from the three methods.



# Notes and Limitations

- When the pipeline is executed, logging output is currently minimal and may not always clearly reflect the progress of each step.
- Improvements to logging and execution reporting are planned for future versions of the pipeline.


# Comments / Questions / Bugs 

- Users are encouraged to report issues, bugs, or unexpected behavior through the GitHub Issues section of this repository.
- Because the pipeline integrates several external resources (Ensembl BioMart, Ensembl Compara API, FEELnc), occasional changes in these services may affect execution.
- Some parts of the pipeline (particularly logging and execution reporting) are still being improved and may evolve in future versions.
- Feedback from users applying the workflow to additional species or custom genome annotations is particularly appreciated.
- Feel free to contact us directly by email for any questions/sugestions.


# Authors

- Fabien Degalez
- Sandrine Lagarrigue


# License

MIT License
