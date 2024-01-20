# eDNA scientific analysis

This document contains scripts and notes for my part of the eDNA analyses

- [Main eDNA analyses](#main-eDNA-analyses)
  - [Phylogenetic/Genetic diversity](#phylogenetic/genetic-diversity)
- [Main comparisons with OBIS data](#main-comparisons-with-OBIS-data)
  - [Fraction of species (all, fish, threatened) detected](#Fraction-of-species-all-fish-threatened-detected)
  - [Functional traits detected](#Functional-traits-detected)
  - [Diversity detected](#Diversity-detected)
  - [Number of new species added and in which taxonomic groups](#Number-of-new-species-added-and-in-which-taxonomic-groups)

## Main eDNA analyses
### Phylogenetic/Genetic diversity

* Loci - 12S(region 1), 12S (region 2), 16S, COI

* Use genetic distance (# nucleotide changes) as a measure of genetic diversity in a species at a site for a locus

* AFD calculates a measure of genetic difference for a species between sites

* Could calculate genetic distance for a locus for a whole site - but this doesn't take into account # of species, obviously the more species, the more diversity


* Haloptypes = ASVs
* Species = OTUs
* Loci
* Sites

* Phylodiversity
1. Reduce to concatenated species phylogeny
2. One sequence per species
3. For all data
4. Calculate PD per site
5. Divide by species to compare sites

Genetic Diversity
1. Align loci
2. Calculate Shannons Index (frequency of haplotypes taken into account) or pairwise distance
3. Divide by # of species (because more species will equal more haplotypes)
4. Average across the 4 loci to get a metric of diversity per site


* Could take sequence and calculate frequency of that sequence across site to create a table


References:

[Environmental DNA reveals the genetic diversity and population structure of an invasive species in the Laurentian Great Lakes](https://doi.org/10.1073/pnas.2307345120)

[Can metabarcoding resolve intraspecific genetic diversity changes to environmental stressors? A test case using river macrozoobenthos](https://doi.org/10.3897/mbmg.4.51925)

[Opportunities and inherent limits of using environmental DNA for population genetics](https://doi.org/10.1002/edn3.448)

[Metabarcoding reveals high-resolution biogeographical and metaphylogeographical patterns through marine barriers](https://doi.org/10.1111/jbi.14548)

## Main comparisons with OBIS data
### Fraction of species (all, fish, threatened) detected
### Functional traits detected
### Diversity detected
### Number of new species added and in which taxonomic groups

