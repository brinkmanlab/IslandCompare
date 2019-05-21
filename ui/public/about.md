# IslandCompare

IslandCompare has been developed to facilitate the analysis of microbial population datasets. While [IslandViewer](http://www.pathogenomics.sfu.ca/islandviewer/browse/) continues to serve as an excellent resource for visual analysis of genomic islands for individual genomes, IslandCompare enables the user to submit sets of genomes to be analysed and visualized together. The interactive visual output is designed to allow users to visually distinguish genomic island and alignment differences between their genomes.

The IslandCompare analysis pipeline is implemented in Galaxy and we have made the tool open source to allow users to install the pipeline independently and modify it as they desire. For more details, please NOLAN FILL IN HERE – Make sure to mention the open source licence and run the content for downloading by Fiona.

This software and website are developed and maintained by the [Brinkman lab](http://www.brinkman.mbb.sfu.ca) at Simon Fraser University, Canada. Please [contact us](#/contact) with any questions or comments.

# Analysis Overview

IslandCompare performs a phylogenetic analysis and alignment of the submitted genomes in order to facilitate the visual presentation of data and help users interpret the genomic island results. Similar to IslandViewer, genomic islands are predicted in IslandCompare by two integrated genomic island prediction software – IslandPath-DIMOB and Sigi-HMM. Unique to IslandCompare is the cross-genome clustering of genomic islands. Each cluster will be consistently coloured across genomes to enable the comparison of genomic island content. Finally, antimicrobial resistance genes are also predicted in the genomes and annotated in the visual output.

# Phylogenomics

**Parnsp** is a tool from the Harvest suite of bioinformatics software [Treangen et al., 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0524-x) that is designed to compute phylogenies for highly similar genome sequences. Parsnp evaluates regions present across all genomes (core genome regions) and identifies single nucleotide polymorphisms (SNPs), as well as other variants in these regions in order to determine phylogenetic relatedness of the isolates under analysis.

# Alignment

**Mauve** implements a seed-and-extend hashing method to identify aligned blocks within a set of sequences [Darling et al., 2010](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011147). Mauve can efficiently compute whole genome alignments. For IslandCompare, these alignments are displayed in the interactive visual to support the investigation of putative genomic

islands. Having the alignment presented alongside the genomic island annotations helps to visually evaluate whether islands annotated in a subset of genomes align elsewhere in adjacent genomes on the phylogeny or seem unique to that subset.

# Genomic Island Prediction

**IslandPath-DIMOB** predicts genomic islands based on dinucleotide bias and the presence of mobility genes. A sliding window encompassing six contiguous open reading frames at a time moves along the genome. For each set of open reading frames, the dinucleotide usage relative to that of the entire genome is measured. An average of each of the 16 dinucleotides is computed and this is taken to be the dinucleotide bias value. Regions with dinucleotide bias that differs significantly from the whole genome are retained and merged with neighboring regions of bias in accordance with a number of cut-offs. Regions of dinucleotide bias that do not contain at least one mobility gene are filtered out. Mobility genes are identified by both a text search and homology analysis using HMMer. For our most recent IslandPath-DIMOB publication, please see [Bertelli et al, 2018](https://academic.oup.com/bioinformatics/article/34/13/2161/4904263).

**Sigi-HMM** (see [Waack et al., 2006](http://www.biomedcentral.com/1471-2105/7/142)) is a sequence composition method that is part of the Columbo software package. This method uses a Hidden Markov Model (HMM) and measures codon usage to identify possible genomic islands. Sigi-HMM has been shown to be one of the highest precision genomic island prediction tools ([Bertelli et al., 2018](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby042/5032564)).

# Clustering

**Mash** is used as part of a two-step clustering to approximate the similarity between genomic island sequences. Mash applies MinHashing to sequence k-mers in order to determine sequence similarity, as measured by an approximation of the Jaccard index [Ondov et al.,2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x).

**MCL (Matrix Cluster Algorithm)** is an approach to resolving similarity matrices into clusters. The algorithm first converts the matrix into a weighted graph and then performs iterative graph traversals, trimming underutilized edges. Denser regions of the graph will ultimately be grouped into clusters [van Dongen, 2000](https://micans.org/mcl/index.html?sec_thesisetc).

# Antimicrobial Resistance Gene Annotations

**Resistance Gene Identifier (RGI)** is a tool developed to make antimicrobial resistance gene predictions using the curated [Comprehensive Antibiotic Resistance Database (CARD)](https://card.mcmaster.ca). CARD contains entries for both protein variant genes (mutations conferring resistance) and protein homolog genes involved in resistance with curated detection cut-offs. For a full description of CARD and RGI, please see [Jia et al., 2017](https://academic.oup.com/nar/article/45/D1/D566/2333912).

# Support for Draft Genomes

Draft genomes are processed in much the same way as IslandViewer. Upon selecting a set of draft genomes to analyse, the contigs in a single genome/file will need to be concatenated in order to perform the analysis in IslandCompare. If draft genomes are included in a job submission, you will need to select either “OPTION-FOR-STITCHING-IN-ORDER” to concatenate the contigs in the order submitted, or select a reference genome against which all genomes in the analysis will be aligned and concatenated. If the option to align against a reference is chosen, contigs unique to the custom genome or contigs that could be placed in several position according to the reference genome (such as identical transposases that could not be solved by short read assembly software) will remain unaligned and placed at the end of the pseudochromosome. These contigs that could not be ordered are shown in ???. Contig gaps are indicated by ???. Contigs placed in this unaligned region should be evaluated with extra caution.

# Analysis Considerations

IslandCompare integrates two sequence composition genomic island prediction methods: Sigi-HMM and IslandPath-DIMOB. Genomic islands that have been more anciently incorporated into a genome can ameliorate into the given genome, meaning their sequence composition will come to more closely resemble that of the host genome. Genomic islands transferred between closely related species will also have a more similar compositional signature. Sequence composition genomic island prediction methods may have difficulty detecting ancient genomic islands due to amelioration and other genomic islands with less distinguishable sequence bias. Additionally, although both methods have high precision, they can make false predictions due to the normal variation in sequence composition that can occur in bacterial genomes. 
