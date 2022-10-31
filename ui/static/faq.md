<a href="/analysis?tour=tour" target="_self" class="button-icon"><i class="icon icon-tutorial"></i>Tutorial</a>

# FAQ

How does IslandCompare remember who I am between sessions?
:   A token is stored in your browsers cache that persists between sessions and uniquely identifies you. If you clear your cache, change browsers, or change computers a new token will be issued and you will not be able to access your previous analysis or uploads.
    This token is also added to the url displayed in the browser address bar for the Analysis or Job History pages. If you bookmark this address, you can retain your token until your account is purged due to inactivity (3 months).  

How do I cite IslandCompare?
:   Please see the [publications](/publications) page. 

Why does IslandCompare not include IslandPick/Islander/virulence factor annotations like IslandViewer does?
:   IslandCompare is currently under active development and this is the first version released for public use. We are working to incorporate these features and plan to make them available in a future version of IslandCompare.

    Furthermore, implementing IslandPick for population datasets will be non-trivial. An approach to selecting reference genomes for a set of genomes under analysis has not yet been developed and would need to carefully be thought out and evaluated before being offered.

How can I view pre-computed complete genomes?

:   For viewing individual, pre-computed genomes, we recommend using [IslandViewer 4](http://www.pathogenomics.sfu.ca/islandviewer/browse/). All complete bacterial and archaeal genomes available for download on NCBI are pre-computed through IslandViewer.

How are the phylogeny/alignment/AMR annotations determined?

:   For detailed information on the tools used to perform each stage of the IslandCompare analysis, please see our [About](/about) page. Briefly, phylogenies are computed based on snps in the core genome using [parsnp](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0524-x). The genomes are aligned with [mauve](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011147). Antimicrobial resistance genes are predicted with the [Resistance Gene Identifier](https://academic.oup.com/nar/article/45/D1/D566/2333912) (uses curated genes in the Comprehensive Antibiotic Resistance Database).

Can I set IslandCompare up independently so that I can customize the workflow for my analysis?

:   Please see the [download page](/download).

Why does my input file fail to upload or run?

:   Most often, problems arise due to incorrectly formatted input files. Please ensure your input file follows GenBank or EMBL formats with all necessary fields, including annotated coding sequences (CDSs). See examples for more details: [GenBank](http://www.pseudomonas.com/downloads/pseudomonas/pgd_r_18_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.gbk), [EMBL](http://www.pseudomonas.com/downloads/pseudomonas/pgd_r_18_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.embl).

    In some cases, SIGI-HMM will fail to run due to issues beyond our control (SIGI-HMM software was written by others), but IslandPath-DIMOB results will still be available. If you notice your jobs are taking longer than a few hours to complete, please don't hesitate to contact us and we would be happy to help identify the problem.

Why do I have to annotate genes in my genomes prior to running IslandCompare?

:   Both IslandPath-DIMOB and Sigi-HMM rely on the gene annotations provided to predict genomic islands (please see our [About](/about) page for more detailed information). As such, the gene annotations provided are essential for making accurate genomic island predictions in IslandCompare. For an excellent, general purpose annotation tool for microbial genomes, we would recommend [prokka](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517), which can be downloaded from [GitHub](https://github.com/tseemann/prokka) or installed as a conda environment.

How should I be annotating my genomes prior to submitting to IslandCompare?
 
:   While we do not endorse a specific annotation tool at this time, an important point we urge all users to keep in mind is to be consistent in your gene annotations. The integrated GI prediction software are dependant on the gene annotations. Therefore, you should be sure to use the same annotation software and version for all of the genomes in your analysis for the most consistent results. Further, when annotating your files, ensure that the annotations are in order and that your file is formatted in adherence with the standards for GenBank and EMBL files.

How do the genomic island prediction tools integrated into IslandCompare work and why were they chosen?

:   IslandCompare integrates predictions from two methods for user genomes: IslandPath-DIMOB, and SIGI-HMM. Both methods have a high precision (>85%) and hence make few but some false positive predictions and we encourage users to carefully check the results. Both tools are sequence composition based tools. IslandPath-DIMOB identifies genomic islands by searching for sets of genes with dinucleotide usage bias and at least one mobility gene, while Sigi-HMM identifies genomic islands by looking for regions with codon usage bias. Please see our [About](/about) page or the most recent publications for [IslandPath-DIMOB](https://academic.oup.com/bioinformatics/article/34/13/2161/4904263) and [Sigi-HMM](http://www.biomedcentral.com/1471-2105/7/142) for details.

Why is an expected genomic island missed from IslandCompare predictions?

:   Genomic islands that have been more anciently incorporated into a genome can ameliorate into the given genome, meaning their sequence composition will come to more closely resemble that of the host genome. Genomic islands transferred between closely related species will also have a weaker compositional signature. This can complicate predictions by sequence-composition based tools like SIGI-HMM and IslandPath-DIMOB. Additionally, these tools were selected with the intent of prioritizing precision, potentially sacrificing recall to some extent.

What are the issues with running an incomplete genome through IslandCompare?

:   Incomplete genomes are first reordered against a user-selected reference genome. The quality of contig reordering will depend on the sequence similarity between the two organisms and the quality of the draft genomes. Contigs unique to the custom genome or contigs that could be placed in several position according to the reference genome (such as identical transposases that could not be solved by short read assembly software will remain unaligned and placed at the end of the pseudochromosome. These contigs that could not be ordered are shown in ???. Contig gaps are indicated by ???.

    Due to the pitfalls of short read sequencing and the unknown quality of contig reordering against a reference, predictions in IslandCompare by the integrated genomic island prediction tools could falsely predict genomic islands, and could miss real genomics islands. A proper assessment of the accuracy of genomic island prediction in incomplete genomes is being performed, and until such assessment is complete, all genomic island predictions in incomplete genomes through IslandCompare should be carefully evaluated for validity.

What if my microorganism has several replicons?

:   For users wishing to analyse a genome of a microorganisms with multiple replicons we recommend submitting a single genbank/embl file with the replicons in the order you wish them to appear in the visual (ensure that homologous replicons are in the same order for each genome submitted if taking this approach. Then select “OPTION FOR THE IN-ORDER-SUBMITTED STITCHER” when submitting the analysis. The replicons will be concatenated for the presentation of results and output of genomic island borders. Alternatively, if you wish to view each replicon separately, or if the replicons are not closed and require rearrangement against a reference in IslandCompare, it may be best to upload each replicon separately and analyse sets of homologous replicons (eg. Analyse Chromosome I from a set of Vibrio genomes as one submission and Chromosome II as a separate submission.

What is the accuracy of IslandCompare genomic island predictions?

:   IslandCompare incorporates two of the most accurate genomic island prediction methods – IslandPath-DIMOB has an accuracy of 77%, while Sigi-HMM has an accuracy of 92%. Both methods have >85% precision. In some cases, you may want higher recall, for which we suggest using Alien_Hunter, but the number of false predictions will also increase greatly. IslandPick is the most precise/specific, if comparison genomes are available. A number of other genomic island prediction tools are available, including GIHunter, which offers high precision and recall. Please see our most recent review for a more comprehensive overview of computational genomic island prediction software [ Bertelli et al., 2018](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby042/5032564).

What is the source of antimicrobial resistance gene annotations? Which resources should I cite if I make use of these annotations?

:   Antimicrobial resistance genes are annotated with the Resistance Gene Identifier, which draws from the curated antimicrobial resistance genes in the [Comprehensive Antibiotic Resistance Database (CARD)](https://card.mcmaster.ca). See [Jia et al., 2017](https://academic.oup.com/nar/article/45/D1/D566/2333912) or our [About](/about) page for more details.

How can I check if a genomic island might be a pathogenicity or resistance island?

:   Antimicrobial resistance genes, as predicted by the Resistance Gene Identifier, are annotated with pink flags appended to the genome in the visual IslandCompare output. You can also explore AMR genes associated with a given genomic island by clicking on the island, which will pull up a separate page with the islands in that cluster and their gene annotations. Keep in mind that the annotations available here are by no means complete and are not available for every genome so if you do not see any annotated islands, there still may be a resistance island in your genome of interest. Further investigation of the types of genes within the genomic island predictions will be crucial for classifying these types of genomic islands.

    Unfortunately, we do not specifically flag virulence factors associated with pathogenicity islands at this time and evaluating these features will need to be performed more manually for the time being. Please check back when updated versions become available for this feature though!

Where can I find more information about performing genomic island prediction using computational methods?

:   If you are interested in learning more about the computational prediction of genomic islands, please see the review article by [Bertelli et al., 2018](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby042/5032564) in Briefings in Bioinformatics and also note this earlier review as well [Langille et al., 2010](http://www.nature.com/nrmicro/journal/v8/n5/full/nrmicro2350.html).

How do I get curated GI predictions included for my genomes?

:   IslandCompare currently supports the detection curated, known Salmonella enterica genomic islands and aims to expand to support curated island detection in additional taxa in the future. In the visualization you will see these islands listed with the "Curated Genomic Islands" predictor. Additionally, if you hover your mouse over the island, it will show the genomic island name detected. If you are analyzing Salmonella enterica genomes in IslandCompare and want to see curated GI predictions for your genomes, make sure that the "Organism" field in your submitted file is filled out and includes the text "Salmonella enterica".
