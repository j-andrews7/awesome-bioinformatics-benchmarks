# Awesome Bioinformatics Benchmarks
A curated list of bioinformatics bench-marking papers and resources.
The credit for this format goes to Sean Davis for his [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) repository and Ming Tang for his [ChIP-seq-analysis](https://github.com/crazyhottommy/ChIP-seq-analysis) repository. 
If you have a benchmarking study that is not yet included on this list, please make a [Pull Request](https://github.com/j-andrews7/Awesome-Bioinformatics-Benchmarks/pulls) to add it.

## Rules for Included Papers
 - Papers must be objective comparisons of 3 or more tools/methods.
 - Papers should **not** be from authors showing why their tool/method is better than others.
 - Benchmarking data should be publicly available or simulation code/methods should be obviously stated.
 
Additional guidelines/rules may be added as necessary.

# Tool/Method Sections
Additional sections/sub-sections can be added as needed.


## ChIP-seq

### Peak Callers

## RNA-seq

### Normalisation methods
 - [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](https://www.ncbi.nlm.nih.gov/pubmed/22988256). 2013. `Total Count (TC)` `Upper Quartile (UQ)` `Median (Med)` `DESeq` `edgeR` `Quantile (Q)` `RPKM`


### Differential Gene Expression

## RNA / cDNA Microarrays

## Variant Callers

### Germline SNP/Indel Callers

- [Systematic comparison of germline variant calling pipelines cross multiple next-generation sequencers](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6597787/). 2019. `GATK` `Strelka2` `Samtools-Varscan2`

- [Comparison of three variant callers for human whole genome sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6294778/). 2018. `DeepVariant` `GATK` `SpeedSeq`

 - [A Comparison of Variant Calling Pipelines Using Genome in a Bottle as a Reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4619817/). 2015. `FreeBayes` `GATK HaplotypeCaller` `GATK UnifiedGenotyper` `SAMtools mpileup` `SNPSVM`


### Somatic SNV/Indel callers

 - [Evaluation of Nine Somatic Variant Callers for Detection of Somatic Mutations in Exome and Targeted Deep Sequencing Data
](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4803342/). 2016. `EBCall` `Mutect` `Seurat` `Shimmer` `Indelocator` `SomaticSniper` `Strelka` `VarScan2` `Virmid`
 - [Comparison of somatic mutation calling methods in amplicon and whole exome sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3986649/). 2014. `GATK UnifiedGenotyper followed by subtraction` `MuTect` `Strelka` `SomaticSniper` `VarScan2`

### CNV Callers

 - [Benchmark of tools for CNV detection from NGS panel data in a genetic diagnostics context](https://www.biorxiv.org/content/10.1101/850958v1). 2019 (bioRxiv). `DECoN` `CoNVaDING` `panelcn.MOPS` `ExomeDepth` `CODEX2`

 - [An evaluation of copy number variation detection tools for cancer using whole exome sequencing data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5452530/). 2017. `ADTEx` `CONTRA` `cn.MOPS` `ExomeCNV` `VarScan2` `CoNVEX`


### SV callers

 - [Comprehensive evaluation and characterisation of short read general-purpose structural variant calling software](https://www.nature.com/articles/s41467-019-11146-4). 2019. `BreakDancer `cortex` `CREST` `DELLY` `GRIDSS` `Hydra` `LUMPY` `manta` `Pindel` `SOCRATES`
 
 - [Comprehensive evaluation of structural variation detection algorithms for whole genome sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5). 2019. `1-2-3-SV` `AS-GENESENG` `BASIL-ANISE` `BatVI` `BICseq2` `BreakDancer` `BreakSeek` `BreakSeq2` `Breakway` `CLEVER` `CNVnator` `Control-FREEC` `CREST` `DELLY` `DINUMT` `ERDS` `FermiKit` `forestSV` `GASVPro` `GenomeSTRiP` `GRIDSS` `HGT-ID` `Hydra-sv` `iCopyDAV` `inGAP-sv` `ITIS` `laSV` `Lumpy` `Manta` `MATCHCLIP` `Meerkat` `MELT` `MELT-numt` `MetaSV` `MindTheGap` `Mobster` `Mobster-numt` `Mobster-vei` `OncoSNP-SEQ` `Pamir` `PBHoney` `PBHoney-NGM` `pbsv` `PennCNV-Seq` `Pindel` `PopIns` `PRISM` `RAPTR` `readDepth` `RetroSeq` `Sniffles` `Socrates` `SoftSearch` `SoftSV` `SoloDel` `Sprites` `SvABA` `SVDetect` `Svelter` `SVfinder` `SVseq2` `Tangram` `Tangram-numt` `Tangram-vei` `Tea` `TEMP` `TIDDIT` `Ulysses` `VariationHunter` `VirusFinder` `VirusSeq` `Wham`
