# Awesome Bioinformatics Benchmarks
[![Build Status](https://travis-ci.org/j-andrews7/Awesome-Bioinformatics-Benchmarks.svg?branch=master)](https://travis-ci.org/j-andrews7/Awesome-Bioinformatics-Benchmarks)

A curated list of bioinformatics benchmarking papers and resources.

The credit for this format goes to Sean Davis for his [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) repository and Ming Tang for his [ChIP-seq-analysis](https://github.com/crazyhottommy/ChIP-seq-analysis) repository. 

If you have a benchmarking study that is not yet included on this list, please make a [Pull Request](https://github.com/j-andrews7/Awesome-Bioinformatics-Benchmarks/pulls).

## Rules for Included Papers
 - Papers must be objective comparisons of 3 or more tools/methods.
 - Papers should **not** be from authors showing why their tool/method is better than others.
 - Benchmarking data should be publicly available or simulation code/methods must be well-documented and reproducible.
 
Additional guidelines/rules may be added as necessary.

## Format
Please include the following information when adding papers.

**Title:**

**Authors:**

**Journal Info:**

**Description:**

**Tools/methods compared:**

**Recommendation(s):**

**Additional links (optional):**

# Tool/Method Sections
Additional sections/sub-sections can be added as needed.


## DNase & ChIP-seq

### Peak Callers

**Title:** [A Comparison of Peak Callers Used for DNase-Seq Data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0096303)

**Authors:** Hashem Koohy, et al.

**Journal Info:** PLoS ONE, May 2014

**Description:** This paper compares four peak callers specificty and sensitivity on DNase-seq data from two publications composed of three cell types, using ENCODE data for the same cell types as a benchmark. 
The authors tested multiple parameters for each caller to determine the best settings for DNase-seq data for each.

**Tools/methods compared:** `F-seq`, `Hotspot`, `MACS2`, `ZINBA`.

**Recommendation(s):** [F-seq](https://github.com/aboyle/F-seq) was the most sensitive, though [MACS2](https://github.com/taoliu/MACS) and [Hotspot](https://github.com/rthurman/hotspot) both performed competitively as well. 
ZINBA was the least performant by a massive margin, requiring much more time to run, and was also the least sensitive.

---

**Title:** [Features that define the best ChIP-seq peak calling algorithms](https://academic.oup.com/bib/article/18/3/441/2453291)

**Authors:** Reuben Thomas, et al.

**Journal Info:** Briefings in Bioinformatics, May 2017

**Description:** This paper compared six peak calling methods on 300 simulated and three real ChIP-seq data sets across a range of significance values. 
Methods were scored by sensitivity, precision, and F-score.

**Tools/methods compared:** `GEM`, `MACS2`, `MUSIC`, `BCP`, `Threshold-based method (TM)`, `ZINBA`.

**Recommendation(s):** Varies. [BCP](http://ranger.sourceforge.net/manual1.18.html) and [MACS2](https://github.com/taoliu/MACS) performed the best across all metrics on the simulated data. 
For Tbx5 ChIP-seq, [GEM](http://groups.csail.mit.edu/cgs/gem/) performed the best, with BCP also scoring highly. 
For histone H3K36me3 and H3K4me3 data, all methods performed relatively comparably with the exception of ZINBA, which the authors could not get to run properly. 
[MUSIC](https://github.com/gersteinlab/MUSIC) and BCP had a slight edge over the others for the histone data. 
More generally, they found that methods that utilize variable window sizes and Poisson test to rank peaks are more powerful than those that use a Binomial test. 

## RNA-seq

### Normalisation Methods

**Title:** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](https://academic.oup.com/bib/article/14/6/671/189645). 

**Authors:** Marie-Agnès Dillies, et al.

**Journal Info:** Briefings in Bioinformatics, November 2013

**Description:** This paper compared seven RNA-seq normalization methods in the context of differential expression analysis on four real datasets and thousands of simulations.

**Tools/methods compared:** `Total Count (TC)`, `Upper Quartile (UQ)`, `Median (Med),` `DESeq`, `edgeR`, `Quantile (Q)`, `RPKM`.

**Recommendation(s):** The authors recommend [DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html) ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) now available as well) or [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), as those methods are robust to the presence of different library sizes and compositions, whereas the (still common) Total Count and RPKM methods are ineffective and should be abandoned.

### Differential Gene Expression


### Cell-Type Deconvolution

**Title:** [Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology](https://academic.oup.com/bioinformatics/article/35/14/i436/5529146)

**Authors:** Markus List\*, Tatsiana Aneichyk\*, et al.

**Journal Info:** Bioinformatics, July 2019

**Description:** This paper benchmarks and compares seven methods for computational deconvolution of cell-type abundance in bulk RNA-seq samples. Each method was tested on both simulated and true bulk RNA-seq samples validated by FACS.

**Tools/methods compared:** `quanTIseq`, `TIMER`, `CIBERSORT`, `CIBERSORT abs. mode`, `MCPCounter`, `xCell`, `EPIC`.

**Recommendation(s):** Varies. In general, the authors recommend [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) and [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) due to their overall robustness and absolute (rather than relative) scoring, though [xCell](http://xcell.ucsf.edu/) is recommended for binary presence/absence of cell types and [MCPcounter](https://github.com/ebecht/MCPcounter) was their recommended relative scoring method.

**Additional links:** The authors created an [R package called immunedeconv](https://github.com/icbi-lab/immunedeconv) for easy installation and use of all these methods. For developers, they have made available their [benchmarking pipeline](https://github.com/icbi-lab/immune_deconvolution_benchmark) so that others can reproduce/extend it to test their own tools/methods.

## RNA/cDNA Microarrays


## Variant Callers

### Germline SNP/Indel Callers

**Title:** [Systematic comparison of germline variant calling pipelines cross multiple next-generation sequencers](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6597787/)

**Authors:** Jiayun Chen, et al.

**Journal Info:** Scientific Reports, June 2019

**Description:** This paper compared three variant callers for WGS and WES samples from NA12878 across five next-gen sequencing platforms

**Tools/methods compared:** `GATK`, `Strelka2`, `Samtools-Varscan2`.

**Recommendation(s):** Though all methods tested generally scored well, [Strelka2](https://github.com/Illumina/strelka) had the highest F-scores for both SNP and indel calling in addition to being the most computationally performant.

---

**Title:** [Comparison of three variant callers for human whole genome sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6294778/) 

**Authors:** Anna Supernat, et al.

**Journal Info:** Scientific Reports, December 2018

**Description:** The paper compared three variant callers for WGS samples from NA12878 at 10x, 15x, and 30x coverage.

**Tools/methods compared:** `DeepVariant`, `GATK`, `SpeedSeq`.

**Recommendation(s):** All methods had similar F-scores, precision, and recall for SNP calling, but [DeepVariant](https://github.com/google/deepvariant) scored higher across all metrics for indels at all coverages.

---

**Title:** [A Comparison of Variant Calling Pipelines Using Genome in a Bottle as a Reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4619817/)

**Authors:** Adam Cornish, et al.

**Journal Info:** BioMed Research International, October 2015

**Description:** This paper compared 30 variant calling pipelines composed of six different variant callers and five different aligners on NA12878 WES data from the "Genome in a Bottle" consortium.

**Tools/methods compared:** 
 - Variant callers: `FreeBayes`, `GATK-HaplotypeCaller`, `GATK-UnifiedGenotyper`, `SAMtools mpileup`, `SNPSVM`
 - Aligners: `bowtie2`, `BWA-mem`, `BWA-sampe`, `CUSHAW3`, `MOSAIK`, `Novoalign`.

**Recommendation(s):** [Novoalign](http://www.novocraft.com/products/novoalign/) with [GATK-UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) exhibited the highest sensitivity while producing few false positives. 
In general, [BWA-mem](https://github.com/lh3/bwa) was the most consistent aligner, and `GATK-UnifiedGenotyper` performed well across the top aligners (BWA, bowtie2, and Novoalign).



### Somatic SNV/Indel callers

**Title:** [Evaluation of Nine Somatic Variant Callers for Detection of Somatic Mutations in Exome and Targeted Deep Sequencing Data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4803342/)

**Authors:** Anne Bruun Krøigård, et al.

**Journal Info:** PLoS ONE, March 2016

**Description:** This paper performed comparisons between nine somatic variant callers on five paired tumor-normal samples from breast cancer patients subjected to WES and targeted deep sequencing.

**Tools/methods compared:** `EBCall`, `Mutect`, `Seurat`, `Shimmer`, `Indelocator`, `SomaticSniper`, `Strelka`, `VarScan2`, `Virmid`.

**Recommendation(s):** [EBCall](https://github.com/friend1ws/EBCall), [Mutect](https://github.com/broadinstitute/mutect), [Virmid](https://sourceforge.net/p/virmid/wiki/Home/), and [Strelka](https://github.com/Illumina/strelka) (now Strelka2) were most reliable for both WES and targeted deep sequencing. 
[EBCall](https://github.com/friend1ws/EBCall) was superior for indel calling due to high sensitivity and robustness to changes in sequencing depths.

---

**Title:** [Comparison of somatic mutation calling methods in amplicon and whole exome sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3986649/)

**Authors:** Huilei Xu, et al.

**Journal Info:** BMC Genomics, March 2014

**Description:** Using the "Genome in a Bottle" gold standard variant set, this paper compared five somatic mutation calling methods on matched tumor-normal amplicon and WES data.

**Tools/methods compared:**  `GATK-UnifiedGenotyper followed by subtraction`, `MuTect`, `Strelka`, `SomaticSniper`, `VarScan2`.

**Recommendation(s):** [MuTect](https://github.com/broadinstitute/mutect) and [Strelka](https://github.com/Illumina/strelka) (now Strelka2) had the highest sensitivity, particularly at low frequency alleles, in addition to the highest specificity.  


### CNV Callers

**Title:** [Benchmark of tools for CNV detection from NGS panel data in a genetic diagnostics context](https://www.biorxiv.org/content/10.1101/850958v1)

**Authors:** José Marcos Moreno-Cabrera, et al.

**Journal Info:** bioRxiv, November 2019.

**Description:**

**Tools/methods compared:** `DECoN`, `CoNVaDING`, `panelcn.MOPS`, `ExomeDepth`, `CODEX2`.

**Recommendation(s):**

---

**Title:** [An evaluation of copy number variation detection tools for cancer using whole exome sequencing data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5452530/)

**Authors:** Fatima Zare, et al.

**Journal Info:** BMC Bioinformatics, May 2017

**Description:**

**Tools/methods compared:** `ADTEx`, `CONTRA`, `cn.MOPS`, `ExomeCNV`, `VarScan2`, `CoNVEX`.

**Recommendation(s):**



### SV callers

**Title:** [Comprehensive evaluation and characterisation of short read general-purpose structural variant calling software](https://www.nature.com/articles/s41467-019-11146-4)

**Authors:** Daniel L. Cameron, et al.

**Journal Info:** Nature Communications, July 2019

**Description:**

**Tools/methods compared:** `BreakDancer`, `cortex`, `CREST`, `DELLY`, `GRIDSS`, `Hydra`, `LUMPY`, `manta`, `Pindel`, `SOCRATES`.

**Recommendation(s):**
 
 ---
 
 **Title:** [Comprehensive evaluation of structural variation detection algorithms for whole genome sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5)

**Authors:** Shunichi Kosugi, et al.

**Journal Info:** Genome Biology, June 2019

**Description:**

**Tools/methods compared:** `1-2-3-SV`, `AS-GENESENG`, `BASIL-ANISE,` `BatVI`, `BICseq2`, `BreakDancer`, `BreakSeek`, `BreakSeq2`, `Breakway`, `CLEVER`, `CNVnator`,
`Control-FREEC`, `CREST`, `DELLY`, `DINUMT`, `ERDS`, `FermiKit`, `forestSV`, `GASVPro`, `GenomeSTRiP`, `GRIDSS`, `HGT-ID`, `Hydra-sv`, `iCopyDAV`, `inGAP-sv`, `ITIS`,
`laSV`, `Lumpy`, `Manta`, `MATCHCLIP`, `Meerkat`, `MELT`, `MELT-numt`, `MetaSV`, `MindTheGap`, `Mobster`, `Mobster-numt`, `Mobster-vei`, `OncoSNP-SEQ`, `Pamir`, `PBHoney`,
`PBHoney-NGM`, `pbsv`, `PennCNV-Seq`, `Pindel`, `PopIns`, `PRISM`, `RAPTR`, `readDepth`, `RetroSeq`, `Sniffles`, `Socrates`, `SoftSearch`, `SoftSV`, `SoloDel`, `Sprites`,
`SvABA`, `SVDetect`, `Svelter`, `SVfinder`, `SVseq2`, `Tangram`, `Tangram-numt`, `Tangram-vei`, `Tea`, `TEMP`, `TIDDIT`, `Ulysses`, `VariationHunter`, `VirusFinder`, `VirusSeq`, `Wham`.

**Recommendation(s):**


## Single Cell

### Trajectory Inference

**Title:** [A comparison of single-cell trajectory inference methods](https://doi.org/10.1038/s41587-019-0071-9)

**Authors:** Wouter Saelens\*, Robrecht Cannoodt\*, et al.

**Journal Info:** Nat Biotech, April 2019

**Description:** A comprehensive evaluation of 45 trajectory inference methods, this paper provides an unmatched comparison of the rapidly evolving field of single-cell trajectory inference. Each method was scored on accuracy, scalability, stability, and usability. Should be considered a gold-standard for other benchmarking studies.

**Tools/methods compared:** `PAGA`, `RaceID/StemID`, `SLICER`, `Slingshot`, `PAGA Tree`, `MST`, `pCreode`, `SCUBA`, `Monocle DDRTree`, `Monocle ICA`, `cellTree maptpx`, `SLICE`, `cellTree VEM`, `EIPiGraph`, `Sincell`, `URD`, `CellTrails`, `Mpath`, `CellRouter`, `STEMNET`, `FateID`, `MFA`, `GPfates`, `DPT`, `Wishbone`, `SCORPIUS`, `Component 1`, `Embeddr`, `MATCHER`, `TSCAN`, `Wanderlust`, `PhenoPath`, `topslam`, `Waterfall`, `EIPiGraph linear`, `ouijaflow`, `FORKS`, `Angle`, `EIPiGraph cycle`, `reCAT`.

**Recommendation(s):** Varies depending on dataset and expected trajectory type, though [PAGA, PAGA Tree](https://scanpy.readthedocs.io/en/latest/examples.html#trajectory-inference), [SCORPIUS](https://github.com/rcannood/SCORPIUS), and [Slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html) all scored highly across all metrics. Authors wrote an [interactive Shiny app](https://dynverse.org/users/3-user-guide/2-guidelines/) to help users choose the best methods for their data.

**Additional links:** The [dynverse site](https://dynverse.org/) contains numerous packages for users to run and compare results from different trajectory methods on their own data without installing each individually by using Docker. Additionally, they provide several tools for developers to wrap and benchmark their own method against those included in the study. 

### Integration/Batch Correction

### Cell Annotation/Inference

# Contributors

 - Jared Andrews ([@j-andrews7](https://github.com/j-andrews7/))
 - Kevin Blighe ([@kevinblighe](https://github.com/kevinblighe), [biostars](https://www.biostars.org/u/41557/))

