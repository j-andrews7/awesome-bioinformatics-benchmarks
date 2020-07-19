# Awesome Bioinformatics Benchmarks [![Build Status](https://travis-ci.org/j-andrews7/Awesome-Bioinformatics-Benchmarks.svg?branch=master)](https://travis-ci.org/j-andrews7/Awesome-Bioinformatics-Benchmarks)

A curated list of bioinformatics benchmarking papers and resources.

The credit for this format goes to Sean Davis for his [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) repository and Ming Tang for his [ChIP-seq-analysis](https://github.com/crazyhottommy/ChIP-seq-analysis) repository. 

If you have a benchmarking study that is not yet included on this list, please make a [Pull Request](https://github.com/j-andrews7/Awesome-Bioinformatics-Benchmarks/pulls).

## Contents
   - [Rules for Included Papers](#rules-for-included-papers)
   - [Format &amp; Organization](#format--organization)
   - [Tool/Method Sections](#toolmethod-sections)
      - [DNase &amp; ChIP-seq](#dnase--chip-seq)
         - [Peak Callers](#peak-callers)
      - [RNA-seq](#rna-seq)
         - [Alignment/Quantification Methods](#alignmentquantification-methods)
         - [Normalisation Methods](#normalisation-methods)
         - [Differential Gene Expression](#differential-gene-expression)
         - [Gene Set Enrichment Analysis](#gene-set-enrichment-analysis)
         - [Differential Splicing](#differential-splicing)
         - [Transcript Assembly and Quantification](#transcript-assembly-and-quantification)
         - [Cell-Type Deconvolution](#cell-type-deconvolution)
      - [RNA/cDNA Microarrays](#rnacdna-microarrays)
      - [Variant Callers](#variant-callers)
         - [Germline SNP/Indel Callers](#germline-snpindel-callers)
         - [Somatic SNV/Indel callers](#somatic-snvindel-callers)
         - [CNV Callers](#cnv-callers)
         - [SV callers](#sv-callers)
      - [Single Cell](#single-cell)
         - [scRNA Sequencing Protocols](#scrna-sequencing-protocols)
         - [scRNA Analysis Pipelines](#scrna-analysis-pipelines)
         - [scRNA Imputation Methods](#scrna-imputation-methods)
         - [scRNA Differential Gene Expression](#scrna-differential-gene-expression)
         - [Trajectory Inference](#trajectory-inference)
         - [Gene Regulatory Network Inference](#gene-regulatory-network-inference)
         - [Integration/Batch Correction](#integrationbatch-correction)
         - [Dimensionality Reduction](#dimensionality-reduction)
         - [Cell Annotation/Inference](#cell-annotationinference)
         - [Variant Calling](#variant-calling)
         - [ATAC-seq](#atac-seq)
      - [Statistics](#statistics)
         - [False Discovery Rates](#false-discovery-rates)
      - [Microbiome](#microbiome)
         - [Diversity analysis](#diversity-analysis)
   - [Contributors](#contributors)


## Rules for Included Papers
 - Papers must be objective comparisons of 3 or more tools/methods.
 - Papers must be **awesome**. This list isn't meant to chronicle every benchmarking study ever performed, only those that are particularly expansive, well done, and/or provide unique insights.
 - Papers should generally **not** be from authors showing why their tool/method is better than others.
 - Benchmarking data should be publicly available or simulation code/methods must be well-documented and reproducible.

Additional guidelines/rules may be added as necessary.

## Format & Organization
Please include the following information when adding papers.

**Title:**

**Authors:**

**Journal Info:**

**Description:**

**Tools/methods compared:**

**Recommendation(s):**

**Additional links (optional):**

Papers within each section should be ordered by publication date, with more recent papers listed first.

# Tool/Method Sections
Additional sections/sub-sections can be added as needed.


## DNase & ChIP-seq

### Peak Callers

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

---

**Title:** [A Comparison of Peak Callers Used for DNase-Seq Data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0096303)

**Authors:** Hashem Koohy, et al.

**Journal Info:** PLoS ONE, May 2014

**Description:** This paper compares four peak callers specificity and sensitivity on DNase-seq data from two publications composed of three cell types, using ENCODE data for the same cell types as a benchmark. 
The authors tested multiple parameters for each caller to determine the best settings for DNase-seq data for each.

**Tools/methods compared:** `F-seq`, `Hotspot`, `MACS2`, `ZINBA`.

**Recommendation(s):** [F-seq](https://github.com/aboyle/F-seq) was the most sensitive, though [MACS2](https://github.com/taoliu/MACS) and [Hotspot](https://github.com/rthurman/hotspot) both performed competitively as well. 
ZINBA was the least performant by a massive margin, requiring much more time to run, and was also the least sensitive.

## ATAC-Seq

### Normalisation Methods

**Title:** [ATAC-seq normalization method can significantly affect differential accessibility analysis and interpretation](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00342-y)

**Authors:** Jake J. Reske, et al.

**Journal Info:** Epigenetics & Chromatin, April 2020

**Description:** This paper compares the effect of normalization method during differential ATAC-seq analysis. 

**Tools/methods compared:**  `MACS2`, `DiffBind`, `csaw`, `voom`, `DEseq2`, `edgeR`, `limma`

**Recommendation(s):**  This paper compares 8 analytical approaches to calculate ATAC-seq differential accessibility (the description of different combination can be seen in paper Table1). The authors found different analytical approaches can produce very differential chromatin accessibility results using MA-plots. 

The authors also proposed a generalized workflow for differential accessibility analysis, which can be found in [Github](https://github.com/reskejak/ATAC-seq)

**Additional links (optional):**  For ATAC-Seq data anlysis, there is another paper: [From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)

## RNA-seq

### Alignment/Quantification Methods

**Title:** [Alignment and mapping methodology influence transcript abundance estimation](https://www.biorxiv.org/content/10.1101/657874v2)

**Authors:** Avi Srivastava\* &  Laraib Malik\*, et al.

**Journal Info:** bioRXiv, October 2019

**Description:** This paper compares the influence of mapping and alignment on the accuracy of transcript quantification in both simulated and experimental data, as well as the effect on subsequent differential expression analysis.

**Tools/methods compared:** `Bowtie2`, `STAR`, `quasi-mapping`, `Selective Alignment`, `RSEM`, `Salmon`.

**Recommendation(s):** When trying to choose an approach, a choice can be made by the user performing the analysis based on any time-accuracy tradeoff they wish to make. In terms of speed, quasi-mapping is the fastest approach, followed by [Selective Alignment](https://github.com/COMBINE-lab/salmon) (SA) then [STAR](https://github.com/alexdobin/STAR). [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) was considerably slower than all three of these approaches. However, in terms of accuracy, SA yielded the best results, followed by alignment to the genome (with subsequent transcriptomic projection) using STAR and SA (using carefully selected decoy sequences). Bowtie2 generally performed similarly to SA, but without the benefit of decoy sequences, seemed to admit more spurious mappings. Finally, lightweight mapping of sequencing reads to the transcriptome showed the lowest overall consistency with quantifications derived from the oracle alignments. Note: Both Selective Alignment and quasi-mapping are part of the [salmon](https://combine-lab.github.io/salmon/) codebase.

---

### Normalisation Methods

**Title:** [Statistical models for RNA-seq data derived from a two-condition 48-replicate experiment](https://academic.oup.com/bioinformatics/article/31/22/3625/240923)

**Authors:** Marek Gierliński\*, Christian Cole\*, Pietá Schofield\*, Nicholas J. Schurch\*, et al.

**Journal Info:** Bioinformatics, November 2015

**Description:** This paper compares the effect of normal, log-normal, and negative binomial distribution assumptions on RNA-seq gene read-counts from 48 RNA-seq replicates.

**Tools/methods compared:** `normal`, `log-normal`, `negative binomial`.

**Recommendation(s):** Assuming a normal distribution leads to a large number of false positives during differential gene expression. A log-normal distribution model works well unless a sample contains zero counts. Use tools that assume a negative binomial distribution (`edgeR`, `DESeq`, `DESeq2`, etc).

---

**Title:** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](https://academic.oup.com/bib/article/14/6/671/189645). 

**Authors:** Marie-Agnès Dillies, et al.

**Journal Info:** Briefings in Bioinformatics, November 2013

**Description:** This paper compared seven RNA-seq normalization methods in the context of differential expression analysis on four real datasets and thousands of simulations.

**Tools/methods compared:** `Total Count (TC)`, `Upper Quartile (UQ)`, `Median (Med),` `DESeq`, `edgeR`, `Quantile (Q)`, `RPKM`.

**Recommendation(s):** The authors recommend [DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html) ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) now available as well) or [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), as those methods are robust to the presence of different library sizes and compositions, whereas the (still common) Total Count and RPKM methods are ineffective and should be abandoned.

### Differential Gene Expression

**Title:** [How well do RNA-Seq differential gene expression tools perform in a complex eukaryote? A case study in Arabidopsis thaliana.](https://www.ncbi.nlm.nih.gov/pubmed/30726870)

**Authors:** Kimon Froussios\*, Nick J Schurch\*, et al.

**Journal Info:** Bioinformatics, February 2019

**Description:** This paper compared nine differential gene expression tools (and their underlying model distributions) in 17 RNA-seq replicates of *Arabidopsis thaliana*.
Handling of inter-replicate variability and false positive fraction were the benchmarking metrics used.

**Tools/methods compared:** `baySeq`, `DEGseq`, `DESeq`, `DESeq2`, `EBSeq`, `edgeR`, `limma`, `Poisson-Seq`, `SAM-Seq`.

**Recommendation(s):** Six of the tools that utilize negative binomial or log-normal distributions ([edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html), [baySeq](https://bioconductor.org/packages/release/bioc/html/baySeq.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), and [EBseq](https://bioconductor.org/packages/release/bioc/html/EBSeq.html) control their identification of false positives well.

**Additional links:** The authors released their benchmarking scripts on [Github](https://github.com/bartongroup/KF_arabidopsis-GRNA).

---

**Title:** [How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?](https://www.ncbi.nlm.nih.gov/pubmed/27022035)

**Authors:** Nicholas J. Schurch\*, Pietá Schofield\*, Marek Gierliński\*, Christian Cole\*, Alexander Sherstnev\*, et al. 

**Journal Info:** RNA, March 2016

**Description:** This paper compared 11 differential expression tools on varying numbers of RNA-seq biological replicates (3-42) between two conditions.
Each tool was compared against itself as a standard (using all replicates) and against the other tools.

**Tools/methods compared:** `baySeq`, `cuffdiff`, `DEGSeq`, `DESeq`, `DESeq2`, `EBSeq`, `edgeR (exact and glm modes)`, `limma`, `NOISeq`, `PoissonSeq`, `SAMSeq`. 

**Recommendation(s):** With fewer than 12 biological replicates, [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) were the top performers. 
As replicates increased, [DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html) did a better job minimizing false positives than other tools.

Additionally, the authors recommend at least six biological replicates should be used, rising to at least 12 if users want to identify all significantly differentially expressed genes no matter the fold change magnitude.

---

**Title:** [Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r95)

**Authors:** Franck Rapaport, et al.

**Journal Info:** Genome Biology, September 2013

**Description:** This paper compared six differential expression methods on three cell line data sets from ENCODE (GM12878, H1-hESC, and MCF-7) and two samples from the SEQC study, which had a large fraction of differentially expressed genes validated by qRT-PCR.
Specificity, sensitivity, and false positive rate were the main benchmarking metrics used. 

**Tools/methods compared:** `Cuffdiff`, `edgeR`, `DESeq`, `PoissonSeq`, `baySeq`, `limma`.

**Recommendation(s):** Though no method emerged as favorable in all conditions, those that used negative binomial modeling ([DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [baySeq](https://bioconductor.org/packages/release/bioc/html/baySeq.html)) generally performed best.

The more replicates, the better. Replicate numbers (both biological and technical) have a greater impact on differential detection accuracy than sequencing depth.

### Gene Set Enrichment Analysis

**Title:** [Toward a gold standard for benchmarking gene set enrichment analysis](https://doi.org/10.1093/bib/bbz158)

**Authors:** Ludwig Geistlinger, et al.

**Journal Info:** Briefings in Bioinformatics, February 2020

**Description:** This paper developed a Bioconductor package for reproducible GSEA benchmarking, and used the package to assess 10 widely used enrichment methods with regard to runtime and applicability to RNA-seq data, fraction of enriched gene sets depending on the null hypothesis tested, and recovery of biologically relevant gene sets. The framework can be extended to additional methods, datasets, and benchmark criteria, and should serve as a gold standard for future GSEA benchmarking studies.

**Tools/methods compared:** The paper quantitatively asseses the performance of 10 enrichment methods (`ORA`, `GSEA`, `GSA`, `PADOG`, `SAFE`, `CAMERA`, `ROAST`, `GSVA`, `GLOBALTEST`, `SAMGS`). The paper also compares 10 frequently used enrichment tools implementing these methods (`DAVID`, `ENRICHR`, `CLUSTER-PROFILER`, `GOSTATS`, `WEBGESTALT`, `G:PROFILER`, `GENETRAIL`, `GORILLA`, `TOPPGENE`, `PANTHER`).


**Recommendation(s):**
- `ORA` for the exploratory analysis of simple gene lists, pre-ranked `GSEA` or pre-ranked `CAMERA` for the analysis of pre-ranked gene lists accompanied by gene scores such as fold changes, </br> 
- For enrichment analysis on the full expression matrix (genes x samples), the paper recommends to provide normalized log2 intensities for microarray data and logTPMs (or logRPKMs/logFPKMs) for RNA-seq data; when given raw read counts, the paper recommends to apply a variance-stabilizing transformation such as `voom` to arrive at library-size normalized logCPMs, </br>
- `ROAST` (sample group comparisons) or `GSVA` (single sample) if the question of interest is to test for association of any gene in the set with the phenotype (self-contained null hypothesis), </br>
- `PADOG` (simple experimental designs) or `SAFE` (extended experimental designs) if the question of interest is to test for excess of differential expression in a gene set relative to genes outside the set (competitive null hypothesis).

**Additional links (optional):** http://bioconductor.org/packages/GSEABenchmarkeR


### Differential Splicing

**Title:** [A survey of software for genome-wide discovery of differential splicing in RNA-Seq data](https://humgenomics.biomedcentral.com/articles/10.1186/1479-7364-8-3)

**Authors:** Joan E Hooper

**Journal Info:** Human Genomics, January 2014

**Description:** This paper compares the methodologies, advantages, and disadvantages of eight differential splicing analysis tools, detailing use-cases and features for each.

**Tools/methods compared:** `Cuffdiff2`, `MISO`, `DEXSeq`, `DSGseq`, `MATS`, `DiffSplice`, `Splicing compass`, `AltAnalyze`.

**Recommendation(s):** This is a true breakdown of each tools' advantages and disadvantages. 
The author makes no recommendation due to the performance reliance on experimental setup, data type (e.g. `AltAnalyze` works best on junction + exon microarrays), and user objectives.
Table 1 provides a good comparison of the features and methodology of each method.

**Title:** [A benchmarking of workflows for detecting differential splicing and differential expression at isoform level in human RNA-seq studies](https://academic.oup.com/bib/article/20/2/471/4524048)

---

**Authors:** Gabriela A Merino

**Journal Info:** Briefings in Bioinformatics, March 2019

**Description:** This paper compares nine most commonly used workflows to detect differential isoform expression and splicing.

**Tools/methods compared:** `EBSeq`, `DESeq2`, `NOISeq`, `Limma`, `LimmaDS`, `DEXSeq`, `Cufflinks`, `CufflinksDS`, `SplicingCompass`.

**Recommendation(s):** DESeq2, Limma and NOISeq for differential isoform expression(DIE) analysis and DEXSeq and LimmaDS for differential splicing (DS) testing.

### Transcript Assembly and Quantification

**Title:** [Benchmark analysis of algorithms for determining and quantifying full-length mRNA splice forms from RNA-seq data](https://doi.org/10.1093/bioinformatics/btv488)

**Authors:** Katharina E. Hayer et al.

**Journal Info:** Bioinformatics, Dec 2015

**Description:** This paper compared both guided and de novo transcript reconstruction algorithms using simulated and in vitro transcription (IVT) generated libraries. Precision/recall metrics were obtained by comparing the reconstructed transcripts to their true models.

**Tools/methods compared:** `Cufflinks`, `CLASS`, `FlipFlop`, `IReckon`, `IsoLasso`, `MiTie`, `StringTie`, `StringTie-SR`, `AUGUSTUS`, `Trinity`, `SOAP`, `Trans-ABySS`.

**Recommendation(s):** 
All tools measured produced less than ideal precision-recall (both <90%) when using imperfect simulated or IVT data and genes producing mulitple isoforms. [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks) and [StringTie](https://github.com/gpertea/stringtie) are among the best performers.

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

**Title:** [Evaluation of Nine Somatic Variant Callers for Detection of Somatic Mutations in Exome and Targeted Deep Sequencing Data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0151664)

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

**Description:** This paper compared five germline copy number variation callers against four genetic diagnostics datasets (495 samples, 231 CNVs validated by MLPA) using both default and optimized parameters.
Sensitivity, specificity, positive predictive value, negative predictive value, F1 score, and various correlation coefficients were used as benchmarking metrics.

**Tools/methods compared:** `DECoN`, `CoNVaDING`, `panelcn.MOPS`, `ExomeDepth`, `CODEX2`.

**Recommendation(s):** Most tools performed well, but varied based on datasets. 
The authors felt [DECoN](https://www.imm.ox.ac.uk/research/units-and-centres/mrc-wimm-centre-for-computational-biology/groups/lunter-group/lunter-group/decon-detection-of-exon-copy-number) and [panelcn.MOPS](https://bioconductor.org/packages/release/bioc/html/panelcn.mops.html) with optimized parameters were sensitive enough to be used as screening methods in genetic dianostics.

**Additional links:** The authors have made their benchmarking code ([CNVbenchmarkeR](https://github.com/TranslationalBioinformaticsIGTP/CNVbenchmarkeR)) available, which can be run to determine optimal parameters for each algorithm for a given user's data.

---

**Title:** [An evaluation of copy number variation detection tools for cancer using whole exome sequencing data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5452530/)

**Authors:** Fatima Zare, et al.

**Journal Info:** BMC Bioinformatics, May 2017

**Description:** This paper compared six copy number variation callers on ten TCGA breast cancer tumor-normal pair WES datasets in addition to simulated datasets from VarSimLab.
Sensitivity, specificity, and false-discovery rate were used as the benchmarking metrics.

**Tools/methods compared:** `ADTEx`, `CONTRA`, `cn.MOPS`, `ExomeCNV`, `VarScan2`, `CoNVEX`.

**Recommendation(s):** All tools suffered from high FDRs (\~30-60%), but [ExomeCNV]https://github.com/cran/ExomeCNV) (a now defunct R package) had the highest overall sensitivity. 
[VarScan2](http://dkoboldt.github.io/varscan/) had moderate sensitivity and specificity for both amplifications and deletions.



### SV callers

**Title:** [Comprehensive evaluation and characterisation of short read general-purpose structural variant calling software](https://www.nature.com/articles/s41467-019-11146-4)

**Authors:** Daniel L. Cameron, et al.

**Journal Info:** Nature Communications, July 2019

**Description:** This paper compared 10 structural variant callers on four cell line WGS datasets (NA12878, HG002, CHM1, and CHM13) with orthogonal validation data.
Precision and recall were the benchmarking metrics used.

**Tools/methods compared:** `BreakDancer`, `cortex`, `CREST`, `DELLY`, `GRIDSS`, `Hydra`, `LUMPY`, `manta`, `Pindel`, `SOCRATES`.

**Recommendation(s):** The authors found [GRIDSS](https://github.com/PapenfussLab/gridss) and [manta](https://github.com/Illumina/manta) consistently performed well, but also provide more general guidelines for both users and developers.

 - Use a caller that utilizes multiple sources of evidence and assembly.
 - Use a caller that can call all events you care about.
 - Ensemble calling is not a cure-all and generally don't outperform the best individual callers (at least on these datasets).
 - Do not use callers that rely only on paired-end data.
 - Calls with high read counts are typically artefacts.
 - Simulations aren't real - benchmarking solely on simulations is a bad idea.
 - Developers - be wary of incomplete trust sets and the potential for overfitting. Test tools on multiple datasets.
 - Developers - make your tool easy to use with basic sanity checks to protect against invalid inputs. Use standard file formats.
 - Developers - use all available evidence and produce meaningful quality scores.

---

**Title:** [Comprehensive evaluation of structural variation detection algorithms for whole genome sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5)

**Authors:** Shunichi Kosugi, et al.

**Journal Info:** Genome Biology, June 2019

**Description:** This study compared 69 structural variation callers on simulated and real (NA12878, HG002, and HG00514) datasets.
F-scores, precision, and recall were the main benchmarking metrics. 

**Tools/methods compared:** `1-2-3-SV`, `AS-GENESENG`, `BASIL-ANISE,` `BatVI`, `BICseq2`, `BreakDancer`, `BreakSeek`, `BreakSeq2`, `Breakway`, `CLEVER`, `CNVnator`,
`Control-FREEC`, `CREST`, `DELLY`, `DINUMT`, `ERDS`, `FermiKit`, `forestSV`, `GASVPro`, `GenomeSTRiP`, `GRIDSS`, `HGT-ID`, `Hydra-sv`, `iCopyDAV`, `inGAP-sv`, `ITIS`,
`laSV`, `Lumpy`, `Manta`, `MATCHCLIP`, `Meerkat`, `MELT`, `MELT-numt`, `MetaSV`, `MindTheGap`, `Mobster`, `Mobster-numt`, `Mobster-vei`, `OncoSNP-SEQ`, `Pamir`, `PBHoney`,
`PBHoney-NGM`, `pbsv`, `PennCNV-Seq`, `Pindel`, `PopIns`, `PRISM`, `RAPTR`, `readDepth`, `RetroSeq`, `Sniffles`, `Socrates`, `SoftSearch`, `SoftSV`, `SoloDel`, `Sprites`,
`SvABA`, `SVDetect`, `Svelter`, `SVfinder`, `SVseq2`, `Tangram`, `Tangram-numt`, `Tangram-vei`, `Tea`, `TEMP`, `TIDDIT`, `Ulysses`, `VariationHunter`, `VirusFinder`, `VirusSeq`, `Wham`.

**Recommendation(s):** Varies greatly depending on type and size of the structural variant in addition to read length. 
`GRIDSS`, `Lumpy`, `SVseq2`, `SoftSV`, and `Manta` performed well calling deletions of diverse sizes.
`TIDDIT`, `forestSV`, `ERDS`, and `CNVnator` called large deletions well, while `pbsv`, `Sniffles`, and `PBHoney` were the best performers for small deletions.
For duplications, good choices included `Wham`, `SoftSV`, `MATCHCLIP`, and `GRIDSS`, while `CNVnator`, `ERDS`, and `iCopyDAV` excelled calling large duplications.
For insertions, `MELT`, `Mobster`, `inGAP-sv`, and methods using long read data were most effective.

---

**Title:** [Evaluating nanopore sequencing data processing pipelines for structural variation identification](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1858-1)

**Authors:** Anbo Zhou, et al.

**Journal Info:** Genome Biology, November 2019

**Description:** This paper evaluated four alignment tools and three SV detection tools on four nanopore datasets (both simulated and real).

**Tools/methods compared:** *aligners* - `minimap2`, `NGMLR`, `GraphMap`, `LAST`. SV Callers - `Sniffles`, `NanoSV`, `Picky`.

**Recommendation(s):** The authors recommend using the [minimap2](https://github.com/lh3/minimap2) aligner in combination with the SV caller [Sniffles](https://github.com/fritzsedlazeck/Sniffles) because of their speed and relatively balanced performance.

**Additional links (optional):** The authors provide [all code used in the study](https://github.com/JXing-Lab/nanopore-sv-evaluation) as well as a singularity package containing pre-installed programs and all seven pipeline.

## Single Cell

### scRNA Sequencing Protocols

**Title:** [Benchmarking single-cell RNA-sequencing protocols for cell atlas projects](https://www.nature.com/articles/s41587-020-0469-4)

**Authors:** Elisabetta Mereu\*, Atefeh Lafzi\*, et al.

**Journal Info:** Nature Biotechnology, April 2020

**Description:** This paper evaluated 13 single cell/nuclei RNA-seq protocols to evaluate their aptitude for use in cell atlas-like projects. Using a single cell, multi-species mixture, the authors measured each protocol's ability to capture cell markers, gene detection power, clusterability (with and without integration with other protocols), mappability, and mixability.

**Tools/methods compared:** `Quartz-seq2`, `Chromium`, `Smart-seq2`, `CEL-seq2`, `C1HT-medium`, `C1HT-small`, `ddSEQ`, `Chromium (single nuclei)`, `Drop-seq`, `inDrop`, `ICELL8`, `MARS-seq`, `gmcSCRB-seq`.

**Recommendation(s):** See figure 6 for a summary of benchmarking results for each method. `Quartz-seq2` was the overall best performing, yielding superior results for gene detection and marker expression over other methods, though `Chromium`, `Smart-seq2`, and `CEL-seq2` were also strong performers.

**Additional links:** The authors provide benchmarking code and analysis code in two different Github repositories - [here](https://github.com/ati-lz/HCA_Benchmarking) and [here](https://github.com/elimereu/matchSCore2).

---

**Title:** [Systematic comparison of single-cell and single-nucleus RNA-sequencing methods](https://www.nature.com/articles/s41587-020-0465-8)

**Authors:** Jiarui Ding, et al.

**Journal Info:** Nature Biotechnology, April 2020

**Description:** This study evaluated seven methods for single-cell and/or single-nucleus RNA-sequencing on three types of samples: cell lines, PBMCs, and brain tissue. Evaluation metrics included the structure and alignment of reads, number of multiplets and detection sensitivity, and ability to recover known biological information.

**Tools/methods compared:** `Smart-seq2`, `CEL-Seq2`, `3' 10X Chromium`, `Drop-Seq`, `Seq-Well`, `inDrops`, `sci-RNA-seq`.

**Recommendation(s):** Overall, the authors found `3' 10X Chromium` to have the strongest consistent performance among the high-throughput methods, yielding the highest sensitivity, though it did not perform any better for cell type classification. When greater sensitivity is required, the authors recommend `Smart-seq2` or `CEL-Seq2`, which both performed similarly. Supplementary table 7 includes an overview of each method's relative merits.

**Additional links:** The authors made their unified analysis pipeline (scumi) available as a [python package](https://bitbucket.org/jerry00/scumi-dev/src/master/), the repo of which also includes their R scripts used for cell filtering and cell type assignment.

### scRNA Analysis Pipelines

**Title:** [A systematic evaluation of single cell RNA-seq analysis pipelines](https://www.nature.com/articles/s41467-019-12266-7)

**Authors:** Beate Veith, et al.

**Journal Info:** Nature Communications, October 2019

**Description:** This study evaluated \~3000 pipeline combinations based on three mapping, three annotation, four imputation, seven normalization, and four differential expression testing approaches with five scRNA-seq library protocols on simulated data.

**Tools/methods compared:** scRNA-seq library prep protocols - `SCRB-seq`, `Smart-seq2`, `CEL-seq2`, `Drop-seq`, `10X Genomics`. Mapping - `bwa`, `STAR`, `kallisto`. Annotation - `gencode`, `refseq`, `vega`. Imputation - `filtering`, `DrImpute`, `scone`, `SAVER`. Normalization - `scran`, `SCnorm`, `Linnorm`, `Census`, `MR`, `TMM`. Differential testing - `edgeR-zingeR`, `limma`, `MAST`, `T-test`.  

**Recommendation(s):** **Figure 5F** contains a flowchart with the authors' recommendations. For alignment, `STAR` with Gencode annotations generally had the highest mapping and assignment rates. All mappers performed best with Gencode annotations. For normalization, `scran` was found to best handle potential assymetric differential expression and large numbers of differentially expressed genes. They also note that normalization is overall the most influential step, particularly if asymmetric DE is present (**Figure 5**). For Smart-seq2 data without spike-ins, the authors suggest `Census` may be the best choice. The authors found little benefit to imputation in most scenarios, particularly if one of the better normalization methods (e.g. `scran`) was used. The authors found library prep and normalization strategies to have a stronger effect on pipeline performance than the choice of differential expression tool, but generally found `limma-trend` to have the most robust performance. 

**Additional links:** The authors made their simulation tool (powsimR) available on [Github](https://github.com/bvieth/powsimR) along with their [pipeline scripts](https://github.com/bvieth/scRNA-seq-pipelines) to reproduce their analyses.

---

**Title:** [Comparison of high-throughput single-cell RNA sequencing data processing pipelines](https://doi.org/10.1093/bib/bbaa116)

**Authors:** Mingxuan Gao, et al.

**Journal Info:** Briefings in Bioinformatics, July 2020

**Description:** This study evaluated 7 scRNA-seq pipelines on 8 data sets.

**Tools/methods compared:** `Drop-seq-tools version-2.3.0`, `Cell Ranger version-3.0.2`, `scPipe version-1.4.1`, `zUMIs version-2.4.5b`, `UMI-tools version-1.0.0`, `umis version-1.0.3`, `dropEst version-0.8.6`

**Recommendation(s):** Cell Ranger shows the highest algorithm complexity and parallelization, whereas scPipe, umis and zUMIs have lower complexity that is suitable for large-scale scRNA-seq integration analysis. UMI-tools show the highest transcript quantification accuracy on ERCC datasets from three scRNA-seq platforms. Integration of expression matrices from different pipelines will introduce confounding factors akin to batch effect. zUMIs and dropEst have higher sensitivity to detectmore genes for single cells, which may also bring unwanted factors. For most downstream analysis, Drop-seq-tools, Cell Ranger and UMI-tools show high consistency, whereas umis and zUMIs show inconsistent results compared with the other pipelines.

### scRNA Imputation Methods

**Title:** [A Systematic Evaluation of Single-cell RNA-sequencing Imputation Methods](https://www.biorxiv.org/content/10.1101/2020.01.29.925974v1.full)

**Authors:** Wenpin Hou, et al.

**Journal Info:** bioRxiv, January 2020

**Description:** This paper evaluated 18 scRNA-seq imputation methods using seven datasets containing cell line and tissue data from several experimental protocols. The authors assessed the similarity of imputed cell profiles to bulk samples and investigated whether imputation improves signal recovery or introduces noise in three downstream applications - differential expression, unsupervised clustering, and trajectory inference.

**Tools/methods compared:** `scVI`, `DCA`, `MAGIC`, `scImpute`, `kNN-smoothing`, `mcImpute`, `SAUCIE`, `DrImpute`, `PBLR`, `SAVER`, `VIPER`, `SAVERX`, `DeepImpute`, `scRecover`, `ALRA`, `bayNorm`, `AutoImpute`, `scScope`.

**Recommendation(s):** Figure 6 provides a performance summary of the tested methods. In general, the authors recommend caution using any of these methods, as they can introduce significant variability and noise into downstream analyses. Of the methods tested, `MAGIC`, `kNN-smoothing`, and `SAVER` outperformed the other methods most consistently, though this varied widely across evaluation criteria, protocols, datasets, and downstream analysis. Many methods show no clear improvement over no imputation, and in some cases, perform significantly worse in downstream analyses.

**Additional links (optional):** The authors placed all of their benchmaking code on [Github](https://github.com/Winnie09/imputationBenchmark).

### scRNA Differential Gene Expression

**Title:** [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612)

**Authors:** Charlotte Soneson\* & Mark D Robinson\*

**Journal Info:** Nature Methods, February 2018

**Description:** This paper evaluated 36 approaches for determining differential gene expression from both synthetic and 36 real scRNA-seq datasets. The authors assess type I error control, FDR control and power, computational efficiency, and consistency.

**Tools/methods compared:** `edgeRQLFDetRate`, `MASTcpmDetRate`, `limmatrend`, `MASTtpmDetRate`, `edgeRQLF`, `ttest`, `voomlimma`, `Wilcoxon`, `MASTcpm`, `MASTtpm`, `SAMseq`, `D3E`, `edgeRLRT`, `metagenomeSeq`, `edgeRLRTcensus`, `edgeRLRTdeconv`, `monoclecensus`, `ROTStpm`, `ROTSvoom`, `DESeq2betapFALSE`, `edgeRLRTrobust`, `monoclecount`, `DESeq2`, `DESeq2nofilt`, `ROTScpm`, `SeuratTobit`, `NODES`, `DESeq2census`, `scDD`, `BPSC`, `SCDE`, `DEsingle`, `monocle`, `SeuratBimodnofilt`, `SeuratBimodlsExpr2`, `SeuratBimod`.

**Recommendation(s):** In general, the authors found that gene prefiltering was essential for good, robust performance from many methods. They note high variability between methods and summarize general performance across all metrics in Figure 5. They do not make recommendations as to a specific method/tool. Of note is that Seurat switched to using the wilcoxon test by default after this study was released, as it performed much better than their previously available methods.

**Additional links:** The authors make their benchmarking pipeline, [conquer](https://github.com/markrobinsonuzh/conquer), available on Github. Their processed data and associated reports have also been [made available](http://imlspenticton.uzh.ch:3838/conquer/) for additional comparisons.

### Trajectory Inference

**Title:** [A comparison of single-cell trajectory inference methods](https://doi.org/10.1038/s41587-019-0071-9)

**Authors:** Wouter Saelens\*, Robrecht Cannoodt\*, et al.

**Journal Info:** Nat Biotech, April 2019

**Description:** A comprehensive evaluation of 45 trajectory inference methods, this paper provides an unmatched comparison of the rapidly evolving field of single-cell trajectory inference. Each method was scored on accuracy, scalability, stability, and usability. Should be considered a gold-standard for other benchmarking studies.

**Tools/methods compared:** `PAGA`, `RaceID/StemID`, `SLICER`, `Slingshot`, `PAGA Tree`, `MST`, `pCreode`, `SCUBA`, `Monocle DDRTree`, `Monocle ICA`, `cellTree maptpx`, `SLICE`, `cellTree VEM`, `EIPiGraph`, `Sincell`, `URD`, `CellTrails`, `Mpath`, `CellRouter`, `STEMNET`, `FateID`, `MFA`, `GPfates`, `DPT`, `Wishbone`, `SCORPIUS`, `Component 1`, `Embeddr`, `MATCHER`, `TSCAN`, `Wanderlust`, `PhenoPath`, `topslam`, `Waterfall`, `EIPiGraph linear`, `ouijaflow`, `FORKS`, `Angle`, `EIPiGraph cycle`, `reCAT`.

**Recommendation(s):** Varies depending on dataset and expected trajectory type, though [PAGA, PAGA Tree](https://icb-scanpy.readthedocs-hosted.com/en/stable/tutorials.html#trajectory-inference), [SCORPIUS](https://github.com/rcannood/SCORPIUS), and [Slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html) all scored highly across all metrics. 

Authors wrote an [interactive Shiny app](https://dynverse.org/users/3-user-guide/2-guidelines/) to help users choose the best methods for their data.

**Additional links:** The [dynverse site](https://dynverse.org/) contains numerous packages for users to run and compare results from different trajectory methods on their own data without installing each individually by using Docker. Additionally, they provide several tools for developers to wrap and benchmark their own method against those included in the study. 

### Gene Regulatory Network Inference

**Title:** [Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data](https://www.nature.com/articles/s41592-019-0690-6)

**Authors:** Aditya Pratapa, et al.

**Journal Info:** Nature Methods, January 2020

**Description:** The authors compared 12 gene regulatory network (GRN) inference techniques to assess the accuracy, robustness, and efficiency of each method on simulated data from synthetic networks, simulated data from curated models, and real scRNA-seq datasets.

**Tools/methods compared:** `GENIE3`, `PPCOR`, `LEAP`, `SCODE`, `PIDC`, `SINCERITIES`, `SCNS`, `GRNVBEM`, `SCRIBE`, `GRNBoost2`, `GRISLI`, `SINGE`.

**Recommendation(s):** In general, the authors recommend `PIDC`, `GENIE3`, or `GRNBoost2`, as they were the leading and consistent performers for both curated models and experimental datasets in terms of accuracy as well as having multithreaded implementations. Both the `GENIE3` and `GRNBoost2` methods can be used from the [SCENIC](https://aertslab.org/#scenic) workflows in either R or the (much faster) python implementation.

**Additional links:** The authors provide their benchmarking framework, [BEELINE](https://github.com/murali-group/BEELINE) on Github, which also provides an easy-to-use and uniform interface to each method in the form of a Docker image.

---

**Title:** [Evaluating methods of inferring gene regulatory networks highlights their lack of performance for single cell gene expression data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2217-z)

**Authors:** Shounan Chen\* & Jessica C. Mar\*

**Journal Info:** BMC Bioinformatics, June 2018

**Description:** This study compared 8 gene regulatory network inference methods (5 bulk RNA-seq, 3 specific to scRNA-seq) for scRNA-seq data for precision and recall.

**Tools/methods compared:** `Partial correlation (Pcorr)`, `Bayesian networks`, `GENIE3`, `ARACNE`, `CLR`, `SCENIC`, `SCODE`, `PIDC`.

**Recommendation(s):** Generally, the authors found relatively poor performance across all methods both for simulated and real data. The results from each method had few similarities with other methods and high false positive rates, and the authors recommend caution when interpreting the networks reconstructed with these methods. The authors showed that many of these methods were dramatically affected by dropout events.

### Integration/Batch Correction

**Title:** [A benchmark of batch-effect correction methods for single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9)

**Authors:** Hoa Thi Nhu Tran et al.

**Journal Info:** Genome Biology, January 2020

**Description:** The authors compared 14 methods in terms of computational runtime, the ability to handle large datasets, and batch-effect correction efficacy while preserving cell type purity.

**Tools/methods compared:**
`Seurat2`, `Seurat3`, `Harmony`, `fastMNN`, `MNN Correct`, `ComBat`, `Limma`, `scGen`, `Scanorama`, `MMD-ResNet`, `ZINB-WaVe`, `scMerge`, `LIGER`, `BBKNN`

**Recommendation(s):** Based on the benchmarking results authors suggest [Harmony](https://github.com/immunogenomics/harmony), [LIGER](https://github.com/MacoskoLab/liger), and [Seurat3](https://satijalab.org/seurat/) as best methods for batch integration. 

### Dimensionality Reduction

**Title:** [Accuracy, Robustness and Scalability of Dimensionality Reduction Methods for Single Cell RNAseq Analysis](https://www.biorxiv.org/content/10.1101/641142v2)

**Authors:** Shiquan Sun, et al.

**Journal Info:** BioRxiv, October 2019

**Description:** A mammoth comparison of 18 different dimension reduction methods on 30 publicly available scRNAseq data sets in addition to 2 simulated datasets for a variety of purposes ranging from cell clustering to trajectory inference to neighborhood preservation. 

**Tools/methods compared:**
`factor analysis (FA)`, `principal component analysis (PCA)`, `independent component analysis (ICA)`, `Diffusion Map`, `nonnegative matrix factorization (NMF)`, `Poisson NMF`, `zero-inflated factor analysis (ZIFA)`, `zero-inflated negative binomial based wanted variation extraction (ZINB-WaVE)`, `probabilistic count matrix factorization (pCMF)`, `deep count autoencoder network (DCA)`, `scScope`, `generalized linear model principal component analysis (GLMPCA)`, `multidimensional scaling (MDS)`, `locally linear embedding (LLE)`, `local tangent space alignment (LTSA)`, `Isomap`, `uniform manifold approximation and projection (UMAP)`, `t-distributed stochastic neighbor embedding (tSNE)`.

**Recommendation(s):** Varies depending on use case. 
[Factor Analysis](https://cran.r-project.org/web/packages/psych/index.html) and [principal component analysis](https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/prcomp) performed well for most use cases. 
See figure 5 for pratical guidelines.

**Additional links:** The authors have made their benchmarking code available [on Github](https://github.com/xzhoulab/DRComparison).

---

**Title:** [Benchmarking principal component analysis for large-scale single-cell RNA-sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1900-3)

**Authors:** Koki Tsuyuzaki, et al.

**Journal Info:** Genome Biology, January 2020

**Description:** This study compared 21 implementations of 10 algorithms across Python, R, and Julia for principal component analysis for scRNA-seq data, measuring scalability, computational efficiency, outlier robustness, t-SNE/UMAP replication, ease of use, and more using both synthetic and real datasets.

**Tools/methods compared:** `PCA (sklearn, full)`, `fit (MultiVariateStats.jl)`, `Downsampling`, `IncrementalPCA (sklearn)`, `irlba (irlba)`, `svds (RSpectra)`, `propack.svd (svd)`, `PCA (sklearn, arpack)`, `irlb (Cell Ranger)`, `svds (Arpack.jl)`, `orthiter (OnlinePCA.jl)`, `gd (OnlinePCA.jl)`, `sgd (OnlinePCA.jl)`, `rsvd (rsvd)`, `oocPCA_CSV (oocRPCA)`, `PCA (sklearn, randomized)`, `randomized_svd (sklearn)`, `PCA (dask-ml)`, `halko (OnlinePCA.jl)`, `algorithm971 (OnlinePCA.jl)`.

**Recommendation(s):** Author recommendations vary based on the language being used and matrix size. See figure 8 for recommendations along with recommended parameter settings.

**Additional links:** The authors published their benchmarking scripts [on Github](https://github.com/rikenbit/onlinePCA-experiments).

### Cell Annotation/Inference

**Title:** [A comparison of automatic cell identification methods for single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z)

**Authors:** Tamim Abdelaal\*, Lieke Michielsen\*, et al.

**Journal Info:** Genome Biology, September 2019

**Description:** The authors benchmarked 22 classification methods that automatically assign cell identities including *single-cell-specific* and *general-purpose* classifiers across 27 publicly available single-cell RNA sequencing datasets of different sizes, technologies, species, and levels of complexity. Two types of experimental setups were used evaluate the performance of each method for within dataset predictions (intra-dataset) and across datasets (inter-dataset) based on accuracy, percentage of unclassified cells, and computation time. 

**Tools/methods compared:** `Garnett`, `Moana`, `DigitalCellSorter`, `SCINA`, `scVI`, `Cell-BLAST`, `ACTINN`, `LAmbDA`, `scmapcluster`, `scmapcell`, `scPred`, `CHETAH`, `CaSTLe`, `SingleR`, `scID`, `singleCellNet`, `LDA`, `NMC`, `RF`, `SVM`, `SVM<sub>rejection</sub>`, `kNN`

**Recommendation(s):** All classifiers performed well. The authors recommended SVM<sub>rejection</sub> classifier (with a linear kernel). Other classifiers include SVM , singleCellNet, scmapcell, and scPred were also of high performances.

 In their experiments, incorporating prior knowledge in the form of marker genes does not improve the performance.

**Additional links:** A Snakemake workflow, [scRNAseq_Benchmark](https://github.com/tabdelaal/scRNAseq_Benchmark/), was provided to automate the benchmarking analyses.

### Variant Calling

**Title:** [Systematic comparative analysis of single-nucleotide variant detection methods from single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4)

**Authors:** Fenglin Liu\*, Yuanyuan Zhang\*, et al.

**Journal Info:** Genome Biology, November 2019

**Description:** This paper compared seven variant callers using both simulation and real scRNA-seq datasets and identified several elements influencing their performance, including read depth, variant allele frequency, and specific genomic contexts. Sensitivity and specificity were the benchmarking metrics used.

**Tools/methods compared:** `SAMtools`, `GATK`, `CTAT`, `FreeBayes`, `MuTect2`, `Strelka2`, `VarScan2`.

**Recommendation(s):** Varies, see figure 7 for a flowchart breakdown. Generally, [SAMtools](http://www.htslib.org/) (most sensitive, lower specificity in intronic or high-identity regions), [Strelka2](https://github.com/Illumina/strelka) (good performance when read depth >5), [FreeBayes](https://github.com/ekg/freebayes) (good specificity/sensitivity in cases with high variant allele frequencies), and CTAT (no alignment step necessary) were top performers.

**Additional links:** The authors made their benchmarking code available [on Github](https://github.com/fenglin0/benchmarking_variant_callers).

### ATAC-seq

**Title:** [Assessment of computational methods for the analysis of single-cell ATAC-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5)

**Authors:** Caleb Lareau\*, Tommaso Andreani\*, Micheal E. Vinyard\*, et al.

**Journal Info:** Genome Biology, November 2019

**Description:** This study compares 10 methods for scATAC-seq processing and featurizing using 13 synthetic and real datasets from diverse tissues and organisms.

**Tools/methods compared:** `BROCKMAN`, `chromVAR`, `cisTopic`, `Cicero`, `Gene Scoring`, `Cusanovich2018`, `scABC`, `Scasat`, `SCRAT`, `SnapATAC`.

**Recommendation(s):** [SnapATAC](https://github.com/r3fang/SnapATAC), [Cusanovich2018](https://www.ncbi.nlm.nih.gov/pubmed/30078704), and [cisTopic](https://github.com/aertslab/cisTopic) were the top performers for separating cell populations of different coverages and noise levels. SnapATAC was the only method capable of analyzing a large dataset (>80k cells).

**Additional links:** The authors have made their benchmarking code available [on Github](https://github.com/pinellolab/scATAC-benchmarking/).

## Statistics 

### False Discovery Rates 

**Title:** [A practical guide to methods controlling false discoveries in computational biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1716-1)

**Authors:** Keegan Korthauer\*, Patrick K. Kimes\*, et al.

**Journal Info:** Genome Biology, June 2019

**Description:** An benchmark comparison of the accuracy, applicability, and ease of use of two classic and six modern methods that control for the false discovery rate (FDR). Used simulation studies as well as six case studies in computational biology (specifically differential expression testing in bulk RNA-seq, differential expression testing in single-cell RNA-seq, differential abundance testing and correlation analysis in 16S microbiome data, differential binding testing in ChIP-seq, genome-wide association testing, and gene set analysis). 

**Tools/methods compared:** Benjamini-Hochberg, Storey’s q-value, conditional local FDR (LFDR), FDR regression (FDRreg), independent hypothesis weighting (IHW), adaptive shrinkage (ASH), Boca and Leek’s FDR regression (BL), and adaptive p-value thresholding (AdaPT).  

**Recommendation(s):** Modern FDR methods that use an informative covariate (as opposed to only _p_-values) leads to more power while controlling the FDR over classic methods. The improvement of the modern FDR methods over the classic methods increases with the informativeness of the covariate, total number of hypothesis tests, and proportion of truly non-null hypotheses. 

**Additional links:** Full analyses of the in silico experiments, simulations, and case studies are provided in Additional files 2–41 at https://pkimes.github.io/benchmark-fdr-html/. The source code to reproduce all results in the manuscript and additional files, as well as all figures, is [available on GitHub](https://github.com/pkimes/benchmark-fdr). An `ExperimentHub` package containing the full set of results objects is available through the Bioconductor project, and a Shiny application for interactive exploration of these results is also [available on GitHub](https://github.com/kdkorthauer/benchmarkfdr-shiny). The source code, ExperimentHub package, and Shiny application are all made available under the MIT license.

## Microbiome 

### Diversity analysis

**Title:** [Evaluating Bioinformatic Pipeline Performance for Forensic Microbiome Analysis](https://doi.org/10.1111/1556-4029.14213)

**Authors:** Sierra F. Kaszubinski\*, Jennifer L. Pechal\*, et al.

**Journal Info:** Journal of Forensic Sciences, 2019

**Description:** Sequence reads from postmortem microbiome samples were analyzed with [mothur v1.39.5](https://mothur.org), [QIIME2 v2018.11](https://qiime.org), and [MG-RAST v4.0.3](https://mg-rast.org). For postmortem data, MG-RAST had a much smaller effect size than mothur and QIIME2 due to the twofold reduction in samples. QIIME2 and Mothur returned similar results, with Mothur showing inflated richness due to unclassified taxa. Adjusting minimum library size had significant effects on microbial community structure, sample size less so except for low abundant taxa.

**Tools/methods compared:** `mothur`, `QIIME2`, `MG-RAST`

**Recommendation(s):** QIIME2 was deemed the most appropriate choice for forensic analysis in this study.

**Additional links:** Sequence data are archived through the European Bioinformatics Institute European Nucleotide Archive (www.ebi.ac.uk/ena) under accession number: PRJEB22642. Pipeline parameters and microbial community analyses are available on GitHub (https://github.com/sierrakasz/postmortem-analysis).

# Contributors

 - Jared Andrews ([@j-andrews7](https://github.com/j-andrews7/))
 - Kevin Blighe ([@kevinblighe](https://github.com/kevinblighe), [biostars](https://www.biostars.org/u/41557/))
 - Ludwig Geistlinger ([@lgeistlinger](https://github.com/lgeistlinger))
 - Jeremy Leipzig ([@leipzig](https://github.com/leipzig))
 - Avi Srivastava ([@k3yavi](https://github.com/k3yavi))
 - Stephanie Hicks ([@stephaniehicks](https://github.com/stephaniehicks))
 - Sridhar N Srivatsan ([@sridhar0605](https://github.com/sridhar0605))
 - Qingzhou Zhang ([@zqzneptune](https://github.com/zqzneptune))

