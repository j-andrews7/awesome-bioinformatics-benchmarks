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
Please follow the format below when adding papers.

# Tool/Method Sections
Additional sections/sub-sections can be added as needed.


## ChIP-seq

### Peak Callers

## RNA-seq

### Differential Gene Expression

### Cell-Type Deconvolution

**Title:** [Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology](https://academic.oup.com/bioinformatics/article/35/14/i436/5529146)

**Authors:** Markus List\*, Tatsiana Aneichyk\*, et al.

**Journal Info:** Bioinformatics, July 2019

**Description:** This paper benchmarks and compares seven methods for computational deconvolution of cell-type abundance in bulk RNA-seq samples. Each method was tested on both simulated and true bulk RNA-seq samples validated by FACS.

**Tools/methods compared:** `quanTIseq`, `TIMER`, `CIBERSORT`, `CIBERSORT abs. mode`, `MCPCounter`, `xCell`, `EPIC`.

**Recommendation(s):** Varies. In general, the authors recommend [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) and [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) due to their overall robustness and absolute (rather than relative) scoring, though [xCell](http://xcell.ucsf.edu/) is recommended for binary presence/absence of cell types and [MCPcounter](https://github.com/ebecht/MCPcounter) was their recommended relative scoring method.

**Additional links:** The authors created an [R package called immunedeconv](https://github.com/icbi-lab/immunedeconv) for easy installation of use of all these methods. For developers, they have made available their [benchmarking pipeline](https://github.com/icbi-lab/immune_deconvolution_benchmark) so that others can reproduce/extend it to test their own tools/methods.

## RNA/cDNA Microarrays

## Variant Callers

### SNP/Indel Callers

### CNV Callers

## Single Cell

### Trajectory Inference

**Title:** [A comparison of single-cell trajectory inference methods](https://doi.org/10.1038/s41587-019-0071-9)

**Authors:** Wouter Saelens\*, Robrecht Cannoodt\*, et al.

**Journal Info:** Nat Biotech, April 2019

**Description:** A comprehensive evaluation of 45 trajectory inference methods, this paper provides an unmatched comparison of the rapidly evolving field of single-cell trajectory inference. Each method was scored on accuracy, scalability, stability, and usability. Should be considered a gold-standard for other benchmarking studies.

**Tools/methods compared:** `PAGA`, `RaceID/StemID`, `SLICER`, `Slingshot`, `PAGA Tree`, `MST`, `pCreode`, `SCUBA`, `Monocle DDRTree`, `Monocle ICA`, `cellTree maptpx`, `SLICE`, `cellTree VEM`, `EIPiGraph`, `Sincell`, `URD`, `CellTrails`, `Mpath`, `CellRouter`, `STEMNET`, `FateID`, `MFA`, `GPfates`, `DPT`, `Wishbone`, `SCORPIUS`, `Component 1`, `Embeddr`, `MATCHER`, `TSCAN`, `Wanderlust`, `PhenoPath`, `topslam`, `Waterfall`, `EIPiGraph linear`, `ouijaflow`, `FORKS`, `Angle`, `EIPiGraph cycle`, `reCAT`.

**Recommendation(s):** Varies depending on dataset and expected trajectory type, though [PAGA, PAGA Tree](https://scanpy.readthedocs.io/en/latest/examples.html#trajectory-inference), [SCORPIUS](https://github.com/rcannood/SCORPIUS), and [Slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html) all scored highly across all metrics. Authors wrote an [interactive Shiny app](https://dynverse.org/users/3-user-guide/2-guidelines/) to help users choose the best methods for their data.

**Additional links:** The [dynverse site](https://dynverse.org/) contains numerous packages for users to run and compare results from different trajectory methods on their own data without installing each individually by using Docker. Additionally, they provide several tools for developers to wrap and benchmark their own method against those included in the study. 

---


# Contributors

 - Jared Andrews ([@j-andrews7](https://github.com/j-andrews7/))

