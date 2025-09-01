# What is {{sc}}?

{{sc}} is a Snakemake pipeline designed to automate the laborious steps of processing single-cell and single-nuclei RNA sequencing data. This comprehensive solution effectively handles data generated using 10x Genomics technology. Starting from raw sequencing files (as well as a metadata file and a config file), {{sc}} performs a series of standard analysis steps, delivering outputs suitable for a variety of downstream analyses. The pipeline supports analyses across multiple samples, enabling researchers to gain a complete understanding of their datasets while ensuring reproducibility and scalability in their workflows.

---

## Steps

### Overview
{{sc}} consists of a series of core steps which can be performed in a single execution of the workflow. Additional optional steps are available to extend the core analyses as needed. The diagram below outlines all steps, with detailed descriptions provided in the following sections.

![workflow_chart](assets/images/scprocess_workflow_diagram_white_bg.png#only-light)
![workflow_chart](assets/images/scprocess_workflow_diagram_black_bg.png#only-dark)

---
<div class="img-caption"> The figure illustrates all steps in {{sc}}. Specific software packages used are listed for individual steps. Several steps generate HTML reports with diagnostic plots, enabling users to inspect the results at key points in the workflow. </div>

### Core pipeline steps

* #### Read alignment and quantification
    
    {{sc}} starts by mapping sequencing reads to the reference genome and performing transcript quantification using `simpleaf`[^1]. This open-source alternative to [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is designed for speed and memory efficiency, while also offering the advantage of reporting counts for spliced and unspliced transcripts separately.

* #### Ambient RNA removal (optional)

    In droplet-based assays, cells or nuclei are encapsulated in droplets, but some freely-floating RNA can also be captured, which is referred to as ambient RNA. Ambient RNA contamination is particularly common in single-nucleus assays due to residual cytoplasmic material and harsh isolation protocols that can cause nuclei to rupture. Since ambient RNA can interfere with results of downstream analysis, it is beneficial to remove this contamination *in silico*. In {{sc}}, users can select from `CellBender`[^2] and `DecontX`[^3] to remove ambient RNA. 

* #### Cell calling

    Single cell data is based on droplets with unique barcodes attached to them. This means that an important step is distinguishing barcodes/droplets that contain cells/nuclei from those that are empty. {{sc}} employs several options for detecting barcodes that correspond to cell-containing droplets. When `CellBender` is used for ambient RNA removal, its built-in method generates a filtered counts matrix. If `DecontX` is selected or ambient RNA removal is skipped, users can choose between `barcodeRanks` and `emptyDrops` methods from the `DropletUtils`[^4] R package. `barcodeRanks` separates cell-containing and empty droplet populations by detecting key transition points on the barcode-rank curve while `emptyDrops` tests each barcode for significant deviations from the ambient profile.

* #### Doublet detection
    
    Doublets, formed when two cells are captured in the same droplet, can distort single-cell RNA-seq results. {{sc}} addresses this by using `scDblFinder`[^5] for doublet detection. Additionally, {{sc}} integrates `scDblFinder`-flagged doublets with cells that pass the QC filtering step. This allows users to detect doublet-enriched clusters which can be removed from downstream analysis. 

* #### QC filtering

    In addition to removing doublets, {{sc}} filters out cells with low library size, low feature counts, high mitochondrial read proportions, and high spliced read proportions using user-defined thresholds. The spliced proportion is a particularly informative metric in single-nuclei RNA-seq, as elevated levels may indicate residual cytoplasmic material (see Montserrat-Ayuso and Esteve-Codina[^6]).

* #### Generating pseudobulks from cells and empty droplets and ambient gene detection

    {{sc}} aggregates counts from both cell-containing and empty droplets into pseudobulk profiles to identify genes enriched in empty droplets—i.e., ambient genes. These ambient genes likely represent residual contamination rather than true biological signals. Identifying them supports cleaner and more accurate downstream analyses and helps prioritize genes less affected by contamination.

* #### Highly variable gene detection

    Highly variable gene detection in {{sc}} is performed using the `Seurat VST`[^7] method in a chunk-wise or sample-wise manner, enabling the computation of ranking metrics for all genes without the need to load the entire dataset into memory. The efficient generation of a reduced matrix containing only highly variable genes ensures optimal performance in downstream analyses and facilitates the processing of larger datasets.

* #### Integration
    
    After identifying highly variable genes, {{sc}} proceeds with dimentionality reduction and clustering. In multi-sample analysis, various factors can introduce batch effects that obscure true biological signals. To mitigate this, {{sc}} offers the option to compute batch-corrected PCA embeddings using `Harmony`[^8]. This ensures that clustering and downstream analyses reflect true biological relationships rather than technical variation.

* #### Marker gene identification
    
    Assigning meaningful labels to clusters in single-cell data is essential for interpretation of single cell data. This is commonly achieved by examining marker genes for each cluster, identified by comparing the expression profile of each cluster against all others. In {{sc}}, transcript counts are aggregated per cluster within each sample to generate "pseudobulk" values, which are then compared using `edgeR`[^9]. This approach avoids the assumption that individual cells from the same sample are independent, thereby enhancing the statistical reliability of the results. {{sc}} also performs gene set enrichment analysis on all marker genes and includes visualizations of user-defined gene sets in the HTML report.


### Optional steps

* #### Processing multiplexed samples
    
    Multiplexing strategies are commonly used to scale up single-cell experiments by enabling the analysis of multiple samples in a single run. Common approaches include labeling cells in individual samples with lipid-bound or antibody-conjugated oligonucleotides (hashtag oligos, or HTOs) prior to pooling. Alternatively, sample labels can be derived based on differences in genetic backgrounds. {{sc}} supports the analysis of multiplexed samples by quantifying HTOs and demultiplexing samples using the `Seurat HTODemux`[^10] function. It also accommodates outputs from external demultiplexing tools, enabling seamless processing of samples regardless of the multiplexing strategy employed.

* #### Cell type labelling

    {{sc}} provides automated cell type annotation of human and mouse brain datasets using `XGBoost` classifiers Classifiers are trained on the following datasets:

    - Human Brain Classifier: [Transcriptomic diversity of cell types across the adult human brain](https://www.science.org/doi/10.1126/science.add7046) [has to be updated]

    - Mouse Brain Classifier: [The molecular cytoarchitecture of the adult mouse brain](https://www.nature.com/articles/s41586-023-06818-7) and [A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain](https://www.nature.com/articles/s41586-023-06812-z) [work in progress]

* #### Subclustering

    {{sc}} offers a subclustering feature that enables users to perform a second round of analysis on a specific subset of cells and includes the following steps: generating pseudobulks from selected cells and detecting ambient genes, indentifying higly variable genes, data integration and marker gene identification. Cell subsets can be defined based on user-provided cell type labels, clusters identified during the primary round of {{sc}}, or cell type labels assigned by the XGBoost classifier. This functionality is particularly valuable when a primary cluster or cell type contains diverse cell states, developmental stages, or activation states that warrant more detailed exploration.

<!-- citations -->
<!-- [Link to paper](). [Link to package](). -->

[^1]: He D, Patro R. _Simpleaf: A simple, flexible, and scalable framework for single-cell data processing using alevin-fry_. Bioinformatics. 2023;39:btad614. [Link to paper](https://academic.oup.com/bioinformatics/article/39/10/btad614/7295550). [Link to package](https://simpleaf.readthedocs.io/en/latest/).

[^2]: Fleming SJ, Chaffin MD, Arduini A, Akkad A-D, Banks E, Marioni JC, et al. _Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender_. Nat Methods. 2023;20:1323–35. [Link to paper](https://www.nature.com/articles/s41592-023-01943-7). [Link to package](https://cellbender.readthedocs.io/en/latest/).

[^3]: Yang S, Corbett SE, Koga Y, Wang Z, Johnson WE, Yajima M, et al. _Decontamination of ambient RNA in single-cell RNA-seq with DecontX_. Genome Biol. 2020;21:57. [Link to paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6). [Link to package](https://bioconductor.org/packages/release/bioc/html/decontX.html).

[^4]: Lun ATL, Riesenfeld S, Andrews T, Dao TP, Gomes T, participants in the 1st Human Cell Atlas Jamboree, et al. _EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data_. Genome Biol. 2019;20:63. [Link to paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y). [Link to package](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html).

[^5]: Germain P-L, Lun A, Macnair W, Robinson MD. _Doublet identification in single-cell sequencing data using scDblFinder._ F1000Res. 2021;10:979. [Link to paper](https://f1000research.com/articles/10-979/v2). [Link to package](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html).

[^6]: Montserrat-Ayuso T, Esteve-Codina A. _High content of nuclei-free low-quality cells in reference single-cell atlases: a call for more stringent quality control using nuclear fraction_. BMC Genomics. 2024;25:1124. [Link to paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-11015-5).

[^7]: Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A, et al. _Integrated analysis of multimodal single-cell data_. Cell. 2021;184:3573-3587.e29. [Link to paper](https://www.sciencedirect.com/science/article/pii/S0092867421005833). [Link to package](https://satijalab.org/seurat/reference/findvariablefeatures).

[^8]: Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, et al. _Fast, sensitive and accurate integration of single-cell data with Harmony_. Nat Methods. 2019;16:1289–96. [Link to paper](https://www.nature.com/articles/s41592-019-0619-0). [Link to package](https://cran.r-project.org/web/packages/harmony/index.html).

[^9]: Robinson MD, McCarthy DJ, Smyth GK. _edgeR: a Bioconductor package for differential expression analysis of digital gene expression data_. Bioinformatics. 2010;26:139–40. [Link to paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/). [Link to package](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

[^10]: Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Mauck WM 3rd, et al. _Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics_. Genome Biol. 2018;19:224. [Link to paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1). [Link to package]().
