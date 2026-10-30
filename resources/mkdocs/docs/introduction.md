# What is {{sc}}?

{{sc}} is a Snakemake-based pipeline for processing single cell and single nuclei RNA sequencing data generated using the 10x Genomics droplet based sequencing technology, supporting both poly(A)-capture and probe-based (“Flex”) assays. Starting from raw sequencing files, {{sc}} performs a series of standard analysis steps, delivering outputs suitable for a variety of downstream analyses. The pipeline supports analyses across multiple samples, enabling researchers to gain a complete understanding of their datasets while ensuring reproducibility and scalability in their workflows.

---

## Steps

### Overview

{{sc}} consists of a series of core steps which can be performed in a single execution of the workflow. Additional functionality is available to extend the core analyses as needed. Detailed descriptions of all steps are provided below.

### Core pipeline steps

![workflow_chart](assets/images/scprocess_workflow_diagram_white_bg.png#only-light)
![workflow_chart](assets/images/scprocess_workflow_diagram_black_bg.png#only-dark)

---
<div class="img-caption"> The figure illustrates all core steps in {{sc}}. Specific software packages used are listed for individual steps. Several steps generate HTML reports with diagnostic plots, enabling users to inspect the results at key points in the workflow. </div>

* #### Read alignment and quantification
    
    {{sc}} starts by mapping sequencing reads to the reference genome (poly(A)-capture assays) or a probe reference (probe-based assays) and performing unique molecular identifier (UMI) quantification using `simpleaf`[@He2023-fx]. This open-source alternative to [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is designed for speed and memory efficiency, while also offering the advantage of reporting counts for spliced and unspliced transcripts separately.

* #### Ambient RNA removal (optional)

    In droplet-based assays, cells or nuclei are encapsulated in droplets, but some freely-floating RNA can also be captured, which is referred to as ambient RNA. Ambient RNA contamination is particularly common in single-nucleus assays due to residual cytoplasmic material and harsh isolation protocols that can cause nuclei to rupture. Since ambient RNA can interfere with results of downstream analysis, it can be beneficial to remove this contamination *in silico*. In {{sc}}, users can select from `CellBender`[@Fleming2023-cx] and `DecontX`[@Yang2020-zz] to remove ambient RNA. 

* #### Cell calling

    An important step in single cell data analysis is distinguishing barcodes/droplets that contain cells/nuclei from those that are empty. {{sc}} employs several options for detecting barcodes that correspond to cell-containing droplets. When `CellBender` is used for ambient RNA removal, its built-in method generates a filtered counts matrix. If `DecontX` is selected or ambient RNA removal is skipped, users can choose between `barcodeRanks` and `emptyDrops` methods from the `DropletUtils`[@Lun2019-mm],[@Griffiths2018-sk] R package. `barcodeRanks` separates cell-containing and empty droplet populations by detecting key transition points on the barcode-rank curve while `emptyDrops` tests each barcode for significant deviations from the ambient profile.

* #### Doublet detection
    
    Doublets, formed when two cells are captured in the same droplet, can distort single-cell RNA-seq results. {{sc}} addresses this by using `scDblFinder`[@Germain2021-iv] for doublet detection. Additionally, {{sc}} integrates `scDblFinder`-flagged doublets with cells that pass the QC filtering step. This allows users to detect doublet-enriched clusters which can be removed from downstream analysis. 

* #### QC filtering

    In addition to removing doublets, {{sc}} filters out cells based on library size, feature counts, mitochondrial read proportions and spliced read proportions based on user-defined thresholds. The spliced proportion is a particularly informative metric in singlenuclei RNA-seq, as elevated levels may indicate residual cytoplasmic material [@Montserrat-Ayuso2024-bm].

* #### Generating pseudobulks from cells and empty droplets and ambient gene detection

    {{sc}} aggregates counts from both cell-containing and empty droplets into pseudobulk profiles to identify genes enriched in empty droplets—i.e., ambient genes. These ambient genes likely represent residual contamination rather than true biological signals. Identifying them supports cleaner and more accurate downstream analyses and helps prioritize genes less affected by contamination.

* #### Highly variable gene detection

    Highly variable gene (HVG) detection in {{sc}} is performed using the `Seurat VST`[@Stuart2019-pv] method in a chunk-wise or sample-wise manner, enabling the computation of ranking metrics for all genes without the need to load the entire dataset into memory. The efficient generation of a reduced matrix containing only highly variable genes ensures optimal performance in downstream analyses and facilitates the processing of larger datasets.

* #### Cell cycle annotation
    
    The cell cycle program is highly conserved across distinct cell lineages and can overpower cell-type-specific signals, causing different cell types to cluster together. {{sc}} offers the option to estimate and then regress out this cell cycle signal to ensure clustering is driven primarily by cell identity. First, {{sc}} uses tricycle [@Zheng2022-gd] to infer a continuous cell cycle position (theta) for each cell. It then fits a harmonic (sine/cosine) ridge regression on theta to remove phase-dependent variation from HVG expression before calculating the PCA embedding.

* #### Integration
    
    After identifying highly variable genes, {{sc}} proceeds with dimentionality reduction and clustering. In multi-sample analysis, various factors can introduce batch effects that obscure true biological signals. To mitigate this, {{sc}} offers the option to compute batch-corrected PCA embeddings using `Harmony`[@Korsunsky2019-rk]. This ensures that clustering and downstream analyses reflect true biological relationships rather than technical variation. To perform the integration step users can choose between a standard `Scanpy`-based workflow [@Wolf2018-lf] or a workflow with equivalent functionality using `RAPIDS-singlecell`. The latter laverages GPU acceleration to achieve massive performance boosts, particularly for clustering and UMAP [@Dicks2023-gx],[@Nolet2022-be]. Additionally, users can optionally run partition-based graph abstraction (PAGA)[@Wolf2019-lw] within the integration step.

* #### Marker gene identification
    
    Assigning meaningful labels to clusters in single-cell data is essential for interpretation of single cell data. This is commonly achieved by examining marker genes for each cluster, identified by comparing the expression profile of each cluster against all others. In {{sc}}, transcript counts are aggregated per cluster within each sample to generate "pseudobulk" values, which are then compared using `edger.pb` (implementation of the `edgeR`[@Robinson2010-qz],[@Chen2025-jo]package using `BPcells`[@Parks2025-cl] as a backend). This approach avoids the assumption that individual cells from the same sample are independent, thereby enhancing the statistical reliability of the results. 

* #### Gene set enrichment analysis

    {{sc}} includes an option to perform gene set enrichment analysis (GSEA)[@Subramanian2005-bd] on the set of identified marker genes using the `fgsea`[@Korotkevich2016-od] R package. By identifying biological processes and pathways unique to each cluster, GSEA can provide additional evidence for characterizing cellular identity. This feature is available for human and mouse datasets.

* #### Cell type labelling

    {{sc}} enables automated cell type annotation of human brain datasets using an `XGBoost` classifier trained on the adult human whole-brain dataset[@Siletti2023-rb]. {{sc}} alo supports cell type annotation using pre-trained models available through `CellTypist`[@Xu2023-al],[@Dominguez_Conde2022_kd].

* #### Processing multiplexed samples

    Multiplexing strategies are commonly used to scale up single cell experiments, allowing for simultaneous analysis of multiple samples within a single reaction. A common approach includes tagging cells or nuclei in individual samples with lipid-bound or antibody-conjugated oligonucleotides (hashtag oligos, or HTOs) prior to pooling. Prior labelling is not required for on-chip multiplexing (OCM) where up to four samples are physically separated inside the microfluidic chip, each receiving a unique set of cell barcodes. For scaling up to 384 samples, the Flex technology leverages gene-specific probes tagged with sample barcodes that hybridize directly within fixed cells prior to sample pooling. Alternatively, in studies involving genetically distinct individuals, sample identities can be recovered based on natural genetic variation. {{sc}} supports the analysis of multiplexed samples by leveraging `simpleaf` to generate HTO count matrices from associated cDNA libraries, which are subsequently processed via the `Seurat HTODemux`[@Stoeckius2018-os] function. Furthermore, {{sc}} accommodates multiplexed Flex data, OCM datasets and outputs from external demultiplexing tools, enabling seamless processing of samples regardless of the multiplexing strategy employed.

### Extended functionality

* #### Subclustering

    {{sc}} offers a subclustering feature that enables users to perform a second round of analysis on a specific subset of cells and includes the following steps: optional reevaluation of QC thresholds, generation of pseudobulks from selected cells and detecting ambient genes, identification of higly variable genes, data integration and marker gene identification. Cell subsets can be defined based on user-provided cell type labels, clusters identified during the primary round of {{sc}}, or cell type labels assigned with a selected classifier. This functionality is particularly valuable when a primary cluster or cell type contains diverse cell states, developmental stages, or activation states that warrant more detailed exploration.

* #### Cross-project analysis

    {{sc}} provides a join module to integrate outputs across multiple projects (full-dataset runs, subclustering results, or a combination of both), requiring only a YAML file that lists individual project-level configuration files. Integration begins with the selection of consensus highly variable genes (HVGs) across projects, followed by marker gene identification.

* #### Interactive visualization of results

    Any {{sc}} output can be turned into a deployable `Shiny` application for interactive data exploration across four modules. These include “explore genes” for assessing gene expression, “explore clusters” for inspecting cluster markers, “explore gene sets” for comparing gene groups (e.g., transcription factors, pathways, or custom sets) across clusters, and “explore prevalence” for evaluating how cluster composition varies across metadata variables.   
    
* #### Training and reusing custom XGBoost classifiers

    Pre-trained classifiers for cell type annotations are not available for certain tissue and species combinations. To avoid repeated time-consuming manual annotation, {{sc}} allows users to train custom `XGBoost` classifiers using their own annotated reference datasets. Such classifiers can be registered alongside built-in models and reused in subsequent labelling runs.


