# What is {{ software_name }}?

{{ software_name }} is a Snakemake pipeline designed to automate the analysis of single-cell and single-nuclei RNA sequencing data. This comprehensive solution effectively handles data generated using 10x Genomics technology. Starting from raw sequencing files, {{ software_name }} performs a series of standard analysis steps, delivering outputs suitable for a variety of downstream analyses. The pipeline supports analyses across multiple samples, enabling researchers to gain a complete understanding of their datasets while ensuring reproducibility and scalability in their workflows.

---

## Steps

{{ software_name }} consists of a series of core steps, organized under the main Snakemake rule - `rule all`. Additional optional steps are available to extend the core analyses as needed. The diagram below outlines all steps, with detailed descriptions provided in the following sections.

![workflow_chart](assets/images/scprocess_workflow_diagram_white_bg.png#only-light)
![workflow_chart](assets/images/scprocess_workflow_diagram_black_bg.png#only-dark)

---
<div class="img-caption"> The figure illustrates all steps in {{ software_name }}, including those grouped under <code>rule all</code> and optional steps. {{ software_name }} requires three input file types: raw single-cell FASTQ files generated with 10x Genomics technology, sample metadata, and a configuration file specifying analysis parameters. Specific software packages used are listed for individual steps. Some steps process samples independently, while others operate on a combined dataset with multiple samples. Several steps also generate HTML reports with diagnostic plots, enabling users to inspect the results at key points in the workflow. </div>

### `rule all` steps

* #### Read alignment and quantification
    
    {{ software_name }} starts by mapping sequencing reads to the reference genome and performing transcript quantification using [simpleaf](https://academic.oup.com/bioinformatics/article/39/10/btad614/7295550). This open-source alternative to [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is designed for speed and memory efficiency, while also offering the advantage of reporting counts for spliced and unspliced transcripts separately.

* #### Ambient RNA removal (optional)

    In droplet-based assays, cells or nuclei are encapsulated in droplets, but some freely floating RNA can also be captured. This RNA is referred to as ambient RNA. Ambient RNA contamination is particularly common in single-nucleus assays due to residual cytoplasmic material and harsh isolation protocols that can cause nuclei to rupture. Since ambient RNA can interfere with results of downstream analysis, it is beneficial to remove this contamination *in silico*. In {{ software_name }}, users can select from [`Cellbender`](https://www.nature.com/articles/s41592-023-01943-7) and [`DecontX`](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6) to remove ambient RNA. 

* #### Cell calling

    {{ software_name }} employs several options for detecting barcodes that correspond to cell-containing droplets. When `Cellbender` is used for ambient RNA removal, its built-in method generates a filtered counts matrix. Alternatively, if `DecontX` is selected or ambient RNA removal is skipped, users can choose between `barcodeRanks` and `emptyDrops` methods from the [`DropletUtils`](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package. `barcodeRanks` separates cell-containing and empty droplet populations by detecting key transition points on the barcode-rank curve while `emptyDrops` tests each barcode for significant deviations from the ambient profile.

* #### Doublet detection
    
    Doublets, formed when two cells are captured in the same droplet, can distort single-cell RNA-seq results. {{ software_name }} addresses this by using [`scDblFinder`](https://f1000research.com/articles/10-979/v2) for doublet detection. Additionally, {{ software_name }} integrates `scDblFinder`-flagged doublets with cells that pass the QC filtering step. This allows users to detect doublet-enriched clusters which can be removed from downstream analysis. 

* #### QC filtering

    In addition to removing doublets, {{ software_name }} filters out cells with low library size, low feature counts, high mitochondrial read proportions, and high spliced read proportions using predefined thresholds. The spliced proportion is a particularly informative metric in single-nuclei RNA-seq, as elevated levels may indicate residual cytoplasmic material. 

* #### Integration (Batch correction)

    In multi-sample analyses, various factors can introduce batch effects that obscure true biological signals. {{ software_name }} uses [`Harmony`](https://www.nature.com/articles/s41592-019-0619-0) to address this by aligning cells across batches based on shared expression profiles, ensuring that clustering and downstream analyses reflect true biological relationships rather than technical variation.

* #### Marker gene identification

    Cluster identities can be defined based on marker genes, typically identified by comparing the expression profile of each cluster against those of all other clusters. In {{ software_name }}, transcript counts for each cluster are aggregated within each sample, and these pseudobulk values are compared using [`edgeR`](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/). This approach avoids the assumption that cells within the same sample are independent, thereby enhancing the statistical reliability of the results. {{ software_name }} also performs gene set enrichment analysis on all marker genes. Additionally, users have the option to visualize the expression of canonical marker genes specific to different tissue types in the HTML report.

### Optional steps

* #### Cell type labelling

    {{ software_name }} provides automated cell type annotation using XGBoost classifiers for various tissue types, including human and mouse brain, as well as human and mouse peripheral blood mononuclear cells (PBMCs). Classifiers for cell type annotation are trained on the following datasets:

    - Human Brain Classifier: [Transcriptomic diversity of cell types across the adult human brain](https://www.science.org/doi/10.1126/science.add7046) [has to be updated]

    - Mouse Brain Classifier: [The molecular cytoarchitecture of the adult mouse brain](https://www.nature.com/articles/s41586-023-06818-7) [work in progress]

    - Human PBMC Classifier: [Multidimensional single-cell analysis of human peripheral blood reveals characteristic features of the immune system landscape in aging and frailty](https://www.nature.com/articles/s43587-022-00198-9#data-availability) [work in progress]

    - Mouse PBMC Classifier: [work in progress]

    In addition to the pre-trained classifiers, {{ software_name }} allows users to provide a custom table with cell type annotations which can be used as input for other optional steps (metacells and pseudobulks)

* #### Subclustering

    {{ software_name }} provides a subclustering feature that allows users to delve deeper into selected clusters. This approach is particularly useful when a primary cluster encompasses diverse cell states, developmental stages, or activation states that warrant closer examination. By selecting one or multiple clusters of interest, users can initiate a secondary round of integration and clustering within those specific groups.

* #### Merging cells into metacells

    In large datasets, high cell counts can severely hinder or even block downstream analyses. To address this, {{ software_name }} uses the [`SuperCell`](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04861-1) package to merge transcriptionally similar cells into "metacells." By reducing the overall number of data points, this approach accelerates computations and makes complex datasets more manageable.

* #### Generating pseudobulks from cells and empty droplets

    After cell type labeling, {{ software_name }} allows users to aggregate counts from identified cell types and empty droplets into pseudobulks. Additionally {{ software_name }} identifies genes enriched in empty droplets. This step supports cleaner, more accurate downstream analyses and helps prioritize non-contaminating genes.


