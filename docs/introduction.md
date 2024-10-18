# What is scprocess?

{{ software_name }} is a Snakemake pipeline designed to streamline and automate critical analysis steps in single-cell and single-nuclei RNA sequencing data processing. This end-to-end, scalable solution efficiently manages data generated with 10x Genomics technology, starting from raw sequencing files and performing key tasks such as read alignment and quantification, quality control filtering, batch correction, marker gene identification, and cell type labeling. {{ software_name }} generates outputs ready for various downstream analyses and produces comprehansive reports. 

---


### Steps

{{ software_name }} performs a series of steps which include:

* Read alignment and quantification: This step utilizes [simpleaf](https://academic.oup.com/bioinformatics/article/39/10/btad614/7295550), a powerful and flexible alternative for 10x's Cellranger. Simpleaf separately reports count for spliced and unspliced transcripts.
* Ambient RNA removal (optional): Users can select between [Cellbender](https://www.nature.com/articles/s41592-023-01943-7) and [Decontx](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6) to remove ambient RNA
* Doublet identification: {{ software_name }} offers standard doublet detection with [scDblFinder](https://f1000research.com/articles/10-979/v2.) Additionally scprocess uses `scDblFinder` estimates to detect clusters that are enrichhed in doublets. 
* QC
* Integration using `harmony`
* Marker gene identification using pseudobulk values 
* Human brain cell type annotation with a classifier trained on human single cell atlas from [Siletti et.al](https://www.science.org/doi/10.1126/science.add7046)
* Subclustering
* Metacells
