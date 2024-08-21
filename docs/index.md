# What is scprocess?

Scprocess is a Snakemake pipeline designed to simplify and automate common analysis steps involved in processing single-cell and single-nuclei RNA sequencing data.

## Scope and purpose
Scprocess offers a comprehansive, automated solution for efficiently processing single cell RNA sequencing data. Starting from raw sequencing files, this pipeline performs crucial tasks such as read alignment, quantification, and extensive quality control,  data integration, batch correction, marker gene identification, and cell type labelling. Scprocess generates publication-ready results, including detailed HTML reports. It's particularly well-suited for handling large datasets. Scprocess is optimized for single cell gene expression data generated with using 10x technology. 

Scprocess performs a series of steps which include:

### Steps
* Read alignment and quantification: This step utilizes [simpleaf](https://academic.oup.com/bioinformatics/article/39/10/btad614/7295550), a powerful and flexible alternative for 10x's Cellranger. Simpleaf separately reports count for spliced and unspliced transcripts.

* Ambient RNA removal (optional): Users can select between [Cellbender](https://www.nature.com/articles/s41592-023-01943-7) and [Decontx](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6) to remove ambient RNA

* Doublet identification: scproces offers standard doublet detection with [scDblFinder](https://f1000research.com/articles/10-979/v2. Additionally scprocess uses scDblFinder estimates to detect clusters that are enrichhed in doublets. 

* QC (SampleQC)
* Integration using harmony
* Marker gene identification using pseudobulk values
* Optimised selection of highly variable genes 
* Human brain celltype annotation with a classifier trained on human single cell atlas from [Siletti et.al](https://www.science.org/doi/10.1126/science.add7046)
* Subclustering
* Metacells?


## Code Annotation Examples


### Codeblocks

Some `code` goes here.

### Plain codeblock

A plain codeblock:

```
Some code here
def myfunction()
// some comment
```

#### Code for a specific language

Some more code with the `py` at the start:

``` py
import tensorflow as tf
def whatever()
```

#### With a title

``` py title="bubble_sort.py"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

#### With line numbers

``` py linenums="1"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

#### Highlighting lines

``` py hl_lines="2 3"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```


Trying out some boxes

!!! note "Title of NOTE goes here"
    This is a test note to see if I can make a nice note in a box

!!! abstract "Title of ABSTRACT goes here"
    This is a test note to see if I can make a nice note in a box


This block should be collapsible (use ???+ to render the block expanded)

??? info "Title of INFO goes here"
    This is a test note to see if I can make a nice note in a box


!!! success "Title of SUCCESS goes here"
    This is a test note to see if I can make a nice note in a box

!!! question "Title of QUESTION goes here"
    This is a test note to see if I can make a nice note in a box

!!! warning "Title of WARNING goes here"
    This is a test note to see if I can make a nice note in a box

!!! failure "Title of FAILURE goes here"
    This is a test note to see if I can make a nice note in a box

!!! danger "Title of DANGER goes here"
    This is a test note to see if I can make a nice note in a box

!!! bug "Title of BUG goes here"
    This is a test note to see if I can make a nice note in a box


!!! example "Title of EXAMPLE goes here"
    This is a test note to see if I can make a nice note in a box


!!! quote "Title of QUOTE goes here"
    This is a test note to see if I can make a nice note in a box

