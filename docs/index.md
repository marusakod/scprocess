# What is scprocess? <img src="assets/images/scprocess_logo.png" alt="image title" align="right" width="360" height="180">

Scprocess is a Snakemake pipeline designed to streamline and automate critical analysis steps in single-cell and single-nuclei RNA sequencing data processing. This end-to-end, scalable solution efficiently manages data generated with 10x Genomics technology, starting from raw sequencing files and performing key tasks such as read alignment and quantification, quality control filtering, batch correction, marker gene identification, and cell type labeling. Scprocess generates outputs ready for various downstream analyses and produces comprehansive reports. 

---


### Steps

Scprocess performs a series of steps which include:

* Read alignment and quantification: This step utilizes [simpleaf](https://academic.oup.com/bioinformatics/article/39/10/btad614/7295550), a powerful and flexible alternative for 10x's Cellranger. Simpleaf separately reports count for spliced and unspliced transcripts.

* Ambient RNA removal (optional): Users can select between [Cellbender](https://www.nature.com/articles/s41592-023-01943-7) and [Decontx](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6) to remove ambient RNA

* Doublet identification: scproces offers standard doublet detection with [scDblFinder](https://f1000research.com/articles/10-979/v2.) Additionally scprocess uses scDblFinder estimates to detect clusters that are enrichhed in doublets. 

* QC
* Integration using harmony
* Marker gene identification using pseudobulk values
* Optimised selection of highly variable genes 
* Human brain cell type annotation with a classifier trained on human single cell atlas from [Siletti et.al](https://www.science.org/doi/10.1126/science.add7046)
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

