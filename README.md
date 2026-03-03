# CorrieCaller
Call potential Corrie domains straight from DiffBind outputs!

## Introduction

**CorrieCaller** is a script that allows for the calling of Corrie domains straight from DiffBind output. Corrie domains are a type of heterochromatic domain evident in mutants with deficient histone recycling. They are characterized by being localized between a local minima and maxima of H3K9me3 change and are associated with transcriptional upregulation and AT-rich transposable elements (TE). For more information read: (insert DOI here).

**Running time:** The runtime was <10min or *Caenorhabditis elegans* (ce11) and >1h for *Mus musculus* (mm39) using 18 CPUs. Total run time will depend on genome size and parameters used.

**Requirements:** This script was tested on Python 3.8.10.

**Note:** This script is likely to produce artefacts. Please verify potential Corrie domains through other methods.

## Algorithm pipeline

![CorrieCaller Algorithm steps](/assets/Slide1.PNG "CorrieCaller algorithm")

## Installation

1. Create a CorrieCaller directory.
> $ mkdir CorrieCaller
2. Download "CorrieCaller.py" and "requirements.txt" and save them in the CorrieCaller directory.
3. Create a Python environment and install the required libraries.
> $ mkdir ./CorrieCaller_env
> 
> $ python3 -m venv ./CorrieCaller_env
> 
> $ source ./CorrieCaller_env/bin/activate
> 
> $ pip install -r requirements.txt
4. CorrieCaller is ready to run!

## Inputs and parameters
### Generating the DiffBind output file

For information on running [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) please consult its documentation. For best CorrieCaller results, I recommend generating the output without filtering as shown bellow:
> res=dba.report(dObj,th=1)
> 
> write.table(res,paste(out_dir,PATH/TO/OUTPUT_FILE.csv, row.names=FALSE,quote=FALSE,sep='\t')

