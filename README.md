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
> $ mkdir .CorrieCaller_env
> 
> $ python3 -m venv .CorrieCaller_env
> 
> $ source .CorrieCaller_env/bin/activate
> 
> $ pip install -r requirements.txt
4. CorrieCaller is ready to run!

## Inputs and parameters
### Generating the DiffBind output file

For information on running [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) please consult its documentation. For best CorrieCaller results, I recommend generating the output without filtering as shown bellow:
> res=dba.report(dObj,th=1)
> 
> write.table(res,paste(out_dir,PATH/TO/OUTPUT_FILE.csv, row.names=FALSE,quote=FALSE,sep='\t')

### Parameters
- output_dir: Directory to save the results to (string).
- test: Identifier of the test group (string). This can be found in the DiffBind output file as everything after "Conc_" in one of the column headers.
- control: Identifier of the test group (string). Generally "wild_type". This can be found in the DiffBind output file as everything after "Conc_" in another of the column headers.
- diffBind_path: Path to the DiffBind output file created in the orevious step (string).
- sWindow_size: Size of each sliding window bin to use (bp, integer). Default: 10000.
- sWindow_shift: Distance between the start (and end) of each sliding window bin (bp, integer). Default: 100.
- modification_baseLine: Threshold of average WT H3K9me3 for a bin to be considered H3K9me3-rich (float). Default: 2.
- bins_toCompare: Number of bins up/downstream used in the comparisson to predict local minima and maxima (integer). Deafult: 250.
- cpus: Number of CPUs to dedicate for running this algorithm (integer). Default: 2.

## Running CorrieCaller

1. Edit the CorrieCaller.py file and add your parameters in the block marked with "###" (see below).

![CorrieCaller parameters](/assets/CorrieCaller_parameters.png "CorrieCaller parameters")

2. Activate the Python environment.
> $ source .CorrieCaller_env/bin/activate

3. Run the script!
> $ python3 CorrieCaller.py

## Outputs
### SlidingWindow.tsv
Bins their asigned values in terms of WT signal and Log2FoldChange in the mutants.
![SlidingWindow](/assets/SlidingWindow.png "SlidingWindow")

**Columns:**

- chr: chromosome.
- start: start of the sliding window bin.
- end: end  of the sliding window bin.
- avg_l2fc: Weighted mean of the log2FoldChange in the mutant.
- std_l2fc: Weighted standard deviation of the log2FoldChange in the mutant.
- avg_control: Weighted mean of the WT signal.
- std_control: Weighted standard of the WT signal.

### enriched_regs.tsv
Regions identified as H3K9me3-rich in WT.
<img width="3760" height="217" alt="image" src="https://github.com/user-attachments/assets/9c41222e-6e3b-4cb1-afd8-561c5f3ac3a5" />

**Columns:**

- chr: chromosome.
- start: start of the H3K9me3-rich region.
- end: end of the H3K9me3-rich region.
- width: width of the H3K9me3-rich region.

### predicted_local_maxima_minima.tsv
Regions predicted as likely to contain H3K9me3 FC local maxima and minima.
<img width="3762" height="280" alt="image" src="https://github.com/user-attachments/assets/cdc0f9f9-a87e-4a80-8fa3-8045b5e1d03f" />

**Columns:**

- chr: chromosome.
- start: start of the region.
- end: end of the region.
- type: "maximum" or "minimum"

### corries.tsv
Genomic coordinates of potential Corrie domains and their direction.
<img width="3752" height="275" alt="image" src="https://github.com/user-attachments/assets/beb76ff7-61b8-4d05-a800-da6772d93c9a" />

**Columns:**

- chr: chromosome.
- start_extended and end_extended: outer edge of the local minimum/maximum at the Corrie domain borders +/- half a sliding window (See diagram bellow).
- start/end: midpoint of the local minimum/maximum at the Corrie domain border (See diagram bellow).
- start_strcit/end_strcit: inner edge of the local minimum/maximum at the Corrie domain borders (See diagram bellow).
- strand: direction of the Corrie domain. "+" means that the minimum comes first, "-" means that the maximum comes firts.
- width: width of the Corrie domain.

