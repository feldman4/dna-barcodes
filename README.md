# dna-barcodes

Generate DNA barcode sets with guaranteed Levenshtein edit distance. 
These barcodes can be used to detect and correct sequencing errors when embedded in known sequence, such as barcoded primers or vectors.
Similar to [DNABarcodes](https://www.bioconductor.org/packages/release/bioc/html/DNABarcodes.html), but may yield larger barcode sets. 

## Install
Requires python 3.6 or higher. The example below uses conda, available here (https://docs.conda.io/en/latest/miniconda.html).

```
conda create -n dnabarcodes python=3.6
conda activate dnabarcodes
conda install -y --file requirements.txt
```
## Run

```
python barcode_design.py -h
python barcode_design.py
```

[example](example.png)
