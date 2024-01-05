# SeqTools {🧬:🛠️}

A small collection of tools for working with sequence files in a Python-based environment.

# Tools:
* `seq_gen`: Generate DNA sequences - either pure sequences, or fasta/fastq formats
* `seq_sub`: Subsample reads files
* `tools` (misc.)
    - k-mers: tools for generating k-mers from a string, as well as counting non-unique k-mers
    - Visualisation: plot histograms for read length and PHRED-score

# Requirements
* **The Big 3 (+1)**
    - `pandas`
    - `numpy`
    - `matplotlib`
    - `seaborn`
* **`Biopython`** - tools for loading, manipulating and exporting read data
* `tqdm` - used for pretty progress bars (currently only in rare cases)

For convenience, a `conda` environment file has been included.
```shell
$ conda env create -f environment.yml
```