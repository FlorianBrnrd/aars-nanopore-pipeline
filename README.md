# aaRS Nanopore Pipeline
Python scripts for identifying and extracting trans-spliced Aminoacyl-tRNA synthetases (aaRSs) in kinetoplastids sequencing libraries using nanopore direct-cDNA sequencing.

___

### Prerequisites

The scripts were developped using:
- python 3.6 
- [pysam](https://pysam.readthedocs.io/en/latest/index.html) 0.15.3
- [parasail-python](https://github.com/jeffdaily/parasail-python) 1.1.12

____


### Extracting aaRS sequences

This script is used for finding and extracting reads mapped to any aaRS transcript after alignment onto a transcriptome of reference.

```
extract_aars.py [-h] -i INPUT -o OUTPUT

optional arguments:
  -h, --help                      Show this help message and exit
Main options:
  -i, --input [SAM/BAM File]      Transcriptome alignments file
  -o, --output [Fasta File]       Output file name
```

___

### Finding splice leader sequence (SL) in reads
This script allows to extract fasta records containing a given motif. By default, it extracts reads containing the *Leishmania tarentolae* splice leader (SL) sequence.
```
search_motif.py [-h] -i INPUT -o OUTPUT [-m MOTIF] [-s SENSITIVITY]

optional arguments:
  -h, --help                      Show this help message and exit
Main options:
  -i, --input [Fasta file..]      Fasta file
  -o, --output [Fasta file..]     Output file 
  -m, --motif [Sequence]          Sequence to be searched (default: SL sequence)
  -s, --sensitivity [0-1]         Sensitivity value (default:0.7)
```

___

### License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/FlorianBrnrd/aars-nanopore-pipeline/blob/main/LICENSE) file for details

