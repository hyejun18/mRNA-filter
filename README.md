# mRNA-filter
Make a set of reference mRNAs: Major isoforms of non-redundant, non-NMD candidates with no wrong ORF

### Quick Start
``` Bash
python3 mRNA-filter.py --refflat REFFLAT --fasta-directory FA_DIR --outfile OUT_FILE
```

### Input
+ NCBI RefFlat file, downloaded from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
+ Directory containing FASTA files of assembly sequence in one file per chromosome, downloaded from [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/#:~:text=13%3A44%20%20%2027K-,hg38.chromFa.tar.gz,-2014%2D01%2D23)

### Output
Tab-separated file (RefFlat format) in extended columns: tx_Size, 5'UTR_Size, ORF_Size, 3'UTR_Size

### Dependencies
+ Python3
+ Python package [tqdm](https://github.com/tqdm/tqdm)

### Warning
mRNA-filter is on-going project. This repository can be edited or removed without any kind of notice. It was only tested by hg38 RefSeq.
