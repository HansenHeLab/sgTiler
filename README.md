# sgTiler
An ultra-fast tool to design tiling sgRNAs for any genomic region

sgTiler predicts efficient sgRNA spacer sequences in any input DNA sequence and designs optimally distributed spacer sequences covering the maximum sequence while minimizing the number of spacer sequences.

### Installation
sgTiler requires python 2.7+ and bowtie installed in the system. 
Please refer to http://bowtie-bio.sourceforge.net/index.shtml for bowtie installation.
path to bowtie must be in the environment, i.e., the command `bowtie` must be executable from any location.

### Running
Download sgTyler.py and run in any command line environment.

### Input files
sgTiler requires four input files:
1. bowtie index file for the desired genome
2. Input fasta file
3. Genome exon regions in bed format
4. Open chromatin or any histone mark in bed format

### Command options
Two major input parameters are `--sg-expected` and `--sg-flex` which determines the distance between two sgRNAs and strictness of the distance, respectively. For example, `--sg-expected 10` indicates that 10 spacer sequences are desired in each 100bp input sequence. Accordingly, sgTiler will detect positions in the input sequences that are evenly disparsed and 10bp away from each other - 5,15,25,35..95. But not always there is a candidate spacer sequence at these exact positions. For this, sgTiler makes windows spanning these positions and picks the best sgRNA in each window. The size of these windows are determined by `--sg-flex`. For example, `sg-expected 10 --sg-flex 3` will makes windows of 2-8, 12-18, 22-28 and so on. The tool implements a scoring method to pick the best sgRNA within each of these windows. An user can play with these two parameters to decide how densely or disparsed the tiling should be. The tool outputs figures of distribution of spacer sequence for each input sequence. The user can check the figures to choose the optimum numbers.

Please type `python sgTiler.py -h` from your command line to see other parameters.

### Example run command
```python
python sgTiler.py -i input.fa --bowtie-index hg19.ebwt/hg19 --dhs dhs.bed --gtf allExons.bed --sg-expected 8 --verbose --dir output_boxplots --output sgTiler_output
```

### Output
The tool output four text files and two pdf files:
1. .all.txt - list of all candidate sgRNAs
2. .sgRNAs.txt - list of filtered sgRNAs
3. .stats.txt - list of sgRNA details for the input sequences 
4. .report.txt - a summary report with important statistics
5. .sgrna_count.pdf - graphical summaries of no. of sgRNAs and
6. .bp_coverage.pdf - sequence coverage per input region.

Additionally, SgTiler generates graphical representation of distribution of sgRNA for each individual input region. Combining the overall statistics, user can predict the success of the screening.

#### Additional help
Pre-built bowtie index for several genomes are available here: ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/
