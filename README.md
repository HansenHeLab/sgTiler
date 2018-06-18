# sgTiler
An ultra-fast tool to design tiling sgRNAs for any genomic region

sgTiler predicts efficient sgRNA spacer sequences and distribute them optimally in any input DNA sequence. sgTiler is a great tool to design tiling sgRNA library that aims to minimalize the number of sgRNAs needed to maximally cover the input DNA sequence. This tool provides great flexibility to users to design the library in MUCH greater speed than any other sgRNA desinging tool currently available. sgTiler is best suited for designing tiling sgRNAs targeting regulatory regions including promoters and enhancers, however, it can also be used to target exons or any other part of the genome.

### Installation
sgTiler requires python 2.7+ and bowtie installed in the system.

1. Download sgTiler.py

2. Install bowtie 1.x. Please refer to http://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie for bowtie installation. **Note:**  path to bowtie must be in the PATH variable, i.e., the command `bowtie` must be executable from any location. Easiest way to do that is to download the appropriate version of bowtie from https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/ to where you downloaded sgTiler.py and make the file executable. You should download [bowtie-1.2.2-mingw-x86_64.zip](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-mingw-x86_64.zip/download) for Windows computer, [bowtie-1.2.2-macos-x86_64.zip](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-macos-x86_64.zip/download) for Mac and [bowtie-1.2.2-linux-x86_64.zip](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download) for Linux.

3. Download bowtie genome index file from [here](<ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/>). E.g., for Hg19, [download this file](ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip) and unzip in the current directory.

3. Download [exons.zip](https://github.com/HansenHeLab/sgTiler/blob/master/exons.zip) in the current directory and unzip.

4. Optionally, downlaod [wgEncodeRegDnaseClusteredV3.consensus.simplified.bed](https://github.com/HansenHeLab/sgTiler/blob/master/wgEncodeRegDnaseClusteredV3.consensus.simplified.bed) in the current directory.
 

### Running
Run `sgTiler.py` in any command line environment.

### Input files
sgTiler requires four input files:

1. bowtie index file for the desired genome. Pre-built bowtie index for several genomes are available here: ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/

2. Input fasta file

3. Genome exon regions in bed format (provided in the github page or you can use your own). The provided exons regions are curanted from GENCODE v19 for Hg19. You can provide your own bed file.

4. Open chromatin or any histone mark in bed format. A list of highly consensus DNase hypersensitive (DHS) regions common in at least 113 cell lines is provided in the github page (wgEncodeRegDnaseClusteredV3.consensus.simplified.bed). This file is downloaded from [ENCODE project](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz) and curated in a way that any DHS site present in less than 113 cell lines are removed. It is recommended to use your own open chromatin region, e.g., DNase sensitive sites or H3K27ac marks or H3K4m1 marks or transcription factor binding sites relevant to the cell line or system you are performing the screening in. However, if you do not have such marks, you can use the file provided which will only likely overestimate the off-target potential score.

**(Please note that the provided files are in Hg19)**

### Example run command
```bash
python sgTiler.py -i input.fa --bowtie-index hg19.ebwt/hg19 --dhs -wgEncodeRegDnaseClusteredV3.consensus.simplified.bed -gtf allExons.sorted.merged.gencodev19.hg19.bed --verbose --dir output_boxplots --output sgTiler_output
```

### Command options
The sgTiler.py has the following command options:
```
-i              Required. input fasta file
--bowtie-index  Required. Path to bowtie index file
--dhs           Required. Path to regulatory regions bed file 
--gtf           Required. Path to exon bed file 
--output        Required. Output prefix
--pam           PAM sequence. Default: NGG
--gc-min        Minimum GC content in percentage. Default: 20
--gc-max        Maximum GC content in percentage. Default: 80
--nthreads      No. of threads for parallel processing. Default: 4
--strand        Strand to find the sgRNAs. Options: positive, negative or both. Default: both
--length        Length of sgRNA without the PAM. Default: 19
--missmatch     Minimum missmatch allowed in offtargets. Default: 2
--sg-expected   Expected approximate # of sgRNAs per 100bp. Default: 7
--sg-flex       Room of flexibility for evenness in distribution. Default: 3
--dir           Directory to store plots
--plot-off      Do not generate plots
--optimize-off  Do not perform optimization
--distribution-off Do not filter for distribution
--save-tmp      Save all temporary files
-v              Turn on verbosity
-h              Show command help 
```
**It is recommended to leave all options to their default values when possible.** However, two major input parameters which highly infleunce the number of sgRNAs are `--sg-expected` and `--sg-flex`. `--sg-expected` and `--sg-flex` determines the distance between two sgRNAs and strictness of the distance, respectively. For example, `--sg-expected 10` indicates that 10 spacer sequences are desired in each 100bp input sequence. Accordingly, sgTiler will detect positions in the input sequences that are evenly disparsed and 10bp away from each other - `5,15,25,35..95`. But not always there is a candidate spacer sequence at these exact positions. For this, sgTiler makes windows spanning these positions and picks the best sgRNA in each window. The size of these windows are determined by `--sg-flex`. For example, `sg-expected 10 --sg-flex 3` will makes windows of `2-8, 12-18, 22-28` and so on. The tool implements a scoring method to pick the best sgRNA within each of these windows. An user can play with these two parameters to decide how densely or disparsed the tiling should be. The tool outputs figures of distribution of spacer sequence for each input sequence. The user can check the figures to choose the optimum numbers.

sgTiler automatically chooses the best sgRNAs combining their efficiency score and off-target potential. The user only has to decide how dense or dispersed the library should be.


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

Please email musaddeque.ahmed@gmail.com for any further help.
