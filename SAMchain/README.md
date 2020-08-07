# SAMChain1.0
SAMchain code organized by owner, sequencer, and clinical node 
This repository contains files to (1) push genomic data from a BAM file to a multichain blockchain, and (2)analyze data stored in the chain.

Please refer to the SAMChain1.0-network-ecosystem file for details on how to exchange data in a SAMChain network. 

### dependencies
To use this code you will need the following:

python (https://www.python.org/downloads/)

multichain (https://www.multichain.com/download-community/)

pysam (https://pysam.readthedocs.io/en/latest/installation.html)

pandas (https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

To develop and test this code, Python 2.7.16 and MultiChain 2.0.3 Daemon (Community Edition, latest protocol 20011) were used on Ubuntu 19.04 (GNU/Linux 5.0.0-37-generic x86_64).

### buildChain.py
Use this file to initialize the chain and initialize streams within the chain
One would do all of this on the command line through the following line:

``` python buildChain.py -cn <my-chain-name> -dr <datadir> ```
  
This will create a multichain with name <my-chain-name> in the directory <datadir>. 
  
Many other arguments can be passed to buildChain.py, including chromosome length and number, stream length, and read length. These are set to default values at the end of the file, and can be altered according to the needs of the user. 

After building a chain, multichain keeps a directory ```<datadir>/<my-chain-name>``` to store the data from the chain. It also continuously runs a process for the chain.
  

Given a chain <my-chain-name>, one can query a position in the genome, perform depth analysis of a position, or perform pileup analysis of a position. One can also rebuild a BAM file from the read data in the chain.


### insertData.py
Use this file to push data from an input BAM file to the chain.
One would do all of this on the command line through the following line:

``` python insertData.py <mydata>.bam -cn <my-chain-name> -dr <datadir> ```
  
This will insert BAM data to a subscribed multichain with name <my-chain-name> in the directory <datadir>.   

### queryReads.py

Use this file to query a position in the genome via the following line:

``` python queryReads.py <my-chain-name> chr<x>:<pos1>-<pos2> -dir <datadir> -ref <ref-genome-filename> ```

This will return a list of all reads in BAM format in the chain that fall in the input position range. 

Example:
``` python queryReads.py chain1 chr1:1880913 -dir ./myChains -ref GRCh38_no_alt_analysis_set_GCA_000001405.15.fa```

might return this read

```FC30V1RHM_20081223:8:80:961:561  0 chr1  1880913 14M1D14M  CAGGCTACGAAGACAGAGTGGGGTAAAA <OPTIONAL TAGS>```


### queryDepth.py


Use this file to query a position in the genome for depth analysis via the following line:

``` python queryDepth.py <my-chain-name> chr<x>:<pos1>-<pos2> -dir <datadir> -ref <ref-genome-filename> ```

This will return a list of the depth statistics for all reads in the chain that fall in the input position range.

Example:

``` python queryDepth.py chain1 chr1:1880913-1880917 -dir ./myChains -ref GRCh38_no_alt_analysis_set_GCA_000001405.15.fa ```

might return these depth statistics:

```
chr1 1880913 1
chr1 1880914 1
chr1 1880915 1
chr1 1880916 1
chr1 1880917 1 
```
### pileup.py


Use this file to query a position in the genome for pileup analysis via the following line:

``` python pileup.py <my-chain-name> chr<x>:<pos1>-<pos2> -dir <datadir> -ref <ref-genome-filename> ```

This will return a list of the pileup statistics for all reads in the chain that fall in the input position range.

Example:

``` python pileup.py chain1 chr1:1880913-1880917 -dir ./myChains -ref GRCh38_no_alt_analysis_set_GCA_000001405.15.fa ```

might return these pileup statistics:

```
chr1 1880913 C
chr1 1880914 A
chr1 1880915 G
chr1 1880916 G
chr1 1880917 C 
```

### buildBAM.py


Use this file to build a BAM file from the read data stored in an input chain via the following line:

``` python buildBAM.py <my-chain-name> -dr <datadir> -bam <filename>.bam -ref <ref-genome-filename> ```

This will make a BAM file <filename>.bam in directory containing all read data stored in chain <my-chain-name>.
