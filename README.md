# RNACache

RNACache is a transcriptomic mapping tool that maps RNA-seq reads directly to a transcriptomes using MinHashing, a locality sensitive hashing technique based on pseudo-random k-mer subampling, and several filters utilizing both per-read and read-global information.

## Installation

#### Requirements
RNACache itself should compile on any platform for which a C++14 conforming compiler is available.

Dependencies for the default mode are included in the repository as source.
For the BAM-format output htslib is required, which is included only for linux x86-64 as a static library.

#### Compile
Run 'make' in the directory containing the Makefile.
This will compile RNACache with the default data type settings which support databases with up to 2^32-1 reference sequences (targets) and k-mer sizes up to 16.

Using the following compilation options you can compile RNACache with support for more or fewer reference sequences and greater k-mer lengths.

##### BAM support
To compile RNACache with support for the BAM output format, htslib is required, which is only included as a static library for linux (x86-64), and must otherwise be installed manually and the Makefile adjusted accordingly. (This will change in the future.)

* To compile with BAM support, start make with the RC_BAM=TRUE environment variable as in:
  ```
  RC_BAM=TRUE make
  ```

##### number of referece sequences (targets)

* support for up to 65,535 reference sequences (needs less memory):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint16_t"
  ```

* support for up to 4,294,967,295 reference sequences (default):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint32_t"
  ```

* support for more than 4,294,967,295 reference sequences (needs more memory)
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint64_t"
  ```

##### reference sequence lenghts

* support for targets up to a length of 65,535 windows (default)
  with default settings (window length, k-mer size) no sequence length must exceed 7.4 million nucleotides
  ```
  make MACROS="-DMC_WINDOW_ID_TYPE=uint16_t"
  ```

* support for targets up to a length of 4,294,967,295 windows (nees more memory)
  with default settings (window length, k-mer size) no sequence length must exceed 485.3 billion nucleotides
  ```
  make MACROS="-DMC_WINDOW_ID_TYPE=uint32_t"
  ```


##### kmer lengths
* support for kmer lengths up to 16 (default):
  ```
  make MACROS="-DMC_KMER_TYPE=uint32_t"
  ```

* support for kmer lengths up to 32 (needs more memory):
  ```
  make MACROS="-DMC_KMER_TYPE=uint64_t"
  ```

You can of course combine these options (don't forget the surrounding quotes):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint32_t -DMC_WINDOW_ID_TYPE=uint32_t"
  ```

**Note that a database can only be queried with the same variant of RNACache (regarding data type sizes) that it was built with.**

In rare cases databases built on one platform might not work with RNACache on other platforms due to bit-endianness and data type width differences. Especially mixing RNACache executables compiled with 32-bit and 64-bit compilers might be probelematic.


## Building Databases
   
### [Building Custom Databases...](docs/building.md)

## Mapping
Once a database is built you can map reads.
* a single FASTA file containing some reads:
  ```
  ./rnacache query refseq my_reads.fa -out results.txt
  ```
* an entire directory containing FASTA/FASTQ files:
  ```
  ./rnacache query refseq my_folder -out results.txt
  ```
* paired-end reads in separate files:
  ```
  ./rnacache query refseq my_reads1.fa my_reads2.fa -pairfiles -out results.txt
  ```
* paired-end reads in one file (a1,a2,b1,b2,...):
  ```
  ./rnacache query refseq my_paired_reads.fa -pairseq -out results.txt
  ```


### [Output Interpretation, Analysis & Formatting Options...](docs/output.md)


## Documentation of Command Line Parameters

* [for mode `build`](docs/mode_build.txt): build database from reference transcriptome sequences
* [for mode `query`](docs/mode_query.txt): query reads against database
* [for mode `modify`](docs/mode_modify.txt): add reference sequences to database
* [for mode `info`](docs/mode_info.txt): obtain information about a database


#### View options documentation from the command line
List available modes:
```
./rnacache help
```
or jump directly to a mode's man page with:
```
./rnacache help build
./rnacache help query
...
```