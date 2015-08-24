# STRBait 

## SUMMARY
STRBait is a collection of three programs to automate the discovery of exact short tandem repeats (and its flanking region) from Illumina fastq sequences. Not yet ready for public consumption.

## REQUIREMENTS
STRBait should work on any standard 64 bit Linux environment with 

- GCC
- Python version >= 2.6
- Google Sparsehash (http://code.google.com/p/sparsehash)
- zlib (http://zlib.net)

## ACKNOWLEDGEMENTS
STRBait uses 128 bit version of MurmurHash as the hash function. MurmurHash
was created by Austin Appleby in 2008. Seungyoung Kim has ported it's cannonical
implementation to C language in 2012 and published it as a part of qLibc
component (dedicated to the public domain by the authors). The version used in 
STRBait is heavily borrowed from that.

## INSTALLATION
Make sure you have installed Google Sparsehash before moving to the next step.
zlib, python, gcc are normally pre-installed on most Linux systems. If not,
please go ahead and install them as well.

Modify the src/Makefile to add the Google Sparsehash header files to line 51
after `-I`. Once that is done please do the following in the top level directory
of the distribution:
```
make && make install
```
This should create a `bin` folder and copy all the required binaries to it.
Please make sure that following files are in the `bin` folder:
```
extend_STR_reads
fastq.py
merge_STR_reads
select_STR_reads
```
## DESCRIPTION
STRBait consists of 3 modules that should be run sequentially to find exact STR's
along with their flanking bases. Here we describe the 3 modules and the general
idea behind them.

### select_STR_reads
Select and annotate reads that harbor a short tandem repeat.
```
    usage:
        select_STR_reads [options] reads1.fq reads2.fq ...

    where the options are:
        -h,--help : print usage and quit.
        -n,--numcopies: require at least this number of copies. [2]
        -f,--flank: require >= these many bps around the STR on each side. [31]
        -v,--version: print version and quit.
        -2,--2mer : print reads that harbor only 2mers.
        -3,--3mer : print reads that harbor only 3mers.
        -4,--4mer : print reads that harbor only 4mers.
        -5,--5mer : print reads that harbor only 5mers.
        -6,--6mer : print reads that harbor only 6mers.
        -r,--removehm: ignore STR that harbor homopolymer runs.
```

This script uses the `re` module in python to find all STR's of length 2-6.
    
#### Notes:
- We throw away the read if that longest motif has an `N` in it.
- We print the motif, number of copies of the motif, 0 based start position
  of the motif, and end position of the motif in the name of the
  read. We also print the information about the motif if the read was 
  reverse complemented and analyzed.
- We include homopolymer runs in the output. For example `AA` is considered
  a valid motif/pattern. This can be disabled using the `-r` option.
- The user should always select a flanking distance that is longer than 
  the kmer length he/she intends to use for the subsequent analyses.
- We print Sanger quality values for the reads that are selected.
- In case where we find an STR that could be the result of repetition of
  two or more mers, then we always select the lower'mer. For example, 
  ATATATAT can be a result of repetition of AT or ATAT. We select AT in 
  this example.
- reads1.fq, reads2.fq, ... are all files in FASTQ format which contain the
  reads from a single sample.
- The read files (reads1.fq, reads2.fq, ...) can be zipped. This module
  can handle them safely if the file names end with a `.gz`.
- The user should use --2mer if they are only interested in finding reads
  that harbor a STR formed as a result of repetition of a 2mer. --2mer
  --3mer will result in finding of reads that harbor either a 2mer STR or a
  3mer STR and so on. By default reads that harbor [2,6]mer all found.

### merge_STR_reads
Merge reads that support the same STR.

```
    usage:
        merge_STR_reads [options] klength reads.str.fq

    where the options are:
        -h,--help : print usage and quit
        -d,--debug: print debug information
        -v,--version: print version and quit
        -x,--min_threshold: discard blocks with < min_threshold reads [3]
        -y,--max_threshold: discard blocks with > max_threshold reads.[10000]
```

- klength refers to the kmer length to be used.
- reads.str.fq refers to the fastq file with reads that have STR's (in most
  cases this is the output from select_STR_reads.

#### Notes:
- This script assumes that the fastq quality values in the input file are
  in Sanger format.
- The kmer length in this module should be <= the size of the flanks in the
  input file reads.str.fq (which would be at least equal to the value of
  the argument -f for select_STR_reads)
- The arguments min_threshold and max_threshold should be set to ignore
  reads that have a lot of errors or are in regions which are large repeats 
  and would be difficult to untangle.
    
### extend_STR_reads
Extend fastq reads based on the kmer structure from Illumina reads.

```
    usage:
        extend_STR_reads [options] gs klen str.reads.fq reads1.fq reads2.fq ...
    
    options:
        help: print this string and quit.[--nohelp]
        debug: print extra debug information in this run.[--nodebug]
        min_threshold: discard kmers that are observed < min_threshold times. 
                       [--min_threshold=2]
        max_threshold: discard kmers that are observed > max_threshold times. 
                       [--max_threshold=255]
        progress: print progress every so many sequences[--progress=1000000]
        flanks: the maximum size of flanks extension[--flanks=1024]
        ploidy: the ploidy of the genome[--ploidy=2]
        heterozygosity: fraction of nucleotides that differ between inherited
                        chromosomes[--heterozygosity=0.001]
        errorrate: expected error rate in sequencing[--errorrate=0.01]
```

- gs is the expected genome size of the sample.
- klen refers to the kmer length that should be used. It should be <= the size
  of the flanks in the input file reads.str.fq (which would be at least equal 
  to the value of the argument -f for select_STR_reads).
  str.reads.fq is the output from merge_STR_reads.
- reads1.fq, reads2.fq, ... are all files in FASTQ format which contain the
  reads from a single sample.
- flanks represents the maximum extension of the reads on either side. By
  defaults, the reads are extended at most 1024 bases on both sides.
- ploidy represents the number of sets of chromosomes in the nucleus of the
  cell. By default we assume a value of 2 to represent mammalian genomes.
- heterozygosity refers to the fraction of nucleotides that differ between
  inherited chromosomes in the species. By default a value of 0.001 is used.
- errorrate is the approximate sequencing error rate, which is assumed to be
  1% by default to represent Illumina short sequences.

#### Notes:
- This module can handle zipped FASTQ files as well, just like
  select_STR_reads and merge_STR_reads.
- This module uses a bloom filter to throw out kmers that are observed less
  than min_threshold times. We iterate the sequences in reads1.fq, 
  reads2.fq... twice to calculate the correct kmer counts and then
  ignore kmers that are either observed less than min_threshold times, or 
  observed greater than max_threshold times.
- We extend the flanks of the STR regions up to 1024 bases on both sides.
  The code can be changed relatively easily to handle larger values, but
  since the idea is to have flanks for PCR amplification, that would not be
  of much value.
- The output of this module is mini-contigs, which contain the STR regions.
- The names of the contigs contain a unique identifier, 0-based [start,end)
  coordinates for the STR as well.
- The sample used in this case is assumed to be diploid. If this module
  finds that the sample has two different copy numbers for the same motif 
  and flanks, then that is specified in the name of the contig. For example
  the name of a contig can be:

 ```   
        >Block2 TC      5,7     191     205
 ```

 This means that we observed the provided flanks around a STR with either
 5 or 7 copies in the sample. The STR in this contig can be found at 
 [191,205), and reflects the STR with 7 copies ((205-191) / length(TC))
- The program ensures that a flanking region is only selected at most
  once in the dataset. This is done by mantaining a flag which is stored
  for each non-erroneous kmer. The program does not extend a read if it
  encounters a kmers during the extension that has been used in a previous
  extension.

## TEST-DATASET
A test dataset is provided with the distribution in the `test_data` folder.
The test_data folder has the following files:

- 2 sequence files (Illumina_100_500_1.fq & Illumina_100_500_2.fq) in the FASTQ format where the qualities are encoded in the Illumina format (q + 64). The sequence reads were created using pIRS (http://www.ncbi.nlm.nih.gov/pubmed/22508794) for a artificial genome of 4000 bps.
- Makefile which shows the various steps that need to be run for this pipeline.

So, let us run the steps, one at a time and examine the pipeline.
```
make select_strs
```
This runs the following command: 
```
../bin/select_STR_reads -i -n 3 -f 29  Illumina_100_500_1.fq Illumina_100_500_2.fq > reads.str.fq
```
This module extracts the reads that harbor an STR which has at least 3 copies
and at least 29 bps flanking the STR on both sides. The option `-i` is to inform
the script that the quality values are encoded in the Illumina format. The
output of the above step is a fastq file reads.str.fq, which contains the above
mentioned reads along with some additional information about the STR.
```
make merge_strs
```
This runs the following command
```
../bin/merge_STR_reads 27 reads.str.fq > merged.reads.str.fq
```
This module analyzes the results of the previous step and attempts to merge
reads that support the same STR. The merged reads are now stored in a fastq file
merged.reads.str.fq 
```
make extend_strs
```
This runs the following command
```
../bin/extend_STR_reads 4000 27 merged.reads.str.fq Illumina_100_500_1.fq Illumina_100_500_2.fq > contigs.str.fa
```
This module reads in all the reads in Illumina_100_500_1.fq and
Illumina_100_500_2.fq. It removes all kmers that are only seen once. It then
reads the sequences in merged.reads.str.fq and attempts to extend them on both
sides. The resulting output has one contig for each such STR in
merged.reads.str.fq. The names of the contigs contain a unique identifier, the
motif that is repeated in the STR, the  0-based [start,end) coordinates for the
STR in the contig.

The expected output should contain 3 contigs, each one containing the flanking
region around a motif in the dataset.
