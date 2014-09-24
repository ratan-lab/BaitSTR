#ifndef FASTQ_H_
#define FASTQ_H_

#include <zlib.h>

#include "utilities.h"

typedef struct FastqSequence_st {
    gzFile fd;

    char* name;          // name of the read
    size_t name_length;  // maximum length of the name seen so far  

    char* bases;          // the actual bases in the sequence   
    size_t bases_length;  // the maximum length of the sequence seen so far

    char* quals;          // the quality values encoded
    size_t quals_length;  // the maximum length of the sequence seen so far

    size_t slen;  // length of the sequence and quals after trimming and
                  // everything else has been done
    Bool is_illumina_encoded;  // TRUE if the quality encoding is q + 64
    Bool do_trim; // should I trim the sequences at the 3' end to ignore
                  // stretches of quality value 2...
}FastqSequence;

// open the fastq file and return the first FastqSequence
FastqSequence* ReadFastqSequence(const char* const file,
                                 const Bool is_quality_illumina_encoded,
                                 const Bool do_trim_reads);

// print the fastq FastqSequence
void PrintFastqSequence(const FastqSequence* const sp);

// get the next fastq FastqSequence to the FastqSequence sp
FastqSequence* GetNextSequence(FastqSequence* const sp);

// free all the resources used by this fastq FastqSequence
void CloseFastqSequence(FastqSequence* sp);

#endif  // FASTQ_H_
