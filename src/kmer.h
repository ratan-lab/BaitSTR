#ifndef Kmer_H_
#define Kmer_H_

#include <stdlib.h>
#include <stdio.h>   // io functions
#include <stdarg.h>  // va_list and other structures
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>

#include "utilities.h"

// datatype to store the Kmer
#ifdef Large
    typedef __uint128_t Kmer;
#else
    typedef uint64_t Kmer;
#endif

// datatype to store the count of the Kmer as well as the information whether
// the Kmer has been seen before or not
typedef struct Kcount_st {
    uint8_t count;
    uint8_t flag;
} Kcount;

// return the first Kmer from this sequence
Kmer BuildIndex(const char* const s, const uint klen);

// update and return the next Kmer in the sequence
Kmer GetNextKmer(const Kmer word,
                 const char* const s,
                 const uint klen,
                 const int i);

// return the reverse complement of the given Kmer of length klen
Kmer ReverseComplementKmer(const Kmer word, const uint klen);

// print the Kmer in human-readable format. The Kmer is written into the buffer
// provided to this function. It is the duty of the caller to check to make sure
// that sufficient space has been allocated for the buffer. It should be of size
// buffer + 1 if you would like to print it as a string later
void ConvertKmerToString(const Kmer word, const uint klen, char** buffer);

// convert the human readable format into the binary representation of the kmer.
Kmer ConvertStringToKmer(const char* const kmer_string, const uint klen);

#endif  // KMER_H_
