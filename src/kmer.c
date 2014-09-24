#include "kmer.h"

// return the first Kmer from this sequence
Kmer BuildIndex(const char* const s, const uint klen) {
#ifdef Large
    assert(klen < 64);
    if (klen < 32) {
        PrintThenDie("You should use a different build to handle kmers of length less than 32");
    }
#else
    assert(klen < 32);
#endif
    Kmer word = 0;
    uint i;
    for (i = 0; i < (klen - 1); i++) {
            word = (word << 2);
            word += fasta_encoding[(int)s[i]];
    }

    return word;
}

// update and return the next Kmer in the sequence
Kmer GetNextKmer(const Kmer word,
                 const char* const s,
                 const uint klen,
                 const int i) {
    Kmer next = (word << ((8 * sizeof(Kmer)) - 2*klen + 2));
    next = (next >> ((8 * sizeof(Kmer)) - 2*klen));
    next += fasta_encoding[(int)s[i+klen-1]];
    return next;
}

// return the reverse complement of the given Kmer of length klen
Kmer ReverseComplementKmer(const Kmer word, const uint klen) {
    Kmer copy = word;
    Kmer rc = 0;

    uint i = 0;
    for (i = 0; i < klen; i++) {
            rc = (rc << 2);
            rc += 3 - (copy & 3);
            copy = (copy >> 2);
    }

    return rc;
}

// print the Kmer in human-readable format. The Kmer is written into the buffer
// provided to this function. It is the duty of the caller to check & make sure
// that sufficient space has been allocated for the buffer. It should be of the
// size buffer + 1 if you would like to print it as a string later
void ConvertKmerToString(const Kmer name, const uint klen, char** pbuffer) {
    char* buffer = *pbuffer;
    Kmer foo = name;
    int i = klen;
    while (i > 0) {
            buffer[--i] = bit_encoding[foo & 3];
            foo = foo >> 2;
    }
}

// convert the human readable format into the binary representation of the kmer.
Kmer ConvertStringToKmer(const char* const kmer_string, const uint klen) {
    Kmer kmer;
    kmer = BuildIndex(kmer_string, klen);
    kmer = GetNextKmer(kmer, kmer_string, klen, 0);
    return kmer;
}
