#ifndef BLOOM_FILTER_H_
#define BLOOM_FILTER_H_

#include <inttypes.h>
#include <math.h>
#include <time.h>

#include "utilities.h"
#include "murmur_hash.h"
#include "kmer.h"

#define MB_to_Bits 8388608

typedef struct Bitset_st
{
    uint64_t size;
    uint64_t byte_size;
    uchar* bits;
}Bitset;

typedef struct BloomFilter_st
{
    struct BloomFilter_st* next;
    uint seed;
    float false_positive_rate;
    uint num_hash_functions;
    uint64_t num_bits;
    uint64_t num_set_bits;
    uint64_t num_entries_added;
    Bitset* bs;
}BloomFilter;

BloomFilter* NewBloomFilter(const float false_positive_rate,
                            const uint64_t num_expected_entries,
                            const uint64_t memory_available);

void AddKmerToBloomFilter(BloomFilter* const bf, const Kmer word);

Bool CheckKmerInBloomFilter(BloomFilter* const bf, const Kmer word);

void PrintStatsForBloomFilter(const BloomFilter* const bf);

void FreeBloomFilter(BloomFilter** pbf);

#endif
