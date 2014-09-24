#include "bloom_filter.h"

extern Bool debug_flag;

static Bitset* NewBitset(const uint64_t size) {
    uint64_t byte_size = size / 8;
    uint bit_index = size % 8;
    if(bit_index != 0){
        byte_size += 1;
    }

    Bitset* bs = CkalloczOrDie(sizeof(Bitset));
    bs->size = size;
    bs->byte_size = byte_size;
    bs->bits = CkalloczOrDie(byte_size);
    
    return bs;
}

static Bool SetBit(Bitset* const bs, const uint64_t index) {
    uint64_t byte_index = index / 8;
    uint bit_index = index % 8;
    
    // the index should be within the boundaries of the bitset
    ForceAssert(index < bs->size);
    
    uchar* byte = bs->bits + byte_index;
    uchar mask = '1' - '0';
    mask = mask << 7;
    mask = mask >> bit_index;
    
    Bool is_set = FALSE;
    if(((*byte) & mask) != 0){
        is_set = TRUE;
    }

    *byte = (*byte) | mask;
    return is_set;
}

static Bool CheckBit(Bitset* const bs, const uint64_t index) {
    uint64_t byte_index = index / 8;
    uint bit_index = index % 8;
    
    // the index should be within the boundaries of the bitset
    ForceAssert(index < bs->size);
    
    uchar* byte = bs->bits + byte_index;
    uchar mask = '1' - '0';
    mask = mask << 7;
    mask = mask >> bit_index;
 
    return (*byte) & mask;  
}

//static const char* byte_to_binary(char x) {
//    static char b[9];
//    b[0] = '\0';
//
//    int z;
//    for (z = 128; z > 0; z >>= 1)
//    {
//        strcat(b, ((x & z) == z) ? "1" : "0");
//    }
//
//    return b;
//}

//static void PrintBitset(const Bitset* const bs) {
//    printf("BitSet: ");
//    
//    uint i;
//    for(i = 0; i < bs->byte_size; i++){
//        const char* binaryRep = byte_to_binary(bs->bits[i]);
//        printf("%s", binaryRep);
//    }
//    printf("\n");    
//}

BloomFilter* NewBloomFilter(const float false_positive_rate,
                            const uint64_t num_expected_entries,
                            const uint64_t memory_available) {
    pre(false_positive_rate < 1.0);
    static uint index = 1;

    // number of bits that I will require assuming an optimal k (which is 7 for
    // m, in most cases)
    uint64_t num_bits;

    // if the user specifies the amount of memory available for this bloom
    // filter, then I will respect that. If not then I will create a bloom
    // filter using the expected false-positive rate
    if (memory_available != 0) {
        // the memory is expected to be given in MB (megabytes)
        num_bits = MB_to_Bits * memory_available;
    } else {
        num_bits = (num_expected_entries * -1.0 * log(false_positive_rate)) 
                 / (log(2) * log(2));
    }
    num_bits += 1;
    uint num_hash_functions = num_bits * log(2) / num_expected_entries;

    // the false-positive rate
    float fdr;     
    if (memory_available != 0) {
        fdr = exp((-1.0 * num_bits * log(2) * log(2)) / num_expected_entries);
    } else {
        fdr = false_positive_rate;
    }

    // lets create the bloom filter
    srand((unsigned)time(NULL));
    
    BloomFilter* bf = CkalloczOrDie(sizeof(BloomFilter));
    bf->seed = rand();
    bf->false_positive_rate = fdr;
    bf->num_hash_functions = num_hash_functions;
    bf->num_bits = num_bits;
    bf->bs = NewBitset(bf->num_bits);

    PrintDebugMessage("Bloom filter%u:", index);
    PrintDebugMessage("\tNumber of expected entries: %"PRIu64, 
    num_expected_entries);
    PrintDebugMessage("\tFalse positive rate: %2.6f", fdr);
    PrintDebugMessage("\tNumber of bits used: %"PRIu64, num_bits);
    PrintDebugMessage("\tNumber of hash functions used: %d\n", 
    num_hash_functions);
    return bf;
}

void AddKmerToBloomFilter(BloomFilter* const bf, const Kmer word) {
    uint128_t last_hash = bf->seed;
    uint idx;
    uint128_t hash;
    Bool is_added = FALSE;
    for (idx = 0; idx < bf->num_hash_functions; idx++){
        if (MurmurHash3_128((char*)&word, 
                            sizeof(Kmer), 
                            last_hash, 
                            &hash) == FALSE) {
            PrintThenDie("could not ascertain hash for kmer");
        }
        if(SetBit(bf->bs, hash % bf->num_bits) == FALSE){
            assert(CheckBit(bf->bs, hash % bf->num_bits) > 0);
            bf->num_set_bits += 1;
            is_added = TRUE;
        }
        last_hash = hash;
    }

    if (is_added == TRUE) {
        bf->num_entries_added += 1;
    }
}

Bool CheckKmerInBloomFilter(BloomFilter* const bf, const Kmer word) {
    uint128_t last_hash = bf->seed;
    uint idx;
    uint128_t hash;
    for (idx = 0; idx < bf->num_hash_functions; idx++) {
        if (MurmurHash3_128((char*)&word, 
                            sizeof(Kmer), 
                            last_hash, 
                            &hash) == FALSE) {
            PrintThenDie("could not ascertain hash for kmer");
        }
        if (CheckBit(bf->bs, hash % bf->num_bits) == 0) { 
            return FALSE;
        }
        last_hash = hash;
    }
    return TRUE;    
}

void PrintStatsForBloomFilter(const BloomFilter* const bf) {
    fprintf(stderr, "\nBloom filter stats:\n");
    fprintf(stderr, "\tFalse positive rate: %2.6f\n", bf->false_positive_rate);
    fprintf(stderr, "\tNumber of bits used: %"PRIu64"\n", bf->num_bits);
    fprintf(stderr, "\tNumber of bits set: %"PRIu64"(%2.2f%%)\n", 
    bf->num_set_bits, bf->num_set_bits * 100.0 / bf->num_bits);
    fprintf(stderr, "\tNumber of hash functions used: %d\n", 
    bf->num_hash_functions);
    fprintf(stderr, "\n");
}

void FreeBloomFilter(BloomFilter** pbf) {
    BloomFilter* bf = *pbf;
    Ckfree(bf->bs->bits);
    Ckfree(bf->bs);
    Ckfree(bf);
}
