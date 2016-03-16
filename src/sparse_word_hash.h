#ifndef SPARSE_WORD_HASH_H_
#define SPARSE_WORD_HASH_H_

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <functional>
#include <sparsehash/sparse_hash_map>

#include <iostream>

extern "C" {
#include "murmur_hash.h"
}

using google::sparse_hash_map;      // namespace where class lives by default
using std::cout;
using std::endl;

typedef struct Copies_st
{
    struct Copies_st* next;
    uint16_t copies;
    uint16_t nsupport;
} Copies;

typedef struct Block_st
{
    struct Block_st* next;
    uint16_t zstart;
    uint16_t end;
    uint16_t slen;
    uint16_t support;
    Copies* supports;
    char* seq;
    char* qual;
} Block;

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
  }
};

struct MurmurHasher {
    size_t operator()(const char* t) const {
        uint128_t hashVal;
        if (MurmurHash3_128(t, strlen(t), 0, &hashVal) == false)
            PrintThenDie("error in determining the hash");
        return hashVal;
    }
};


typedef sparse_hash_map<const char* , Block*, MurmurHasher, eqstr> SparseWordHashMap;

Bool CheckInSparseWordHashMap(SparseWordHashMap& blocks,
                         const char* key) 
{
    SparseWordHashMap::iterator it;
    it = blocks.find(key);

    if (it == blocks.end()) {
        return FALSE;
    }
    return TRUE; 
}

#endif  // SPARSE_WORD_HASH_H_
