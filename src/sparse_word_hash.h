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

typedef struct Block_st
{
    uint zstart;
    uint end;
    uint slen;
    uint support;
    char* seq;
    char* qual;
    int copies[3];
    int nsupport[3];
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
