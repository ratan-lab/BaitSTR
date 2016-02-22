#ifndef SPARSE_KMER_HASH_H_
#define SPARSE_KMER_HASH_H_

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <sparsehash/sparse_hash_map>

#include <iostream>

extern "C" {
#include "kmer.h"
#include "murmur_hash.h"
}

using google::sparse_hash_map;      // namespace where class lives by default
using std::cout;
using std::endl;
//using __gnu_cxx::hash;
//using tr1::hash
//using ext::hash; 

// simple hash adapter for types without pointers
template<typename T>
struct SparseMurmurHasher {
    size_t operator()(const T& t) const{
        uint128_t hashVal;
        if (MurmurHash3_128((char*)&t, sizeof(t), 0, &hashVal) == false)
            PrintThenDie("error in determining the hash");
        return hashVal;
    }
};

struct SparseEqKmer {
    Bool operator()(const Kmer s1, const Kmer s2) const{
        return s1 == s2;
    }
};

// A file serialized with this should start with
// 24687531...
// when I use
// xxd -p file | head
struct SparseKmerSerializer {
    Bool operator()(FILE* fp, const std::pair<const Kmer, 
                                              Kcount> value) const {
        // write the key 
        if (fwrite(&value.first, sizeof(Kmer), 1, fp) != 1) {
            return FALSE;
        }

        // write the value
        if (fwrite(&value.second, sizeof(Kcount), 1, fp) != 1) {
            return FALSE;
        }
        return TRUE;
    }

    Bool operator()(FILE* fp, std::pair<const Kmer,
                                        Kcount>* value) const {
        // read the key.
        if (fread(const_cast<Kmer*>(&value->first), sizeof(Kmer), 1, fp) != 1) {
            return FALSE;
        }
        // at this point *((Kmer*)value->first) should be equal to the kmer
        // value

        // read the value
        Kcount cnt;
        if (fread(&cnt, sizeof(Kcount), 1, fp) != 1) {
            return FALSE;
        }
        value->second = cnt;

        return TRUE;
    }
};

typedef sparse_hash_map<Kmer,Kcount,SparseMurmurHasher<Kmer>, SparseEqKmer> SparseHashMap;

Bool CheckKmerInSparseHashMap(SparseHashMap& kmers,
                              const Kmer kmer) {
    SparseHashMap::iterator it;
    it = kmers.find(kmer);

    if (it == kmers.end()) {
        return FALSE;
    }
    return TRUE;
}

#endif  // SPARSE_KMER_HASH_H_
