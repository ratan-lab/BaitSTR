#include "murmur_hash.h"

// MurmurHash3 was created by Austin Appleby  in 2008. The cannonical
// implementations are in C++ and placed in the public.
//
//    https://sites.google.com/site/murmurhash/
//
//  Seungyoung Kim has ported it's cannonical implementation to C language
//  in 2012 and published it as a part of qLibc component. This implementation
//  is heavily borrowed from that.

Bool MurmurHash3_32(const void *data,
                    const size_t nbytes,
                    const uint32_t seed,
                    void* retbuf) {
    if (data == NULL || nbytes == 0) return FALSE;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    const int nblocks = nbytes / 4;
    const uint32_t *blocks = (const uint32_t *)(data);
    const uint8_t *tail = (const uint8_t *)(data + (nblocks * 4));

    uint32_t h = seed;

    int i;
    uint32_t k;
    for (i = 0; i < nblocks; i++) {
        k = blocks[i];

        k *= c1;
        k = (k << 15) | (k >> (32 - 15));
        k *= c2;

        h ^= k;
        h = (h << 13) | (h >> (32 - 13));
        h = (h * 5) + 0xe6546b64;
    }

    k = 0;
    switch (nbytes & 3) {
        case 3:
            k ^= tail[2] << 16;
        case 2:
            k ^= tail[1] << 8;
        case 1:
            k ^= tail[0];
            k *= c1;
            k = (k << 13) | (k >> (32 - 15));
            k *= c2;
            h ^= k;
    }

    h ^= nbytes;

    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;

    ((uint32_t *)retbuf)[0] = h;

    return TRUE;
}

Bool MurmurHash3_128(const void *data,
                     const size_t nbytes,
                     const uint32_t seed,
                     void *retbuf) {
    if (data == NULL || nbytes == 0) return FALSE;

    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;

    const int nblocks = nbytes / 16;
    const uint64_t *blocks = (const uint64_t *)(data);
    const uint8_t *tail = (const uint8_t *)(data + (nblocks * 16));

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    int i;
    uint64_t k1, k2;
    for (i = 0; i < nblocks; i++) {
        k1 = blocks[i*2+0];
        k2 = blocks[i*2+1];

        k1 *= c1;
        k1  = (k1 << 31) | (k1 >> (64 - 31));
        k1 *= c2;
        h1 ^= k1;

        h1 = (h1 << 27) | (h1 >> (64 - 27));
        h1 += h2;
        h1 = h1*5+0x52dce729;

        k2 *= c2;
        k2  = (k2 << 33) | (k2 >> (64 - 33));
        k2 *= c1;
        h2 ^= k2;

        h2 = (h2 << 31) | (h2 >> (64 - 31));
        h2 += h1;
        h2 = h2*5+0x38495ab5;
    }

    k1 = k2 = 0;
    switch (nbytes & 15) {
        case 15:
            k2 ^= (uint64_t)(tail[14]) << 48;
        case 14:
            k2 ^= (uint64_t)(tail[13]) << 40;
        case 13:
            k2 ^= (uint64_t)(tail[12]) << 32;
        case 12:
            k2 ^= (uint64_t)(tail[11]) << 24;
        case 11:
            k2 ^= (uint64_t)(tail[10]) << 16;
        case 10:
            k2 ^= (uint64_t)(tail[ 9]) << 8;
        case  9:
            k2 ^= (uint64_t)(tail[ 8]) << 0;
            k2 *= c2;
            k2  = (k2 << 33) | (k2 >> (64 - 33));
            k2 *= c1;
            h2 ^= k2;
        case  8:
            k1 ^= (uint64_t)(tail[ 7]) << 56;
        case  7:
            k1 ^= (uint64_t)(tail[ 6]) << 48;
        case  6:
            k1 ^= (uint64_t)(tail[ 5]) << 40;
        case  5:
            k1 ^= (uint64_t)(tail[ 4]) << 32;
        case  4:
            k1 ^= (uint64_t)(tail[ 3]) << 24;
        case  3:
            k1 ^= (uint64_t)(tail[ 2]) << 16;
        case  2:
            k1 ^= (uint64_t)(tail[ 1]) << 8;
        case  1:
            k1 ^= (uint64_t)(tail[ 0]) << 0;
            k1 *= c1;
            k1  = (k1 << 31) | (k1 >> (64 - 31));
            k1 *= c2;
            h1 ^= k1;
    }

    //----------
    // finalization

    h1 ^= nbytes;
    h2 ^= nbytes;

    h1 += h2;
    h2 += h1;

    h1 ^= h1 >> 33;
    h1 *= 0xff51afd7ed558ccdULL;
    h1 ^= h1 >> 33;
    h1 *= 0xc4ceb9fe1a85ec53ULL;
    h1 ^= h1 >> 33;

    h2 ^= h2 >> 33;
    h2 *= 0xff51afd7ed558ccdULL;
    h2 ^= h2 >> 33;
    h2 *= 0xc4ceb9fe1a85ec53ULL;
    h2 ^= h2 >> 33;

    h1 += h2;
    h2 += h1;

    ((uint64_t *)retbuf)[0] = h1;
    ((uint64_t *)retbuf)[1] = h2;

    return TRUE;
}
