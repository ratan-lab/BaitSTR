#ifndef MURMURHASH_H_
#define MURMURHASH_H_

#include <stdint.h>
#include <stdlib.h>

#include "utilities.h"

Bool MurmurHash3_32(const void *data,
                    const size_t nbytes,
                    const uint32_t seed,
                    void* retbuf);

Bool MurmurHash3_128(const void *data,
                     const size_t nbytes,
                     const uint32_t seed,
                     void *retbuf);

#endif  // MURMURHASH_H_
