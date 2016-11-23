extern "C" {
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include "inttypes.h"

#include "utilities.h"
#include "clparsing.h"
#include "kmer.h"
#include "fastq_seq.h"
#include "bloom_filter.h"
}

#include "sparse_kmer_hash.h"

char* program_version       = "";
char* program_name          = "extend_STR_reads";
char* program_description   = 
    "Extend fastq reads based on the kmer structure";
char* program_use           = 
    "extend_STR_reads [options] gs cov klen str.reads.fq reads1.fq reads2.fq ...";

Bool debug_flag;

// generic buffer to write Kmers
char* kmer_buffer = NULL;

// the maximum extension
uint flank_chunk = 1024;

static void ReadAndCountNonSingletonKmers(SparseHashMap& kmers,
                                          const uint64_t num_expected_kmers,
                                          const uint kmer_length,
                                          char** const argv,
                                          const uint nameidx,
                                          const uint progress_chunk,
                                          const uint min_threshold,
                                          const uint max_threshold) {

    // all the singleton kmers shall be stored here.
    BloomFilter* singletons = NewBloomFilter(0.1, num_expected_kmers, 0);  

    // read the kmers the first time and identify kmers that might be present
    // more than once.
    int idx;
    static uint64_t num_hash_entries = 0;
    ReportMemoryUsage();
    for (idx = 4; idx < nameidx; idx++) {
        // read the kmers from this fastq file. 
        FastqSequence* sequence = ReadFastqSequence(argv[idx], FALSE, FALSE);
        Kmer word, antiword, stored;
        SparseHashMap::iterator it;
        Bool already_in_hash;
    
        uint64_t num_kmers_added = 0;
        uint64_t num_sequence_processed = 0;
    
        while (sequence) {
            if (sequence->slen >= kmer_length) {
                if (debug_flag == TRUE) {
                    PrintDebugMessage("1. Processing %s", sequence->name + 1);
                }
                // print progress
                num_sequence_processed += 1;
                if ((num_sequence_processed - 1) % progress_chunk == 0) {
                    PrintDebugMessage("1. Processing read number %"PRIu64": %s",
                    num_sequence_processed, sequence->name + 1);
                }
    
                // a load factor greater than 0.7-0.8 is a sign that the user did
                // not select the expected number of kmers judiciously. Lets warn
                // the user, as increasing the size of the hashtable can be very
                // slow.
                if (kmers.load_factor() > 0.8) {
                    PrintWarning("Current load factor: %2.6f",  
                    kmers.load_factor());
                    PrintWarning(
                    "Try increasing expected number of kmers from %"PRIu64, 
                    num_expected_kmers);
                }
    
                // let account for all the kmers in this sequence
                word = BuildIndex(sequence->bases, kmer_length);
                uint num_kmers = strlen(sequence->bases) - kmer_length + 1;

                for (uint i = 0; i < num_kmers; i++) {
                    word = GetNextKmer(word, sequence->bases, kmer_length, i);
                    
                    antiword = ReverseComplementKmer(word, kmer_length);
                    stored = word < antiword ? word : antiword;
    
                    already_in_hash = CheckKmerInSparseHashMap(kmers, stored);

                    if (already_in_hash == FALSE) {
                        if (CheckKmerInBloomFilter(singletons, stored) == TRUE) {
                            // this kmer has already been seen once, so add K 
                            // to the hashtable
                            kmers[stored].count = 0;
                            kmers[stored].flag = 0;
                            if (debug_flag == TRUE) {
                                ConvertKmerToString(word, 
                                                    kmer_length, 
                                                    &kmer_buffer);
                                PrintDebugMessage("[[ %d ]] 1. Adding kmer %s",
                                              num_kmers_added, kmer_buffer);
                            }
                        } else {
                            // add it only to the bloom filter
                            AddKmerToBloomFilter(singletons, stored);
                        }
                    }
                }
            }
    
            sequence = GetNextSequence(sequence);
        }
        CloseFastqSequence(sequence);        
        PrintDebugMessage("1. Done with all the sequences in %s", argv[idx]);
        PrintDebugMessage("1. Counted %zu different kmers", kmers.size());
        ReportMemoryUsage();
    }       

    // I am done with the bloom filter.
    PrintStatsForBloomFilter(singletons);
    FreeBloomFilter(&singletons);

    // lets iterate through the kmers once more and remove the false positives
    for (idx = 4; idx < nameidx; idx++) {
        // read the kmers from this fastq file. 
        FastqSequence* sequence = ReadFastqSequence(argv[idx], FALSE, FALSE);
        Kmer word, antiword, stored;
        SparseHashMap::iterator it;
        Bool already_in_hash;
    
        uint64_t num_kmers_added = 0;
        uint64_t num_sequence_processed = 0;
    
        while (sequence) {
            if (sequence->slen >= kmer_length) {
                if (debug_flag == TRUE) {
                    PrintDebugMessage("2. Processing %s", sequence->name + 1);
                }
                // print progress
                num_sequence_processed += 1;
                if ((num_sequence_processed - 1) % progress_chunk == 0) {
                    PrintDebugMessage("2. Processing read number %"PRIu64": %s",
                    num_sequence_processed, sequence->name + 1);
                }
    
                // let account for all the kmers in this sequence
                word = BuildIndex(sequence->bases, kmer_length);
                uint num_kmers = strlen(sequence->bases) - kmer_length + 1;

                for (uint i = 0; i < num_kmers; i++) {
                    word = GetNextKmer(word, sequence->bases, kmer_length, i);
                    antiword = ReverseComplementKmer(word, kmer_length);
                    stored = word < antiword ? word : antiword;
    
                    already_in_hash = CheckKmerInSparseHashMap(kmers, stored);

                    if (already_in_hash == TRUE) {
                            uint8_t old_kcnt = kmers[stored].count;
                            if (old_kcnt <= (umaxof(uint8_t) - 1)) {
                                kmers[stored].count += 1;
                            } else {
                                kmers[stored].count = umaxof(Kcount);
                            }
                            if (debug_flag == TRUE) {
                                ConvertKmerToString(word, 
                                                    kmer_length, 
                                                    &kmer_buffer);
                                PrintDebugMessage("[[ %d ]] 2. Incrementing kmer %s count to %d", num_kmers_added, kmer_buffer, kmers[stored]);
                            }
                        }
                    }
                }
    
            sequence = GetNextSequence(sequence);
        }
        CloseFastqSequence(sequence);        
        PrintDebugMessage("2. Done with all the sequences in %s", argv[idx]);
        ReportMemoryUsage();
    }   

    // go through and mark kmers as deleted if they occur less than a number of
    // times 
    SparseHashMap::iterator it;
    for (it = kmers.begin(); it != kmers.end(); it++) {
        if (((*it).second.count < min_threshold) || ((*it).second.count > max_threshold))  {
            kmers.erase(it);
        }
    }   
    kmers.resize(0);
}

static void RvKmers(Kmer word, 
                    Kmer*& rvs, 
                    const uint kmer_length) {
    Kmer indx = 0;
    for (indx = 0; indx < 4; indx++) {
        Kmer right_shifted = word >> 2;
        Kmer to_add = (indx <<  (2*kmer_length - 2));
        rvs[indx] = right_shifted + to_add;
    }    
}

static void FwKmers(Kmer word,
                    Kmer*& fws,
                    const uint kmer_length) {
    Kmer indx = 0;
    for (indx = 0; indx < 4; indx++) {
        Kmer left_shifted = (word << ((8 * sizeof(Kmer)) - 2*kmer_length + 2));
        Kmer right_shifted = (left_shifted >>((8 * sizeof(Kmer))-2*kmer_length));
        right_shifted += indx;
        fws[indx] = right_shifted;
    }
}

static Bool CheckForSNPBackwards(const Kmer kmer,
                                 SparseHashMap& kmers,
                                 const uint kmer_length) 
{
    uint indx;
    Bool issnp = FALSE;

    Kmer* rvs = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));
    Kmer* fws = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));

    Kmer curr = kmer;
    Kmer rev  = ReverseComplementKmer(curr, kmer_length);

    RvKmers(curr, rvs, kmer_length);
    FwKmers(rev, fws, kmer_length);
    //ConvertKmerToString(curr, kmer_length, &kmer_buffer); 
    //printf("CheckForSNP: %s\n", kmer_buffer);

    uint num_extensions = 0;
    Kmer* extensions = (Kmer*)CkalloczOrDie(3 * sizeof(Kmer));

    for (indx = 0; indx < 4; indx++) {
        if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
            if (num_extensions < 3) {
                extensions[num_extensions++] = rvs[indx];
            }
            //ConvertKmerToString(rvs[indx], kmer_length, &kmer_buffer);
            //printf("Possible extension: %s\n", kmer_buffer);
        }
        if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
            if (num_extensions < 3) {
                extensions[num_extensions++] = ReverseComplementKmer(fws[indx], kmer_length);
            }
            //ConvertKmerToString(fws[indx], kmer_length, &kmer_buffer);
            //printf("Possible extension: %s\n", kmer_buffer);
        }
    }
    ForceAssert(num_extensions == 2);

    uint numsteps;
    Kmer extension, extension1, extension2;

    // extend the first candidate
    numsteps = 0;
    curr = extensions[0];
    while (numsteps <= kmer_length) {
        rev = ReverseComplementKmer(curr, kmer_length);
        //ConvertKmerToString(curr, kmer_length, &kmer_buffer);
        //printf("Extension 1: %s\n", kmer_buffer);
        RvKmers(curr, rvs, kmer_length);
        FwKmers(rev, fws, kmer_length);

        num_extensions = 0;
        for (indx = 0; indx < 4; indx++) {
            if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
                num_extensions++;
                extension = rvs[indx];
            }
            if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
                num_extensions++;
                extension = ReverseComplementKmer(fws[indx], kmer_length);
            }
        }  
        if (num_extensions != 1) break;
        curr = extension;
        numsteps += 1;
    }
    if (numsteps != (kmer_length+1)) goto backwardsdecide;
    extension1 = curr;
    
    // extend the second candidate
    numsteps = 0;
    curr = extensions[1];
    while (numsteps <= kmer_length) {
        rev = ReverseComplementKmer(curr, kmer_length);
        //ConvertKmerToString(curr, kmer_length, &kmer_buffer);
        //printf("Extension 2: %s\n", kmer_buffer);
        RvKmers(curr, rvs, kmer_length);
        FwKmers(rev, fws, kmer_length);

        num_extensions = 0;
        for (indx = 0; indx < 4; indx++) {
            if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
                num_extensions++;
                extension = rvs[indx];
            }
            if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
                num_extensions++;
                extension = ReverseComplementKmer(fws[indx], kmer_length);
            }
        }  
        if (num_extensions != 1) break;
        curr = extension;
        numsteps += 1;
    }
    if (numsteps != (kmer_length+1)) goto backwardsdecide;
    extension2 = curr;
 
    if (extension1 == extension2) issnp = TRUE;    

backwardsdecide:
    Ckfree(rvs);
    Ckfree(fws);
    Ckfree(extensions);

    return issnp;
}

static Bool CheckForSNPForwards(const Kmer kmer,
                                SparseHashMap& kmers,
                                const uint kmer_length) 
{
    uint indx;
    Bool issnp = FALSE;

    Kmer* rvs = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));
    Kmer* fws = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));

    Kmer curr = kmer;
    Kmer rev  = ReverseComplementKmer(curr, kmer_length);

    FwKmers(curr, fws, kmer_length);
    RvKmers(rev, rvs, kmer_length);
    //ConvertKmerToString(curr, kmer_length, &kmer_buffer); 
    //printf("CheckForSNP: %s\n", kmer_buffer);

    uint num_extensions = 0;
    Kmer* extensions = (Kmer*)CkalloczOrDie(3 * sizeof(Kmer));

    for (indx = 0; indx < 4; indx++) {
        if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
            if (num_extensions < 3) {
                extensions[num_extensions++] = fws[indx];
            }
            //ConvertKmerToString(extensions[num_extensions-1], kmer_length, &kmer_buffer);
            //printf("Possible extension: %s\n", kmer_buffer);
        }
        if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
            if (num_extensions < 3) {
                extensions[num_extensions++] = ReverseComplementKmer(rvs[indx],
kmer_length);
            }
            //ConvertKmerToString(extensions[num_extensions-1], kmer_length, &kmer_buffer);
            //printf("Possible extension: %s\n", kmer_buffer);
        }
    }
    ForceAssert(num_extensions == 2);

    uint numsteps;
    Kmer extension, extension1, extension2;

    // extend the first candidate
    numsteps = 0;
    curr = extensions[0];
    while (numsteps <= kmer_length) {
        rev = ReverseComplementKmer(curr, kmer_length);
        //ConvertKmerToString(curr, kmer_length, &kmer_buffer);
        //printf("Extension 1: %s\n", kmer_buffer);
        FwKmers(curr, fws, kmer_length);
        RvKmers(rev, rvs, kmer_length);

        num_extensions = 0;
        for (indx = 0; indx < 4; indx++) {
            if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
                num_extensions++;
                extension = ReverseComplementKmer(rvs[indx], kmer_length);
            }
            if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
                num_extensions++;
                extension = fws[indx];
            }
        }  
        if (num_extensions != 1) {
            //printf("Number of extensions possible: %d\n", num_extensions);
            break;
        }
        curr = extension;
        numsteps += 1;
    }
    if (numsteps != (kmer_length+1)) goto forwardsdecide;
    extension1 = curr;
 
    // extend the second candidate   
    numsteps = 0;
    curr = extensions[1];
    while (numsteps <= kmer_length) {
        rev = ReverseComplementKmer(curr, kmer_length);
        //ConvertKmerToString(curr, kmer_length, &kmer_buffer);
        //printf("Extension 2: %s\n", kmer_buffer);
        FwKmers(curr, fws, kmer_length);
        RvKmers(rev, rvs, kmer_length);

        num_extensions = 0;
        for (indx = 0; indx < 4; indx++) {
            if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
                num_extensions++;
                extension = ReverseComplementKmer(rvs[indx], kmer_length);
            }
            if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
                num_extensions++;
                extension = fws[indx];
            }
        }  
        if (num_extensions != 1) break;
        curr = extension;
        numsteps += 1;
    }
    if (numsteps != (kmer_length+1)) goto forwardsdecide;
    extension2 = curr;
 
    if (extension1 == extension2) issnp = TRUE;    

forwardsdecide:
    Ckfree(rvs);
    Ckfree(fws);
    Ckfree(extensions);

    return issnp;
}

/*
 * Simple Algorithm 
 * ----------------
 *  rvkmers = [kmer]
 *
 *  curr = kmer
 *  while True:
 *      rev = ReverseComplement(curr)
 *
 *      extensions = [x for x in RvKmer(curr) if x in kmersdict]
 *      extensions.extend([x for x in FwKmer(rev) if x in kmersdict])
 *
 *      # Quit if there is more than one extension possible
 *      if len(extensions) != 1:
 *          break
 *
 *      # Quit if the extension has already been seen in this contig.
 *      rev_extension = ReverseComplement(extensions[0])
 *      if extensions[0] in rvkmers or rev_extension in rvkmers:
 *          break
 *
 *      if extensions[0] in RvKmer(curr):
 *          curr = extensions[0]
 *      else:
 *          assert extensions[0] in FwKmer(rev)
 *          curr = rev_extension
 *  
 *      rvkmers.append(curr)
 *  return rvkmers
 * 
 * We return at most flank_chunk bases on each end of the STR. 
 */
static char* ExtendBackward(SparseHashMap& kmers, 
                            const char* const bases, 
                            const uint motif_zstart, 
                            const uint kmer_length) {
    uint rvkmers_used = 0;
    uint rvkmers_allocated = flank_chunk;

    #ifdef INFEXPAND
    Kmer* rvkmers = (Kmer*)CkallocOrDie(rvkmers_allocated * sizeof(Kmer));
    #else
    Kmer rvkmers[flank_chunk];
    #endif

    Kmer kmer = ConvertStringToKmer(bases+motif_zstart-kmer_length, kmer_length);
    rvkmers[rvkmers_used++] = kmer;
    
    Bool rvflag;
    Bool already_seen = FALSE;
    Kmer curr, rev, extension, rev_extension;
    uint indx;
    Kmer* rvs = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));
    Kmer* fws = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));

    curr = kmer;
    while (TRUE) {
        rev = ReverseComplementKmer(curr, kmer_length);

        RvKmers(curr, rvs, kmer_length);
        FwKmers(rev, fws, kmer_length);

        uint num_extensions = 0;
        for (indx = 0; indx < 4; indx++) {
            if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
                extension = rvs[indx];
                num_extensions += 1;
                rvflag = TRUE;
            }
            if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
                extension = fws[indx];
                num_extensions += 1;
                rvflag = FALSE;
            }
        }

        // quit, if there is more than one extension possible
        if (num_extensions == 2) {
            //fprintf(stderr,"Backward: More than one extension possible\n");
            // this could be rescued if this is only a substitution polymorphism
            if (CheckForSNPBackwards(curr, kmers, kmer_length) == FALSE) {
                break;
            }
        } else if (num_extensions == 1) {
        } else if (num_extensions == 0) {
            //fprintf(stderr,"Backward: No extensions possible\n");
            break;
        } else {
            break;
        }
        rev_extension = ReverseComplementKmer(extension, kmer_length);
        
        // quit, if the extension has already been seen in this contig.
        for (indx = 0; indx < rvkmers_used; indx++) {
            if ((extension == rvkmers[indx]) || 
                (rev_extension == rvkmers[indx])) {
                already_seen = TRUE;
                break;
            }
        }
        if (already_seen == TRUE) break;

        // quit, if we have used this kmer for another STR
        ForceAssert(CheckKmerInSparseHashMap(kmers, extension) == TRUE);
        // if (kmers[extension].flag != 0) {
        //     Ckfree(rvs);
        //     Ckfree(fws);
        //     return NULL;
        // }
        // kmers[extension].flag = 1;

        if (rvflag == TRUE) {
            curr = extension;
        } else {
            curr = rev_extension;
        }

        rvkmers[rvkmers_used++] = curr;

        #ifdef INFEXPAND
        if (rvkmers_used == rvkmers_allocated) {
            rvkmers_allocated += flank_chunk;
            rvkmers = (Kmer*)CkreallocOrDie(rvkmers, rvkmers_allocated*sizeof(Kmer));
        }
        #else
        if (rvkmers_used == rvkmers_allocated) break;    
        #endif
    }

    char* lflank = (char*)CkalloczOrDie(kmer_length + rvkmers_used + 1);
    
    for (indx = 0; indx < (rvkmers_used - 1); indx++) {
        ConvertKmerToString(rvkmers[indx], kmer_length, &kmer_buffer);
        lflank[rvkmers_used+kmer_length-indx-2] = kmer_buffer[kmer_length-1];
    }
    ConvertKmerToString(rvkmers[indx], kmer_length, &kmer_buffer);
    for (indx = 0; indx < kmer_length; indx++) {
        lflank[indx] = kmer_buffer[indx];
    }
    
    Ckfree(rvs);
    Ckfree(fws);
    #ifdef INFEXPAND
    Ckfree(rvkmers);
    #endif
        

    return lflank;
}

/*   
 * Simple Algorithm:
 * ----------------   
 *    fwkmers = [kmer]
 *
 *    curr = kmer
 *    while True:
 *        rev = ReverseComplement(curr)
 *
 *        extensions = [x for x in FwKmer(curr) if x in kmersdict]
 *        extensions.extend([x for x in RvKmer(rev) if x in kmersdict])
 *         
 *        # Quit if there is more than one extension possible
 *        if len(extensions) != 1:
 *            break
 *         
 *         # Quit if the extension has already been seen in this contig.
 *        rev_extension = ReverseComplement(extensions[0])
 *        if extensions[0] in fwkmers or rev_extension in fwkmers:
 *            break
 *
 *        if extensions[0] in FwKmer(curr):
 *            curr = extensions[0]
 *        else:
 *            assert extensions[0] in RvKmer(rev)
 *            curr = rev_extension
 *
 *        fwkmers.append(curr)
 *    return fwkmers
 *
 * We only return at most flank_chunk base pairs on either side of the STR.
 */
static char* ExtendForward(SparseHashMap& kmers, 
                           const char* const bases, 
                           const uint motif_end, 
                           const uint kmer_length) {
    uint fwkmers_used = 0;
    uint fwkmers_allocated = flank_chunk;

    #ifdef INFEXPAND
    Kmer* fwkmers = (Kmer*)CkallocOrDie(fwkmers_allocated * sizeof(Kmer));
    #else
    Kmer fwkmers[flank_chunk];
    #endif

    Kmer kmer = ConvertStringToKmer(bases+motif_end, kmer_length);
    fwkmers[fwkmers_used++] = kmer;
    
    Bool rvflag;
    Bool already_seen = FALSE;
    Kmer curr, rev, extension, rev_extension;
    uint indx;
    Kmer* rvs = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));
    Kmer* fws = (Kmer*)CkalloczOrDie(4 * sizeof(Kmer));

    curr = kmer;
    while (TRUE) {
        rev = ReverseComplementKmer(curr, kmer_length);

        FwKmers(curr, fws, kmer_length);
        RvKmers(rev, rvs, kmer_length);

        uint num_extensions = 0;
        for (indx = 0; indx < 4; indx++) {
            if (CheckKmerInSparseHashMap(kmers, rvs[indx]) == TRUE) {
                extension = rvs[indx];
                num_extensions += 1;
                rvflag = TRUE;
            }
            if (CheckKmerInSparseHashMap(kmers, fws[indx]) == TRUE) {
                extension = fws[indx];
                num_extensions += 1;
                rvflag = FALSE;
            }
        }

        // quit, if there is more than one extension possible
        if (num_extensions == 2) {
            //fprintf(stderr,"Forward: More than one extension possible\n");
            // this could be rescued if this is only a substitution polymorphism
            if (CheckForSNPForwards(curr, kmers, kmer_length) == FALSE) {
                break;
            }
        } else if (num_extensions == 1) { 
        } else if (num_extensions == 0) {
            //fprintf(stderr,"Forward: No extensions possible\n");
            break;
        } else {
            break;
        }
        rev_extension = ReverseComplementKmer(extension, kmer_length);
        
        // quit, if the extension has already been seen in this contig.
        for (indx = 0; indx < fwkmers_used; indx++) {
            if ((extension == fwkmers[indx]) || 
                (rev_extension == fwkmers[indx])) {
                already_seen = TRUE;
                break;
            }
        }
        if (already_seen == TRUE) break;

        // quit, if this extension has been used in some other STR
        ForceAssert(CheckKmerInSparseHashMap(kmers, extension) == TRUE);
        // if (kmers[extension].flag != 0) {
        //     Ckfree(rvs);
        //     Ckfree(fws);
        //     return NULL;
        // }
        // kmers[extension].flag = 1; 

        if (rvflag == FALSE) {
            curr = extension;
        } else {
            curr = rev_extension;
        }

        fwkmers[fwkmers_used++] = curr;
        #ifdef INFEXPAND
        if (fwkmers_used == fwkmers_allocated) {
            fwkmers_allocated += flank_chunk;
            fwkmers = (Kmer*)CkreallocOrDie(fwkmers, fwkmers_allocated*sizeof(Kmer));
        }
        #else
        if (fwkmers_used == fwkmers_allocated) break;
        #endif
    }

    char* rflank = (char*)CkalloczOrDie(kmer_length + fwkmers_used + 1);
    
    for (indx = 0; indx < fwkmers_used; indx++) {
        ConvertKmerToString(fwkmers[indx], kmer_length, &kmer_buffer);
        rflank[indx] = kmer_buffer[0];
    }
    for(uint i = 1; i < kmer_length; i++) {
        rflank[indx++] = kmer_buffer[i];
    }
    
    Ckfree(rvs);
    Ckfree(fws);
    #ifdef INFEXPAND
    Ckfree(fwkmers);
    #endif

    return rflank;
}

/* Calculate the percent identity between sequence1, and 
 * a) sequence2[index:] if do_reverse == FALSE
 * b) sequence2[:index] if do_reverse == TRUE
 */
static float PercentIdentity(const char* const sequence1, const char* const
sequence2, const int index, const Bool do_reverse) {
    int indx = 0;
    uint matches = 0;
    uint mismatches = 0;

    if (do_reverse == FALSE) {
        for (indx = 0; 
             indx < MIN(strlen(sequence2)-index, strlen(sequence1)); 
             indx++) {
            if (sequence1[indx] == sequence2[index + indx]) {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    } else {
        int indx1 = 0, indx2 = 0;

        for (indx1 = strlen(sequence1)-1, indx2 = index-1;
             (indx1 >= 0) && (indx2 >= 0);
             indx1--, indx2--) {
            if (sequence1[indx1] == sequence2[indx2]) {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }

    float pid = (matches * 100.00) / (matches + mismatches);
    return pid;
}

/* Print the extended contig.
 */
static void PrintContig(const char* const name,
                        const char* const motif,
                        const char* const copies,
                        const char* const lflank, 
                        const char* const bases, 
                        const uint read_zstart,
                        const uint read_end,
                        const char* const rflank,
                        const uint motif_zstart, 
                        const uint motif_end,
                        const uint kmer_length) {
    uint indx;
    if ((lflank == NULL) || (rflank == NULL)) return;

    printf(">%s\t%s:%s:%d:%d\n", 
           name+1, motif, copies, 
           strlen(lflank) + motif_zstart - read_zstart, 
           strlen(lflank) + motif_end - read_zstart);

    if (lflank != NULL) {
        printf("%s", lflank);
    }

    for (indx = read_zstart; indx < read_end; indx++) {
        printf("%c", bases[indx]);
    }

    if (rflank != NULL) {
        printf("%s", rflank);
    }
    printf("\n");
}

static uint FindFirstGoodKmer(SparseHashMap& kmers, 
                              const char* const bases,
                              const uint num_kmers,
                              const Bool return_on_first,
                              const uint kmer_length) {
    uint result = 0;    

    if (return_on_first == TRUE) {
        Kmer word = BuildIndex(bases, kmer_length);
        Kmer antiword, stored;
        int i = 0;
    
        for (i = 0; i < num_kmers; i++) {
            word = GetNextKmer(word, bases, kmer_length, i);
            antiword = ReverseComplementKmer(word, kmer_length);
            stored = word < antiword ? word : antiword;

            if (CheckKmerInSparseHashMap(kmers, stored) == TRUE) break; 
        }
    
        result += (i + kmer_length);
    } else {
        Kmer word = BuildIndex(bases, kmer_length);
        Kmer antiword, stored;
        int i = 0, j = 0;
    
        for (i = 0; i < num_kmers; i++) {
            word = GetNextKmer(word, bases, kmer_length, i);
            antiword = ReverseComplementKmer(word, kmer_length);
            stored = word < antiword ? word : antiword;

            if (CheckKmerInSparseHashMap(kmers, stored) == TRUE) j = i;
        }

        result += j;
    }

    return result;
}

static void IssueWarningAboutKmerLength(const uint kmer_length)
{
    fprintf(stderr, "\n==========================================================================\n");
    fprintf(stderr, "Please note:\n");
    fprintf(stderr, "Kmer length %u for extension is greater than the flank requirement used for \n", kmer_length);
    fprintf(stderr, "STR discovery in select_STR_reads. This can lead to some blocks being ignored\n");
    fprintf(stderr, "during extension.\n");
    fprintf(stderr, "============================================================================\n\n");
}

/*
 * Assumptions and notes:
 *  a) We use a combination of a hash table and a bloom filter to select kmers 
       in the read dataset (reads1.fq, reads2.fq, ...) that are seen more than
       once. 
    b) We assume that the quality values in the files are in Sanger format. I 
       do not test to confirm this as of now. It would be prudent to add some 
       test for this.
    c) We do trim the low quality 3' end of the reads. Maybe I should provide 
       this as a user option.
    b) We iterate through the reads with putative STR sequences, and try
       to extend them on both ends till we reach a point where more than one 
       extension is possible. 
    c) In this version of the code, we quit extending if we hit a heterozygous 
       location (assuming that the species is diploid). An improved version 
       should be able to extend past the single-nucleotide polymorphisms by 
       identifying them as bubbles in the de-bruijn graph.
    d) We should add [start,stop] of the motif in the output as per request from
       Logan.
 */
static void ExtendShortTandemRepeatReads(const uint64_t haploid_genome_size,
                                         const uint kmer_length,
                                         const char* const fqname,
                                         char** const argv,
                                         const uint nameidx,
                                         const uint progress_chunk,
                                         const uint min_threshold,
                                         const uint max_threshold,
                                         const uint ploidy,
                                         const double heterozygosity,
                                         const uint expected_coverage,
                                         const double error_rate) {
    uint64_t genome_size = haploid_genome_size * (1 + heterozygosity * (ploidy - 1) * kmer_length);    
    uint64_t num_expected_kmers = genome_size * (1 + (expected_coverage * (1 - pow((1-error_rate),kmer_length))));
    PrintDebugMessage("Expecting %"PRIu64" kmers in this dataset with haploid genome size %"PRIu64" bps.\n", num_expected_kmers, haploid_genome_size);

    // all the non-singleton kmers shall be stored here.
    SparseHashMap kmers;
    kmers.rehash(genome_size);
    #ifdef Large
    kmers.set_deleted_key(0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF);
    #else
    kmers.set_deleted_key(0xFFFFFFFFFFFFFFFF);
    #endif

    // read and count the non-singleton kmers in the dataset.
    ReadAndCountNonSingletonKmers(kmers, 
                                  num_expected_kmers,
                                  kmer_length, 
                                  argv,
                                  nameidx, 
                                  progress_chunk,
                                  min_threshold,
                                  max_threshold);
    PrintDebugMessage("Read %zu kmers that are observed at least 2 times.", kmers.size());

    // traverse the reads with the STR's and try to extend them on both ends
    FastqSequence* sequence = ReadFastqSequence(fqname, FALSE, FALSE);

    uint indx1, indx2;
    char name[1024];
    char motif[7];
    char copies[1024];
    int zstart, end;
    char* lflank = NULL;
    char* rflank = NULL;
    Bool extensionWarningSet = FALSE;
    uint64_t numNotExtended = 0;

    while (sequence) {
        lflank = NULL;
        rflank = NULL;

        // parse the name of the read.
        if (sscanf(sequence->name, 
                   "%s\t%s\t%s\t%d\t%d\n", 
                   name, motif, copies, &zstart, &end) != 5) {
            PrintMessageThenDie("Error in parsing read name %s",sequence->name);
        }
        if (debug_flag) PrintDebugMessage("------------ %s ------------", name);
        

        // lets see if we can extend this read towards the 5' end.
        if (kmer_length > zstart) {
            // this block cannot be extended on this side
            indx1 = 0;
            numNotExtended++;
            if (extensionWarningSet == FALSE) {
                IssueWarningAboutKmerLength(kmer_length);
                extensionWarningSet = TRUE;
            }
        } else {
            indx1 = FindFirstGoodKmer(kmers, 
                                      sequence->bases, 
                                      zstart - kmer_length + 1, 
                                      TRUE,
                                      kmer_length);
            lflank = ExtendBackward(kmers, sequence->bases, indx1, kmer_length);
            
            // does this extension look correct?
            if ((lflank == NULL) || 
                (PercentIdentity(lflank, sequence->bases, indx1, TRUE) < 95.00)){
                goto next;
            }
        }
    
        // Lets see if we can extend this reads towards the 3' end.
        if ((sequence->slen - end) < kmer_length) {
            indx2 = 0;
            numNotExtended++;
        } else {
            indx2 = FindFirstGoodKmer(kmers, 
                                 sequence->bases + end, 
                                 sequence->slen - end - kmer_length + 1, 
                                 FALSE,
                                 kmer_length);
        }
        indx2 += end;
        rflank = ExtendForward(kmers, sequence->bases, indx2, kmer_length);
        
        // does this extension look correct?
        if ((rflank == NULL) || (PercentIdentity(rflank, sequence->bases, indx2, FALSE) < 95.00)) {
            goto next;
        }
               
        // print this contig after the extension.
        PrintContig(name, motif, copies, 
                    lflank, 
                    sequence->bases, indx1, indx2, 
                    rflank,
                    zstart, end, kmer_length);

    next:
        if (lflank) Ckfree(lflank);
        if (rflank) Ckfree(rflank);
        sequence = GetNextSequence(sequence);
    }
    
    CloseFastqSequence(sequence);
    PrintDebugMessage("Done with all extensions in %s", fqname);

    if (extensionWarningSet == TRUE) {
        fprintf(stderr, "\n==========================================================================\n");
        fprintf(stderr, "%"PRIu64" ends of merged reads were not extended as the kmer length for extension \n", numNotExtended);
        fprintf(stderr, "is greater than the flank requirement used for STR discovery in select_STR_reads\n");
        fprintf(stderr, "\n==========================================================================\n");
    }
}

int main(int argc, char** argv) {
    // start time management
    t0 = time(0);

    // set the version number
    program_version = VERSION;

    // parse the command line
    CommandLineArguments* cl_options = NewCommandLineArguments();

    // these are the valid options for the various commands
    AddOption(&cl_options, "min_threshold", "2", TRUE, TRUE, 
    "Discard kmers that are observed < min_threshold", NULL);
    AddOption(&cl_options, "max_threshold", "255", TRUE, TRUE,
    "Discard kmers that are observed > max_threshold", NULL);
    AddOption(&cl_options, "progress", "1000000", TRUE, TRUE,
    "print progress every so many sequences", NULL);
    AddOption(&cl_options, "flanks", "1024", TRUE, TRUE,
    "the maximum size of flanks extension", NULL);
    AddOption(&cl_options, "ploidy", "2", TRUE, TRUE,
    "the ploidy of the genome", NULL);
    AddOption(&cl_options, "heterozygosity", "0.001", TRUE, TRUE,
    "fraction of nucleotides that differ between inherited chromosomes", NULL);
    AddOption(&cl_options, "errorrate", "0.01", TRUE, TRUE,
    "expected error rate in sequencing", NULL);

    ParseOptions(&cl_options, &argc, &argv);

    // does the user just want some help
    Bool print_help = GetOptionBoolValueOrDie(cl_options, "help");
    if (print_help == TRUE) {
        PrintSimpleUsageString(cl_options);
        return EXIT_SUCCESS;
    }

    // does the user know what he/she is doing?
    if (argc < 7){
        PrintSimpleUsageString(cl_options);
        return EXIT_FAILURE;
    }

    uint64_t genome_size;
    if (sscanf(argv[1], "%"PRIu64, &genome_size) != 1) {
        PrintMessageThenDie("Expected genome size should be an integer: %s",
        argv[1]);
    }

    uint expected_coverage;
    if (sscanf(argv[2], "%u", &expected_coverage) != 1) {
        PrintMessageThenDie("Expected coverage should be an integer: %s",
        argv[2]);
    } 


    uint kmer_length;
    if (sscanf(argv[3], "%u", &kmer_length) != 1) {
        #ifdef Large
        PrintMessageThenDie("Kmer length should be an odd integer < 64: %s",
        argv[3]);
        #else
        PrintMessageThenDie("Kmer length should be an odd integer < 32: %s",
        argv[3]);
        #endif
    }
    if (kmer_length % 2 == 0) {
        PrintWarning("Kmer length should be an odd integer, using %u",
        --kmer_length);
    }
    kmer_buffer = (char*)CkalloczOrDie(kmer_length + 1);

    char* str_reads_name = argv[4];

    // kmers seen less than these many times should be ignored.
    uint min_threshold = GetOptionUintValueOrDie(cl_options, "min_threshold");
    uint max_threshold = GetOptionUintValueOrDie(cl_options, "max_threshold");

    // the maximum extension on both sides
    flank_chunk = GetOptionUintValueOrDie(cl_options, "flanks");

    // do I need additional debug info
    debug_flag = GetOptionBoolValueOrDie(cl_options, "debug");

    // how often should I print progress?
    uint progress_chunk = GetOptionUintValueOrDie(cl_options, "progress");

    // the expected ploidy
    uint ploidy = GetOptionUintValueOrDie(cl_options, "ploidy");

    // the heterozygosity in the genome
    double heterozygosity = GetOptionDoubleValueOrDie(cl_options,"heterozygosity");

    // error rate 
    double error_rate = GetOptionDoubleValueOrDie(cl_options,"errorrate");

    ExtendShortTandemRepeatReads(genome_size, 
                                 kmer_length, 
                                 str_reads_name,
                                 argv, 
                                 argc, 
                                 progress_chunk,
                                 min_threshold,
                                 max_threshold,
                                 ploidy,
                                 heterozygosity,
                                 expected_coverage,
                                 error_rate);

    Ckfree(kmer_buffer);
    FreeParseOptions(&cl_options, &argv);      
    return EXIT_SUCCESS;
}

