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

#include "sparse_word_hash.h"

char* program_version_major = "0";
char* program_version_minor = "30";
char* program_revision_date = "20160218";
char* program_name          = "merge_STR_reads";
char* program_description   = 
    "Merge reads that support the same STR";
char* program_use           = 
    "merge_STR_reads [options] klength reads.str.fq";

Bool debug_flag;

// A percent identity threshold. 
const double pid_threshold = 90.0;

static void Reverse(char* const seq, const uint len)
{
    char* s = seq;
    char* p = seq + len - 1;
    
    while (s <= p) {
        uchar c;
        c = *s;
        *s = *p;
        *p = c;
        ++s; --p;
    }
}

// align the two sequences, and calculate the consensus sequence
static void Align(const char* const seq1, const char* const qual1, 
                  const uint zstart1, const uint end1,
                  const char* const seq2, const char* const qual2,
                  const uint zstart2, const uint end2,
                  uint* gaps, float* pid, char** pseq, char** pqual, 
                  const uint slen, const Bool is_right_gapped)
{
    char* seq = *pseq;
    char* qual = *pqual;

    uint len1 = end1 - zstart1 + 1;
    uint len2 = end2 - zstart2 + 1;

    const char* t1 = seq1 + zstart1;
    const char* q1 = qual1 + zstart1;
    const char* t2 = seq2 + zstart2;
    const char* q2 = qual2 + zstart2;

    // counters in loops
    uint i, j;

    // initialize the scores
    int match    = 1;
    int mismatch = -1;
    int gap  = -3;

    int** A = (int**)CkalloczOrDie(len1 * sizeof(int*));
    for(i = 0; i < len1; i++) {
        A[i] = (int*)CkalloczOrDie(len2 * sizeof(int));
    }
    int** B = (int**)CkalloczOrDie(len1 * sizeof(int*));
    for(i = 0; i < len1; i++) {
        B[i] = (int*)CkalloczOrDie(len2 * sizeof(int));
    }

    int best = 0, optlox = 0, optloy = 0;
    int score1, score2, score3;

    for (i = 1; i < len1; i++) {
        for (j = 1; j < len2; j++) {
            score1 = A[i][j-1] + gap;

            score2 = A[i-1][j] + gap;

            if (t1[i-1] == t2[j-1]) {
                score3 = A[i-1][j-1] + match;
            } else {
                score3 = A[i-1][j-1] + mismatch;
            }

            A[i][j] = MAX(MAX(score1,score2),MAX(score3,0));

            if (A[i][j] == score3) {
                B[i][j] = 0;
            } else if (A[i][j] == score2) {
                B[i][j] = 1;
            } else if (A[i][j] == score1) {
                B[i][j] = 2;
            }
    
            if (A[i][j] >= best) {
                best = A[i][j];
                optlox = i;
                optloy = j;
            }
        }
    }

    // trace backwards to find the best location 
    int matches = 0, mismatches = 0;
    int num_gaps = 0;

    i = optlox;
    j = optloy;
    int max_score = best;
    int sindx = 0;

    if (is_right_gapped) {
        uint a = 0;
        if (len1 > len2) {
            int a = len1-1;
            while (a > j) {
                seq[sindx] = t1[a-1];
                qual[sindx++] = q1[a-1];
                a--;
            }
        } else {
            int a = len2-1;
            while (a > j) {
                seq[sindx] = t2[a-1];
                qual[sindx++] = q2[a-1];
                a--;
            }
        }
    } else {
        num_gaps = len2 - optloy - 1;
    }

    while ((max_score > 0) && (i >= 1) && (j >= 1)) {
        if (B[i][j] == 0) {
            if (q1[i-1] > q2[j-1]) {
                seq[sindx] = t1[i-1];
                qual[sindx] = q1[i-1];
            } else {
                seq[sindx] = t2[j-1];
                qual[sindx] = q2[j-1];
            }
            sindx++;
            if (t1[i-1] != t2[j-1]) {
                mismatches++;
            } else {
                matches++;
            }
            i--; j--;
        } else if (B[i][j] == 1) {
            if (q1[i-1] > '5') {
                seq[sindx] = t1[i-1];
                qual[sindx++] = q1[i-1];
            }
            i--;
            num_gaps++;
        } else if (B[i][j] == 2) {
            if (q2[j-1] > '5') {
                seq[sindx] = t2[j-1];
                qual[sindx++] = q2[j-1];
            }
            j--;
            num_gaps++;
        }
        
        max_score = A[i][j];
    }

    if (is_right_gapped) {
        num_gaps = MAX(i,j) - MIN(i,j);
    } else {
        uint a = 0;
        if (len1 > len2) {
            int a = i;
            while (a > 0) {
                seq[sindx] = t1[a-1];
                qual[sindx++] = q1[a-1];
                a--;
            }
        } else {
            int a = j;
            while (a > 0) {
                seq[sindx] = t2[a-1];
                qual[sindx++] = q2[a-1];
                a--;
            }
        }
    }

    Reverse(seq, sindx);
    Reverse(qual, sindx);
    seq[sindx] = 0;
    qual[sindx] = 0;

    *gaps = num_gaps;
    *pid = matches * 100.0 / (matches + mismatches);

    for(i = 0; i < len1; i++) {
        Ckfree(A[i]);
    }
    Ckfree(A);
    for(i = 0; i < len1; i++) {
        Ckfree(B[i]);
    }
    Ckfree(B);
}

static Bool AlignFlanks(Block* const block,
                        const FastqSequence* const seq,
                        const char* const motif,
                        const int copies,
                        const int zstart,
                        const int end,
                        const uint max_threshold)
{
    float pid;    
    uint gaps = 0;
    
    int slen = MAX(block->slen,seq->slen);
    
    // Align the sequences to the left of the STR
    char* lseq   = (char*)CkalloczOrDie(slen);
    char* lqual  = (char*)CkalloczOrDie(slen);

    Align(block->seq, block->qual, 0, block->zstart, 
    seq->bases, seq->quals, 0, zstart, &gaps, &pid, &lseq, &lqual, slen, FALSE);

    if ((pid < pid_threshold) || (gaps > 2)) {
        if (debug_flag) {
            PrintDebugMessage(
            "Low pid (%2.2f) or too many gaps (%d) for the left flank.", 
            pid, gaps);
        }
        Ckfree(lseq);
        Ckfree(lqual);
        return NULL;
    }

    // Align the sequences to the right of the STR
    char* rseq   = (char*)CkalloczOrDie(slen);
    char* rqual  = (char*)CkalloczOrDie(slen);

    Align(block->seq, block->qual, block->end, block->slen, 
    seq->bases, seq->quals, end, seq->slen, &gaps, &pid,&rseq,&rqual,slen, TRUE);

    if ((pid < pid_threshold) || (gaps > 2)) {
        if (debug_flag) {
            PrintDebugMessage(
            "Low pid (%2.2f) or too many gaps (%d) for the right flank.",
            pid, gaps);
        }
        Ckfree(lseq);
        Ckfree(lqual);
        Ckfree(rseq);
        Ckfree(rqual);
        return FALSE;
    }

    // If I am here then these this block and this read seem like they support
    // the same STR. Lets combine them to create a new block
    block->zstart = strlen(lseq);
    block->support++;
    Bool found = FALSE;
    uint maxcopies = 0, i = 0, j = 0;
    for (i = 0, j = 0; i < 3; i++) {
        if (block->copies[i] != 0) {
            if (block->copies[i] > maxcopies) { 
                maxcopies = block->copies[i];
            }
            j++;
            if (block->copies[i] == copies) {
                found = TRUE;
            }
        }     
    }
    if ((found == FALSE) && (j <= 2)) {
        block->copies[j] = copies;
    }
    block->end = block->zstart + strlen(motif) * maxcopies;
    block->slen = block->end + strlen(rseq);
    Ckfree(block->seq);
    block->seq = (char*)CkalloczOrDie(block->slen+1);
    memcpy(block->seq, lseq, block->zstart);
    for (int i = 0; i < maxcopies; i++) {
        sprintf(block->seq + block->zstart + (i * strlen(motif)), "%s", motif);
    }
    memcpy(block->seq + block->end, rseq, block->slen - block->end);
    Ckfree(block->qual);
    block->qual = (char*)CkalloczOrDie(block->slen+1);
    memcpy(block->qual, lqual, block->zstart);
    for (int i = 0; i <= (strlen(motif) * maxcopies); i++) {
        sprintf(block->qual + block->zstart + i, "!");
    }
    memcpy(block->qual + block->end, rqual, block->slen - block->end);

    Ckfree(lseq);
    Ckfree(lqual);
    Ckfree(rseq);
    Ckfree(rqual);

    return TRUE;
}

static void MergeShortTandemRepeatReads(const uint klength, 
                                        const char* const fqname,
                                        char** const argv, 
                                        const uint nameidx, 
                                        const uint progress_chunk,
                                        const uint min_threshold,
                                        const uint max_threshold)
{
    uint64_t num_sequence_processed = 0;

    // we store the blocks here
    SparseWordHashMap blocks;

    // buffer to write the keys into
    char* buffer = (char*)CkalloczOrDie(1024);

    FastqSequence* sequence = ReadFastqSequence(fqname, FALSE, FALSE);
    
    char name[1024];
    char fmotif[7], rmotif[7];
    int fcopies, fzstart, fend, rcopies, rzstart, rend;
    char lflank[1024], rflank[1024];

    Block* block = NULL;
    Bool merged_block = FALSE;

    while (sequence) {
        if (debug_flag == TRUE) {
            PrintDebugMessage("Processing %s", sequence->name + 1);
        } else {
            // print progress
            num_sequence_processed += 1;
            if ((num_sequence_processed - 1) % progress_chunk == 0) {
                PrintDebugMessage("Processing read number %"PRIu64": %s",
                num_sequence_processed, sequence->name + 1);
            }
        }
        
        // parse the name of the read.
        if (sscanf(sequence->name,
                   "%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n",
                   name, fmotif, &fcopies, &fzstart, &fend, 
                   rmotif, &rcopies, &rzstart, &rend) != 9) {
            PrintMessageThenDie("Error in parsing read name %s",sequence->name);
        }
        
        // Check if the sequence aligns to some reads that we have already
        // processed
        merged_block = FALSE;
          
        memcpy(lflank, sequence->bases + (fzstart - klength), klength);
        lflank[klength] = '\0';
        memcpy(rflank, sequence->bases + fend, klength);
        rflank[klength] = '\0';
        sprintf(buffer, "%s %s %s", fmotif, lflank, rflank);

        if (CheckInSparseWordHashMap(blocks, buffer) == TRUE) {
            block = blocks[buffer];
            merged_block = AlignFlanks(block, sequence, fmotif, fcopies, 
                                       fzstart, fend, max_threshold);
        }            

        if (merged_block == FALSE) {
            ReverseComplementSequence(sequence);
            memcpy(lflank, sequence->bases + (rzstart - klength), klength);
            lflank[klength] = '\0';
            memcpy(rflank, sequence->bases + rend, klength);
            rflank[klength] = '\0';
            sprintf(buffer, "%s %s %s", rmotif, lflank, rflank);

            if (CheckInSparseWordHashMap(blocks, buffer) == TRUE) {
                block = blocks[buffer];
                merged_block = AlignFlanks(block, sequence, rmotif, rcopies, 
                                           rzstart, rend, max_threshold);
            }
        }

        if (merged_block == FALSE) {
            block = (Block*)CkalloczOrDie(sizeof(Block));
            block->zstart = rzstart;
            block->end = rend;
            block->slen = sequence->slen;
            block->support = 1;
            block->seq = CopyString(sequence->bases);
            block->qual = CopyString(sequence->quals);
            block->copies[0] = rcopies;
            blocks[CopyString(buffer)] = block;
        }

        if (debug_flag) 
            fprintf(stderr, "-----------------------------------------------\n");
        sequence = GetNextSequence(sequence);
    }

    Ckfree(buffer);

    // lets print the merged blocks
    SparseWordHashMap::iterator it;
    uint bindex = 1;
    for (it = blocks.begin(); it != blocks.end(); it++) {
        if (sscanf((*it).first, "%s %*s %*s", fmotif) != 1) {
            PrintMessageThenDie("Error in parsing key : %s", (*it).first);
        }
        Block* block = (*it).second;

        if ((block->support >= min_threshold) && 
            (block->support <= max_threshold) && 
            (block->copies[2] == 0) && 
            (block->copies[1] != 0)) {
            printf("@Block%d\t%s\t", bindex++, fmotif);
            for (int i = 0; i < 3; i++) {
                if (block->copies[i] != 0) {
                    if (i != 0) printf (",");
                    printf("%d",block->copies[i]);
                }
            }
            printf("\t%d\t%d\n", block->zstart, block->end);
            printf("%s\n", block->seq);
            printf("+\n");
            printf("%s\n", block->qual);

            Ckfree(block->seq);
            Ckfree(block->qual);
            Ckfree(block);
        }
    }
}


int main(int argc, char** argv) {
    // start time management
    t0 = time(0);

    // parse the command line
    CommandLineArguments* cl_options = NewCommandLineArguments();

    // these are the valid options for the various commands
    AddOption(&cl_options, "min_threshold", "3", TRUE, TRUE, 
    "Discard blocks that include < min_threshold reads", NULL);
    AddOption(&cl_options, "max_threshold", "10000", TRUE, TRUE,
    "Discard blocks that include > max_threshold reads", NULL);
    AddOption(&cl_options, "progress", "1000000", TRUE, TRUE,
    "print progress every so many sequences", NULL);

    ParseOptions(&cl_options, &argc, &argv);

    // does the user just want some help
    Bool print_help = GetOptionBoolValueOrDie(cl_options, "help");
    if (print_help == TRUE) {
        PrintSimpleUsageString(cl_options);
        return EXIT_SUCCESS;
    }

    // does the user know what he/she is doing?
    if (argc < 2){
        PrintSimpleUsageString(cl_options);
        return EXIT_FAILURE;
    }

    uint kmer_length;
    if (sscanf(argv[1], "%u", &kmer_length) != 1) {
        #ifdef Large
        PrintMessageThenDie("Kmer length should be an odd integer < 64: %s",
        argv[1]);
        #else
        PrintMessageThenDie("Kmer length should be an odd integer < 32: %s",
        argv[1]);
        #endif
    }
    if (kmer_length % 2 == 0) {
        PrintWarning("Kmer length should be an odd integer, using %u",
        --kmer_length);
    }

    char* str_reads_name = argv[2];

    // kmers seen less than these many times should be ignored.
    uint min_threshold = GetOptionUintValueOrDie(cl_options, "min_threshold");
    uint max_threshold = GetOptionUintValueOrDie(cl_options, "max_threshold");

    // do I need additional debug info
    debug_flag = GetOptionBoolValueOrDie(cl_options, "debug");

    // how often should I print progress?
    uint progress_chunk = GetOptionUintValueOrDie(cl_options, "progress");

    MergeShortTandemRepeatReads(kmer_length, 
                                str_reads_name,
                                argv, 
                                argc, 
                                progress_chunk,
                                min_threshold,
                                max_threshold);

    FreeParseOptions(&cl_options, &argv);      
    return EXIT_SUCCESS;
}

