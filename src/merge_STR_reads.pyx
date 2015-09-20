from sys import argv, stderr, stdin, exit, stdout
from getopt import getopt, GetoptError
from profile import run
from time import time

from fastq import *

cdef float PercentIdentity(char* seq1, char* seq2):
    cdef:
        int matches = 0
        int mismatches = 0

    for a,b in zip(seq1,seq2):
        if a == b:
            matches += 1
        else:
            mismatches += 1

    cdef float fraction = matches * 1.0 / (matches + mismatches)
    return fraction * 100.0

def Consensus(flank1, flank2, flip = False):
    seq1,qual1 = flank1
    seq2,qual2 = flank2

    seq = ""
    qual = ""

    if flip:
        seq1 = seq1[::-1]
        qual1 =  qual1[::-1]
        seq2 = seq2[::-1]
        qual2 = qual2[::-1]

        for i,(s1,q1,s2,q2) in enumerate(zip(seq1,qual1,seq2,qual2)):
            if q1 > q2:
                seq += s1
                qual += q1
            else:
                seq += s2
                qual += q2

        if len(seq1) > (i + 1):
            assert len(seq2) == (i + 1)
            for s,q in zip(seq1[i+1:],qual1[i+1:]):
                seq += s
                qual += q
        elif len(seq2) > (i + 1):
            assert len(seq1) == (i + 1)
            for s,q in zip(seq2[i+1:],qual2[i+1:]):
                seq += s
                qual += q

        seq = seq[::-1]
        qual = qual[::-1]
    else:

        for i,(s1,q1,s2,q2) in enumerate(zip(seq1,qual1,seq2,qual2)):
            if q1 > q2:
                seq += s1
                qual += q1
            else:
                seq += s2
                qual += q2  

        # Only one of the two sequences should be longer than (i+1)
        if len(seq1) > (i + 1): 
            assert len(seq2) == (i + 1)
            for s,q in zip(seq1[i+1:],qual1[i+1:]):
                seq += s
                qual += q
        elif len(seq2) > (i + 1):
            assert len(seq1) == (i + 1)
            for s,q in zip(seq2[i+1:],qual2[i+1:]):
                seq += s
                qual += q

    return seq,qual

class ScoreParam:
    def __init__(self, int gap, int match, int mismatch):
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

def make_matrix(int sizex, int sizey):
    return [[0]*sizey for i in xrange(sizex)]

cdef float Align(char* x, char* y, float pid_threshold, 
                 int* pgaps, score = ScoreParam(-3,1,-1)):
    # lets do something quick to eliminate the ones which definitely will not
    # align
    cdef float pid = PercentIdentity(x,y)
    if pid >= pid_threshold:
        pgaps[0] = 0
        return pid

    cdef:
        int xmax = len(x) + 1
        int ymax = len(y) + 1

    A = make_matrix(xmax, ymax)
    B = make_matrix(xmax, ymax)

    cdef:
        int best = 0
        int optlox = 0
        int optloy = 0

        int score1
        int score2 
        int score3

        int i
        int j

    # fill in A in the right order
    for i in xrange(1, xmax): 
        for j in xrange (1, ymax):
            score1 = A[i][j-1] + score.gap
            score2 = A[i-1][j] + score.gap
            if x[i-1] == y[j-1]:
                score3 = A[i-1][j-1] + score.match 
            else:
                score3 = A[i-1][j-1] + score.mismatch

            A[i][j] = max(score1,score2,score3,0)
    
            if A[i][j] == score3:
                B[i][j] = 0
            elif A[i][j] == score2:
                B[i][j] = 1
            elif A[i][j] == score1:
                B[i][j] = 2

            # track the cell with the largest score
            if A[i][j] >= best: 
                best = A[i][j] 
                optlox = i
                optloy = j

    # lets trace backwards from the best location to find the number of
    # differences between the two in the region of the alignments
    cdef int max_score = best
    i = optlox
    j = optloy

    cdef: 
        int num_mismatches = 0
        int num_matches = 0
        int num_gaps = 0

    while max_score > 0:
        if B[i][j] == 0:
            #print x[i-1],y[j-1]
            i -= 1
            j -= 1
            if x[i-1] != y[j-1]:
                num_mismatches += 1
            else:
                num_matches += 1
        elif B[i][j] == 1:
            #print x[i-1],"-"
            i -= 1
            num_gaps += 1
        elif B[i][j] == 2:
            #print "-",y[j-1]
            j -= 1
            num_gaps += 1
        max_score = A[i][j]
        
    if i != 0 and j != 0:
        pgaps[0] = 0
        return 0   

    cdef float identity = num_matches * 100.0 / (num_matches + num_mismatches)

    # return the opt score and the best location
    pgaps[0] = num_gaps
    return identity

def AlignFlanks(block, s, motif, copies, zstart, end, max_threshold, 
                pid_threshold, debug_flag):
    assert s.seq[zstart:end] == motif * copies

    cdef:
        int zstart1
        int end1
        int support1
        char* seq1
        char* qual1

    # The first sequence is from a block of reads
    zstart1,end1,support1,seq1,qual1,names = block

    # This is a fastq sequence I am looking at the first time
    cdef int zstart2 = zstart
    cdef int end2 = end

    cdef char* seq2 = s.seq
    cdef char* qual2 = s.qual

    if debug_flag:
        print >> stderr, "Aligning against %s" % ",".join(names)
        print >> stderr, "Block sequence: %s" % seq1
    # Lets look at the matches on the left of the motif
    lflank1 = (seq1[:zstart1], qual1[:zstart1])
    lflank2 = (seq2[:zstart2], qual2[:zstart2])

    cdef int gaps

    pid = Align(lflank1[0][::-1], lflank2[0][::-1], pid_threshold, &gaps)
    if pid < pid_threshold or gaps > 2:
        if debug_flag:
            print >> stderr, "Low percent identity (%2.2f) for the left flank."\
            % pid
        return None
   
    rflank1 = (seq1[end1:], qual1[end1:])
    rflank2 = (seq2[end2:], qual2[end2:])

    pid = Align(rflank1[0], rflank2[0], pid_threshold, &gaps)
    if pid < pid_threshold or gaps > 2: 
        if debug_flag:
            print >> stderr, "Low percent identity (%2.2f) for the right flank."\
            % pid
        return None

    mflank1 = (seq1[zstart1:end1], qual1[zstart1:end1])
    mflank2 = (seq2[zstart2:end2], qual2[zstart2:end2])

    if len(block[5]) > max_threshold:
        return block

    # If I am here then these this block and this read seem like they support
    # the same STR. Lets combine them to create a new block
    lseq,lqual = Consensus(lflank1, lflank2, True)
    mseq,mqual = Consensus(mflank1, mflank2)
    rseq,rqual = Consensus(rflank1, rflank2)
    
    names.append(s.name)
    return (len(lseq),
            len(lseq) + len(mseq),
            support1 + 1,
            lseq + mseq + rseq,
            lqual + mqual + rqual,
            names)


#def FastAlignFlanks(block, s, motif, copies, zstart, end, pid_threshold):
#    assert s.seq[zstart:end] == motif * copies
#
#    # The first sequence is from a block of reads
#    zstart1,end1,support1,seq1,qual1,names = block
#
#    # This is a fastq sequence I am looking at the first time
#    zstart2 = zstart
#    end2 = end
#
#    seq2 = s.seq
#    qual2 = s.qual
#
#    if debug_flag:
#        print >> stderr, "Aligning against %s" % ",".join(names)
#        print >> stderr, "Block sequence: %s" % seq1
#    # Lets look at the matches on the left of the motif
#    lflank1 = (seq1[:zstart1], qual1[:zstart1])
#    lflank2 = (seq2[:zstart2], qual2[:zstart2])
#
#    pid = PercentIdentity(lflank1[0][::-1], lflank2[0][::-1])
#    if pid < pid_threshold:
#        if debug_flag:
#            print >> stderr, "Low percent identity (%2.2f) for the left flank."\
#            % pid
#        return None
#   
#    rflank1 = (seq1[end1:], qual1[end1:])
#    rflank2 = (seq2[end2:], qual2[end2:])
#
#    pid = PercentIdentity(rflank1[0], rflank2[0])
#    if pid < pid_threshold: 
#        if debug_flag:
#            print >> stderr, "Low percent identity (%2.2f) for the right flank."\
#            % pid
#        return None
#
#    mflank1 = (seq1[zstart1:end1], qual1[zstart1:end1])
#    mflank2 = (seq2[zstart2:end2], qual2[zstart2:end2])
#
#    # If I am here then these this block and this read seem like they support
#    # the same STR. Lets combine them to create a new block
#    lseq,lqual = Consensus(lflank1, lflank2, True)
#    mseq,mqual = Consensus(mflank1, mflank2)
#    rseq,rqual = Consensus(rflank1, rflank2)
#    
#    names.append(s.name)
#    return (len(lseq),
#            len(lseq) + len(mseq),
#            support1 + 1,
#            lseq + mseq + rseq,
#            lqual + mqual + rqual,
#            names)

def merge_str_reads(int klength, char* filename, int min_threshold, 
                    int max_threshold,
                    float pid_threshold, int debug_flag):
    records = fastq(filename)

    # Store all the merged reads here. 
    blocks = {}
    start_time = time()

    cdef int indx = 0
    for r in records:
        s = r.fastqsequence
        indx += 1
        if indx % 1000000 == 0:
            print >> stderr, "Processing read %d [%d sec. elapsed]" % (indx, time() - start_time)        

        if debug_flag:
            print >> stderr, "-"*79
            print >> stderr, s.name
            print >> stderr, s.seq
            print >> stderr, ""

        name,m1,c1,s1,e1,m2,c2,s2,e2 = s.name.strip().split("\t")
        c1 = int(c1)
        s1 = int(s1)
        e1 = int(e1)
        c2 = int(c2)
        s2 = int(s2)
        e2 = int(e2)

        # Check if the sequence aligns to some reads that we have already
        # processed.
        does_align = False
        aligned_block = None

        # We will keep only one of the two motifs (the motif on the forward
        # strand or the motif on the reverse strand)
        motif  = m1
        copies = c1
        zstart = s1
        end    = e1

        lflank = s.seq[zstart - klength : zstart]
        rflank = s.seq[end : end + klength]

        if (motif,lflank,rflank) in blocks:
            for block in blocks[(motif,lflank,rflank)]:
                merged_block = AlignFlanks(block,s,motif,copies,zstart,end,
                                          max_threshold,pid_threshold,debug_flag)
                if debug_flag: print >> stderr, ""
                if merged_block != None:
                    aligned_block = block
                    does_align = True
                    break

        if does_align == False:
            motif  = m2
            copies = c2
            zstart = s2
            end    = e2
            s.reverse_complement()
        
            lflank = s.seq[zstart - klength : zstart]
            rflank = s.seq[end : end + klength]

            if (motif,lflank,rflank) in blocks:
                for block in blocks[(motif,lflank,rflank)]:
                    merged_block=AlignFlanks(block,s,motif,copies,zstart,end,
                                          max_threshold,pid_threshold,debug_flag)
                    if debug_flag: print >> stderr, ""
                    if merged_block != None:
                        aligned_block = block
                        does_align = True
                        break

        if does_align == True:
            blocks[(motif,lflank,rflank)].remove(aligned_block)
            blocks[(motif,lflank,rflank)].append(merged_block)
            continue

        if (motif,lflank,rflank) not in blocks: 
            blocks[(motif,lflank,rflank)] = []
        blocks[(motif,lflank,rflank)].append((zstart,end,1,s.seq,s.qual,[s.name]))

    records.close()
    print >> stderr, "Processed %d reads [%d sec. elapsed]" % (indx, time() - start_time)        

    # Lets print out the merged reads
    indx = 1
    for (motif,lflank,rflank),reads in blocks.items():
        for zstart,end,num_members,sequence,qual,names in reads:
            copies = [str(x.split("\t")[2]) for x in names]
            copies = list(set(copies))

            if len(names) >= min_threshold and \
               len(names) <= max_threshold and \
               len(copies) <= 2:
                print "@Block%d\t%s\t%s\t%d\t%d" % \
                    (indx,motif,",".join(copies),zstart,end)
                print sequence
                print "+"
                print qual
    
                if debug_flag:
                    print >> stderr, "Block%d\t%s\t%s\t%d\t%d" % \
                    (indx,motif,",".join(copies),zstart,end)
                    print >> stderr, sequence
                    for name in names:
                        print >> stderr, "\t%s" % name

                indx += 1  