from sys import argv, stderr, stdin, exit, stdout
from getopt import getopt, GetoptError
from itertools import izip_longest
from string import maketrans

import re

from fastq import *

class MotifDetails:
    """Simple class to store additional information about an STR.
    """
    pass

def ReadSTRList(filename):
    """Read and save the STR information from a file.

    The primary reason to do this is when we are debugging or running simulations
    to understand why certain motif:copies are not detected by our method. We
    will save the motif, the number of copies of the motif in the STR, 20 bp 
    flank sequence 5' to the STR, 20 bps sequence 3' to the STR. If less than
    20 bps are present on any side of the STR then we just store those for now.
    """
    strlist = {}

    # If we are not looking to check any STRs, then just return an empty list.
    if filename == None: return strlist

    file = open(filename, "r")

    for line in file:
        motif, copies, lflank, rflank = line.strip().split()       
        copies = int(copies)

        if (motif,copies) not in strlist:
            strlist[(motif,copies)] = []

        motiflist = strlist[(motif,copies)]

        md = MotifDetails()
        md.lflank = lflank
        md.rflank = rflank
        md.found  = 0

        motiflist.append(md)
    file.close()
    return strlist

def IsHomopolymer(motif):
    """Return True if the input string is a homopolymer.
    """
    if len(list(set(list(motif)))) == 1:
        return True
    return False

#def MatchPattersAndReportBest(string, patterns, remove_homopolymers,
#                              flanking_distance, debug_flag):
#    """Find and return the best motif:copies in this string.
#
#    In cases where the STR could be from 2mer and 4mer for example, I will always
#    report the 2mer. This is expected behavior as per Logan.
#    """
#    matches = []
#    for pattern in patterns:
#        for match in re.finditer(pattern, string, flags = re.IGNORECASE):
#            STR = match.group(1)
#            motif = match.group(2)
#            copies = len(STR)/len(motif)
#            start = match.start()
#            end = match.end()
#
#            # Special case: motif's with Ns are not desirable. Lets throw
#            # those away
#            if motif.find("N") != -1: 
#                if debug_flag: print >> stderr, "N in motif."
#                continue
#
#            # Throw away matches that do not satisfy the flank requr.
#            if start < flanking_distance: 
#                if debug_flag: print >> stderr, "Insufficient 5' flank."
#                continue
#            if (len(string) - end) < flanking_distance: 
#                if debug_flag: print >> stderr, "Insufficient 3' flank."
#                continue
#                
#            # ignore homopolymer runs if that is requested by the user
#            if remove_homopolymers == True and IsHomopolymer(motif):
#                continue
#
#            # This match could be interesting. Lets save it.
#            matches.append((len(STR), copies, STR, motif, start, end))
#
#    if len(matches) == 0: return None
#
#    # Sort so that the longest stretch of STR's is the first match. If 
#    # two matches are equal in length, then I want the one with more
#    # number of copies to be the first match. The second member of the
#    # stored set assures of that.
#    matches.sort(reverse = True)
#    match = matches[0]
#
#    return match

def MatchPattersAndReportBest(string, patterns, remove_homopolymers,
                              flanking_distance, debug_flag):
    """Find and return the best motif:copies in this string.

    In cases where the STR could be from 2mer and 4mer for example, I will always
    report the 2mer. This is expected behavior as per Logan.
    """
    matches = []
    for pattern in patterns:
        for match in re.finditer(pattern, string, flags = re.IGNORECASE):
            STR = match.group(1)
            motif = match.group(2)
            copies = len(STR)/len(motif)
            start = match.start()
            end = start + len(STR)

            # Special case: motif's with Ns are not desirable. Lets throw
            # those away
            if motif.find("N") != -1: 
                if debug_flag: print >> stderr, "N in motif."
                continue

            # This match could be interesting. Lets save it.
            matches.append((len(STR), copies, STR, motif, start, end))

    if len(matches) == 0: return None

    # Sort so that the longest stretch of STR's is the first match. If 
    # two matches are equal in length, then I want the one with more
    # number of copies to be the first match. The second member of the
    # stored set assures of that.
    matches.sort(reverse = True)
    match = matches[0]
    motif = match[3]
    start = match[4]
    end   = match[5]
    # Throw away matches that do not satisfy the flank requr.
    if start < flanking_distance: 
        if debug_flag: print >> stderr, "Insufficient 5' flank.: %s" % start
        return None
    if (len(string) - end) < flanking_distance: 
        if debug_flag: 
            print >> stderr, "Insufficient 3' flank.: %s" % (len(string)-end)
        return None
                
    # ignore homopolymer runs if that is requested by the user
    if remove_homopolymers == True and IsHomopolymer(motif):
        return None

    return match

    #eligible_match = None
    #
    #for match in matches:
    #    motif = match[3]
    #    start = match[4]
    #    end   = match[5]
    #    # Throw away matches that do not satisfy the flank requr.
    #    if start < flanking_distance: 
    #        if debug_flag: print >> stderr, "Insufficient 5' flank.: %s" % start
    #        continue
    #    if (len(string) - end) < flanking_distance: 
    #        if debug_flag: 
    #            print >> stderr, "Insufficient 3' flank.: %s" % (len(string)-end)
    #        continue
    #                
    #    # ignore homopolymer runs if that is requested by the user
    #    if remove_homopolymers == True and IsHomopolymer(motif):
    #        continue

    #    eligible_match = match
    #    break

    #return eligible_match


def PrintSequence(s, match, rc_match, illumina_quals):
    """Print this fastq sequence with some additional information about the STR.
    """
    print "@%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d" % \
        (s.name, 
         match[3], match[1], match[4], match[5], 
         rc_match[3], rc_match[1], rc_match[4], rc_match[5])
    print "%s" % s.seq
    print "+"
    if illumina_quals:
        qual = [chr(ord(x)-64+33) for x in s.qual]
        qual = "".join(qual)
    else:
        qual = s.qual
    print qual

def ReverseComplement(seq):
    complement = maketrans('atcgnATCGN', 'tagcnTAGCN')
    return seq.translate(complement)[::-1]

def PercentIdentity(md, lflank, rflank):
    s = 0 # same in flanks
    d = 0 # different in flanks

    for r,q in izip_longest(md.lflank[::-1],lflank[::-1]):
        if r != None and q != None: 
            if r.upper() == q.upper():
                s += 1
            else:
                d += 1
    for r,q in izip_longest(md.rflank,rflank):
        if r != None and q != None: 
            if r.upper() == q.upper():
                s += 1
            else:
                d += 1
    return ((s * 100.0) / (s + d)) 

def RegisterSTR(string, match, rc_match, to_check):
    max_pid = 0
    motif = match[3]
    rc_motif = ReverseComplement(motif) 

    e1 = match[4]
    s1 = e1 - 20
    if s1 < 0: s1 == 0

    s2 = match[5]
    e2 = s2 + 20
    if e2 > len(string): e2 = len(string)

    lflank = string[s1:e1]
    rflank = string[s2:e2]   

    if (motif,match[1]) in to_check:
        mds = to_check[(motif,match[1])]

        for md in mds:
            pid = PercentIdentity(md, lflank, rflank) 
            if pid > 90:
                md.found += 1
                return
            if pid > max_pid: max_pid = pid

    if (rc_motif,match[1]) in to_check:
        mds = to_check[(rc_motif,match[1])]
        n_lflank = ReverseComplement(rflank)
        n_rflank = ReverseComplement(lflank)
        
        for md in mds:
            pid = PercentIdentity(md, n_lflank, n_rflank)
            if pid > 90:
                md.found += 1
                return
            if pid > max_pid: max_pid = pid

    # Lets look at the match for the reverse complement
    rc_string = ReverseComplement(string)
    motif = rc_match[3]
    rc_motif = ReverseComplement(motif) 

    e1 = rc_match[4]
    s1 = e1 - 20
    if s1 < 0: s1 == 0

    s2 = rc_match[5]
    e2 = s2 + 20
    if e2 > len(rc_string): e2 = len(rc_string)

    lflank = rc_string[s1:e1]
    rflank = rc_string[s2:e2]   

    if (motif,rc_match[1]) in to_check:
        mds = to_check[(motif,rc_match[1])]

        for md in mds:
            pid = PercentIdentity(md, lflank, rflank)
            if pid > 90:
                md.found += 1
                return
            if pid > max_pid: max_pid = pid

    if (rc_motif,rc_match[1]) in to_check:
        mds = to_check[(rc_motif,rc_match[1])]
        n_lflank = ReverseComplement(rflank)
        n_rflank = ReverseComplement(lflank)
        
        for md in mds:
            pid = PercentIdentity(md, n_lflank, n_rflank)
            if pid > 90:
                md.found += 1
                return
            if pid > max_pid: max_pid = pid

    print >> stderr,"No match for %s:%d:%s:%s. Best pid:%2.2f" % \
        (motif,match[1],lflank,rflank, max_pid)

def PrintPerformanceLog(to_check):
    for (motif,copies),mds in to_check.items():
        for md in mds:
            print >> stderr, "Found %d reads supporting %s:%d:%s:%s" % \
                (md.found, motif, copies, md.lflank, md.rflank)
   
def select_str_reads(filenames, num_minimum_copies, flanking_distance,
         illumina_quals, check_strs, only_k_mers, remove_homopolymers,
         debug_flag):
    """Find the reads that harbor STR and satisfy certain conditions.
    """
    # STRs are short sequences of DNA, normally of length 2-6 base pairs, that 
    # are repeated numerous times. The following python pattern can accurately 
    # describe STRs.
    patterns = []
    for k in only_k_mers:
        patterns.append(r'((?=(.{%d}))\2{%d,})' % (k,num_minimum_copies))

    # Is this a run where we want to check the reason why certain STRs were not
    # found? Read the file to save the STRs I am to look out for.
    to_check = ReadSTRList(check_strs)

    for filename in filenames:
        records = fastq(filename)

        for r in records:
            s = r.fastqsequence
            if debug_flag: 
                print >> stderr, "-"*79
                print >> stderr, "%s" % s.name

            # Does this sequence have a short tandem repeat in it?
            match = MatchPattersAndReportBest(s.seq, 
                    patterns, remove_homopolymers, flanking_distance, debug_flag)
            if match == None:
                if debug_flag: print >> stderr, "0 valid matches."
                continue
            if debug_flag: print >> stderr, "Found %s:%d." % (motif,copies)  

            # What about the matches when the read is reverse complemented?
            rc_match = MatchPattersAndReportBest(ReverseComplement(s.seq),
                       patterns, remove_homopolymers, flanking_distance, 
                       debug_flag)
            if rc_match == None:
                if debug_flag: print >> stderr, "0 valid matches."
                continue

            # Lets register the fact that we found this STR associated read.
            if check_strs: RegisterSTR(s.seq, match, rc_match, to_check)

            # Print the read details along with the motif, position of the STR
            PrintSequence(s, match, rc_match, illumina_quals)

        records.close()
        print >> stderr, "Done processing %s" % filename
                    
    if check_strs:
        PrintPerformanceLog(to_check)



