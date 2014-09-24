#!/usr/bin/env python

"""
Read and write fastq format files.
"""

from sys import argv

import gzip
import string

__author__ = "Aakrosh Ratan"
__email__  = "ratan@bx.psu.edu"

class fastqsequence:
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq  = seq
        self.qual = qual

    def fq2seqfa(self):
        str  = ">%s\n" % self.name
        str += "%s\n" % self.seq
        return str[:-1]

    def fq2qfa(self, offset):
        str  = ">%s\n" % self.name
        for q in self.qual:
            str += "%d " % (ord(q) - offset)
        return str[:-1]

    def __str__(self):
        str = "@%s\n" % self.name
        str += "%s\n" % self.seq
        str += "+\n"
        str += "%s\n" % self.qual
        return str[:-1]

    def __len__(self):
        assert len(self.seq) == len(self.qual)
        return len(self.seq)

    def reverse_complement(self):
        complement = string.maketrans('atcgnATCGN', 'tagcnTAGCN')
        self.seq = self.seq.translate(complement)[::-1]
        self.qual = self.qual[::-1]

class fastq:
    def __init__(self, filename):
        if filename[-3:] == ".gz":
            self.file = gzip.open(filename, "r")
        else:
            self.file = open(filename, "r")
        self.fastqsequence = None

    def __del__(self):
        self.file.close()

    def __iter__(self):
        return self

    @staticmethod
    def checkbases(line):
        for base in line[:-1]:
            if base not in ["a","c","g","t","n","A","C","G","T","N","."]:
                return False
        return True

    def close(self):
        self.file.close()

    def next(self):
        line = self.file.readline()
        if not line:
            self.file.close()
            raise StopIteration
        assert line[0] == "@", "header should begin with a @ in a fastq file"
        name = line[1:-1]
        
        line = self.file.readline()
        assert fastq.checkbases(line) == True, \
        "read %s should only have ACGTN" % (self.name)
        seq = line[:-1].upper()
        
        line = self.file.readline()
        assert line[0] == "+", \
        "separator line in fastq file should begin with a +"

        line = self.file.readline()
        qual = line[:-1]

        self.fastqsequence = fastqsequence(name, seq, qual)
    
        return self

if __name__ == "__main__": 
    records = fastq(argv[1])
        
    for r in records:
        s = r.fastqsequence
        print s

    records.close()
