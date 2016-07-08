CC  = gcc
CC1 = g++

CFLAGS  = -W -Wformat -Wimplicit -Wreturn-type -Wall -Werror \
          -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
          -Wcast-align -Winline -Wextra 
CPFLAGS = -W -Wformat -Wreturn-type -Wall -w \
          -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
          -Wcast-align -Winline -Wextra

CFLAGS += -g -ggdb
CFLAGS += -O3

CPFLAGS += -g -ggdb -w
CPFLAGS += -O3 

# compile the version allowing upto 32mers
all: compile

# compile the version allowing upto 64mers
large: CFLAGS  += -DLarge 
large: CPFLAGS += -DLarge 
large: compile

.PHONY: clean archive

clean:
	cd sparsehash && $(MAKE) clean
	-rm *.o
	-rm -rf *.dSYM
	-rm *.pyc
	-rm -rf ../bin
	-rm *.so
	-rm select_STR_reads.c
	-rm fastq.c
	-rm merge_STR_reads
	-rm extend_STR_reads

compile: fastq.pyx \
		 select_STR_reads select_STR_reads.pyx \
		 utilities.h utilities.c \
		 sllist.h sllist.c \
		 clparsing.h clparsing.c \
    	 kmer.h kmer.c \
    	 murmur_hash.h murmur_hash.c \
    	 bloom_filter.h bloom_filter.c \
    	 fastq_seq.h fastq_seq.c \
		 sparse_word_hash.h \
		 sparse_kmer_hash.h \
		 merge_STR_reads.c \
		 extend_STR_reads.c	
	cd sparsehash && ./configure && $(MAKE)
	python setup.py build_ext --inplace
	rm -rf build
	$(CC)  $(CFLAGS) -c utilities.c
	$(CC)  $(CFLAGS) -c sllist.c
	$(CC)  $(CFLAGS) -c clparsing.c
	$(CC)  $(CFLAGS) -c kmer.c
	$(CC)  $(CFLAGS) -c murmur_hash.c
	$(CC)  $(CFLAGS) -c bloom_filter.c
	$(CC)  $(CFLAGS) -c fastq_seq.c
	$(CC1) $(CPFLAGS) -D'VERSION="$(shell cat VERSION .)"' \
		-o merge_STR_reads \
		-Isparsehash/src \
        utilities.o sllist.o clparsing.o kmer.o murmur_hash.o bloom_filter.o \
	    fastq_seq.o \
		merge_STR_reads.c -lz -lm
	$(CC1) $(CPFLAGS) -D'VERSION="$(shell cat VERSION .)"' \
		-o extend_STR_reads \
		-Isparsehash/src \
        utilities.o sllist.o clparsing.o kmer.o murmur_hash.o bloom_filter.o \
	    fastq_seq.o \
		extend_STR_reads.c -lz -lm
	mkdir -p ../bin
	-rm select_STR_reads.c
	-rm fastq.c
	-rm *.o
	-rm -rf *.dSYM
