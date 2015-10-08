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
all: install

# compile the version allowing upto 64mers
large: CFLAGS  += -DLarge 
large: CPFLAGS += -DLarge 
large: install

.PHONY: clean archive

clean:
	-rm *.o
	-rm -rf *.dSYM
	-rm *.pyc
	-rm -rf ../bin
	-rm *.so
	-rm select_STR_reads.c
	-rm merge_STR_reads.c
	-rm fastq.c
	-rm extend_STR_reads

install: fastq.pyx \
		 select_STR_reads select_STR_reads.pyx \
		 merge_STR_reads merge_STR_reads.pyx \
		 utilities.h utilities.c \
		 sllist.h sllist.c \
		 clparsing.h clparsing.c \
    	 kmer.h kmer.c \
    	 murmur_hash.h murmur_hash.c \
    	 bloom_filter.h bloom_filter.c \
    	 fastq_seq.h fastq_seq.c \
		 sparse_hash.h \
		 extend_STR_reads.c	
	python setup.py build_ext --inplace
	rm -rf build
	$(CC)  $(CFLAGS) -c utilities.c
	$(CC)  $(CFLAGS) -c sllist.c
	$(CC)  $(CFLAGS) -c clparsing.c
	$(CC)  $(CFLAGS) -c kmer.c
	$(CC)  $(CFLAGS) -c murmur_hash.c
	$(CC)  $(CFLAGS) -c bloom_filter.c
	$(CC)  $(CFLAGS) -c fastq_seq.c
	$(CC1) $(CPFLAGS)  \
		-o extend_STR_reads \
        utilities.o sllist.o clparsing.o kmer.o murmur_hash.o bloom_filter.o \
	    fastq_seq.o \
		extend_STR_reads.c -lz -lm
	mkdir -p ../bin
	-rm select_STR_reads.c
	-rm merge_STR_reads.c
	-rm fastq.c
	-rm *.o
	-rm -rf *.dSYM
