BINARIES = extend_STR_reads

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
all: $(BINARIES)

# compile the version allowing upto 64mers
large: CFLAGS  += -DLarge 
large: CPFLAGS += -DLarge 
large: $(BINARIES)

.PHONY: clean archive

clean:
	-rm $(BINARIES) *.o
	-rm -rf *.dSYM
	-rm *.pyc

extend_STR_reads: utilities.h utilities.c \
				  sllist.h sllist.c \
				  clparsing.h clparsing.c \
                  kmer.h kmer.c \
                  murmur_hash.h murmur_hash.c \
                  bloom_filter.h bloom_filter.c \
                  fastq.h fastq.c \
				  sparse_hash.h \
				  extend_STR_reads.c	
	$(CC)  $(CFLAGS) -c utilities.c
	$(CC)  $(CFLAGS) -c sllist.c
	$(CC)  $(CFLAGS) -c clparsing.c
	$(CC)  $(CFLAGS) -c kmer.c
	$(CC)  $(CFLAGS) -c murmur_hash.c
	$(CC)  $(CFLAGS) -c bloom_filter.c
	$(CC)  $(CFLAGS) -c fastq.c
	$(CC1) $(CPFLAGS) -I \
		-o extend_STR_reads \
        utilities.o sllist.o clparsing.o kmer.o murmur_hash.o bloom_filter.o \
	    fastq.o \
		extend_STR_reads.c -lz -lm
