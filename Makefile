all:
	cd src && $(MAKE)

install:
	mkdir -p bin
	cp src/fastq.so  bin/
	cp src/select_STR_reads.so  bin/
	cp src/select_STR_reads  bin/
	cp src/merge_STR_reads.so  bin/
	cp src/merge_STR_reads  bin/
	cp src/extend_STR_reads  bin/

clean:
	-rm -rf bin
	cd src && $(MAKE) clean
