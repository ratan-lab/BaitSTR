all:
	cd src && $(MAKE)

install:
	mkdir -p bin
	cp src/fastq.py src/select_STR_reads src/merge_STR_reads src/extend_STR_reads bin/

clean:
	-rm -rf bin
	cd src && $(MAKE) clean
