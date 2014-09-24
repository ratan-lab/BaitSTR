C=../bin

all: select_strs merge_strs extend_strs

select_strs:
	$C/select_STR_reads -i -n 3 -f 29  Illumina_100_500_1.fq Illumina_100_500_2.fq > reads.str.fq 

merge_strs:
	$C/merge_STR_reads 27 reads.str.fq > merged.reads.str.fq

extend_strs:
	$C/extend_STR_reads 4000 20 27 merged.reads.str.fq Illumina_100_500_1.fq Illumina_100_500_2.fq > contigs.str.fa

.PHONY: clean

clean:
	@-rm reads.str.fq merged.reads.str.fq contigs.str.fa 
