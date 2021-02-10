run: 
	gcc -DSAMPLE_NUM=`cat SAMPLE_NUM` CM3600_12_96.c bm.c controlbits.c decap.c decrypt.c encap.c encrypt.c evaluation.c permutation.c little_endian.c field_ordering.c fq.c keygen.c pk_gen.c poly_gen.c randnum.c syndrome.c -o CM3600_12_96 -lkeccak -lcrypto
	./CM3600_12_96