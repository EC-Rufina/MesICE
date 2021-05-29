#!/bin/bash

#you can rename the file of single gene in this way
for filename in coregenes/*.fasta; do
	[ -f "$filename" ] || continue
	mv "$filename" "${filename/_singlegene.txt_coregenes}"
done

#now align them with mafft
for filename in coregenes/*.fasta; do
	name=${filename/coregenes\/}
	mafft --thread 16 --auto "$filename" > coregenes_alinment/al_"$name"
done

#remove gaps
for filename in coregenes_alinment/*.fasta; do
	name=${filename/coregenes_alinment\/}
	goalign clean sites -i "$filename" > coregenes_alinment_cleaned/"$name"
done

#concatenate the alinments 
goalign concat -i coregenes_alinment_cleaned/*.fasta -o clean_concatenated_core_aln.fasta


#build the tree
raxmlHPC -f a -p $RANDOM -x $RANDOM -N 100 -m GTRCATX -T 16 -s clean_concatenated_core_aln.fasta -n mygenusphylogeny 
