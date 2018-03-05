#!/usr/bin/env nextflow
# vipr2t.py -1 denv1/232-N701-S504_S31_1.fastq.gz -2 denv1/232-N701-S504_S31_2.fastq.gz -o denv1_232 -r denv1/new_ref.fasta -p denv1/primer_D1.fasta -n denv1_232 --no-run

# docker run --rm -it -v $PWD:$PWD -w $PWD vipr2 vipr2t.py -1 denv1/232-N701-S504_S31_1.fastq.gz -2 denv1/232-N701-S504_S31_2.fastq.gz -o denv1_232 -r denv1/new_ref.fasta -p denv1/primer_D1.fasta -n denv1_232

nextflow vipr2.nf
