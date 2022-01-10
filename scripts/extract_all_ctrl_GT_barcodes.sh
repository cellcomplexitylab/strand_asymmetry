#!/bin/bash

zcat mismatches/GT_LA_no_ISceI_rep1_seq.txt.gz | cut -f1 | sort -u
zcat mismatches/GT_LA_no_ISceI_rep2_seq.txt.gz | cut -f1 | sort -u
zcat mismatches/GT_LA_no_ISceI_rep3_seq.txt.gz | cut -f1 | sort -u
zcat mismatches/GT_LA_no_ISceI_rep4_seq.txt.gz | cut -f1 | sort -u
