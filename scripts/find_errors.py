#!/usr/bin/env python

import gzip
import os
import sys

from collections import defaultdict 

L = "GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC"
R = "TGCAACGAATTCATTAGTGCGGATGATCTTGTCGGTGAAGATCACGCTGTCCTCGGGGAAGCCGGTGCCCACCACCTTGAAGTCGCCGATCA"

def is_probable_meCpG(seq, REF):
   if len(seq) != len(REF): return None
   for i in range(len(seq)):
      if seq[i] == REF[i]: continue
      if REF[i:i+2] == "CG" and seq[i:i+2] == "TG":
         return True
      elif i > 0 and REF[i-1:i+1] == "CG" and seq[i-1:i+1] == "CA":
         return True
   return False

def highlight(seq, REF):
   if len(seq) != len(REF): return seq.lower()
   return "".join([seq[i].lower() if seq[i] != REF[i] else seq[i] for i in range(len(REF))])

if __name__ == "__main__":
   left_variants = defaultdict(lambda: defaultdict(int))
   right_variants = defaultdict(lambda: defaultdict(int))

   MMcode = os.path.basename(sys.argv[1])[:2]

   with gzip.open(sys.argv[1]) as f:
      for line in f:
         items = line.decode("ascii").split()
         if items[2] != L and int(items[3]) > 2:
            left_variants[items[0]][items[2]] += 1
         if items[5] != R and int(items[6]) > 2:
            left_variants[items[0]][items[5]] += 1

   for brcd in left_variants:
      for variant, count in left_variants[brcd].items():
         if len(variant) == len(L) and count > 0 and is_probable_meCpG(variant, L):
            print MMcode, brcd

   for brcd in right_variants:
      for variant, count in right_variants[brcd].items():
         if len(variant) == len(R) and count > 0 and is_probable_meCpG(variant, R):
            print MMcode, brcd
