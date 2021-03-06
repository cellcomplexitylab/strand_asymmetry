#!/usr/bin/env python

import gzip
import os
import sys

from collections import defaultdict 

L = "GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC"
R = "TGCAACGAATTCATTAGTGCGGATGATCTTGTCGGTGAAGATCACGCTGTCCTCGGGGAAGCCGGTGCCCACCACCTTGAAGTCGCCGATCA"

def has_one_mismatch(S1, S2):
   if (S1 == S2) or (len(S1) != len(S2)): return False
   differences = 0
   for i in range(len(S1)):
      if S1[i] != S2[i]: differences += 1
      if differences > 1: return False
   # Not equal, same length, not more than 1 difference...
   return True


def is_probable_meCpG(seq, REF, fwd):
   for i in range(len(seq)):
      if seq[i] == REF[i]: continue
      if fwd and REF[i:i+2] == "CG" and seq[i:i+2] == "TG":
         return True
      elif not fwd and i > 0 and REF[i-1:i+1] == "CG" and seq[i-1:i+1] == "CA":
         return True
   return False


def highlight(seq, REF):
   if len(seq) != len(REF): return seq.lower()
   return "".join([seq[i].lower() if seq[i] != REF[i] else seq[i] for i in range(len(REF))])


if __name__ == "__main__":
   left_variants = defaultdict(lambda: defaultdict(int))
   right_variants = defaultdict(lambda: defaultdict(int))

   umi_count_per_barcode = defaultdict(int)
   with gzip.open(sys.argv[1]) as f:
      for line in f:
         items = line.decode("ascii").split()
         umi_count_per_barcode[items[0]] += 1
         if has_one_mismatch(items[2], L) and int(items[3]) > 2:
            left_variants[items[0]][items[2]] += 1
         if has_one_mismatch(items[5], R) and int(items[6]) > 2:
            right_variants[items[0]][items[5]] += 1

   # At least two UMIs with the same variant, representing 5% of the total.
   for brcd in left_variants:
      if len(left_variants[brcd].keys()) != 1: continue
      for variant, count in left_variants[brcd].items():
         if is_probable_meCpG(variant, L, True):
            print brcd, highlight(variant, L), umi_count_per_barcode[brcd], count, 1
         else:
            print brcd, highlight(variant, L), umi_count_per_barcode[brcd], count, 0

   for brcd in right_variants:
      if len(right_variants[brcd].keys()) != 1: continue
      for variant, count in right_variants[brcd].items():
         if is_probable_meCpG(variant, R, False):
            print brcd, highlight(variant, R), umi_count_per_barcode[brcd], count, 1
         else:
            print brcd, highlight(variant, R), umi_count_per_barcode[brcd], count, 0
