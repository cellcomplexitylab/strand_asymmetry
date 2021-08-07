#!/usr/bin/env python

'''
python gather_gRNA_counts.py A1.co.gz ... A4.co.gz  B1.co.gz ... B4.co.gz
                          ---------------------  ---------------------
                                guide #1               guide #2
'''

import gzip
import sys

from collections import defaultdict

# Store values for gRNA1 and gRNA2.
ratio_GC = defaultdict(lambda: [None] * 2)

def process(group, idx):
   AT = defaultdict(int)
   GC = defaultdict(int)
   for fname in group:
      with gzip.open(fname) as f:
         _ = next(f) # Discard header.
         for line in f:
            barcode, _FF, _AT, _GC  = line.split()
            if int(_FF) >= int(_AT) and int(_FF) >= int(_GC):
               continue
            AT[barcode] += int(_AT)
            GC[barcode] += int(_GC)
   for brcd in AT:
      ratio_GC[brcd][idx] = \
         float(GC[brcd]) / (float(GC[brcd]) + float(AT[brcd]))


def main():
   for i in range(2):
      process(sys.argv[(4*i+1):(4*i+5)], i)
   for barcode, counts in ratio_GC.items():
      if None in counts: continue
      print("\t".join([barcode] + [str(round(x,3)) for x in counts]))


if __name__ == "__main__":
   main()
