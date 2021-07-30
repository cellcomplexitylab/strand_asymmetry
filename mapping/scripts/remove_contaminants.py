#!/usr/bin/env python

'''
Some reads from one mapping sometimes went to the duplicate
experiment. We can easily attribute the barcode to a given
replicate by comparing the number of reads in support for
a given insertion site.
'''

import gzip
import sys

insertion_records = dict()

def read_other_insertion_file(f):
   for line in f:
      brcd, chrom, _, pos, nreads, _, _ = line.split()
      insertion_records[brcd] = [chrom, int(pos), int(nreads)]

def read_self_insertion_file(f):
   for line in f:
      brcd, chrom, _, pos, nreads, _, _ = line.split()
      # Check if barcode was found in the other file.
      if brcd in insertion_records:
         xchrom, xpos, xnreads = insertion_records[brcd]
         if chrom == xchrom and abs(int(pos) - xpos) < 100:
            # Same barcode at the same site. Keep it only
            # if we have more reads than the other file.
            if int(nreads) < xnreads: continue
      # All checked. Keep the barcode.
      sys.stdout.write(line)

if __name__ == "__main__":
   with gzip.open(sys.argv[2]) as f:
      read_other_insertion_file(f)
   with gzip.open(sys.argv[1]) as f:
      read_self_insertion_file(f)
