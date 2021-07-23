#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from itertools import combinations

def makeset(f):
   # Collect barcodes from open file.
   S = set()
   for line in f:
      # The barcode is the first element.
      S.add(line.split()[0])
   return S

def main(fname_list):
   L = list()
   # Go through the files one by one and
   # Put their barcodes in the list.
   for fname in fname_list:
      with open(fname) as f:
         # Collect their barcodes.
         L.append(makeset(f))

   BL = set()

   #Any barcode that is in more than one
   # experiments is in the black list.
   for S1,S2 in combinations(L,2):
      BL.update(S1 & S2)

   for brcd in BL:
      print brcd

if __name__ == '__main__':
   main(sys.argv[1:])
