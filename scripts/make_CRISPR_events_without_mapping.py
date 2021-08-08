#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import gzip
import sys

def get_black_list(f):
   S = set()
   for line in f:
      S.add(line.rstrip())
   return S


def main(f, gRNA, BL):

   mmcode = "GT"
   tcode = 24
   lacode = "6xPCR"
   ctrl = "test"

   mm = ['FF', 'AT', 'GC']

   rep = None
   # Discard header.
   next(f)
   for line in f: 
      bcd,FF,AT,GC = line.split()
      # Skip barcodes in the black list
      if bcd in BL: continue
      # Skip truncated barcodes.
      if len(bcd) < 18: continue
      scores = [float(a) for a in (FF, AT, GC)]
      winner = max((0,1,2), key=lambda x: scores[x])
      # Winning with only 1 UMI is cheating.
      if winner == 0 or scores[winner] < 2: continue
      ratio = scores[1] / (scores[1] + scores[2])
      print "%s\t%.3f\t%s\t%d\t%s\t%s\t%s\t%s" % \
         (bcd, ratio, mmcode, tcode, lacode, ctrl, gRNA, fname)



if __name__ == '__main__':
   with gzip.open(sys.argv[1]) as f:
      BL = get_black_list(f)

   # debug info.
   sys.stderr.write("processing {}\n".format(sys.argv[3]))
   with gzip.open(sys.argv[3]) as f:
      fname = os.path.basename(sys.argv[3])
      main(f, sys.argv[2], BL)
