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


def main(f, fname, BL):

   mmcode = fname[:2]
   tcode = 48 if '_48h_' in fname else 24
   lacode = 'LA' if '_LA_' in fname.upper() else '6xPCR'
   ctrl = 'ctrl' if '_no_ISceI_' in fname else 'test'

   mm = ['FF', 'AT', 'GC']

   rep = None
   # Discard header.
   next(f)
   for line in f: 
      #if line[0].isspace(): continue
      bcd,FF,AT,GC = line.split()
      # Skip barcodes in the black list.
      if bcd in BL: continue
      scores = [float(a) for a in (FF, AT, GC)]
      # Must have more than one UMI against FF.
      if scores[0] >= max(scores)-1:continue
      winner = max((0,1,2), key=lambda x: scores[x])
      if scores[winner] < 2: continue
      ratio = scores[1] / (scores[1] + scores[2])
      print "%s\t%.3f\t%s\t%d\t%s\t%s\t%s" % \
         (bcd, ratio, mmcode, tcode, lacode, ctrl, fname)



if __name__ == '__main__':
   with gzip.open(sys.argv[1]) as f:
      BL = get_black_list(f)

   # debug info.
   sys.stderr.write("processing {}\n".format(sys.argv[2]))
   with gzip.open(sys.argv[2]) as f:
      fname = os.path.basename(sys.argv[2])
      main(f, fname, BL)
