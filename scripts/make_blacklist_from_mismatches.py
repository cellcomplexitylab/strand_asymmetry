#!/usr/bin/env python
# -*- coding:utf-8 -*-

import gzip
import sys

BL = set()

def main(f):

   # Discard header.
   next(f)
   for line in f: 
      bcd,FF,AT,GC = line.split()
      scores = [float(a) for a in (FF, AT, GC)]
      # Must have more than one UMI against FF.
      if scores[0] >= max(scores)-1:continue
      winner = max((0,1,2), key=lambda x: scores[x])
      if scores[winner] < 2: continue
      BL.add(bcd)


if __name__ == '__main__':
   for fname in sys.argv[1:]:
      # debug info.
      sys.stderr.write("processing {}\n".format(fname))
      with(gzip.open(fname)) as f:
         main(f)
   for bcd in BL:
      print bcd
