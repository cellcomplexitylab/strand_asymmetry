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


def makedict(f):
   D = dict()
   for line in f:
      items = line.split()
      D[items[0]] = items
   return D


def main(f, fname, barcode_dict_1, barcode_dict_2, BL):

   mmcode = 'GT'
   tcode = 24
   lacode = '6xPCR'
   ctrl = 'test'

   mm = ['FF', 'AT', 'GC']

   rep = GC1 = GC2 = None
   # Discard header.
   next(f)
   for line in f: 
      #if line[0].isspace(): continue
      bcd,FF,AT,GC = line.split()
      # Skip barcodes in the black list.
      if bcd in BL: continue
      # Skip truncated barcodes.
      if len(bcd) < 18: continue
      scores = [float(a) for a in (FF, AT, GC)]
      winner = max((0,1,2), key=lambda x: scores[x])
      # Winning with only 1 UMI is cheating.
      if winner == 0 or scores[winner] < 2: continue
      ratio = scores[1] / (scores[1] + scores[2])
      if bcd in barcode_dict_1:
         the_dict = barcode_dict_1
         rep = 1
      elif bcd in barcode_dict_2:
         the_dict = barcode_dict_2
         rep = 2
      else:
         continue
      the_fields = the_dict[bcd]
      chrom  = the_fields[1]
      strand = the_fields[2]
      pos    = the_fields[3]
      GC1    = the_fields[5]
      GC2    = the_fields[6]
      print "%s\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % \
         (bcd, ratio, mmcode, tcode, lacode,
               ctrl, rep, GC1, GC2, chrom, strand, pos, fname)



if __name__ == '__main__':
   with gzip.open(sys.argv[1]) as f:
      BL = get_black_list(f)
   with gzip.open(sys.argv[3]) as f:
      barcode_dict_1 = makedict(f)
   with gzip.open(sys.argv[4]) as f:
      barcode_dict_2 = makedict(f)

   # debug info.
   sys.stderr.write("processing {}\n".format(sys.argv[2]))
   with gzip.open(sys.argv[2]) as f:
      fname = os.path.basename(sys.argv[2])
      main(f, fname, barcode_dict_1, barcode_dict_2, BL)
