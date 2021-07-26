#!/usr/bin/env python
# -*- coding:utf-8 -*-

# Standard library packages.
import gzip
import re
import sys

from itertools import izip

# Others.
import seeq

def trimSuffix(matcher, txt):
   return matcher.matchPrefix(txt, False) or ''

########  Mapping Pipeline ###############################################

def extract_reads_from_PE_fastq(fname_iPCR_PE1, fname_iPCR_PE2):
   """This function takes the 2 pair-end sequencing files and extracts
   the barcode making sure that the other read contains the
   transposon."""

   # This is the scarcode that allows to identify which
   # experiment is sequenced (must be CT).
   matcher = seeq.compile('CGCTAATTAATGGAATCATG', 3)

   outf1 = gzip.open('CT_TCT.fasta.gz', 'w')
   outf2 = gzip.open('CT_ACG.fasta.gz', 'w')

   # There are many errors in the index, especially in the
   # first base. The most frequent errors are hard coded
   # in the dictionary so that the reads are written to the
   # proper file.
   outfiles = {
      'TCT': outf1,
      'GCT': outf1,
      'ACT': outf1,
      'ACG': outf2,
      'AGG': outf2,
      'CCG': outf2,
   }

   with gzip.open(fname_iPCR_PE1) as f, gzip.open(fname_iPCR_PE2) as g:
      for lineno,(line1,line2) in enumerate(izip(f,g)):
         # Take sequence lines of the fastq files.
         if lineno % 4 != 1: continue

         brcd = trimSuffix(matcher, line1)
         # If we find a barcode between 13 and 25 nucleotides
         # then the scarcode must have been the right one.
         if len(brcd) < 13 or len(brcd) > 25: continue

         # Remove first 25 nucleotides.
         suff = line2.rstrip()[25:].split('CATG')[0]
         # Cut genome fragment after the first CATG.
         genome = re.sub(r'CATG.*', 'CATG', suff)

         # Avoid short strings that are unmappable.
         if len(genome) < 20:
            genome = 'gatcctgatgctagtgactgatgagctgctgaagctgga'

         # The first 3 nucleotides of the reverse read are the
         # index. Check that it belongs to the right group.
         idx = line2[:3]
         if idx in outfiles:
            outf = outfiles[idx]
            outf.write('>%s\n%s\n' % (brcd,genome))

if __name__ == '__main__':
   extract_reads_from_PE_fastq(sys.argv[1], sys.argv[2])
