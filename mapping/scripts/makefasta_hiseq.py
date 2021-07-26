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

def extract_reads_from_PE_fastq(fname_iPCR_PE1, fname_iPCR_PE2, flag):
   """This function takes the 2 pair-end sequencing files and extracts the
   barcode making sure that the other read contains the transposon."""

   # Those are the scarcodes that allow to identify which
   # experiment is sequenced (CA, GA, GT or TC mismatch).
   if flag == 'CA_GT_GA':
      matchers = {
         'CA': seeq.compile('GCTAGCAGTCAGGAATCATG', 3),
         'GA': seeq.compile('GCTAGCTCGTTGGAATCATG', 3),
         'GT': seeq.compile('GCTAGCTCCGCAGAATCATG', 3),
      }
   elif flag == 'TC':
      matchers = {
         'TC': seeq.compile('GCTAGCGCGCGTGAATCATG', 3),
      }
   else:
      raise Exception('wrong flag')
   
   if flag == 'CA_GT_GA':
      indexes = {
         'CA': frozenset(['AAC', 'ACA', 'AGG', 'TTC']),
         'GA': frozenset(['ATT', 'CCG', 'TAA', 'TGC']),
         'GT': frozenset(['ACT', 'ATC', 'TGA', 'TGT']),
      }
   elif flag == 'TC':
      indexes = {
         'TC': frozenset(['ACT', 'AAC', 'CCG', 'TTC']),
      }
   else:
      raise Exception('wrong flag')

   # Assign all valid triplets to a single fasta file for
   # the CT mismatch. Other files can be properly demultiplexed.
   if flag == 'CA_GT_GA':
      outfiles = {
         ('CA','AAC'): gzip.open('CA_AAC.fasta.gz', 'w'),
         ('CA','ACA'): gzip.open('CA_ACA.fasta.gz', 'w'),
         ('CA','AGG'): gzip.open('CA_AGG.fasta.gz', 'w'),
         ('CA','TTC'): gzip.open('CA_TTC.fasta.gz', 'w'),

         ('GT','ACT'): gzip.open('GT_ACT.fasta.gz', 'w'),
         ('GT','ATC'): gzip.open('GT_ATC.fasta.gz', 'w'),
         ('GT','TGA'): gzip.open('GT_TGA.fasta.gz', 'w'),
         ('GT','TGT'): gzip.open('GT_TGT.fasta.gz', 'w'),

         ('GA','ATT'): gzip.open('GA_ATT.fasta.gz', 'w'),
         ('GA','CCG'): gzip.open('GA_CCG.fasta.gz', 'w'),
         ('GA','TAA'): gzip.open('GA_TAA.fasta.gz', 'w'),
         ('GA','TGC'): gzip.open('GA_TGC.fasta.gz', 'w'),
      }
   elif flag == 'TC':
      outfiles = {
         ('TC','ACT'): gzip.open('TC_ACT.fasta.gz', 'w'),
         ('TC','AAC'): gzip.open('TC_AAC.fasta.gz', 'w'),
         ('TC','CCG'): gzip.open('TC_CCG.fasta.gz', 'w'),
         ('TC','TTC'): gzip.open('TC_TTC.fasta.gz', 'w'),
      }
   else:
      raise Exception('wrong flag')

   # End of the pT2 transposon sequence.
   pT2 = seeq.compile('AAACTTCCGACTTCAACTGTA', 3)

   try:
      with gzip.open(fname_iPCR_PE1) as f, gzip.open(fname_iPCR_PE2) as g:
         # Aggregate iterator of f,g iterators -> izip(f,g).
         for lineno,(line1,line2) in enumerate(izip(f,g)):
            # Take sequence lines of the fastq file.
            if lineno % 4 != 1: continue

            # Use the scarcode to identify the experiment.
            for MM,matcher in matchers.items():
               brcd = trimSuffix(matcher, line1)
               # If we find a barcode between 13 and 25 nucleotides
               # then the scarcode must have been the right one.
               if len(brcd) < 13 or len(brcd) > 25: continue

               # Find pT2 on the reverse read. Abort if we cannot.
               suff = pT2.matchSuffix(line2.rstrip(), False)
               if suff is None: continue

               # Cut genome fragment after the first CATG.
               genome = re.sub(r'CATG.*', 'CATG', suff)

               # Avoid short strings that are unmappable.
               if len(genome) < 20:
                  genome = 'gatcctgatgctagtgactgatgagctgctgaagctgga'

               # The first 3 nucleotides of the reverse read are the
               # index. Check that it belongs to the right group.
               idx = line2[:3]
               if idx in indexes[MM]:
                  outf = outfiles[(MM,idx)]
                  outf.write('>%s\n%s\n' % (brcd,genome))

               # If the script reaches this point, the sample was
               # already identified (even though nothing may be
               # printed) so there is no need to check the others.
               break

   finally:
      for f in outfiles.values():
         f.close()

if __name__ == '__main__':
   extract_reads_from_PE_fastq(sys.argv[1], sys.argv[2], sys.argv[3])
