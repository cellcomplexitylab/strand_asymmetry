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
   # experiment is sequenced (T:G, A:C, T:C or C:A mismatch).
   if flag == 'AG_AC_TC':
      matchers = {
         'AG': seeq.compile('GCTAGCAGTCAGGAATCATG', 3),
         'TC': seeq.compile('GCTAGCTCCGCAGAATCATG', 3),
         'AC': seeq.compile('GCTAGCTCGTTGGAATCATG', 3),
      }
   elif flag == 'GT':
      matchers = {
         'GT': seeq.compile('GCTAGCGCGCGTGAATCATG', 3),
      }
   else:
      raise Exception('wrong flag')
   
   if flag == 'AG_AC_TC':
      indexes = {
         'AG': frozenset(['AAC', 'ACA', 'AGG', 'TTC']),
         'TC': frozenset(['ACT', 'ATC', 'TGA', 'TGT']),
         'AC': frozenset(['ATT', 'CCG', 'TAA', 'TGC']),
      }
   elif flag == 'GT':
      indexes = {
         'GT': frozenset(['ACT', 'AAC', 'CCG', 'TTC']),
      }
   else:
      raise Exception('wrong flag')

   # Assign all valid triplets to a single fasta file for
   # the CT mismatch. Other files can be properly demultiplexed.
   if flag == 'AG_AC_TC':
      outfiles = {
         ('AG','AAC'): gzip.open('AG_AAC.fasta.gz', 'w'),
         ('AG','ACA'): gzip.open('AG_ACA.fasta.gz', 'w'),
         ('AG','AGG'): gzip.open('AG_AGG.fasta.gz', 'w'),
         ('AG','TTC'): gzip.open('AG_TTC.fasta.gz', 'w'),

         ('TC','ACT'): gzip.open('TC_ACT.fasta.gz', 'w'),
         ('TC','ATC'): gzip.open('TC_ATC.fasta.gz', 'w'),
         ('TC','TGA'): gzip.open('TC_TGA.fasta.gz', 'w'),
         ('TC','TGT'): gzip.open('TC_TGT.fasta.gz', 'w'),

         ('AC','ATT'): gzip.open('AC_ATT.fasta.gz', 'w'),
         ('AC','CCG'): gzip.open('AC_CCG.fasta.gz', 'w'),
         ('AC','TAA'): gzip.open('AC_TAA.fasta.gz', 'w'),
         ('AC','TGC'): gzip.open('AC_TGC.fasta.gz', 'w'),
      }
   elif flag == 'GT':
      outfiles = {
         ('GT','ACT'): gzip.open('GT_ACT.fasta.gz', 'w'),
         ('GT','AAC'): gzip.open('GT_AAC.fasta.gz', 'w'),
         ('GT','CCG'): gzip.open('GT_CCG.fasta.gz', 'w'),
         ('GT','TTC'): gzip.open('GT_TTC.fasta.gz', 'w'),
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
