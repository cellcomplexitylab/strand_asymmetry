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
   if flag == 'TG_AC_TC':
      matchers = {
         'TG': seeq.compile('GCTAGCAGTCAGGAATCATG', 3),
         'AC': seeq.compile('GCTAGCTCCGCAGAATCATG', 3),
         'TC': seeq.compile('GCTAGCTCGTTGGAATCATG', 3),
      }
   elif flag == 'CA':
      matchers = {
         'CA': seeq.compile('GCTAGCGCGCGTGAATCATG', 3),
      }
   else:
      raise Exception('wrong flag')
   
   if flag == 'TG_AC_TC':
      indexes = {
         'TG': frozenset(['AAC', 'ACA', 'AGG', 'TTC']),
         'AC': frozenset(['ACT', 'ATC', 'TGA', 'TGT']),
         'TC': frozenset(['ATT', 'CCG', 'TAA', 'TGC']),
      }
   elif flag == 'CA':
      indexes = {
         'CA': frozenset(['ACT', 'AAC', 'CCG', 'TTC']),
      }
   else:
      raise Exception('wrong flag')

   # Assign all valid triplets to a single fasta file for
   # the CT mismatch. Other files can be properly demultiplexed.
   if flag == 'TG_AC_TC':
      outfiles = {
         ('TG','AAC'): gzip.open('TG_AAC.fasta.gz', 'w'),
         ('TG','ACA'): gzip.open('TG_ACA.fasta.gz', 'w'),
         ('TG','AGG'): gzip.open('TG_AGG.fasta.gz', 'w'),
         ('TG','TTC'): gzip.open('TG_TTC.fasta.gz', 'w'),

         ('AC','ACT'): gzip.open('AC_ACT.fasta.gz', 'w'),
         ('AC','ATC'): gzip.open('AC_ATC.fasta.gz', 'w'),
         ('AC','TGA'): gzip.open('AC_TGA.fasta.gz', 'w'),
         ('AC','TGT'): gzip.open('AC_TGT.fasta.gz', 'w'),

         ('TC','ATT'): gzip.open('TC_ATT.fasta.gz', 'w'),
         ('TC','CCG'): gzip.open('TC_CCG.fasta.gz', 'w'),
         ('TC','TAA'): gzip.open('TC_TAA.fasta.gz', 'w'),
         ('TC','TGC'): gzip.open('TC_TGC.fasta.gz', 'w'),
      }
   elif flag == 'CA':
      outfiles = {
         ('CA','ACT'): gzip.open('CA_ACT.fasta.gz', 'w'),
         ('CA','AAC'): gzip.open('CA_AAC.fasta.gz', 'w'),
         ('CA','CCG'): gzip.open('CA_CCG.fasta.gz', 'w'),
         ('CA','TTC'): gzip.open('CA_TTC.fasta.gz', 'w'),
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
