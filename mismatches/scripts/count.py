#!/usr/bin/env python
# -*- coding:utf-8 -*-

import datetime
import gzip
import random
import seeq
import sys

from collections import defaultdict

MIN_BCD_LEN = 5

class NonCanonicalBarcodeException(Exception):
   pass

class WrongScarcodeException(Exception):
   pass

class TieException(Exception):
   pass


#  ------      IMPORTANT NOTE      ------

# The two-letter codes are the key nucleotides in the reverse
# order of how they appear in the product of gene synthesis.
# The mismatches as they happen during the repair process are
# as follows:
#    CT    G:T
#    CA    G:A
#    GA    C:A
#    GT    C:T
#    TC    A:C

SCARCODES = {
   'CT': 'CGCTAATTAATG',
   'CA': 'GCTAGCAGTCAG',
   'GA': 'GCTAGCTCGTTG',
   'GT': 'GCTAGCTCCGCA',
   'TC': 'GCTAGCGCGCGT',
}


class ScarcodeRemover():
   '''Utility to remove and check scarcodes.
   Note that if the first two nucleotides of the scarcode
   are mutated, they will be added to the barcode.
   
   For instance, if the tag has the following structure

            agatgctgagctaggcc tcCTAATTAATG
                barcode         scarcode

   the recovered barcode will be

            agatgctgagctaggcctc
   '''

   def __init__(self, MM):
      # Assign a fixed matcher upon instantiation.
      self.matcher = {
         'CT': seeq.compile(SCARCODES['CT'], 2),
         'CA': seeq.compile(SCARCODES['CA'], 2),
         'GA': seeq.compile(SCARCODES['GA'], 2),
         'GT': seeq.compile(SCARCODES['GT'], 2),
         'TC': seeq.compile(SCARCODES['TC'], 2),
      }[MM]

   def remove(self, bcd_scar):
      # Get the barcode proper by removing the scarcode
      bcd = self.matcher.matchPrefix(bcd_scar, False)
      if bcd is None:
         # In case the scarcode was not
         # found raise an exception
         raise WrongScarcodeException
      return bcd



class EventCounter:
   '''Basic information about the experiment.'''

   def __init__(self, normalizer, info):

      # Events are recorded in a nested counter.
      # barcode: umi: variant: count.
      self.events = defaultdict(lambda:
            defaultdict(lambda: defaultdict(int)))
      self.seq1 = defaultdict(lambda: defaultdict(list))
      self.seq2 = defaultdict(lambda: defaultdict(list))

      self.info = info

      self.normalizer = normalizer
      self.normalize_tag = normalizer.normalize

      assign = {
         #          FF          AT          GC
         'GA': (('A', 'C'), ('A', 'T'), ('G', 'C')),
         'GT': (('T', 'C'), ('T', 'A'), ('G', 'C')),
         'CA': (('A', 'G'), ('A', 'T'), ('C', 'G')),
         'CT': (('T', 'G'), ('T', 'A'), ('C', 'G')),
         'TC': (('G', 'T'), ('A', 'T'), ('G', 'C')),
      }

      # Use the starcode file of the normalizer in order
      # to retrieve a representative fraction of the
      # scarcode and to infer the mismatch code.
      self.info.get_MMcode(normalizer)

      self.MM = assign[self.info.MMcode]
      self.scar = ScarcodeRemover(self.info.MMcode)


   def consensus(self, seqlist):
      '''Find the consensus sequence among a list of sequences.'''

      # Remove the sequences with a different size.
      L = max([len(seq) for seq in seqlist])
      seqlist = [seq for seq in seqlist if len(seq) == L]

      con = bytearray('N'*L)
      score = 1.0
      for i in range(L):
         counter = defaultdict(float)
         for seq in seqlist: counter[seq[i]] += 1.0
         winner = max(counter, key=counter.get)
         con[i] = winner
         this_score = counter[winner] / sum(counter.values())
         if this_score < score: score = this_score

      return str(con), score


   def outputseq(self, seqoutf, variant, bcd, umi):
      '''Write the consensus sequences to file for given UMI.'''

      decode = {
         ('A', 'C'): 'FF',
         ('T', 'C'): 'FF',
         ('A', 'G'): 'FF',
         ('T', 'G'): 'FF',
         ('A', 'T'): 'AT',
         ('T', 'A'): 'AT',
         ('G', 'C'): 'GC',
         ('C', 'G'): 'GC',
      }

      SEQ1,s1 = self.consensus(self.seq1[bcd][umi])
      SEQ2,s2 = self.consensus(self.seq2[bcd][umi])

      n1 = len(self.seq1[bcd][umi])
      n2 = len(self.seq2[bcd][umi])

      # Output to file.
      if variant in decode:
         seqoutf.write('%s\t%s\t%s\t%d\t%f\t%s\t%d\t%f\n' % \
               (bcd, decode[variant], SEQ1, n1, s1, SEQ2, n2, s2))


   def count(self, f, outf=sys.stdout, seqoutf=None):
      '''Processing function to convert the reads to repair events.'''

      # Create a filter for UMIs used multiple times.
      reverse_lookup = defaultdict(lambda: defaultdict(int))

      for line in f:
         self.info.nreads += 1
         # Assume that preprocessed file is tab-separated.
         tag, V1, V2, SEQ1, SEQ2 = line.split()
         try:
            bcd_scar, umi = self.normalize_tag(tag)
            bcd = self.scar.remove(bcd_scar)
            reverse_lookup[umi][bcd] += 1
         except NonCanonicalBarcodeException:
            continue
         except WrongScarcodeException:
            self.info.wrong_scarcode += 1
            continue
         if len(bcd) < MIN_BCD_LEN:
            self.info.barcode_too_short += 1
            continue
         self.events[bcd][umi][(V1,V2)] += 1
         self.seq1[bcd][umi].append(SEQ1)
         self.seq2[bcd][umi].append(SEQ2)

      # Header of the file.
      outf.write('barcode\tFF\tAT\tGC\n')

      # Run through the barcodes one at a time.
      for bcd,dict_of_umis in self.events.items():
         counter = defaultdict(int)
         for umi,dict_of_variants in dict_of_umis.items():
            # Discard all the events with a single read.
            for key in [a for a,b in dict_of_variants.items() if b < 2]:
               dict_of_variants.pop(key)
               self.info.too_few_reads += 1
            # If nothing is left just move on.
            if len(dict_of_variants) < 1:
               continue
            # Discard UMIs used multiple times (with at least two reads).
            S = [1 for (a,b) in reverse_lookup[umi].items() if b > 1]
            if len(S) > 1:
               self.info.non_unique_UMI += sum(dict_of_variants.values())
               continue
            variant = self.info.normalize_variant(bcd,
                  umi, dict_of_variants)
            counter[variant] += 1
            if seqoutf:
               self.outputseq(seqoutf, variant, bcd, umi)
         # If all UMIs were lost, the counter
         # is empty and there is nothing to show.
         if not counter:
            continue
         # Show counts and sequences (FF, AT, GC)
         counts = '\t'.join(['%d' % counter[a] for a in self.MM])
         outf.write('%s\t%s\n' % (bcd, counts))



class TagNormalizer:
   '''Correct reading errors. This object encapsulates a starcode
   file in order to normalize barcodes and UMIs. It also contains
   an iterator to read the starcode file, so for practical purposes
   a 'TagNormalizer' "is" a starcode file.'''

   def __init__(self, f):
      '''Construct a normalizer from open Starcode file.'''

      self.canonical = dict()
      # Build a dictionary from Starcode file.
      for line in f:
        centroid,discard,others = line.split()
        for bcd in others.split(','):
            self.canonical[bcd] = centroid

   def __iter__(self):
      '''Make this class iterable so that statements of the form
      'for tag in TagNormalizer' iterate over the canonical
      tags.'''
      return iter(set(self.canonical.values()))

   def normalize(self, tag):
      '''Replace the tag by its canonical sequence (where errors
      are reverted). Sperate the barcode from the UMI and returns
      both as a pair. In case of failure, a 'NonCanonicalBarcodeException' 
      is raised.'''

      try:
         barcode, umi = self.canonical[tag].split('ATGCTACG')
      except (KeyError, ValueError):
         # In case the spacer is not present or the tag
         # is not indexed, report aberrant tag.
         raise NonCanonicalBarcodeException

      return barcode, umi


class CountingInfo:
   '''Store general information regarding the couting process for
   quality control and troubleshooting.'''

   # MMcodes and their scarcodes.
   MM = {
      SCARCODES['CT']: 'CT',
      SCARCODES['CA']: 'CA',
      SCARCODES['GA']: 'GA',
      SCARCODES['GT']: 'GT',
      SCARCODES['TC']: 'TC',
   }


   def __init__(self, fname1, fname2):
      self.fname1 = fname1
      self.fname2 = fname2

      self.MMcode = None

      self.nreads = 0
      self.used_reads = 0

      self.vart_conflicts = dict()
      self.wrong_scarcode = 0
      self.barcode_too_short = 0
      self.too_few_reads = 0
      self.non_unique_UMI = 0
      self.minority_report = 0


   def normalize_variant(self, bcd, umi, dict_of_variants):
      '''In several cases the barcode/UMI pair is associated with
      multiple variants. Read errors and template switching during the
      PCR can cause this. In any event, a barcode/UMI pair corresponds
      to a unique molecule and therefore a single repair event that we
      can try to infer.'''

      # Only one variant.
      if len(dict_of_variants) == 1:
         return dict_of_variants.keys()[0]

      # The most frequent variant is assigned to the tag.
      # In case of tie, assign a winner at random.
      maxval = max(dict_of_variants.values())
      winners = [a for a,b in dict_of_variants.items() if b == maxval]
      variant = random.choice(winners)

      # If variants are not unanimous for the barcode/UMI pair,
      # create an exception entry for the records.
      self.vart_conflicts[(bcd,umi)] = dict_of_variants.copy()
      total = sum(dict_of_variants.values())
      kept = dict_of_variants[variant]
      self.minority_report += (total - kept)
      self.used_reads += kept

      return variant


   # Side effect: updates 'self.MMcode'
   def get_MMcode(self, tags):
      '''The barcodes have a scar (the scarcode) that identifies the
      kind of mismatch that is generated during the DNA repair. This
      information can be extracted to specify the variants that are
      expected in a given set of reads.'''

      # Scarcodes and the corresponding mismatches.
      refscars = {
         SCARCODES['CT']: 0,
         SCARCODES['CA']: 0,
         SCARCODES['GA']: 0,
         SCARCODES['GT']: 0,
         SCARCODES['TC']: 0,
      }

      # Note: the class 'TagNormalizer' was made
      # iterable for this 'for' loop to work.
      for tag in tags:
         try:
            bcd,umi = tag.split('ATGCTACG')
            scarcode = bcd[-12:]
         except ValueError:
            # The tag may have been tampered with
            # during the sequence clustering.
            continue
         if scarcode in refscars:
            refscars[scarcode] += 1
         # Count maximum 10000 barcodes.
         if sum(refscars.values()) > 10000:
            break

      winner = max(refscars, key=refscars.get)

      # Modify the object in place.
      self.MMcode = self.MM[winner]

      return self.MMcode


   def write_to_file(self, f):
      '''Log the counting process.'''

      # Time stamp.
      dt = datetime.datetime
      f.write(dt.strftime(dt.now(),
         '%Y-%m-%d %H:%M:%S Counting summary\n'))

      # Processed files.
      f.write('%s\n' % self.fname1)
      f.write('%s\n' % self.fname2)

      # Mismatch type.
      f.write('MM type: %s\n' % self.MMcode)

      # Total and percent used reads.
      percent = 100 * float(self.used_reads) / self.nreads
      f.write('Used reads:\t%d (%.2f%%)\n' % \
            (self.used_reads, percent))

      # Total and percent reads lost.
      f.write('Reads lost to:\n')
      percent = 100 * float(self.wrong_scarcode) / self.nreads
      f.write('  Wrong scarcode:\t%d (%.2f%%)\n' % \
            (self.wrong_scarcode, percent))
      percent = 100 * float(self.barcode_too_short) / self.nreads
      f.write('  Barcode too short:\t%d (%.2f%%)\n' % \
            (self.barcode_too_short, percent))
      percent = 100 * float(self.too_few_reads) / self.nreads
      f.write('  Too few reads:\t%d (%.2f%%)\n' % \
            (self.too_few_reads, percent))
      percent = 100 * float(self.non_unique_UMI) / self.nreads
      f.write('  Non unique UMI:\t%d (%.2f%%)\n' % \
            (self.non_unique_UMI, percent))
      percent = 100 * float(self.minority_report) / self.nreads
      f.write('  Minority report:\t%d (%.2f%%)\n' % \
            (self.minority_report, percent))

      # End of report.
      f.write('---\n')



def main(fname1, fname2, info, seqoutf=None):
   # Instantiate and run tools for the analysis.
   with gzip.open(fname1) as f, gzip.open(fname2) as g:
      # Create a normalizer encapsulating the starcode file.
      normalizer = TagNormalizer(g)
      counter = EventCounter(normalizer, info)
      counter.count(f, seqoutf=seqoutf)


if __name__ == '__main__':
   # Initialize CountingInfo object to monitor the process.
   info = CountingInfo(sys.argv[1], sys.argv[2])
   seqoutf = gzip.open(sys.argv[3], 'w') if len(sys.argv) > 3 else None
   try:
      random.seed(123)
      main(sys.argv[1], sys.argv[2], info, seqoutf)
   finally:
      info.write_to_file(open('counting_logs.txt', 'a'))
      if seqoutf: seqoutf.close()
