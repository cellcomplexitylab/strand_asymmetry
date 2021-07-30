#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest

from StringIO import StringIO
from textwrap import dedent

import preprocess


class TestLaneInfo(unittest.TestCase):

   def test_constructor(self):
      info = preprocess.LaneInfo('fname1', 'fname2')
      self.assertEqual(info.fname1, 'fname1')
      self.assertEqual(info.fname2, 'fname2')
      self.assertEqual(info.ntotal, 0)
      self.assertEqual(info.naberrant, 0)

   def test_write_to_file(self):
      info = preprocess.LaneInfo('fname1', 'fname2')
      info.ntotal = 1000000
      info.naberrant = 100000

      buffer = StringIO()
      info.write_to_file(buffer)

      txt = '''fname1
         fname2
         Reads lost:\t100000 (10.00 %)
         ---'''

      check = '\n'.join(buffer.getvalue().splitlines()[1:])
      self.assertEqual(check, dedent(' '*9 + txt))


class TestExtractor(unittest.TestCase):

   def test_init(self):
      ex = preprocess.Extractor()
      self.assertIsNone(ex.seq_after_tag)
      self.assertIsNone(ex.seq_before_variant)


   def test_Read1Extractor(self):
      ex = preprocess.Read1Extractor()

      # Test case 1.
      seq = 'aaaaGAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGCctga'

      BCD,SNP1,SEQ1 = ex.extract_all(seq)
      self.assertEqual(BCD, 'aaaa')
      self.assertEqual(SNP1, 'c')
      self.assertEqual(SEQ1, 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC')

      # Test case 2 (10 substitutions)
      seq = 'aaaaGAATCATtAACACCaGCATCaAGAAGaACGAGGACGGCGaCGTaCTGCACGTGtGCTTCAcCTtCCGCTACGAGcCCGGCCGCctga'
      #                 ^      ^     ^     ^            ^   ^         ^      ^  ^          ^

      BCD,SNP1,SEQ1 = ex.extract_all(seq)
      self.assertEqual(BCD, 'aaaa')
      self.assertEqual(SNP1, 'c')
      self.assertEqual(SEQ1, 'GAATCATtAACACCaGCATCaAGAAGaACGAGGACGGCGaCGTaCTGCACGTGtGCTTCAcCTtCCGCTACGAGcCCGGCCGC')

      # Test case 3 (10 deletions)
      seq = 'aaaaGAATCATAACACCGCATCAGAAGACGAGGACGGCGCGTCTGCACGTGGCTTCACTCCGCTACGAGCCGGCCGCctga'

      BCD,SNP1,SEQ1 = ex.extract_all(seq)
      self.assertEqual(BCD, 'aaaa')
      self.assertEqual(SNP1, 'c')
      self.assertEqual(SEQ1, 'GAATCATAACACCGCATCAGAAGACGAGGACGGCGCGTCTGCACGTGGCTTCACTCCGCTACGAGCCGGCCGC')

      # Test case 4 (11 substitutions)
      seq = 'aaaaGAATCATtAACACCaGCATCaAGAAGaACGAGcACGGCGaCGTaCTGCACGTGtGCTTCAcCTtCCGCTACGAGcCCGGCCGCctga'
      #                 ^      ^     ^     ^     ^      ^   ^         ^      ^  ^          ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 5 (11 deletions)
      seq = 'aaaaGAATCATAACACCGCATCAGAAGACGAGACGGCGCGTCTGCACGTGGCTTCACTCCGCTACGAGCCGGCCGCctga'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 6 (missing TGA).
      seq = 'aaaaGAATCATAACACCGCATCAGAAGACGAGACGGCGCGTCTGCACGTGGCTTCACTCCGCTACGAGCCGGCCGCctta'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)



   def test_Read2Extractor(self):
      ex = preprocess.Read2Extractor()

      # Test case 1.
      seq = 'aaaaTGCAACGAATTCATTAGTGCGGATGATCTTGTCGGTGAAGATCACGCTGTCCTCGGGGAAGCCGGTGCCCACCACCTTGAAGTCGCCGATCAcgcg'

      UMI,SNP2,SEQ2 = ex.extract_all(seq)
      self.assertEqual(UMI, 'aaaa')
      self.assertEqual(SNP2, 'c')
      self.assertEqual(SEQ2, 'TGCAACGAATTCATTAGTGCGGATGATCTTGTCGGTGAAGATCACGCTGTCCTCGGGGAAGCCGGTGCCCACCACCTTGAAGTCGCCGATCA')

      # Test case 2 (10 substitutions).
      seq = 'aaaaTGaAACGAATTCAaTAGTcCGGAgGATCTTGTCGGTGAcGATCACGCTcTCCTCGGGGAtGCCGGTGgCCACgACCTTGAAGaCGCCGATCAcgcg'
      #            ^          ^    ^    ^              ^         ^          ^       ^    ^         ^

      UMI,SNP2,SEQ2 = ex.extract_all(seq)
      self.assertEqual(UMI, 'aaaa')
      self.assertEqual(SNP2, 'c')
      self.assertEqual(SEQ2, 'TGaAACGAATTCAaTAGTcCGGAgGATCTTGTCGGTGAcGATCACGCTcTCCTCGGGGAtGCCGGTGgCCACgACCTTGAAGaCGCCGATCA')

      # Test case 3 (10 deletions).
      seq = 'aaaaTGAACGAATTCATAGTCGGAGATCTTGTCGGTGAGATCACGCTTCCTCGGGGAGCCGGTGCCACACCTTGAAGCGCCGATCAcgcg'

      UMI,SNP2,SEQ2 = ex.extract_all(seq)
      self.assertEqual(UMI, 'aaaa')
      self.assertEqual(SNP2, 'c')
      self.assertEqual(SEQ2, 'TGAACGAATTCATAGTCGGAGATCTTGTCGGTGAGATCACGCTTCCTCGGGGAGCCGGTGCCACACCTTGAAGCGCCGATCA')

      # Test case 4 (11 substitutions).
      seq = 'aaaaTGaAACGAATTCAaTAGTcCGGAgGATCTTtTCGGTGAcGATCACGCTcTCCTCGGGGAtGCCGGTGgCCACgACCTTGAAGaCGCCGATCAcgcg'
      #            ^          ^    ^    ^      ^       ^         ^          ^       ^    ^         ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 5 (11 deletions).
      seq = 'aaaaTGAACGAATTCATAGTCGGAGATCTGTCGGTGAGATCACGCTTCCTCGGGGAGCCGGTGCCACACCTTGAAGCGCCGATCAcgcg'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 6 (missing GCG).
      seq = 'aaaaTGAACGAATTCATAGTCGGAGATCTGTCGGTGAGATCACGCTTCCTCGGGGAGCCGGTGCCACACCTTGAAGCGCCGATCAcgtg'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)



if __name__ == '__main__':
   unittest.main(verbosity=2)
