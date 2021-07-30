#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest

from StringIO import StringIO
from textwrap import dedent

import count


class TestScarcodeRemover(unittest.TestCase):

   def test_init_fail(self):
      with self.assertRaises(KeyError):
         count.ScarcodeRemover('AA')

   def test_remove(self):

      # List of scarcodes:
      # CT: CGCTAATTAATG
      # CA: GCTAGCAGTCAG
      # GA: GCTAGCTCGTTG
      # GT: GCTAGCTCCGCA

      scar = count.ScarcodeRemover('CT')
      self.assertEqual(scar.remove('aaaCGCTAATTAATG'), 'aaa')
      self.assertEqual(scar.remove('aaaCGCTAtTTAATG'), 'aaa')
      self.assertEqual(scar.remove('aaagGCTAtTTAATG'), 'aaag')
      self.assertEqual(scar.remove('aaaCaCTAtTTAATG'), 'aaa')

      scar = count.ScarcodeRemover('CA')
      self.assertEqual(scar.remove('aaaGCTAGCAGTCAG'), 'aaa')
      self.assertEqual(scar.remove('aaaGCTAGCAaTCAG'), 'aaa')
      self.assertEqual(scar.remove('aaatCTAGCAGTCAG'), 'aaat')
      self.assertEqual(scar.remove('aaaGtTAGCAGTCAG'), 'aaa')

      scar = count.ScarcodeRemover('GA')
      self.assertEqual(scar.remove('aaaGCTAGCTCGTTG'), 'aaa')
      self.assertEqual(scar.remove('aaaGCTAGCTCaTTG'), 'aaa')
      self.assertEqual(scar.remove('aaaaCTAGCTCGTTG'), 'aaaa')

      scar = count.ScarcodeRemover('GT')
      self.assertEqual(scar.remove('aaaGCTAGCTgCGCA'), 'aaa')
      self.assertEqual(scar.remove('aaaGCTAGCTCCGCA'), 'aaa')
      self.assertEqual(scar.remove('aaatCTAGCTCCGCA'), 'aaat')
      self.assertEqual(scar.remove('aaaGaTAGCTCCGCA'), 'aaa')


   def test_remove_fail(self):
      # Test if the remove raises an exception if the scarcode
      # is not found.

      scar = count.ScarcodeRemover('CT')
      with self.assertRaises(count.WrongScarcodeException):
         scar.remove('aaaAAAAAAAAAAAA')
      with self.assertRaises(count.WrongScarcodeException):
         scar.remove('aaaCGCaAgTTtATG')
      with self.assertRaises(count.WrongScarcodeException):
         scar.remove('aaaTAATTAATG')
      with self.assertRaises(count.WrongScarcodeException):
         scar.remove('aaatatTAATTAATG')
      with self.assertRaises(count.WrongScarcodeException):
         scar.remove('aaaGCTAGCAGTCAG')



class TestTagNormalizer(unittest.TestCase):

   def setUp(self):
      # Mini starcode file.
      f = StringIO(
         'AGATGCTACGCG\t3\tACATGCTACGGC,ACATGCTACGCC\n' \
         'CCATGCTACGAA\t9\tCCATGCTACGAC,CCATGCTACGCA'
      )

      self.normalizer = count.TagNormalizer(f)

   def test_init(self):
      # Make sure that the internal dictionary has been
      # updated upon construction.
      self.assertNotEqual(self.normalizer.canonical, dict())

      # Make sure that an exception is raised when the input
      # file is not properly formatted.
      f = StringIO('AGATGCTACGCG')

      with self.assertRaises(ValueError):
         self.normalizer = count.TagNormalizer(f)


   def test_normalize(self):
      # Make sure that the normalizer can normalize tags.
      bcd,umi = self.normalizer.normalize('ACATGCTACGGC')
      self.assertEqual((bcd,umi), ('AG', 'CG'))

      bcd,umi = self.normalizer.normalize('ACATGCTACGCC')
      self.assertEqual((bcd,umi), ('AG', 'CG'))

      bcd,umi = self.normalizer.normalize('CCATGCTACGAC')
      self.assertEqual((bcd,umi), ('CC', 'AA'))

      bcd,umi = self.normalizer.normalize('CCATGCTACGCA')
      self.assertEqual((bcd,umi), ('CC', 'AA'))


   def test_iter(self):

      # Make sure that the normalizer can be called directly
      # in a for statement, and use list comprehension for the test.
      tags = sorted([tag for tag in self.normalizer])
      self.assertEqual(tags, ['AGATGCTACGCG', 'CCATGCTACGAA'])



class TestCountingInfo(unittest.TestCase):

   def test_init(self):
      info = count.CountingInfo('dummy_fname1', 'dummy_fname2')
      self.assertEqual(info.fname1, 'dummy_fname1')
      self.assertEqual(info.fname2, 'dummy_fname2')
      self.assertIsNone(info.MMcode)
      self.assertEqual(info.nreads, 0)
      self.assertEqual(info.used_reads, 0)
      self.assertEqual(info.vart_conflicts, {})
      self.assertEqual(info.wrong_scarcode, 0)
      self.assertEqual(info.barcode_too_short, 0)
      self.assertEqual(info.too_few_reads, 0)
      self.assertEqual(info.non_unique_UMI, 0)
      self.assertEqual(info.minority_report, 0)


   def test_get_MMcode(self):

      # Need a CountingInfo instance.
      info = count.CountingInfo('dummy_fname1', 'dummy_fname2')

      # Test case 1 (GA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCGTTGATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCGTTGATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCGTTGATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'GA')

      # Test case 2 (CT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGCGCTAATTAATGATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGCGCTAATTAATGATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCCGCTAATTAATGAATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'CT')

      # Test case 3 (CA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCAGTCAGATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCAGTCAGATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCAGTCAGATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'CA')

      # Test case 4 (GT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCCGCAATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCCGCAATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCCGCAATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'GT')


   def test_normalize_variant(self):
      # Instantiate info.
      info = count.CountingInfo('dummy_fname1', 'dummy_fname2')

      dict_of_variants = { ('A','T'): 14, ('G','C'): 1, }
      norm = info.normalize_variant('ACTGGACGC', 'AAA', dict_of_variants)

      self.assertEqual(norm, ('A','T'))
      self.assertEqual(info.minority_report, 1)
      self.assertEqual(info.used_reads, 14)
      self.assertEqual(info.vart_conflicts,
         {('ACTGGACGC', 'AAA'): {('G', 'C'): 1, ('A', 'T'): 14}})


   def test_write_to_file(self):
      info = count.CountingInfo('dummy_fname1', 'dummy_fname2')
      info.MMcode = 'GA'
      info.nreads = 100
      info.used_reads = 80
      info.wrong_scarcode = 2
      info.barcode_too_short = 1
      info.too_few_reads = 3
      info.non_unique_UMI = 2
      info.minority_report = 1

      buffer = StringIO()
      info.write_to_file(buffer)

      txt = '''dummy_fname1
         dummy_fname2
         MM type: GA
         Used reads:\t80 (80.00%)
         Reads lost to:
           Wrong scarcode:\t2 (2.00%)
           Barcode too short:\t1 (1.00%)
           Too few reads:\t3 (3.00%)
           Non unique UMI:\t2 (2.00%)
           Minority report:\t1 (1.00%)
         ---'''

      check = '\n'.join(buffer.getvalue().splitlines()[1:])
      self.assertEqual(check, dedent(' '*9 + txt))


class TestEventCounter(unittest.TestCase):

   def test_init(self):

      # Mini starcode file (GA).
      f = StringIO(
         'GCTAGCTCGTTGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCGTTGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCGTTGATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('dummy_fname1', 'dummy_fname2')
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'GA') 
      
      # Mini starcode file (CT).
      f = StringIO(
         'CGCTAATTAATGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'CGCTAATTAATGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'CGCTAATTAATGATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'CT') 

      # Mini starcode file (CA).
      f = StringIO(
         'GCTAGCAGTCAGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCAGTCAGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCAGTCAGATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'CA') 

      # Mini starcode file (GT).
      f = StringIO(
         'GCTAGCTCCGCAATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCCGCAATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCCGCAATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'GT') 


   def test_count(self):
      
      # Mini starcode file (GA).
      f = StringIO(
         'AAGATGCTAGCTCGTTGATGCTACGTAC\t1\tAAGATGCTAGCTCGTTGATGCTACGTAC\n' \
         'AGGATGCTAGCTCGTTGATGCTACGGGG\t1\tAGGATGCTAGCTCGTTGATGCTACGGGG\n' \
         'TTGATGCTAGCTCGTTGATGCTACGAAA\t1\tTTGATGCTAGCTCGTTGATGCTACGAAA\n' \
         'AGGTTTGGTGGTGGAGGATGCTACGTAA\t1\tAGGTTTGGTGGTGGAGGATGCTACGTAA\n' \
         'GATGCTAGCTCGTTGATGCTACGTAC\t1\tGATGCTAGCTCGTTGATGCTACGTAC\n'     \
         'AAGATGCTAGCTCGTTGATGCTACGGGG\t1\tAAGATGCTAGCTCGTTGATGCTACGGGG\n' \
         'TTGATGCTAGCTCGTTGATGCTACGGGG\t1\tTTGATGCTAGCTCGTTGATGCTACGGGG\n'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('dummy_fname1', 'dummy_fname2')
      counter = count.EventCounter(normalizer, info)

      # Mini pps file (GA).
      # barcode GCTAGCTCGTTG ATGCTACG UMI...
      f = StringIO(
         # Read  1 has a wrong scarcode.
         # Reads 2-4 are identical and normal.
         # Reads 5-6 is in conflict with reads 2-4.
         # Reads 7-8 are identical and normal.
         # Read  9 is normal but discarded (too_few_reads).
         # Read  10 is not registered in the starcode file.
         # Read  11 has a too short barcode.
         # Reads 12-15 create a UMI conflict.
         'AGGTTTGGTGGTGGAGGATGCTACGTAA\tA\tC\tGAG\tATA\n' \
         'AAGATGCTAGCTCGTTGATGCTACGTAC\tA\tC\tGAG\tGTA\n' \
         'AAGATGCTAGCTCGTTGATGCTACGTAC\tA\tC\tGAG\tAAA\n' \
         'AAGATGCTAGCTCGTTGATGCTACGTAC\tA\tC\tGAG\tATG\n' \
         'AAGATGCTAGCTCGTTGATGCTACGTAC\tA\tT\tGAG\tATA\n' \
         'AAGATGCTAGCTCGTTGATGCTACGTAC\tA\tT\tGAG\tATA\n' \
         'TTGATGCTAGCTCGTTGATGCTACGAAA\tA\tC\tGAG\tATA\n' \
         'TTGATGCTAGCTCGTTGATGCTACGAAA\tA\tC\tGAG\tATA\n' \
         'AGGATGCTAGCTCGTTGATGCTACGGGG\tA\tC\tGAG\tATA\n' \
         'AACTGAGAGGCGATGAGCGAGTAATAGC\tA\tC\tGAG\tATA\n' \
         'GATGCTAGCTCGTTGATGCTACGTAC\tA\tC\tGAG\tATA\n'   \
         'AAGATGCTAGCTCGTTGATGCTACGGGG\tA\tC\tGAG\tATA\n' \
         'AAGATGCTAGCTCGTTGATGCTACGGGG\tA\tC\tGAG\tATA\n' \
         'TTGATGCTAGCTCGTTGATGCTACGGGG\tA\tC\tGAG\tATA\n' \
         'TTGATGCTAGCTCGTTGATGCTACGGGG\tA\tC\tGAG\tATA\n'
      )

      out = StringIO()
      counter.count(f, out)

      self.assertEqual(out.getvalue(),
         'barcode\tFF\tAT\tGC\n'
         'TTGAT\t1\t0\t0\n'
         'AAGAT\t1\t0\t0\n'
      )

      # Test info gathering.
      self.assertEqual(info.MMcode, 'GA')
      self.assertEqual(info.nreads, 15)
      self.assertEqual(info.wrong_scarcode, 1)
      self.assertEqual(info.barcode_too_short, 1)
      self.assertEqual(info.too_few_reads, 1)
      self.assertEqual(info.non_unique_UMI, 4)
      self.assertEqual(info.minority_report, 2)
      self.assertEqual(len(info.vart_conflicts), 1)


if __name__ == '__main__':
   #suite = unittest.TestSuite()
   #suite.addTest(TestCountingInfo('test_write_to_file'))
   #unittest.TextTestRunner(verbosity=2).run(suite)
   unittest.main(verbosity=2)
