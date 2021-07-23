#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

def main(f):
   for line in f:
      items = line.split()
      barcode = items[0]
      print barcode

if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      main(f)
