#!/usr/bin/env python

import math
import sys

from collections import defaultdict
from itertools import combinations

# MM -> etype/time -> barcode -> [outcomes]
outcomes = defaultdict(lambda: \
           defaultdict(lambda: \
                 defaultdict(list)))

with open(sys.argv[1]) as f:
   for line in f:
      barcode, bias, MM, time, method, etype, expt = line.split()
      ex = "mock" if etype == "ctrl" else time
      outcomes[MM][ex][barcode].append(float(bias))

# Outcome pairs are 0,0, 0,1, 1,0, 1,1.
# MM -> etype/time -> [outcomes_pairs]
outcome_pairs = defaultdict(lambda: \
      defaultdict(lambda: [0,0,0,0]))

for MM in ("AG", "TG", "AC", "TC", "GT"):
   for ex in ("mock", "24", "48"):
      for list_of_outcomes in outcomes[MM][ex].values():
         if len(list_of_outcomes) < 2: continue
         for a,b in combinations(list_of_outcomes, 2):
            idx = int(round(a)) + 2*int(round(b))
            outcome_pairs[MM][ex][idx] += 1

for MM in ("AG", "TG", "AC", "TC", "GT"):
   for ex in ("mock", "24", "48"):
      total = float(sum(outcome_pairs[MM][ex]))
      _00, _01, _10, _11 = [x/total for x in outcome_pairs[MM][ex]]
      __0 = _00 + _10
      __1 = _01 + _11
      _0_ = _00 + _01
      _1_ = _10 + _11
      MI = _00 * math.log( _00 / __0 / _0_, 2) + \
           _01 * math.log( _01 / __1 / _0_, 2) + \
           _10 * math.log( _10 / __0 / _1_, 2) + \
           _11 * math.log( _11 / __1 / _1_, 2)
      sys.stdout.write("{}\t{}\t".format(MM, ex))
      sys.stdout.write("\t".join([str(x) for x in outcome_pairs[MM][ex]]))
      sys.stdout.write("\t{:.3f}\n".format(MI))
