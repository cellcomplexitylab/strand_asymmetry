#!/bin/bash

python scripts/make_blacklist_from_mapping.py mapping/*.ins.gz
python scripts/make_blacklist_from_mismatches.py mismatches/*_no_ISceI_*.co.gz
