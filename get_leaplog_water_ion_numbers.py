#!/usr/bin/env python
import re
import sys

leap_file = sys.argv[1]

with open(leap_file, "r") as f:
    text = f.read()

pat = re.compile("\s+WAT\s+(\d+)")
matches = re.findall(pat, text)
water_num = matches[0]

pat = re.compile("(\d+) Na\+ ion[s]{0,1} required to neutralize.")
matches = re.findall(pat, text)
Na_num = matches[0] if len(matches) > 0 else 0

pat = re.compile("(\d+) Cl- ion[s]{0,1} required to neutralize.")
matches = re.findall(pat, text)
Cl_num = matches[0] if len(matches) > 0 else "0"

print(water_num, Na_num, Cl_num)
