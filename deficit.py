#!/usr/bin/env python

import sys

target= int(sys.argv[1])

for l in sys.stdin:
    bits= l.strip().split(' ')
    if len(bits) == 2:
        c, count= map(int, bits)
        if count < target:
            print c, target - count
