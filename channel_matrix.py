#!/usr/bin/env python

import gzip

import sys

logfile= sys.stdin
Mstep= int(sys.argv[1])
Rbin= int(sys.argv[2])

# Build a histogram

print >>sys.stderr, "Loading results"

counts= {}
samples= None
Mmin= None
Mmax= None
Bmin= None
Bmax= None
#mapM= {}
#averages= {}
counts= {}

for line in logfile:
    if line.strip()[0] == '#':
        #field, value= line.strip().split('\t')
        #if field == "# Samples per modulation": samples= int(value)
        continue

    M,R= map(int,line.strip().split())

    B= R / Rbin

    if Mmin is None or M < Mmin: Mmin= M
    if Mmax is None or M > Mmax: Mmax= M
    if Bmin is None or B < Bmin: Bmin= B
    if Bmax is None or B > Bmax: Bmax= B

    if not (B,M) in counts: counts[B,M]= 0
    counts[B,M]+= 1

    #if not M in averages: averages[M]= 0.0
    #averages[M]+= (B*Rbin + Rbin/2) / samples

    if not M in counts: counts[M]= 0
    counts[M]+= 1

    #if not B in mapM:
        #mapM[B]= M
    #if counts[B,M] > counts[B,mapM[B]]:
        #mapM[B]= M

for M in xrange(Mmin,Mmax+1,Mstep):
    for B in xrange(Bmin,Bmax+1):
        if (B,M) in counts:
            print M, B*Rbin + Rbin/2, float(counts[B,M]) / counts[M]
        else:
            print M, B*Rbin + Rbin/2, 0.0

#print
#print
#for M in xrange(Mmin,Mmax+1,Mstep):
    #print M, averages[M]

#print
#print
#for B in mapM:
    #print B*Rbin + Rbin/2, mapM[B]
