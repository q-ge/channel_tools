#!/usr/bin/env python

from scipy.sparse import coo_matrix, dia_matrix, dok_matrix
from numpy import matrix, ones, log2, exp2, asarray, asscalar, asmatrix
import numpy
import gzip
import sys

from blahut_arimoto import find_bw

logfile= sys.stdin
#logfile= gzip.open(sys.argv[1])
e=       float(sys.argv[1])
#e=       float(sys.argv[2])

# Build a histogram

print >>sys.stderr, "Loading results"

sM=  set()
sR=  set()
sRM= set()

counts= {}

samples= {}

for line in logfile:
    if line.strip()[0] == '#':
        #field, value= line.strip().split('\t')
        #if field == "# Samples per modulation": samples= int(value)
        continue

    M,R= map(int,line.strip().split())

    sM.add(M)
    sR.add(R)
    sRM.add((R,M))

    if not (R,M) in counts: counts[R,M]= 0
    counts[R,M]+= 1

    if not M in samples: samples[M]= 0
    samples[M]+= 1

#print >>sys.stderr, "%d samples per modulation" % samples
print >>sys.stderr, "Modulation range: %d-%d" % (min(sM),max(sM))
print >>sys.stderr, "Result range:     %d-%d" % (min(sR),max(sR))

# Construct the channel matrix

print >>sys.stderr, "Building channel matrix"

nM= len(sM)
nR= len(sR)

lM= list(sM)
lR= list(sR)

# Construct dense indices
iM= {}
for m in xrange(nM):
    iM[lM[m]]= m

iR= {}
for r in xrange(nR):
    iR[lR[r]]= r

# Compress and reindex the histogram.  Permuting symbols has
# no effect on channel bandwidth.
lRM= list(sRM)
Ms= [iM[m] for r,m in lRM]
Rs= [iR[r] for r,m in lRM]
# Scale to probabilities
Ps= [float(counts[r,m]) / samples[m] for r,m in lRM]
#Cs= [counts[r,m] for r,m in lRM]
del sM, sR, sRM, lM, lR, iM, iR, lRM

mCount= coo_matrix((Ps,(Rs,Ms)))
del Ms, Rs, Ps

# Scale to probabilities
Q= mCount.tocsr()# / float(samples)
del mCount

print >>sys.stderr, "Channel matrix size " + str(Q.shape) + \
                    ", %d entries" % Q.getnnz() + \
                    ", density %.6f%%" % \
                        (float(Q.getnnz())/(Q.shape[0]*Q.shape[1]))

# Calculate channel bandwidth using the Blahut-Arimoto algorithm

print >>sys.stderr, \
    "Computing channel bandwith, target error (Iu-Il) is %gb/s" % e

Il, Iu= find_bw(Q, e, verbose=True)

print (Il + Iu) / 2
