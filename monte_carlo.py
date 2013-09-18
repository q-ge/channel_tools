#!/usr/bin/env python

import os.path
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 9:
    print >>sys.stderr, "Usage: %s <data root> <chip> <channel> " \
        "<countermeasure> <timeslice> <divide limit> <size limit> " \
        "<precision>" % \
        sys.argv[0]
    sys.exit(1)

rootdir= os.path.dirname(sys.argv[0])

channel_hist= os.path.join(rootdir, "channel_hist")
channel_matrix= os.path.join(rootdir, "channel_matrix")
capacity= os.path.join(rootdir, "capacity")

dataroot, chip, channel, countermeasure= sys.argv[1:5]
timeslice, divlimit, limit= map(int, sys.argv[5:8])

epsilon= float(sys.argv[8])

print "Monte-carlo simulation of hypothetical empty-channel bandwidth."
print "Chip: %s, Channel: %s, Countermeasure: %s, Timeslice: %d" % \
    (chip, channel, countermeasure, timeslice)
print "Experimental data in %s" % os.path.abspath(dataroot)

run_name= "%s.%s.%s.%d" % (chip, channel, countermeasure, timeslice)

divisors= []
d= 1
for i in xrange(divlimit):
    divisors.append(d)
    d*= 2
sizes= [limit / d for d in divisors]

print "Sampling at sizes [" + ", ".join(map(str,sizes)) + "]"

channel_dir= os.path.join(dataroot, chip, channel)
runs_dir= os.path.join(channel_dir, countermeasure, "TS_%d" % timeslice)

range_str= None
info= file(os.path.join(channel_dir, "info"), "r")
for l in info:
    m= re.match("^modulation range:\s*(\d+)\s*-\s*(\d+)", l)
    if m:
        range_str= m.group(1, 2)
del info
mod_range= (int(range_str[0]), int(range_str[1]))

print "Decompressing and collating samples...",
sys.stdout.flush()
samples_file= tempfile.mkstemp()
subprocess.call("find %s -name \"*.xz\" | xargs xzcat > %s" % \
    (runs_dir, samples_file[1]), shell=True)
print " done."

hist_file= tempfile.mkstemp()
subprocess.call("%s %s %d %d < %s" \
    % (channel_hist, hist_file[1], mod_range[0], mod_range[1], \
       samples_file[1]), shell=True)

print "Calculating bandwidth for partitioned datasets"

partitioned_output= file(run_name + ".part", "w")

for s in sizes:
    print "Building subsampled channel matrices of size %d" % s
    discards= []
    S= 0
    n= 0
    while S < limit:
        discards.append(S)
        S+= s
        n+= 1
    subsampled_channels= [tempfile.mkstemp() for s in sizes]
    pipes= [subprocess.Popen("%s %s %d %d %d %d < %s > /dev/null" % \
        (channel_matrix, subsampled_channels[i][1], \
         mod_range[0], mod_range[1], s, discards[i], samples_file[1]), \
        shell=True, stdout=subprocess.PIPE) for i in xrange(n)]

    for p in pipes:
        sodata, sedata= p.communicate()

    # Calculate capacities
    pipes= [subprocess.Popen("%s %s %f -q" % \
        (capacity, subsampled_channels[i][1], epsilon), \
        shell=True, stdout=subprocess.PIPE) for i in xrange(n)]

    for p in pipes:
        sodata, sedata= p.communicate()
        cap= float(sodata.strip())
        partitioned_output.write("%d %.12f\n" % (s, cap))
         
    for ssc in subsampled_channels:
        os.unlink(ssc[1])
    print " done."

del(partitioned_output)

os.unlink(hist_file[1])
os.unlink(samples_file[1])
