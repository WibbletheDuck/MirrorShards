#!/usr/bin/env python

from optparse import OptionParser
import subprocess, sys, os

parser = OptionParser("Usage: %prog -r <run number> [options]")
parser.add_option("-r","--run",dest="run",type="int",default=0,
    help="Run number [required].")

(opt, args) = parser.parse_args()

if opt.run==0:
    print("Must provide a run number with -r.")
    sys.exit()

# gather distribution
subprocess.call(["matlab","-nodisplay","-r",
    "try; mirror_shards_gather('mshards-r{0}-master.mat'); catch; end; quit".format(opt.run)])

