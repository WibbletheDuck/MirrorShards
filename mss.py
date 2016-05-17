#!/usr/bin/env python

from optparse import OptionParser
import subprocess, sys, os

parser = OptionParser("Usage: %prog -r <run number> [options]")
parser.add_option("-r","--run",dest="run",type="int",default=0,
    help="Run number [required].")
parser.add_option("-c","--cell",dest="cell",type="string",default="j",
    help="Cell to use [default: %default].")
parser.add_option("-s","--shards",dest="shards",type="int",default=20,
    help="Number of shards to break distribution into [%default].")
parser.add_option("-n","--cores",dest="cores",type="int",default=24,
    help="Number of cores per shard [%default].")
parser.add_option("-w","--wall",dest="wall",type="int",default=50,
    help="Wall time per core [%default].")
parser.add_option("-m","--modulo",dest="modulo",type="int",default=0,
    help="Manually specify modulo point.  The first M shards will be given N cores per shard, the remainder will be given N-1.  Set to -1 to disable auto-detect.")

(opt, args) = parser.parse_args()

if opt.run==0:
    print("Must provide a run number with -r.")
    sys.exit()

# break up distribution
subprocess.call(["matlab","-nodisplay","-r",
	"try; mirror_shards_distribute({0},{1}); catch; end; quit".format(opt.run,opt.shards)])

def ms_size(i,t):
    # returns size of an mshard file
    fname = "mshard-r4-{0}of{1}-input.mat".format(i,t)
    return os.path.getsize(fname)

if opt.modulo==0:
    # try to figure out where Distribute put the modulo point
    size = ms_size(1,opt.shards)
    for i in range(2,opt.shards+1):
        if size != ms_size(i,opt.shards):
            opt.modulo = i-1
            break
        if i==opt.shards:
            opt.modulo = opt.shards

elif opt.modulo == -1:
    # auto-detect disabled
    opt.modulo = opt.shards

def mss_files(cores,wall,cell,run,mtag):
    # function to create the PBS scripts
    pbsfn = "mss_PBS-r{0}-m{1}.sh".format(run,mtag)
    pbsfile = open(pbsfn,"w")

    subprocess.call(["sed",
        "s/@@PPN@@/{0}/; s/@@WALL@@/{1}/; s/@@CELL@@/cell{2}/; s/@@RUN@@/{3}/;".format(cores,wall,cell,run),
	    "./mirror_shards.PBStemplate"],
	    stdout=pbsfile)
    pbsfile.close()

# create primary PBS and job submission scripts
mss_files(opt.cores,opt.wall,opt.cell,opt.run,0)

subfile = open("mss_submit-r{0}.sh".format(opt.run),"w")
subfile.write("#!/bin/bash\n\n")
subfile.write("for ((i=1 ; i<={0} ; i++)); do\n".format(opt.modulo))
subfile.write("\tqsub -t $i mss_PBS-r{0}-m{1}.sh\n".format(opt.run,0))
subfile.write("\tsleep 7\n")
subfile.write("done\n")

if opt.modulo != opt.shards:
    # if modulo, create second PBS script, add second part to submission script
    mss_files(opt.cores-1,opt.wall,opt.cell,opt.run,1)

    subfile.write("\nfor ((i={0} ; i<={1} ; i++)); do\n".format(opt.modulo+1,opt.shards))
    subfile.write("\tqsub -t $i mss_PBS-r{0}-m{1}.sh\n".format(opt.run,1))
    subfile.write("\tsleep 7\n")
    subfile.write("done\n")

subfile.close()

# clean up output from Distribute
subprocess.call("rm -rf Job1*",shell=True)
