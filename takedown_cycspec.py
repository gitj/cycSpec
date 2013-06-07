#!/bin/env python
import os
print "Killing all processes"
cmd = 'gpu_run "pkill -f gjones"'
print cmd
os.system(cmd)

print "Stopping pyro nameserver"
cmd = 'pkill -f Pyro'
print cmd
os.system(cmd)