#!/bin/env python
import os
import time
print "Restarting all osf processes"

print "Killing all processes"
cmd = 'gpu_run "pkill -f gjones"'
print cmd
os.system(cmd)

#print "Checking data buffers"
#cmd = 'gpu_run "source /home/gpu/gjones/puppi.sh; /home/gpu/gjones/guppi_daq/bin/check_guppi_databuf -c -i1 -n24 -s32"'
#print cmd
#os.system(cmd)

print "Checking status memory"
cmd = 'gpu_run "source /home/gpu/gjones/puppi.sh; /home/gpu/gjones/guppi_daq/bin/check_guppi_status"'
print cmd
os.system(cmd)

print "Stopping pyro nameserver"
cmd = 'pkill -f Pyro'
print cmd
os.system(cmd)

print "Starting Pyro nameserver"
cmd = 'source /home/gpu/gjones/puppi.sh; nohup python -m Pyro4.naming -n master < /dev/null &> /dev/null &'
print cmd
os.system(cmd)
time.sleep(2)
print "Starting osf servers"

cmd = 'gpu_run "source /home/gpu/gjones/puppi.sh; nohup python /home/gpu/gjones/osfServer.py < /dev/null &> /dev/null &"'
print cmd
os.system(cmd)