#!/usr/bin/env python

'''Take a smi file with NSC compounds and call nscavail on each one.
For some reason every so often the connection hangs and mechanize's timeout
doesn't fire, so we are manually calling the script with a timeout'''

import sys,subprocess32,os,signal

script = "nscavail.py"
if len(sys.argv) > 2:
    script = sys.argv[2]
for line in open(sys.argv[1]):
    (smi,name) = line.split()
    for i in range(3):
      try: # 3 attempts
        p = subprocess32.Popen("%s %s" % (script,name), stdout=subprocess32.PIPE, shell=True,preexec_fn=os.setsid)    
        p.wait(timeout=60)
        if p.stdout.read().startswith("True"):
            print(line.rstrip())
        break
      except subprocess32.TimeoutExpired:
        os.killpg(os.getpgid(p.pid), signal.SIGTERM) 
        sys.stderr.write('Problem with %s\n' % name)

