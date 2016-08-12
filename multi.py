import subprocess

nSessions = 2

for i in xrange(nSessions):
  subprocess.Popen(["./spatiocyte-core"])
