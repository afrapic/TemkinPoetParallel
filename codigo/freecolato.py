#!/usr/bin/python

import string, os
import re
from math import sqrt
from operator import itemgetter


space = re.compile(r'\s+')
cmd = 'qstat -f | grep -v "\-\-" | grep -v "used" '  
#cmd = 'cat queue.txt'
num = 0
mpinodes ={}
for file in os.popen(cmd).readlines():     # run command
	name = file[:-1]                       # strip '\n'
	na = str(name)
        if na.find("local")>0:
		name,qt,nodes,load,arch,stat = space.split(name)
	else:
		break

	#node = name[0]
	#load = name[2]
	#   Select working nodes
	if load.rfind("NA")!=1:
		# get free nodes
		used = nodes.split("/",1)[0]
		avail = 4-int(used)
		if(avail>0) :
			name = name.split(".",1)[0]
			name = name.split("@")[1]
			#print name+":",avail
			mpinodes[name]=avail
			num = num+avail			
#print "Available nodes:",num
mpib = int(sqrt(num))		
#print "Best mpi size:",mpib*mpib
sortedmpi = sorted(mpinodes.items(),key=itemgetter(1), reverse=True)
#print sortedmpi
i=0
#f = open("mpd.hosts", "w")
for n,s in sortedmpi:
        st = n+":"+str(s)
        #f.write(st)
	print st
	i = i+s
	if i>mpib*mpib:
		break
print i
