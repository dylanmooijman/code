#Single molecule FISH probe designer
#Dylan Mooijman 2013

import sys
from string import maketrans

file_name=str(sys.argv[1])
probe_len=int(sys.argv[2])
probe_len_float=float(format(probe_len, '.2f'))

#declare variables
sense='atcg'
a_sense='tagc'
rc=maketrans(sense, a_sense)

spacer=probe_len
oldcounter=-probe_len

good_probes=list()
final_probes=list()


#print 'probe length', probe_len

#def functions
#gc content
def gc(string):
	gc_string=(string.count('c')+string.count('g'))/probe_len_float
	return gc_string

#read mRNA
gene=open(file_name, 'r')
lines = gene.readlines()

#delete fasta header
del lines[0]

#remove newline junk
string_w_n=' '.join(lines)
clean_string=string_w_n.replace('\n', '')
clean_string=clean_string.replace(' ', '')
#reverse complement
rc_string=clean_string[::-1].translate(rc)

#get all robes sequences
probe_seq=[rc_string[x:x+probe_len] for x in range(0, len(rc_string), 1)]

#get mean gc content
import numpy as np
gc_all=[gc(x) for x in probe_seq]
gc_median=np.median(gc_all)
#print gc_median

gc_mean=(float(sum(gc_all)))/(float(len(gc_all)))
#print gc_mean

gc_5=12.1/100
gc_lower=gc_mean-gc_5
gc_upper=gc_mean+gc_5

#test gc_content
import matplotlib.pyplot as plt
plt.plot(gc_all, 'b')
plt.axhline(gc_lower, color='red')
plt.axhline(gc_upper, color='red')
plt.show()


#for probes in gc_all:
#	print probes

#probe location list
probes=[x for x in range(0, len(rc_string), 1) if gc(rc_string[x:x+probe_len])>gc_lower and gc(rc_string[x:x+probe_len])<gc_upper]
#print len(probes)

#probe spacing
probe_space_loc=list()
for probe_loc in probes:
#	print 'prev loc', oldcounter
#	print 'probe_loc', probe_loc
#	print 'prev-current', probe_loc-oldcounter
	if probe_loc-oldcounter>15:
		good_probes.append(probe_loc)
		oldcounter=probe_loc+15
#		print 'move to', oldcounter
		probe_space_loc.append(oldcounter)
#print good_probes
#map probe location back to probes
for x in good_probes:
	final_probes.append(probe_seq[x])
	#print x
	#print probe_seq[x]	
#for x in probe_space_loc:
#	print x

print 'you have found' ,len(final_probes),'probes'

#open file for writing
newfile_name=file_name.replace('.fa', '')
probe_write=open(newfile_name+'.probes', 'w')
write_counter=''
counter=0
for probe in final_probes:
	write_counter='>'+str(counter)
	probe_write.write("%s\n" % write_counter)
	probe_write.write("%s\n" % probe)
	counter+=1
