import os, time, math, tempfile, fileinput
from os.path import join
from glob import glob
import pandas as pd

from MasterScript_v3 import RunFromQueue, HomogenousCrystal


def DefectiveImages(desc_file):
	description = pd.read_csv( desc_file,header=None,sep=': ',engine='python', index_col=0, names=['value'] )
	mlfsom_path = os.path.expanduser('~/Desktop/MLFSOM')
	dir_name = os.path.dirname(desc_file)
	cwd = os.getcwd()
	os.chdir(dir_name)
	img_list = [x for x in os.listdir('.') if x.endswith('.img')]
	frame_list = range(int(description.loc['frames','value']))

	for fnumber in frame_list[:]:
		img = '_'.join([ description.loc['prefix','value'], str(fnumber), '001.img'])
		if os.path.isfile(img) and os.path.getsize(img) > 10**3:
			frame_list.remove(fnumber)
	os.chdir(cwd)
	return frame_list


def UpdateQueueFile(desc_file):
	defectives = DefectiveImages(desc_file)
	dir_name = os.path.dirname(desc_file)
	cwd = os.getcwd()
	os.chdir(dir_name)
	queue_file = glob('*queue.txt')[0]
	for line in fileinput.FileInput(queue_file,inplace=True):
		fnumber = eval(line)[0]
		if fnumber in defectives:
			print line.strip()
	os.chdir(cwd)


# **********************************************************************
pdb_file= "2W3B.pdb"
prefix= "stills_2W3B_2.0A_10fr_0.1deg_350mm"
stills= True
frames= 10
start_mos= 0.05
k_mos= 0.1
k_cell= 0.0004
k_bfactor= 6
resolution= 2.0
solvent_B= 35
osc= 0.1
distance= 350
flux= 1.5e+11
threads= 10

time_initial = time.time()
HomogenousCrystal(pdb_file=pdb_file, prefix=prefix, stills=stills, start_mos=start_mos, k_mos=k_mos,\
	k_cell=k_cell, k_bfactor=k_bfactor, frames=frames, resolution=resolution, solvent_B=solvent_B, osc=osc,\
	 distance=distance, flux=flux, threads=threads)

desc_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-desc.txt")
queue_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-queue.txt")
defectives = DefectiveImages(desc_file)

print '\nFINISHED BATCH NUMBER: 1'
print "\nLIST OF DEFECTIVE IMAGES: %s" %str(defectives)

batch_number = 2
ideal_threads = 10
while defectives != []:
	N_defective = len(DefectiveImages(desc_file))
	if N_defective <= ideal_threads:
		threads = N_defective
	else:
		threads = ideal_threads
	UpdateQueueFile(desc_file)
	RunFromQueue(desc_file,queue_file,threads=threads)
	defectives = DefectiveImages(desc_file)
	print '\nFINISHED BATCH NUMBER: %i' %batch_number
	batch_number += 1

time_final = time.time()
print 'GRAND TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)

#******************************************************************************
pdb_file= "2W3B.pdb"
prefix= "stills_2W3B_2.0A_10fr_0.3deg_350mm"
stills= True
frames= 10
start_mos= 0.05
k_mos= 0.1
k_cell= 0.0004
k_bfactor= 6
resolution= 2.0
solvent_B= 35
osc= 0.3
distance= 350
flux= 1.5e+11
threads= 10

time_initial = time.time()
HomogenousCrystal(pdb_file=pdb_file, prefix=prefix, stills=stills, start_mos=start_mos, k_mos=k_mos,\
	k_cell=k_cell, k_bfactor=k_bfactor, frames=frames, resolution=resolution, solvent_B=solvent_B, osc=osc,\
	 distance=distance, flux=flux, threads=threads)

desc_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-desc.txt")
queue_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-queue.txt")
defectives = DefectiveImages(desc_file)

print '\nFINISHED BATCH NUMBER: 1'
print "\nLIST OF DEFECTIVE IMAGES: %s" %str(defectives)

batch_number = 2
ideal_threads = 10
while defectives != []:
	N_defective = len(DefectiveImages(desc_file))
	if N_defective <= ideal_threads:
		threads = N_defective
	else:
		threads = ideal_threads
	UpdateQueueFile(desc_file)
	RunFromQueue(desc_file,queue_file,threads=threads)
	defectives = DefectiveImages(desc_file)
	print '\nFINISHED BATCH NUMBER: %i' %batch_number
	batch_number += 1

time_final = time.time()
print 'GRAND TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)


# *******************************************************************
pdb_file= "2W3B.pdb"
prefix= "stills_2W3B_2.0A_10fr_0.5deg_350mm"
stills= True
frames= 10
start_mos= 0.05
k_mos= 0.1
k_cell= 0.0004
k_bfactor= 6
resolution= 2.0
solvent_B= 35
osc= 0.5
distance= 350
flux= 1.5e+11
threads= 10

time_initial = time.time()
HomogenousCrystal(pdb_file=pdb_file, prefix=prefix, stills=stills, start_mos=start_mos, k_mos=k_mos,\
	k_cell=k_cell, k_bfactor=k_bfactor, frames=frames, resolution=resolution, solvent_B=solvent_B, osc=osc,\
	 distance=distance, flux=flux, threads=threads)

desc_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-desc.txt")
queue_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-queue.txt")
defectives = DefectiveImages(desc_file)

print '\nFINISHED BATCH NUMBER: 1'
print "\nLIST OF DEFECTIVE IMAGES: %s" %str(defectives)

batch_number = 2
ideal_threads = 10
while defectives != []:
	N_defective = len(DefectiveImages(desc_file))
	if N_defective <= ideal_threads:
		threads = N_defective
	else:
		threads = ideal_threads
	UpdateQueueFile(desc_file)
	RunFromQueue(desc_file,queue_file,threads=threads)
	defectives = DefectiveImages(desc_file)
	print '\nFINISHED BATCH NUMBER: %i' %batch_number
	batch_number += 1

time_final = time.time()
print 'GRAND TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)

# *****************************************************************
pdb_file= "2W3B.pdb"
prefix= "stills_2W3B_2.0A_10fr_1.0deg_350mm"
stills= True
frames= 10
start_mos= 0.05
k_mos= 0.1
k_cell= 0.0004
k_bfactor= 6
resolution= 2.0
solvent_B= 35
osc= 1.0
distance= 350
flux= 1.5e+11
threads= 10

time_initial = time.time()
HomogenousCrystal(pdb_file=pdb_file, prefix=prefix, stills=stills, start_mos=start_mos, k_mos=k_mos,\
	k_cell=k_cell, k_bfactor=k_bfactor, frames=frames, resolution=resolution, solvent_B=solvent_B, osc=osc,\
	 distance=distance, flux=flux, threads=threads)

desc_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-desc.txt")
queue_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-queue.txt")
defectives = DefectiveImages(desc_file)

print '\nFINISHED BATCH NUMBER: 1'
print "\nLIST OF DEFECTIVE IMAGES: %s" %str(defectives)

batch_number = 2
ideal_threads = 10
while defectives != []:
	N_defective = len(DefectiveImages(desc_file))
	if N_defective <= ideal_threads:
		threads = N_defective
	else:
		threads = ideal_threads
	UpdateQueueFile(desc_file)
	RunFromQueue(desc_file,queue_file,threads=threads)
	defectives = DefectiveImages(desc_file)
	print '\nFINISHED BATCH NUMBER: %i' %batch_number
	batch_number += 1

time_final = time.time()
print 'GRAND TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)

# *************************************************************************
pdb_file= "2W3B.pdb"
prefix= "stills_2W3B_2.0A_10fr_2.0deg_350mm"
stills= True
frames= 10
start_mos= 0.05
k_mos= 0.1
k_cell= 0.0004
k_bfactor= 6
resolution= 2.0
solvent_B= 35
osc= 2.0
distance= 350
flux= 1.5e+11
threads= 10

time_initial = time.time()
HomogenousCrystal(pdb_file=pdb_file, prefix=prefix, stills=stills, start_mos=start_mos, k_mos=k_mos,\
	k_cell=k_cell, k_bfactor=k_bfactor, frames=frames, resolution=resolution, solvent_B=solvent_B, osc=osc,\
	 distance=distance, flux=flux, threads=threads)

desc_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-desc.txt")
queue_file = join( os.path.expanduser('~/Desktop/MLFSOM'), 'data_'+prefix, "exp-"+prefix+"-queue.txt")
defectives = DefectiveImages(desc_file)

print '\nFINISHED BATCH NUMBER: 1'
print "\nLIST OF DEFECTIVE IMAGES: %s" %str(defectives)

batch_number = 2
ideal_threads = 10
while defectives != []:
	N_defective = len(DefectiveImages(desc_file))
	if N_defective <= ideal_threads:
		threads = N_defective
	else:
		threads = ideal_threads
	UpdateQueueFile(desc_file)
	RunFromQueue(desc_file,queue_file,threads=threads)
	defectives = DefectiveImages(desc_file)
	print '\nFINISHED BATCH NUMBER: %i' %batch_number
	batch_number += 1

time_final = time.time()
print 'GRAND TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)