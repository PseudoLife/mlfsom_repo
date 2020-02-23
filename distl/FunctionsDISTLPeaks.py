import os, pickle, matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

import img2cbf
reload(img2cbf)  # reload module for new changes(if any) to take effect
from img2cbf import img2cbf


def distl_wrapper(fake_cbf,distl_params):
	'''
	Takes a cbf_file and distl_params txt file
	Returns distl bash command to run
	e.g. distl_wrapper('fake_01_mos_0.10.cbf','distl_params.txt')
	'''
	with open(distl_params) as myfile:
		lines = myfile.readlines()
	lines = [x.strip() for x in lines]
	fixed_params = ['distl.signal_strength', fake_cbf, 'distl.verbose=True']
	return ' '.join(fixed_params+lines)


def ReadDISTLPeakIntensities(image_folder,header_contents,distl_params):
	"""
	Image names must be like: XXX_frameNumber_XXX.img e.g. gaussian_N4_1.5A_3_001.img
	Skips empty/corrupted img files and continues
	Takes a folder of fake detector images with .img extension...
	and header_contents file describing simulation params...
	and distl_params file describing peak detection algorithm
	Returns a data frame with peak loc and intensity vs. frame no 
	1. Converts from img to cbf using img2cbf
	2. Finds peaks and calculates intensities using ReadPeakIntensities	
	3. Data frame is sorted by avg int of peaks over all frames 
	4. Writes the data frame to a pickle file for later use
	"""
	owd = os.getcwd()
	os.chdir(image_folder) # change wd for now, change back to owd at the end

	img_list = [x for x in os.listdir('.') if x.endswith('.img')]
	# sort img_list by numerical value of frame number
	img_list = sorted(img_list, key=lambda x: int(x.split('_')[-2]))
	dic_peaks = {} # stores individual peak intensities
	eps = 5 # two peaks are considered same if within +/- eps in x and y

	for fake_img in img_list:
		fnumber = int(fake_img.split('_')[-2])
		# Skip corrupted img files (file size <10 KB) which cant be converted into cbf
		if os.path.getsize(fake_img) < 10**4:
			print "WARNING: Frame %i is corrupted, will be skipped..." %fnumber
		else:			
			# convert img to cbf
			img2cbf(fake_img,header_contents,keep_original=True)
			fake_cbf = fake_img[0:-3]+'cbf'

			# run DISTL to find peaks using custom params in distl_params file
			# peak finding is highly sensitive to some particular params
			# more info: http://cci.lbl.gov/labelit/html/spotfinder.html
			bash_distl = distl_wrapper(fake_cbf,distl_params)		
			stdout = os.popen(bash_distl).read()
			
			# store individual peak intensities
			stdout_lines = stdout.splitlines()
			for line_num, line in enumerate(stdout_lines):
				if 'Total integrated signal=' in line:
					I = float(line.split()[7][7:])
					x = int(stdout_lines[line_num+2].split('x=')[1][0:4])
					y = int(stdout_lines[line_num+2].split('y=')[1][0:4])

					# check if (x,y) already in dic
					is_match = False
					for (x0,y0) in dic_peaks.keys():
						if x<=x0+eps and x>=x0-eps and y<=y0+eps and y>=y0-eps:
							is_match = True
							x_match = x0; y_match=y0
							break
					
					if is_match == False:
						dic_peaks[(x,y)] = [[fnumber],[I]]
					elif is_match == True:
						dic_peaks[(x_match,y_match)][0].append(fnumber)
						dic_peaks[(x_match,y_match)][1].append(I)

	# write peak intensities to a data frame for convenience
	df_peaks = pd.DataFrame( columns=['x','y'] + range(len(img_list)) )
	for peak in dic_peaks.keys():
		fnumbers,intensities = dic_peaks[peak]
		temp_dic = dict(zip(fnumbers,intensities))
		temp_dic['x'] = peak[0]
		temp_dic['y'] = peak[1]
		df_peaks = df_peaks.append(temp_dic,ignore_index=True)

	# calculate average intensity of each peak over all frames
	df_peaks['I-avg'] = np.nanmean(df_peaks.iloc[:,2::],axis=1)

	# sort peaks by the average intensity of peaks over all frames
	df_peaks.sort_values(by='I-avg', inplace=True, ascending=False)

	# save the data frame to a pickle file
	pickle.dump(df_peaks, open('df_peaks.pickle','wb'))
	os.chdir(owd)

	return df_peaks


def PlotDISTLPeakIntensities(df_distlPeaks,first_N_peaks=None):
	'''
	Takes data frame of peak intensities generated by ReadPeakIntensities
	Plots first_N_peaks at the top of the data frame	
	Sort the data frame by the desired feature before inputting
	By default the data frame is sorted by avg int of peaks over all frames
	
	example1 - process img files from scratch and sort by first frame int: 
	df_peaks = ReadDISTLPeakIntensities(image_folder,header_contents,distl_params)
	df_peaks.sort_values(by=0, inplace=True, ascending=False)
	PlotPeakIntensities(df_peaks,50)

	example2 - load pickle file of the data frame from pre-processed images: 
	df_peaks = pickle.load(open('...path to images/df_peaks.pickle','rb'))
	PlotDISTLPeakIntensities(df_peaks,100)
	'''
	df_peaks = df_distlPeaks.copy()  # make a copy of input dataframe so that it is not modified
	plt.close('all')
	fig, ax = plt.subplots()
	plt.ion()

	if first_N_peaks == None:
		first_N_peaks = df_peaks.shape[0]  # plot all peaks if not provided

	frame_numbers = df_peaks.columns[2:-1]
	# IMPORTANT: Fill NaN intensity values with 0
	df_peaks[frame_numbers] = df_peaks[frame_numbers].fillna(0)	
	for ind in range(first_N_peaks):
		intensities = df_peaks.iloc[ind][2:-1]
		ax.plot(frame_numbers,intensities,marker='o',lw=2)

	ax.set_xlabel('Frame number',fontsize='x-large')
	ax.set_ylabel('DISTL peak intensity\n(pixel value)',fontsize='x-large')
	ax.tick_params(labelsize='large')
	ax.set_ylim(ymin = 0)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	ax.tick_params(axis='both', which='both', length=0)
	ax.grid(which='major',linewidth=0.4)
	ax.grid(which='minor',linewidth=0.2)				
	ax.set_facecolor('0.96')
	for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
	plt.tight_layout()
	plt.show()	


def MlfsomVSdistl(dataframe_mlfsom,dataframe_distl,fnumber):
	"""
	Compares the peak intensities captured from temp files of MLFSOM with
	the distl peak intensities for a desired frame number
	Returns df with mlfsom and distl int. with fractional discrepancy for each peak
	"""
	# Make a copy of input df's to prevent them to be modified
	df_mlfsom = dataframe_mlfsom.copy() 
	df_distl = dataframe_distl.copy()

	# Get the peaks having non-NaN intensity in the fnumber
	sub_peaks = df_mlfsom[~df_mlfsom[fnumber].isna()][['hor','ver','res',fnumber]]
	sub_peaks.sort_values(fnumber,ascending=False,inplace=True)
	sub_distl = df_distl[~df_distl[fnumber].isna()][['y','x',fnumber]]
	sub_distl.sort_values(fnumber,ascending=False,inplace=True)

	eps = 5  # two peaks are considered same if within +/- eps in x and y
	for hkl in sub_peaks.index:
		hor,ver,I_peak = sub_peaks.loc[hkl][['hor','ver',fnumber]]
		match =	sub_distl[sub_distl.y.between(hor-eps,hor+eps) & sub_distl.x.between(ver-eps,ver+eps)]
		if match.shape[0] != 0:
			I_distl = match.iloc[0][fnumber]
			sub_peaks.loc[hkl,'I_distl'] = I_distl

	sub_peaks['disc'] = (1-sub_peaks['I_distl']/sub_peaks[fnumber])
	print 'Frame number: %i' %fnumber
	print 'Matching number of peaks: %i out of %i' %(sum(~sub_peaks.disc.isna()),sub_peaks.shape[0])
	print 'Percent discrepancy (1-I_distl/I_mlfsom) mean: %.1f%%, std: %.1f%%' \
	%(100*np.nanmean(sub_peaks.disc), 100*np.nanstd(sub_peaks.disc))
	return sub_peaks