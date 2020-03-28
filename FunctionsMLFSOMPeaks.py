import pandas as pd
import os, pickle
from os.path import join
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from matplotlib.ticker import AutoMinorLocator

from toPrecision import toPrecision

cell_a = 77.25
cell_b = 77.25
cell_c = 38.66


def ReadSingleXYI(infile):
	"""
	Read peak intensties from a single .XYI file e.g. mlfsom_tempfile1_preds_1.XYI
	Outputs pandas dataframe with index of hkl values and cols: hor-coord, ver-coord and int. of the peaks
	In calculation of resolution, rectangular unit cell is assumed like in lysozyme
	"""
	df = pd.read_csv(infile,sep=' ',header=None, \
		names=['ver', 'hor', 'I', 'psi', 'rsize', 'tsize', 'h', 'k', 'l'])
	df['hkl'] = df.h.apply(lambda x:str(x)) + ' ' + df.k.apply(lambda x:str(x)) + ' ' + df.l.apply(lambda x:str(x))
	df['res'] = 1/np.sqrt((df.h/cell_a)**2 + (df.k/cell_b)**2 + (df.l/cell_c)**2) # assume rect. cell
	df.drop(columns=['psi','rsize','tsize','h','k','l'],inplace=True)
	df_hkl = df.groupby('hkl').agg({'ver':'mean', 'hor':'mean', 'res':'mean', 'I':'sum'})
	df_hkl.sort_values('I',ascending=False,inplace=True)

	# reorder the columns
	cols = df_hkl.columns.tolist()
	cols = cols[1::] + cols[0:1]
	df_hkl = df_hkl[cols]
	return df_hkl


def ReadPeakIntensities(image_folder):
	"""
	Takes an image folder of .XYI files and outputs a data frame with
	hkl, peak intensities, hor and ver coords, resolution etc. 
	Skips missing or corrupted XYI files
	"""
	owd = os.getcwd()
	os.chdir(image_folder) # change wd for now, change back to owd at the end

	df_peaks = pd.DataFrame()
	for XYI_file in glob('mlfsom_tempfile*.XYI'):
		# skip corrupted XYI files whose size is 0 bytes
		if os.path.getsize(XYI_file) != 0:
			fnumber = int(XYI_file.split('_')[1][8::]) # fnumber is a two digit str, format '07' or '26'
			df_hkl = ReadSingleXYI(XYI_file)
			df_hkl.rename(columns={'I':fnumber,'hor':'hor_'+str(fnumber),'ver':'ver_'+str(fnumber),\
				'res':'res_'+str(fnumber)},inplace=True)
			df_peaks = df_peaks.join(df_hkl, how='outer')

	# Add average I, hor, ver and res columns
	intensity_cols = [col for col in df_peaks.columns if type(col)==int]
	df_peaks['I-avg'] = np.nanmean(df_peaks[intensity_cols],axis=1)

	hor_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('hor_')]
	df_peaks['hor'] = np.nanmean(df_peaks[hor_cols],axis=1)

	ver_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('ver_')]
	df_peaks['ver'] = np.nanmean(df_peaks[ver_cols],axis=1)

	res_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('res_')]
	df_peaks['res'] = np.nanmean(df_peaks[res_cols],axis=1)	

	# Return only avg. hor, ver, res, int and int for each frame (drop e.g. hor)
	keep_cols = ['hor','ver','res'] + np.sort([col for col in df_peaks.columns if type(col)==int]).tolist() + ['I-avg']
	df_peaks = df_peaks[keep_cols]
	df_peaks.sort_values('I-avg',ascending=False,inplace=True) # sort by average int. of peaks
	
	os.chdir(owd)
	return df_peaks


def PlotPeakIntensities(df_mlfsomPeaks,first_N_peaks=None,save_img=False):
	"""
	Takes data frame of peak intensities generated by ReadPeakIntensities
	Plots first_N_peaks at the top of the data frame, plots all peaks if None
	Plots either on "linear" or "log" scale
	Sort the data frame by the desired feature before inputting
	By default the data frame is sorted by avg int of peaks over all frames
	
	example1 - process img files from scratch and sort by first frame int: 
	df_peaks = ReadXYIPeakIntensities(image_folder)
	df_peaks.sort_values(by=0, inplace=True, ascending=False)
	PlotPeakIntensities(df_peaks,50)

	example2 - load pickle file of the data frame from pre-processed images: 
	df_peaks = pickle.load(open('...path to images/df_peaks.pickle','rb'))
	PlotPeakIntensities(df_peaks,100)
	"""
	df_peaks = df_mlfsomPeaks.copy()  # make a copy of input dataframe so that it is not modified
	plt.close('all')
	fig, ax = plt.subplots()
	plt.ion()

	if first_N_peaks == None:
		first_N_peaks = df_peaks.shape[0]  # plot all peaks if not provided

	frame_numbers = np.sort([col for col in df_peaks.columns if type(col)==int]).tolist()
	# IMPORTANT: Fill NaN intensity values with 0.1 to avoid numerical problems
	df_peaks[frame_numbers] = df_peaks[frame_numbers].fillna(0.1)
	
	zeros = int(np.log10(df_peaks[frame_numbers].max().max())) # Number of excessive zeros to rescale y axis
	for ind in range(first_N_peaks):
		intensities = df_peaks.iloc[ind].loc[frame_numbers]		
		ax.plot(frame_numbers,intensities/10**zeros,marker='o',lw=2)

	ax.set_ylim(ymin = 0)
	ax.set_xlabel('Frame Number',fontsize='x-large')
	ax.set_ylabel('MLFSOM Peak Intensity\n(pixel value$\\times 10^%i$)' %zeros,fontsize='x-large')
	ax.tick_params(labelsize='large')
	ax.xaxis.set_minor_locator(AutoMinorLocator())	
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	ax.tick_params(axis='both', which='both', length=0)
	ax.grid(which='major',linewidth=0.4)
	ax.grid(which='minor',linewidth=0.2)				
	ax.set_facecolor('0.96')
	for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
	plt.tight_layout()
	if save_img:
		fig.savefig("fig_XYI_PeakIntensities.png",dpi=200)
	plt.show()


def PlotIntegratedIntensities(df_mlfsomPeaks,scale='log',save_img=False):
	"""
	Plots the total integrated intensity vs. frame number
	On "linear" or "log" scale
	"""
	df_peaks = df_mlfsomPeaks.copy()
	plt.close('all')
	fig, ax = plt.subplots()
	plt.ion()

	frame_numbers = np.sort([col for col in df_peaks.columns if type(col)==int]).tolist()	
	int_intensities = df_peaks[frame_numbers].apply(np.nansum)
	if scale == 'linear':	
		zeros = int(np.log10(max(int_intensities)))
		int_intensities = int_intensities/10**zeros
		ax.plot(frame_numbers,int_intensities,marker='o',lw=2)
		ax.set_ylabel('MLFSOM Integrated Intensity\n(pixel value$\\times 10^%i$)' %zeros,fontsize='x-large')
		ax.yaxis.set_minor_locator(AutoMinorLocator())		
	elif scale == 'log':
		ax.plot(frame_numbers,int_intensities,marker='o',lw=2)
		ax.set_yscale('log')
		ax.set_ylabel('MLFSOM Integrated Intensity\n(pixel value)',fontsize='x-large')
		
	ax.set_xlabel('Frame Number',fontsize='x-large')
	ax.tick_params(labelsize='large')
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.tick_params(axis='both', which='both', length=0)
	ax.grid(which='major',linewidth=0.4)
	ax.grid(which='minor',linewidth=0.2)				
	ax.set_facecolor('0.96')
	for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
	plt.tight_layout()
	if save_img:
		fig.savefig("fig_XYI_intIntensities.png",dpi=200)
	plt.show()	


def ShellIntensities(df_mlfsomPeaks,N_shells=10,normalized=False):
	"""
	Takes individual peak intensties and number of shells
	Returns dataframe of shells with int. int., d_min, d_max, n-peaks and D_half in each shell
	"""
	df_peaks = df_mlfsomPeaks.copy()
	df_peaks.sort_values('res',ascending=False,inplace=True)

	# Equal number of peaks in each shell
	peaks_per_shell = int(len(df_peaks)/float(N_shells))
	peaks_each_shell = np.ones(N_shells).astype(int) * peaks_per_shell 

	frame_numbers = np.sort([col for col in df_peaks.columns if type(col)==int]).tolist()
	df_shells = pd.DataFrame(columns=['d_min','d_max','N_peaks','D_half']+frame_numbers)
	# df_shells stores the integrated intensity in each shell with the res. range and n-peaks

	start_index = 0
	for shell_index in range(N_shells):
		peaks_x_shell = peaks_each_shell[shell_index]   # N of peaks in a given shell
		end_index = start_index + peaks_x_shell
		d_max = df_peaks.iloc[start_index]['res']
		d_max = round(d_max,2)
		d_min = df_peaks.iloc[end_index-1]['res']
		d_min = round(d_min,2)
		intensities = df_peaks.iloc[start_index:end_index][frame_numbers].apply(np.nansum)

		# calculate half_dose by fitting exponential to dose curves
		fit_c1, fit_c2 = np.polyfit(frame_numbers, np.log(intensities.replace(0,0.1).to_list()), 1)
		half_dose = round(-np.log(2)/fit_c1,2)

		df_shells.loc[shell_index,['d_min','d_max','N_peaks','D_half']] = [d_min,d_max,peaks_x_shell,half_dose]
		df_shells.loc[shell_index,frame_numbers] = intensities
		start_index += peaks_x_shell

	if normalized:
		cols = [col for col in df_shells.columns if type(col)==int]
		# normalize shell intensities by the first frame
		for ind in df_shells.index:
			df_shells.loc[ind,cols] = df_shells.loc[ind,cols]/df_shells.loc[ind,cols[0]]
	return df_shells


def PlotShellIntensities(df_mlfsomPeaks,N_shells=10,normalized=True,save_img=False):
	"""
	Takes peak intensities and plots peak intensities by shell on log scale
	Returns dataframe of shells with int. int., d_min, d_max, n-peaks and D_half in each shell
	"""
	df_shells = ShellIntensities(df_mlfsomPeaks,N_shells,normalized=normalized)
	plt.close('all')
	fig, ax = plt.subplots()
	plt.ion()

	frame_numbers = np.sort([col for col in df_shells.columns if type(col)==int]).tolist()
	for shell_index in df_shells.index:
		d_min,d_max,half_dose = df_shells.loc[shell_index,['d_min','d_max','D_half']]
		intensities = df_shells.loc[shell_index,frame_numbers]

		ax.plot(frame_numbers,intensities,marker='o',lw=2,label='%s - %s | %s' \
			%(toPrecision(d_max,3),toPrecision(d_min,3),toPrecision(half_dose,3)))

	ax.set_yscale('log')
	leg = ax.legend(fontsize='x-small',markerscale=0)
	leg.set_title('Resolution ($\mathrm{\AA)}$ | $D_{1/2}$',prop={'size':'x-small'})
	ax.set_xlabel('Frame Number',fontsize='x-large')
	ax.set_ylabel('MLFSOM Integrated Intensity\n(pixel value)',fontsize='x-large')
	ax.tick_params(labelsize='large')
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.tick_params(axis='both', which='both', length=0)
	ax.grid(which='major',linewidth=0.4)
	ax.grid(which='minor',linewidth=0.2)
	ax.set_facecolor('0.95')
	for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
	plt.tight_layout()
	if save_img:
		fig.savefig("fig_XYI_shellIntensities.png",dpi=200)
	plt.show()
	return df_shells


def AggregateCustomFramesPeakIntensities(image_folder,start_frame_id,end_frame_id):
	"""
	Gaussian beam only
	Takes an image folder of .XYI files and a custom range of frame ids
	Returns weighted sum of intensties
	Image_folder should have the exp. desc file e.g. exp-gaussian_N4_1.5A-desc.txt...
	which has the weights of each frame id
	Start and end id's inclusive
	Returns None if there is any missing or corrupted XYI file between start and end id
	"""
	owd = os.getcwd()
	os.chdir(image_folder) # change wd for now, change back to owd at the end

	# Get the line where (frame,weight,experiment) starts
	desc_file = glob(join(image_folder,'exp-*desc.txt'))[0]
	with open(desc_file) as myfile:
		line_list = myfile.readlines()
		for line_num, line in enumerate(line_list):
			if "frame" in line and "weight" in line and "experiment_ID" in line:
				start_line_num = line_num
				break
	# Read in ID's and weights of the frames 
	df_weights = pd.read_csv(desc_file,header=None,skiprows=start_line_num+1,\
		names=['macro','weight','ID'],index_col=2)
	df_weights.index = df_weights.index.str.replace(')','').astype(int)
	df_weights.macro = df_weights.macro.str.replace('(','').astype(int)

	df_peaks = pd.DataFrame()
	XYI_file_list = [fi for fi in os.listdir('.') if fi.endswith('.XYI')]
	for fnumber in range(start_frame_id,end_frame_id+1):
		XYI_file = 'mlfsom_tempfile'+str(fnumber)+'_preds_1.XYI'
		# XYI_file is None if there is any missing XYI file between start_id and end_id
		if not os.path.isfile(XYI_file) or os.path.getsize(XYI_file) == 0:
			print "WARNING: At least one XYI file between frames %i-%i missing or corrupted" \
			%(start_frame_id,end_frame_id)
			return None
		else:
			df_hkl = ReadSingleXYI(XYI_file)
			df_hkl.rename(columns={'I':fnumber,'hor':'hor_'+str(fnumber),'ver':'ver_'+str(fnumber),\
				'res':'res_'+str(fnumber)},inplace=True)
			df_peaks = df_peaks.join(df_hkl, how='outer')

	# Add weighted average I, and basic average for hor, ver and res columns
	# PS: I don't bother with weighted average for hor, ver and res columns as... 
	# they don't vary much from frame to frame
	intensity_cols = [col for col in df_peaks.columns if type(col)==int]
	weighted_I_sum = pd.Series([0]*len(df_peaks)) # initialize weighted sum with zeros
	for fnumber in intensity_cols:
		weighted_I_sum = np.nansum(\
			[weighted_I_sum, df_peaks[fnumber]*df_weights.loc[fnumber,'weight']],axis=0)
	df_peaks['I-sum'] = weighted_I_sum

	hor_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('hor_')]
	df_peaks['hor'] = np.nanmean(df_peaks[hor_cols],axis=1)

	ver_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('ver_')]
	df_peaks['ver'] = np.nanmean(df_peaks[ver_cols],axis=1)

	res_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('res_')]
	df_peaks['res'] = np.nanmean(df_peaks[res_cols],axis=1)	

	# Return only avg. hor, ver, res, int and int for each frame
	keep_cols = ['hor','ver','res'] + np.sort([col for col in df_peaks.columns if type(col)==int]).tolist() + ['I-sum']
	df_peaks = df_peaks[keep_cols]
	df_peaks.sort_values('I-sum',ascending=False,inplace=True) # sort by average int. of peaks
	
	os.chdir(owd)
	return df_peaks


def ReadGaussianPeakIntensities(image_folder):
	"""
	Gaussian beam only
	Takes image_folder of a Gaussian beam run
	Returns peaks with combined (summed) intensities from different IDs
	image_folder should have the exp. desc file e.g. exp-gaussian_N4_1.5A-desc.txt
	Skips the macro frames whose at least one XYI file is missing or corrupted
	"""
	owd = os.getcwd()
	os.chdir(image_folder) # change wd for now, change back to owd at the end

	desc_file = glob(join(image_folder,'exp-*desc.txt'))[0]
	desc_df = pd.read_csv(desc_file,sep=': ',header=None,\
		index_col=0,engine='python',names=['param','value'])
	N_grid = int(desc_df.loc['N_grid','value'])
	frames_per_macro = N_grid*(N_grid+2)/8
	N_macros = int(desc_df.loc['frames','value'])
	begin_id = 0
	df_peaks = pd.DataFrame()
 	for macro_num in range(N_macros):
		df_macro = AggregateCustomFramesPeakIntensities(image_folder,begin_id,begin_id+frames_per_macro-1)
		# df_macro is None if any missing XYI files between start and end ids
		if type(df_macro) == pd.core.frame.DataFrame:
			df_macro.drop(columns=[col for col in df_macro.columns if type(col)==int],inplace=True)
			df_macro.rename(columns={'I-sum':macro_num,'hor':'hor_'+str(macro_num),'ver':'ver_'+str(macro_num),\
				'res':'res_'+str(macro_num)},inplace=True)
			df_peaks = df_peaks.join(df_macro, how='outer')
		begin_id += frames_per_macro

	# Add average I, hor, ver and res columns
	intensity_cols = [col for col in df_peaks.columns if type(col)==int]
	df_peaks['I-avg'] = np.nanmean(df_peaks[intensity_cols],axis=1)

	hor_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('hor_')]
	df_peaks['hor'] = np.nanmean(df_peaks[hor_cols],axis=1)

	ver_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('ver_')]
	df_peaks['ver'] = np.nanmean(df_peaks[ver_cols],axis=1)

	res_cols = [col for col in df_peaks.columns if type(col)==str and col.startswith('res_')]
	df_peaks['res'] = np.nanmean(df_peaks[res_cols],axis=1)	

	# Return only avg. hor, ver, res, int and int for each macro-frame
	keep_cols = ['hor','ver','res'] + np.sort([col for col in df_peaks.columns if type(col)==int]).tolist() + ['I-avg']
	df_peaks = df_peaks[keep_cols]
	df_peaks.sort_values('I-avg',ascending=False,inplace=True) # sort by average int. of peaks

	os.chdir(owd)
	return df_peaks