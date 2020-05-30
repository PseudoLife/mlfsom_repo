import os, fileinput, subprocess, pickle
from glob import glob
from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from toPrecision import toPrecision
from CalcResolution import CalcResolution


def ReadSingleDIALSAscii(in_ascii):
	"""
	Takes hkl asci file generated by DIALS in XDS format
	Returns a dataframe with Fhkl and hkl indices 
	FUTURE: Generalize res calc (works for cubic, tetragonal, orthorhombic) to other space groups  
	"""
	asci = pd.read_csv(in_ascii,header=None, sep='\s+', engine='python', \
	names=['h','k','l','Iasci','sigIasci','xd','yd','zd','rlp','peak','cor','psi'],\
	skiprows=34,skipfooter=1)
	asci['hkl'] = asci.h.apply(lambda x:str(x)) + ' ' + asci.k.apply(lambda x:str(x)) +\
 	' ' + asci.l.apply(lambda x:str(x))
	asci.set_index('hkl',inplace=True)
	asci.sort_values(by='Iasci',ascending=False,inplace=True)
	
	# read refined unit cell constants from asci file
	for line in fileinput.FileInput(in_ascii):
		if 'UNIT_CELL_CONSTANTS' in line:
			unit_cell = line.split()[1::]
			break
	unit_cell = tuple([float(x) for x in unit_cell])
	
	# calculate resolution and create a new 'res' column
	hkl_series = asci.apply(lambda x: (x.h, x.k, x.l), axis=1)
	asci['res'] = hkl_series.apply(lambda x: CalcResolution(unit_cell,x))
	return asci


class DIALSAscii:

	def __init__(self,name_template,space_group):
		self.name_template = name_template
		self.renamed = False
		self.dir_name = os.path.dirname(self.name_template)
		self.base_name = os.path.basename(self.name_template)
		self.space_group = space_group
		# read desc file into a df assuming it is inside folder
		if glob(join(self.dir_name,'*desc.txt')) == []:
			print "MyError: Description file is missing in the folder!"
		else:
			self.description = pd.read_csv( glob(join(self.dir_name,'*desc.txt'))[0],\
				header=None,sep=': ',engine='python', index_col=0,names=['value'] )
		
			# read unit cell params and space group from the original pdb file if exists
			pdb_file = join(self.dir_name,self.description.loc['pdb_file','value'])
			if not os.path.isfile(pdb_file):
				print "MyError: Please place the original pdb file in the folder!"
			else:
				with open(pdb_file) as myfile:
					for line in myfile:
						if 'CRYST1' in line:
							unit_cell = tuple(line.split()[1:7])
							self.unit_cell = tuple([float(x) for x in unit_cell])
							break


	def ReadAllAscii(self):
		"""
		Reads multi DIALS.HKL files and merges them into a single dataframe
		"""
		cwd = os.getcwd()
		os.chdir(self.dir_name)
		asci_list = [x for x in os.listdir(self.dir_name) if x.startswith('DIALS_') and x.endswith('.HKL')]
		asci_byframe = pd.DataFrame()
		for asci_file in asci_list:
			fnumber = int(asci_file.split('_')[1][0:3])
			asci = ReadSingleDIALSAscii(asci_file)

			# Drop duplicate hkl values - drop them all
			# when total oscillation angle large e.g. 80 deg, I observed duplicate hkl values
			# Only ~0.2% of peaks were duplicated but they mess up things, so I just discarded them 
			asci = asci[~asci.index.duplicated(keep=False)]
			
			asci.drop(columns=['h','k','l','sigIasci','psi','rlp'],inplace=True)
			asci.rename(columns={'Iasci':'I_'+str(fnumber), 'xd':'xd_'+str(fnumber), 'yd':'yd_'+str(fnumber),\
				'zd':'zd_'+str(fnumber), 'peak':'peak_'+str(fnumber), \
				'cor':'cor_'+str(fnumber), 'res':'res_'+str(fnumber)}, inplace=True)
			asci_byframe = asci_byframe.join(asci, how='outer')

		# calculate avg peak loc, partiality etc. over all frames
		xd_cols = [x for x in asci_byframe.columns if x.startswith('xd_')]
		asci_byframe['xd'] = asci_byframe[xd_cols].apply(np.mean,axis=1)

		yd_cols = [x for x in asci_byframe.columns if x.startswith('yd_')]
		asci_byframe['yd'] = asci_byframe[yd_cols].apply(np.mean,axis=1)

		zd_cols = [x for x in asci_byframe.columns if x.startswith('zd_')]
		asci_byframe['zd'] = asci_byframe[zd_cols].apply(np.mean,axis=1)

		peak_cols = [x for x in asci_byframe.columns if x.startswith('peak_')]
		asci_byframe['peak'] = asci_byframe[peak_cols].apply(np.mean,axis=1)

		cor_cols = [x for x in asci_byframe.columns if x.startswith('cor_')]
		asci_byframe['cor'] = asci_byframe[cor_cols].apply(np.mean,axis=1)

		res_cols = [x for x in asci_byframe.columns if x.startswith('res_')]
		asci_byframe['res'] = asci_byframe[res_cols].apply(np.mean,axis=1)

		# Drop peak loc etc. by frame, keep only avg values calculated above
		drop_cols = [x for x in asci_byframe.columns if x.startswith('xd_') or \
		x.startswith('yd_') or x.startswith('zd_') or x.startswith('peak_') or \
		x.startswith('cor_') or x.startswith('res_')]
		asci_byframe.drop(columns=drop_cols,inplace=True)

		# Sort cols so that I_1, I_2, ..., xd, yd, zd, ... 
		cols_I = [x for x in asci_byframe.columns if x.startswith('I_')]
		cols_I.sort( key=lambda x: int(x.split('_')[1]) )
		cols_others = [x for x in asci_byframe.columns if not x.startswith('I_')]
		asci_byframe = asci_byframe[cols_I + cols_others]

		# Sort cols from low res to high res
		asci_byframe.sort_values(by='res',ascending=False,inplace=True)
		
		os.chdir(cwd)
		self.peaks = asci_byframe
		return asci_byframe


	def PlotPeakIntensities(self,asci_byframe,first_N_peaks=None,save_img=True):
		"""
		Takes data frame of peak intensities generated by ReadPeakIntensities
		Plots first_N_peaks at the top of the data frame, plots all peaks if None
		Plots either on "linear" or "log" scale
		Sort the data frame by the desired feature before inputting
		By default the data frame is sorted by avg int of peaks over all frames
		"""
		df_peaks = asci_byframe.copy()  # make a copy of input dataframe so that it is not modified
		plt.close('all')
		fig, ax = plt.subplots()
		plt.ion()

		if first_N_peaks == None:
			first_N_peaks = df_peaks.shape[0]  # plot all peaks if not provided

		# get the relevant intensity columns
		cols = [x for x in df_peaks.columns if x.startswith('I_')]
		cols.sort(key=lambda x: int(x.split('_')[1]))
		frame_numbers = [int(x.split('_')[1]) for x in cols]

		# IMPORTANT: Fill NaN intensity values with 0.1 to avoid numerical problems
		df_peaks[cols] = df_peaks[cols].fillna(0.1)
		
		zeros = int(np.log10(df_peaks[cols].max().max())) # Number of excessive zeros to rescale y axis
		for ind in range(first_N_peaks):
			intensities = df_peaks.iloc[ind].loc[cols]
			ax.plot(frame_numbers, intensities/10**zeros,marker='o',lw=2)

		ax.set_ylim(ymin = 0)
		ax.set_xlabel('Frame Number',fontsize='x-large')
		ax.set_ylabel('DIALS Peak Intensity\n(indexed-corrected $\\times 10^%i$)' %zeros,fontsize='x-large')
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
			fig.savefig( join(self.dir_name,"fig_DIALS_PeakIntensities.png"), dpi=200 )
		plt.show()


	def PlotIntegratedIntensities(self,asci_byframe,scale='log',save_img=True):
		"""
		Plots the total integrated intensity vs. frame number
		On "linear" or "log" scale
		"""
		df_peaks = asci_byframe.copy()
		plt.close('all')
		fig, ax = plt.subplots()
		plt.ion()

		# get the relevant intensity columns
		cols = [x for x in df_peaks.columns if x.startswith('I_')]
		cols.sort(key=lambda x: int(x.split('_')[1]))
		frame_numbers = [int(x.split('_')[1]) for x in cols]

		int_intensities = df_peaks[cols].apply(np.nansum)
		if scale == 'linear':	
			zeros = int(np.log10(max(int_intensities)))
			int_intensities = int_intensities/10**zeros
			ax.plot(frame_numbers,int_intensities,marker='o',lw=2)
			ax.set_ylabel('DIALS Integrated Intensity\n(indexed-corrected $\\times 10^%i$)' %zeros,fontsize='x-large')
			ax.yaxis.set_minor_locator(AutoMinorLocator())		
		elif scale == 'log':
			ax.plot(frame_numbers,int_intensities,marker='o',lw=2)
			ax.set_yscale('log')
			ax.set_ylabel('DIALS Integrated Intensity \n(indexed-corrected)',fontsize='x-large')
			
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
			fig.savefig( join(self.dir_name,"fig_DIALS_IntIntensities.png"), dpi=200 )
		plt.show()


	def ShellIntensities(self,asci_byframe,N_shells=10,npeak_ratio=1,normalized=False):
		"""
		Takes individual peak intensties and number of shells
		Returns dataframe of shells with int. int., d_min, d_max, n-peaks and D_half in each shell
		npeak_ratio: int, ratio of the N of peaks in highest to lowest res. shell
		npeak_ratio = 1 for equal number of peaks per shell
		"""
		df_peaks = asci_byframe.copy()
		df_peaks.sort_values('res',ascending=False,inplace=True)

		print 'npeak_ratio: %i' %npeak_ratio
		alpha = np.log(npeak_ratio)/N_shells  
		x = np.arange(N_shells)
		peaks_each_shell = np.exp(alpha*x)/sum(np.exp(alpha*x)) * len(df_peaks)
		peaks_each_shell = peaks_each_shell.astype(int)

		# get the relevant intensity columns
		int_cols = [x for x in df_peaks.columns if x.startswith('I_')]
		int_cols.sort(key=lambda x: int(x.split('_')[1]))
		frame_numbers = [int(x.split('_')[1]) for x in int_cols]

		df_shells = pd.DataFrame(columns=['d_min','d_max','N_peaks','D_half']+int_cols)
		# df_shells stores the integrated intensity in each shell with the res. range and n-peaks

		start_index = 0
		for shell_index in range(N_shells):
			peaks_x_shell = peaks_each_shell[shell_index]   # N of peaks in a given shell
			end_index = start_index + peaks_x_shell
			d_max = df_peaks.iloc[start_index]['res']
			d_max = round(d_max,2)
			d_min = df_peaks.iloc[end_index-1]['res']
			d_min = round(d_min,2)
			intensities = df_peaks.iloc[start_index:end_index][int_cols].apply(np.nansum)

			# calculate half_dose by fitting exponential to dose curves
			intensities[intensities<=0] = 0.1  # replace negative values with 0.1 to avoid numerical problems
			fit_c1, fit_c2 = np.polyfit(frame_numbers, np.log(intensities.to_list()), 1)
			half_dose = round(-np.log(2)/fit_c1,2)

			df_shells.loc[shell_index,['d_min','d_max','N_peaks','D_half']] = [d_min,d_max,peaks_x_shell,half_dose]
			df_shells.loc[shell_index,int_cols] = intensities.to_list()
			start_index += peaks_x_shell

		if normalized:
			# normalize shell intensities by the first frame
			for ind in df_shells.index:
				df_shells.loc[ind,int_cols[0]::] = df_shells.loc[ind,int_cols[0]::]/df_shells.loc[ind,int_cols[0]]

		self.shells = df_shells
		return df_shells


	def PlotShellIntensities(self,asci_byframe,N_shells=10,npeak_ratio=1,normalized=True,save_img=True):
		"""
		Takes peak intensities and plots peak intensities by shell on log scale
		Returns dataframe of shells with int. int., d_min, d_max, n-peaks and D_half in each shell
		npeak_ratio: int, ratio of the N of peaks in highest to lowest res. shell
		npeak_ratio = 1 for equal number of peaks per shell
		"""
		df_shells = self.ShellIntensities(asci_byframe,N_shells,npeak_ratio,normalized=normalized)
		plt.close('all')
		fig, ax = plt.subplots()
		plt.ion()

		# get the relevant intensity columns
		cols = [x for x in df_shells.columns if x.startswith('I_')]
		cols.sort(key=lambda x: int(x.split('_')[1]))
		frame_numbers = [int(x.split('_')[1]) for x in cols]

		for shell_index in df_shells.index:
			d_min,d_max,half_dose = df_shells.loc[shell_index,['d_min','d_max','D_half']]
			intensities = df_shells.loc[shell_index,cols]

			#ax.plot(frame_numbers,intensities,marker='o',lw=2,label='%s - %s | %s' \
			#	%(toPrecision(d_max,3),toPrecision(d_min,3),toPrecision(half_dose,3)))

			ax.scatter(frame_numbers,intensities,label='%s - %s | %s' \
				%(toPrecision(d_max,3),toPrecision(d_min,3),toPrecision(half_dose,3)))

		ax.set_yscale('log')
		#leg = ax.legend(fontsize='x-small',markerscale=0)
		leg = ax.legend(fontsize='x-small')
		leg.set_title('Resolution ($\mathrm{\AA)}$ | $D_{1/2}$',prop={'size':'x-small'})
		ax.set_xlabel('Frame Number',fontsize='x-large')
		ax.set_ylabel('DIALS Integrated Intensity \n(indexed-corrected)',fontsize='x-large')
		ax.tick_params(labelsize='large')
		ax.xaxis.set_minor_locator(AutoMinorLocator())
		ax.tick_params(axis='both', which='both', length=0)
		ax.grid(which='major',linewidth=0.4)
		ax.grid(which='minor',linewidth=0.2)
		ax.set_facecolor('0.95')
		for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
		plt.tight_layout()
		if normalized:
			ymin = max(10**-2,df_shells[cols].min().min())  ########################################
			ax.set_ylim(ymin=ymin)  # don't show noise caused by NaN intensities (reassigned to 0.1) 
		if save_img:
			fig.savefig(join(self.dir_name,"fig_DIALS_Shellintensities.png"),dpi=200)
		plt.show()
		self.shells = df_shells
		return df_shells


	def PlotDhalfvsResolution(self,asci_byframe,N_shells=10,npeak_ratio=1,\
		K2=0.2,K1=1.0,weights=(0.75,0.25),save_img=True):
		"""
		Plots Dhalf vs. resolution
		K2: Constant in Dhalf = K * res**2 e.g. 2.0
		K1: Constant in Dhalf = K * res**1
		weights: tuple, weights of d_min, d_max respectively in d_space = d_min*w_min + d_max*w_max
		w_min should be greater than w_max as there are more peaks at higher resolutions
		npeak_ratio: int, ratio of the N of peaks in highest to lowest res. shell 
		npeak_ratio = 1 for equal number of peaks per shell
		"""
		df_shells = self.ShellIntensities(asci_byframe,N_shells,npeak_ratio,normalized=True)
		
		plt.close('all')
		fig, ax = plt.subplots()
		plt.ion()

		# for a given shell define d_space = d_min*w_min + d_max*w_max
		d_space = df_shells.d_min*weights[0] + df_shells.d_max*weights[1]
		D_half = df_shells.D_half
		plt.scatter(d_space,D_half,edgecolors='black',s=150, label=None)
		# fitting line
		x_min = d_space.iloc[-1] * 0.90
		x_max = d_space.iloc[0] * 1.1
		x = np.linspace(x_min, x_max, 100)
		y2 = K2 * x**2
		y1 = K1 * x**1
		ax.plot(x, y2, linestyle = '--', color='tab:red', linewidth=2, label='$D_{1/2} \propto d^{2}$')
		ax.plot(x, y1, linestyle = '--', color='tab:brown', linewidth=2, label='$D_{1/2} \propto d^{1}$')
		ax.legend(fontsize='xx-large')

		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlabel('Resolution ($\mathrm{\\AA}$)',fontsize='x-large')
		ax.set_ylabel('DIALS Half-dose (arbitrary)',fontsize='x-large')
		ax.tick_params(axis='both', which='both', length=0, labelsize='x-large',pad=10)
		ax.grid(which='major',linewidth=0.4)
		ax.grid(which='minor',linewidth=0.2)
		ax.set_facecolor('0.95')
		for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
		plt.tight_layout()
		if save_img:
			fig.savefig(join(self.dir_name,"fig_DIALS_DhalfvsRes.png"),dpi=200)
		plt.show()
		return df_shells

"""
if __name__ == "__main__":
	name_template = ...

"""
