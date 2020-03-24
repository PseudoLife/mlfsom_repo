import os, fileinput, subprocess
from glob import glob
from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from toPrecision import toPrecision
from CalcResolution import CalcResolution


def ReadSingleXDSAscii(in_ascii):
	"""
	Takes hkl asci file generated by XDS-correct
	Returns a dataframe with Fhkl and hkl indices 
	FUTURE: Generalize res calc (works for cubic, tetragonal, orthorhombic) to other space groups  
	"""
	asci = pd.read_csv(in_ascii,header=None, sep='\s+', engine='python', \
	names=['h','k','l','Iasci','sigIasci','xd','yd','zd','rlp','peak','corr','psi'],\
	skiprows=47,skipfooter=1)
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
	for ind in asci.index:
		hkl = tuple(asci.loc[ind,['h','k','l']])
		asci.loc[ind,'res'] = CalcResolution(unit_cell,hkl)
	return asci


class XDSAscii:

	def __init__(self,name_template,space_group):
		self.name_template = name_template
		self.dir_name = os.path.dirname(self.name_template)
		self.base_name = os.path.basename(self.name_template)
		self.space_group = space_group
		# read desc file into a df assuming it is inside folder
		if glob(join(self.dir_name,'*desc.txt')) == []:
			print "MyError: Description file is missing in the current folder!"
		else:
			self.description = pd.read_csv( glob(join(self.dir_name,'*desc.txt'))[0],\
				header=None,sep=': ',engine='python', index_col=0,names=['value'] )
		
			# read unit cell params and space group from the original pdb file if exists
			pdb_file = join(self.dir_name,self.description.loc['pdb_file','value'])
			if not os.path.isfile(pdb_file):
				print "MyError: Please place the original pdb file (%s) in the folder!" %pdb_file
			else:
				with open(pdb_file) as myfile:
					for line in myfile:
						if 'CRYST1' in line:
							unit_cell = tuple(line.split()[1:7])
							self.unit_cell = unit_cell
							break


	def ProcessStills(self):
		"""
		Processes fake img files using XDS (indexing, integrating, correcting)
		Renames and saves CORRECT.LP and XDS_ASCII.HKL files from eachs still frame
		desc file and the original pdb file must be under the folder for this to work
		"""
		cwd = os.getcwd()
		os.chdir(self.dir_name)

		img_list = glob(self.base_name.replace('???','*'))
		for img in img_list:
			fnumber = img.split(self.base_name.split('???')[0])[1][0:3]

			for line in fileinput.FileInput('XDS.INP',inplace=True):
				line = line.strip()
				if 'NAME_TEMPLATE_OF_DATA_FRAMES' in line:
					line = 'NAME_TEMPLATE_OF_DATA_FRAMES= %s' %self.name_template
				if 'OSCILLATION_RANGE' in line:
					line = 'OSCILLATION_RANGE= %s' %self.description.loc['osc','value']
				if 'INCLUDE_RESOLUTION_RANGE' in line:
					line = 'INCLUDE_RESOLUTION_RANGE= 100 %s' %self.description.loc['resolution','value']
				if 'DETECTOR_DISTANCE' in line:
					line = 'DETECTOR_DISTANCE= %s' %self.description.loc['distance','value']
				if 'DATA_RANGE' in line:
					line = 'DATA_RANGE= %i %i' %(int(fnumber),int(fnumber))
				if 'BACKGROUND_RANGE' in line:
					line = 'BACKGROUND_RANGE= %i %i' %(int(fnumber),int(fnumber))
				if 'SPOT_RANGE' in line:
					line = 'SPOT_RANGE= %i %i' %(int(fnumber),int(fnumber))
				if 'UNIT_CELL_CONSTANTS' in line:
					line = 'UNIT_CELL_CONSTANTS= %s' \
					%(' '.join([str(x) for x in self.unit_cell]))
				if 'SPACE_GROUP_NUMBER' in line:
					line = 'SPACE_GROUP_NUMBER= %i' %(self.space_group)
				print line

			subprocess.call(['xds'])
			subprocess.call(['mv', 'XDS_ASCII.HKL', 'XDS_ASCII_%s.HKL' %fnumber])
			subprocess.call(['mv', 'CORRECT.LP', 'CORRECT_%s.LP' %fnumber])
		os.chdir(cwd)


	def ReadAllAscii(self):
		"""
		Reads multi XDS_ASCII.HKL files and merges them into a single dataframe
		"""
		cwd = os.getcwd()
		os.chdir(self.dir_name)
		asci_list = [x for x in os.listdir(self.dir_name) if x.startswith('XDS_ASCII_') and x.endswith('.HKL')]
		asci_byframe = pd.DataFrame()
		for asci_file in asci_list:
			fnumber = int(asci_file.split('_')[2][0:3])
			asci = ReadSingleXDSAscii(asci_file)

			asci.drop(columns=['h','k','l','sigIasci','psi'],inplace=True)
			asci.rename(columns={'Iasci':'I_'+str(fnumber), 'xd':'xd_'+str(fnumber), 'yd':'yd_'+str(fnumber),\
				'zd':'zd_'+str(fnumber), 'rlp':'rlp_'+str(fnumber), 'peak':'peak_'+str(fnumber), \
				'corr':'corr_'+str(fnumber), 'res':'res_'+str(fnumber)}, inplace=True)
			asci_byframe = asci_byframe.join(asci, how='outer')

		os.chdir(cwd)
		return asci_byframe


	def ReadMosCell(self):
		"""
		Reads mosaicity and unit cell-a and unit cell-c from CORRECT.LP files
		Returns dataframe
		"""
		cwd = os.getcwd()
		os.chdir(self.dir_name)

		LP_list = [x for x in os.listdir(self.dir_name) if x.startswith('CORRECT_') and x.endswith('.LP')]
		mosaicity_cell_byframe = pd.DataFrame(columns=['mosaicity','cellA','cellB','cellC'])

		for LP_file in LP_list:
			for line in fileinput.FileInput(LP_file):
				fnumber = int(LP_file.split('_')[1][0:3])
				if 'CRYSTAL MOSAICITY (DEGREES)' in line:
					mosaicity_cell_byframe.loc[fnumber,'mosaicity'] = float(line.split()[-1])
				if 'UNIT CELL PARAMETERS' in line:
					mosaicity_cell_byframe.loc[fnumber,'cellA'] = float(line.split()[3])
					mosaicity_cell_byframe.loc[fnumber,'cellB'] = float(line.split()[4])
					mosaicity_cell_byframe.loc[fnumber,'cellC'] = float(line.split()[5])
		
		mosaicity_cell_byframe.sort_index(inplace=True)
		os.chdir(cwd)
		return mosaicity_cell_byframe


	def PlotPeakIntensities(self,asci_byframe,first_N_peaks=None,save_img=False):
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
		ax.set_ylabel('XDS Peak Intensity\n(indexed-corrected $\\times 10^%i$)' %zeros,fontsize='x-large')
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
			fig.savefig( join(self.dir_name,"fig_peakIntensities.png"), dpi=300 )
		plt.show()


	def PlotIntegratedIntensities(self,asci_byframe,scale='log',save_img=False):
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
			ax.set_ylabel('XDS Integrated Intensity\n(indexed-corrected $\\times 10^%i$)' %zeros,fontsize='x-large')
			ax.yaxis.set_minor_locator(AutoMinorLocator())		
		elif scale == 'log':
			ax.plot(frame_numbers,int_intensities,marker='o',lw=2)
			ax.set_yscale('log')
			ax.set_ylabel('XDS Integrated Intensity \n(indexed-corrected)',fontsize='x-large')
			
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
			fig.savefig( join(self.dir_name,"fig_intIntensities.png"), dpi=300 )
		plt.show()


	def ShellIntensities(self,asci_byframe,N_shells=10,normalized=False):
		"""
		Takes individual peak intensties and number of shells
		Returns dataframe of shells with int. int., d_min, d_max, n-peaks and D_half in each shell
		FUTURE: 1-more peaks in low res shells 2-normalized int
		"""
		df_peaks = asci_byframe.copy()
		res_cols = [x for x in asci_byframe.columns if x.startswith('res_')]
		df_peaks['res'] = df_peaks[res_cols].apply(np.mean,axis=1)
		df_peaks.sort_values('res',ascending=False,inplace=True)

		# Equal number of peaks in each shell
		peaks_per_shell = int(len(df_peaks)/float(N_shells))
		peaks_each_shell = np.ones(N_shells).astype(int) * peaks_per_shell 

		# get the relevant intensity columns
		cols = [x for x in df_peaks.columns if x.startswith('I_')]
		cols.sort(key=lambda x: int(x.split('_')[1]))
		frame_numbers = [int(x.split('_')[1]) for x in cols]

		df_shells = pd.DataFrame(columns=['d_min','d_max','N_peaks','D_half']+cols)
		# df_shells stores the integrated intensity in each shell with the res. range and n-peaks

		start_index = 0
		for shell_index in range(N_shells):
			peaks_x_shell = peaks_each_shell[shell_index]   # N of peaks in a given shell
			end_index = start_index + peaks_x_shell
			d_max = df_peaks.iloc[start_index]['res']
			d_max = round(d_max,2)
			d_min = df_peaks.iloc[end_index-1]['res']
			d_min = round(d_min,2)
			intensities = df_peaks.iloc[start_index:end_index][cols].apply(np.nansum)

			# calculate half_dose by fitting exponential to dose curves
			intensities[intensities<0] = 0.1  # replace negative values with 0.1 to avoid numerical problems 
			fit_c1, fit_c2 = np.polyfit(frame_numbers, np.log(intensities.to_list()), 1)
			half_dose = round(-np.log(2)/fit_c1,2)

			df_shells.loc[shell_index,['d_min','d_max','N_peaks','D_half']] = [d_min,d_max,peaks_x_shell,half_dose]
			df_shells.loc[shell_index,cols] = intensities.to_list()
			start_index += peaks_x_shell

		if normalized:
			# normalize shell intensities by the first frame
			for ind in df_shells.index:
				df_shells.loc[ind,cols[0]::] = df_shells.loc[ind,cols[0]::]/df_shells.loc[ind,cols[0]]

		return df_shells


	def PlotShellIntensities(self,asci_byframe,N_shells=10,normalized=True,save_img=False):
		"""
		Takes peak intensities and plots peak intensities by shell on log scale
		Returns dataframe of shells with int. int., d_min, d_max, n-peaks and D_half in each shell
		"""
		df_shells = self.ShellIntensities(asci_byframe,N_shells,normalized=True)
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

			ax.plot(frame_numbers,intensities,marker='o',lw=2,label='%s - %s | %s' \
				%(toPrecision(d_max,3),toPrecision(d_min,3),toPrecision(half_dose,3)))

		ax.set_yscale('log')
		leg = ax.legend(fontsize='xx-small',markerscale=0)
		leg.set_title('Resolution ($\mathrm{\AA)}$ | $D_{1/2}$',prop={'size':'xx-small'})
		ax.set_xlabel('Frame Number',fontsize='x-large')
		ax.set_ylabel('XDS Integrated Intensity \n(indexed-corrected)',fontsize='x-large')
		ax.tick_params(labelsize='large')
		ax.xaxis.set_minor_locator(AutoMinorLocator())
		ax.tick_params(axis='both', which='both', length=0)
		ax.grid(which='major',linewidth=0.4)
		ax.grid(which='minor',linewidth=0.2)
		ax.set_facecolor('0.95')
		for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
		plt.tight_layout()
		if save_img:
			fig.savefig(join(self.dir_name,"fig_shellIntensities.png"),dpi=300)
		plt.show()
		return df_shells


"""
if __name__ == "__main__":
	name_template = '/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img'
"""