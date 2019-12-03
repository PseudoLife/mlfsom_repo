import os, pickle
import pandas as pd
import numpy as np


def mlfsomVSdistl(dataframe_mlfsom,dataframe_distl,fnumber):
	"""
	Compares the peak intensities captured from temp files of MLFSOM with
	the distl peak intensities for a desired frame number
	Returns df with mlfsom and distl int. with fractional discrepancy for each peak
	"""
	# Make a copy of input df's to prevent them to be modified
	df_mlfsom = dataframe_mlfsom.copy() 
	df_distl = dataframe_distl.copy()

	# Get the peaks having non-NaN intensity in the fnumber
	sub_peaks = df_mlfsom[~df_mlfsom[fnumber].isna()][['hor','ver',fnumber]]
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
	print 'Fractional discrepancy between raw and distl int: %.2f' %np.nanmean(sub_peaks.disc)
	return sub_peaks
