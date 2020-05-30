import os, fileinput, subprocess, pickle, fabio
from glob import glob
from os.path import join
import pandas as pd
import numpy as np

#name_template = '/home/thorne20/Desktop/MLFSOM/data_trial_gau/trial_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_trial_gaussian/trial_gaussian_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_trial_gaussian_1H87_2.5A_Ngrid10_fr50/trial_gaussian_2.5A_fr50_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_trial_gaussian_6K8S/trial_gaussian_6K8S_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_gaussian_1H87_2.0A_3deg_grid10_fr20/gaussian_1H87_2.0A_3deg_grid10_fr20_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_gaussian_1H87_2.0A_3deg_grid10_fr20_highflux/gaussian_1H87_2.0A_3deg_grid10_fr20_highflux_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_gaussian_1H87_2.0A_3deg_grid10_fr20_wild/gaussian_1H87_2.0A_3deg_grid10_fr20_wild_???_001.img'
#name_template = '/home/thorne20/Desktop/MLFSOM/data_gaussian_1H87_2.0A_3deg_grid6_fr20/gaussian_1H87_2.0A_3deg_grid6_fr20_???_001.img'
name_template = '/home/thorne20/Desktop/MLFSOM/data_gaussian_1H87_2.5A_3deg_grid20_fr20_xtal400/gaussian_1H87_2.5A_3deg_grid20_fr20_xtal400_???_001.img'

dir_name = os.path.dirname(name_template)
base_name = os.path.basename(name_template)
if glob(join(dir_name,'*desc.txt')) == []:
	print "MyError: Description file is missing in the folder!"
else:
	description = pd.read_csv( glob(join(dir_name,'*desc.txt'))[0],\
		header=None,sep=': ',engine='python', index_col=0, names=['value'] )

weights_list = list(description.loc['(frame,weight,experiment_ID)'::].index)[1::]
df_weights = pd.DataFrame([eval(x) for x in weights_list], columns=['fnumber','weight','ID'])
df_weights.set_index('ID',inplace=True)

for fnumber in range(int(description.loc['frames','value'])):
	frame_IDs = df_weights.index[df_weights.fnumber == fnumber]
	ID0 = frame_IDs[0]
	im0 = fabio.open( name_template.replace('???',str(ID0)) )
	overlaid_data = im0.data.astype('float') * df_weights.weight[ID0]
	sum_weights = df_weights.weight[ID0]
	for ID in frame_IDs[1::]:
		im = fabio.open( name_template.replace('???',str(ID)) )
		overlaid_data += im.data.astype('float') * df_weights.weight[ID]
		sum_weights += df_weights.weight[ID]
	im0.data = (overlaid_data/sum_weights).astype('uint16')
	#im0.data = overlaid_data
	im0.save( name_template.replace('???','aggregate_'+str(fnumber+1).zfill(3)) )