import os, subprocess
from os.path import join

def RenameImages(loc_idx):
	img_list = [x for x in os.listdir('.') if x.endswith('.img')]
	for img in img_list:
		ori_li = img.split('_')
		mod_li = ori_li[0:loc_idx] + [ori_li[loc_idx].zfill(3)] + ori_li[loc_idx+1::]
		subprocess.call(['mv',img,'_'.join(mod_li)])