import fabio, os

def img2cbf(fake_img,header_contents,keep_original=True):
	'''
	Converts fake images generated by MLFSOM from img to cbf
	Modifies cbf header according to header_contents file...
	so that it can be properly read by DISTL or DIALS
	If keep_original=False, deletes the original file
	For this to work, need to tweak fabio module a little (look for lines "I CHANGED HERE"):
	/Users/matt/Library/Enthought/Canopy/edm/envs/User/lib/python2.7/site-packages/fabio/cbfimage.py
	'''
	owd = os.getcwd()
	dirname = os.path.dirname(os.path.abspath(fake_img))
	basename = os.path.basename(fake_img)
	os.chdir(dirname) # change wd for now, change back to owd at the end

	# Covert from img to cbf ***************************************
	cbf_temp = basename[0:-3]+'cbf'
	bash_2cbf = 'echo %s %s | 2cbf' %(basename,cbf_temp)
	os.system(bash_2cbf)

	# Modify header of cbf, otherwise DISTL won't work *************
	im = fabio.open(cbf_temp)
	with open(header_contents) as myfile:
	    content_lines = myfile.readlines()
	content_lines = [x.strip() for x in content_lines]
	contents = '\r\n'.join(content_lines)
	im.header['_array_data.header_contents'] = contents
	im.update_header()
	cbf_out = cbf_temp  # overwrite temp file
	im.save(cbf_out)

	if not keep_original:
		os.remove(basename)

	os.chdir(owd)
