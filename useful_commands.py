import pandas as pd

def ReadXYI(infile):
	df = pd.read_csv(infile,sep=' ',header=None, \
		names=['ver', 'hor', 'I', 'psi', 'rsize', 'tsize', 'h', 'k', 'l'])
	return df

def ReadXY(infile):
	df = pd.read_csv(infile,sep=' ',header=None, \
		names=['ver', 'hor', 'h', 'k', 'l', 'coeff', 'frac', 'phi', 'psi', 'rsize', 'tsize'])
	return df	

def SliceDF(df,hor0,ver0,hor_width=10,ver_width=10):
	return df[df.ver.between(ver0-ver_width,ver0+ver_width) & df.hor.between(hor0-hor_width,hor0+hor_width)]

# row,col counts of a data frame
df.shape

# isnan function for a data frame or series
df.isna()

# total number of nan's in a series
df['my_column'].isna().sum()

img_list = [x for x in os.listdir('.') if x.endswith('.img')]
var_loc = 6
for img in img_list:
	new_name = '_'.join(img.split('_')[0:var_loc] + [('00'+str(int(img.split('_')[var_loc])+1))[-3::]] + img.split('_')[var_loc+1::])
    os.rename(img,new_name)


#df = pd.read_csv('mlfsom_tempfile0_preds_1.XYI',sep=' ',header=None)
#df.rename(columns={6:'h',7:'k',8:'l',0:'hor',1:'ver'},inplace=True)
#sub = df[df.ver.between(1567,1587) & df.hor.between(955,975)]

# XYI
# format: XDET YDET photons psi Rsize Tsize  h k l
# print X,Y,darwin*I[h,k,l]/lp*frac,psi,Rsize,Tsize,h,k,l

# XY
# format: XDET YDET H K L Lorentz*polar*subspots/yield frac_recorded PHI psi Rsize Tsize