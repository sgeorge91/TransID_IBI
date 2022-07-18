""""""""
Author: SVG
This code takes in all the IBI files with the extension .EVT in a directory and calculates the quantifiers of interest (Mean, Stadard deviation, Higuchi dimension and Multiscale entropy) for each file and writes the results into a csv file.
Please download the libraries imported.
""""""""
import numpy as np
import glob
from astropy.time import Time
import hfda
import neurokit2 as nk2
import scipy.stats as sp
import pandas as pd

#fn is the filename, hl is headerlength, fl is footer length
#ranges takes in a file name associated with R peak times, returns the mean, standard deviation, higuchi dimension and multiscale entropy of a block! 
def ranges(fn,hl=15,fl=1):
	al=np.genfromtxt(fn,skip_header=15, skip_footer=1)
	nl=[]
	flag=0
	for i in al:
		if(i[0]==1):
			flag=1
		if(flag==1):
			nl.append(i[1])
		if(i[0]==2):
			break
	nl2=[]
	for i in range (0,len(nl)-1):
		nl2.append(nl[i+1]-nl[i])
	nl3=np.array(nl2)
	return(np.mean(nl2),np.std(nl2),hfda.measure(nl3, 5),nk2.complexity_mse(nl3, dimension=1))

#Takes in the date-time and converts it into a format that astropy can read (iso)
def t2isot(temp,temp2):
	
	r1=str(temp[0:4]+'-'+temp[4:6]+'-'+temp[6:8]+'T'+temp2[0:2]+":"+temp2[2:4]+":"+temp2[4:6])
	
	return(r1)


flist=glob.glob('*.EVT')
flist.sort()
resl=[]
kmax=5
temp=" "
temp2=" "
tarr=[]

for fn in flist:
	temp=fn.split('_')
	temp2=temp[1].split('.')
	tm=t2isot(temp[0],temp2[0])
	tm2 = Time(tm, format='isot', scale='utc')
	m1,d1,hf1,mse1=ranges(fn)
	resl.append([tm2.mjd,m1,d1,hf1,mse1])#float(temp[0])+(float(temp2[0])*.0001)/24
df = pd.DataFrame(resl, columns=['MJD', 'Mean', 'Dev','HigDim','MSE'])
print(sp.kendalltau(df['MJD'],df['MSE']),'MSE')
print(sp.kendalltau(df['MJD'],df['Mean']),'Mean')
print(sp.kendalltau(df['MJD'],df['Dev']),'SD')
print(sp.kendalltau(df['MJD'],df['HigDim']),'HD')
df.to_csv('Result_1141.dat' ,sep='\t', index=False)
	
