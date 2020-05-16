## Script to make the 2MRS Groups spheres with R200 radii. 

import numpy as np 

infile='Group_Properties.txt'
GID,R200=np.loadtxt(infile,usecols=(0,11),unpack=True,skiprows=2,dtype=str)

infile='2MRS_Groups.txt'
ids,x,y,z = np.loadtxt(infile,usecols=(-1,8,9,10),skiprows=2,unpack=True,dtype=str)

f=open('2MRS_Groups_R200.speck','w')
for i in range(len(GID)):
	group=np.where(ids==GID[i])[0]
	f.write(x[group][0]+' '+y[group][0]+' '+z[group][0]+' ellipsoid -r '+R200[i]+' -c 10 -s wire -n 24 # '+GID[i]+' \n')
f.close()