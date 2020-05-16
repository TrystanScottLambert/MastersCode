#########################################################################
#
# Convert 2MRS GroupGalaxies to a csv for table2 for paper 1
# Trystan Lambert
# 
#########################################################################

import numpy as np 
import IAPUC
infile='2MRS_GroupGalaxies.txt'

ID_2mass,ra,dec,l,b,vcmb,GID = np.loadtxt(infile,usecols=(0,1,2,3,4,6,-1),skiprows=2,dtype=str,unpack=True)
order=ra.astype(float).argsort()
ID_2mass,ra,dec,l,b,vcmb,GID = ID_2mass[order],ra[order],dec[order],l[order],b[order],vcmb[order],GID[order]

vcmb = np.round(vcmb.astype(float)).astype(int)
vcmb = vcmb.astype(str)


infile='helloworld.txt'
sub_ID_2mass,sub_ID = np.loadtxt(infile,dtype=str,unpack=True)

### have to get rid of the subs which are no longer subs. That is to say that there parent ID should no longer exist in the 2MRS catalog. ###

#firs split sub_id intp parent and child

parent=[]
child=[]
for i in range(len(sub_ID)):
	val=sub_ID[i].split('-')
	parent.append(val[0])
parent=np.array(parent)

f=open('2MRS_Galaxies_Table2.csv','w')
f.write('2MASS ID,RA J2000,Dec J2000,$l$,$b$,$V_{cmb}$,Group ID,Sub ID \n')
f.write(' , ,[deg],[deg],[deg],[deg],[km s$^{-1}$], \n')
for i in range(len(ID_2mass)):#np.random.randint(0,len(ID_2mass),30):# range(len(ID_2mass)):
	f.write('$'+ID_2mass[i]+'$'+','+ra[i]+','+dec[i]+','+l[i]+','+b[i]+','+vcmb[i]+',')
	val=np.where(sub_ID_2mass==ID_2mass[i])[0]
	if len(val)== 0:
		f.write(GID[i]+','+'-'+' \n')
	elif len(np.where(GID==parent[val])[0])==0:
		f.write(GID[i]+','+'-'+' \n')
	else:
		f.write(GID[i]+','+sub_ID[val][0]+' \n')
f.close()


######################################################################
#### Writing the IAPUC table for submission to MNRAS ###############

infile ='2MRS_Galaxies_Table2.csv'
MNRAS_2MASS_ID,MNRAS_RA,MNRAS_Dec,MNRAS_l,MNRAS_b,MNRAS_V_cmb,MNRAS_Group_ID,MNRAS_Sub_ID = np.loadtxt(infile,skiprows=2,unpack=True,delimiter=',',dtype=str)

#removing the $$ from the 2MASS ID 

neat_ids=[]
for i in range(len(MNRAS_2MASS_ID)):
	neat_ids.append(MNRAS_2MASS_ID[i].split('$')[1])
neat_ids=np.array(neat_ids)
MNRAS_2MASS_ID = neat_ids

data = [MNRAS_2MASS_ID,MNRAS_RA,MNRAS_Dec,MNRAS_l,MNRAS_b,MNRAS_V_cmb,MNRAS_Group_ID,MNRAS_Sub_ID]
units = ['string','deg','deg','deg','deg','km/s','string','string']
labels=['2MASXJ ID','RA','Dec','l','b','vcmb','Group ID','Sub Group ID']


IAPUC.Generate_IAPUC('MNRAS_Table_3.txt',data,labels,units)


#####################################################################
#####################################################################

#########################################################

# Creating the subgroup csv file. #########################

good = np.intersect1d(parent,GID)
good_sub_ids=[]
for i in range(len(good)):
	val=np.where(parent==good[i])[0]
	for k in range(len(val)):
		good_sub_ids.append(sub_ID[val][k])
good_sub_ids=np.array(good_sub_ids)

good_subs=np.unique(good_sub_ids)



f=open('SubGroup_Table3.csv','w')
f.write('Sub ID, Host Group ID, RA J2000, Dec J2000, $l$, $b$, $V_cmb$, Nmem \n')
f.write('text,text,[deg],[deg],deg],[deg],km s$^{-1}$,number \n')
for i in range(len(good_subs)):
	val=np.where(sub_ID==good_subs[i])[0]
	local_2mass_ids=sub_ID_2mass[val]
	
	pos=[]
	for k in range(len(local_2mass_ids)):
		pos.append(np.where(ID_2mass==local_2mass_ids[k])[0][0])
	pos=np.array(pos)

	f.write(good_subs[i]+','+good_subs[i].split('-')[0]+','+str(round(np.mean(ra.astype(float)[pos]),2))+','+str(round(np.mean(dec.astype(float)[pos]),2))+','+str(round(np.mean(l.astype(float)[pos]),2))+','+str(round(np.mean(b.astype(float)[pos]),2))+','+str(int(round(np.mean(vcmb.astype(float)[pos]))))+','+str(len(pos))+' \n')
f.close()


######################################################################
#### Writing the IAPUC table for submission to MNRAS ###############

infile ='SubGroup_Table3.csv'
MNRAS_Sub_ID, MNRAS_Host_Group_ID, MNRAS_RA, MNRAS_Dec, MNRAS_l, MNRAS_b, MNRAS_V_cmb, MNRAS_Nmem = np.loadtxt(infile,skiprows=2,unpack=True,delimiter=',',dtype=str)

#going to sort by RA for consistency
arg_ra = MNRAS_RA.astype(float)
order = arg_ra.argsort()

MNRAS_Sub_ID, MNRAS_Host_Group_ID, MNRAS_RA, MNRAS_Dec, MNRAS_l, MNRAS_b, MNRAS_V_cmb, MNRAS_Nmem = MNRAS_Sub_ID[order], MNRAS_Host_Group_ID[order], MNRAS_RA[order], MNRAS_Dec[order], MNRAS_l[order], MNRAS_b[order], MNRAS_V_cmb[order], MNRAS_Nmem[order]

data=[MNRAS_Sub_ID, MNRAS_Host_Group_ID, MNRAS_RA, MNRAS_Dec, MNRAS_l, MNRAS_b, MNRAS_V_cmb, MNRAS_Nmem]
units=['string','string','deg','deg','deg','deg','km/s','number']
labels=['Sub Group ID','Host Group ID','RA','Dec','l','b','vcmb','Nmem in sub']

IAPUC.Generate_IAPUC('MNRAS_Table_2.txt',data,labels,units)



#####################################################################
#####################################################################
