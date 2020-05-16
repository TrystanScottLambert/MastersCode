
#################################
# Algorithm testing  function
#MOST RECENT USE THIS TO MAKE THE GENRAL ONE
###############################


import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumtrapz
import datetime    
import time 
import os
from mpl_toolkits.mplot3d import Axes3D
from astropy import units as u
from astropy.coordinates import SkyCoord
import warnings
from tqdm import tqdm
import GM 


#import all the constants that we need for the functions
h,MagLim,vf,red_start,redlim,alpha,M_star,Phi_star=np.loadtxt('Params.txt',usecols=(1))
#Calculated Constants
Phi_star=Phi_star*(h**3)
H0 = 100*h
M_lim = MagLim-25-5*np.log10(vf/H0)
lum_const=0.4*np.log(10)*Phi_star



warnings.filterwarnings('error')  #catch the cos(dec) error that sometimes we get an invoke long astropy method in these instances.


#infile = 'Virgo_Test_data.txt'#VIRGO TEST SET 
infile="2mrs_1175_done.dat"
Ra, Dec, l, b, v, K = np.loadtxt(infile, usecols=(1,2,3,4,24,5),unpack=True)  #readining in numerical data
MASXJIDs,other_names=np.loadtxt(infile,usecols=(0,28),dtype=str,unpack=True) #reading in float data

galaxyid=np.arange(0,len(Ra),1)
cut=np.where((v>red_start)&(v<redlim))
Ra,Dec,l,b,v,K,galaxyid,MASXJIDs,other_names=Ra[cut],Dec[cut],l[cut],b[cut],v[cut],K[cut],galaxyid[cut],MASXJIDs[cut],other_names[cut]

#galaxyid


#infile='CrookGalaxies.dat'
#Ra, Dec, l, b, v, K = np.loadtxt(infile, usecols=(0,1,2,3,4,5),unpack=True)

#l,b,v=np.loadtxt('2MRS2012.dat',usecols=(0,1,2),unpack=True)
#c=SkyCoord(l=l*u.degree,b=b*u.degree,frame='galactic')
#Ra=c.icrs.ra.value
#Dec=c.icrs.dec.value
#K=np.ones(len(Ra))*999

'''h=0.73
MagLim=11.75
vf=1000.
redlim =np.max(v)#12800.0
red_start=300.
prefix_string='2MRS2018_'
lim_rad=1.
galaxyid=np.arange(0,len(Ra),1)
H0 = 100*h
M_lim = MagLim-25-5*np.log10(vf/H0)
group_vlim=1000.  #pm 
MagLim=11.75
M_lim = MagLim-25-5*np.log10(vf/H0)
cut=np.where((v>red_start)&(v<redlim))
Ra,Dec,l,b,v,K=Ra[cut],Dec[cut],l[cut],b[cut],v[cut],K[cut]

#Shecter Parameters
alpha=-1.02
M_star=-24.2
Phi_star=1.08*(10**(-2))*(h**3)
lum_const=0.4*np.log(10)*Phi_star
'''

def WrapMean(array):  
	if (np.max(array)-np.min(array)>=180) and len(np.where((array>90) & (array<270))[0])==0:
		left=[]
		right=[]
		for k in range(len(array)):
			if array[k]<180:
				right.append(array[k])
			else:
				left.append(array[k])
		left_avg=np.mean(left)-360
		right_avg=np.mean(right)
		avg=np.mean([left_avg,right_avg])
		if avg<0:
			avg+=360
	else:
		avg=np.mean(array)
	return avg



def WrapMedian(array):  
	if (np.max(array)-np.min(array)>=180) and len(np.where((array>90) & (array<270))[0])==0:
		left=[]
		right=[]
		for k in range(len(array)):
			if array[k]<180:
				right.append(array[k])
			else:
				left.append(array[k])
		left_avg=np.median(left)-360
		right_avg=np.median(right)
		avg=np.median([left_avg,right_avg])
		if avg<0:
			avg+=360
	else:
		avg=np.median(array)
	return avg



def M12(v_avg):
	return MagLim-25-5*np.log10(v_avg/H0)


def LuminosityFunction(M):  #Used in all the integrals 
	t2=10**(0.4*(alpha+1)*(M_star-M))
	t3=np.exp(-10**(0.4*(M_star-M)))
	return lum_const*t2*t3


def integrate(lower_bound,upper_bound,function,dt):
	x=np.arange(lower_bound,upper_bound,dt)
	y=function(x)                                   
	ysum=np.sum(y)             
	val=ysum*dt
	return val	

'''D_0=0.56#MPc
integral1, err1 = quad(LuminosityFunction,-300,M_lim)
V_0=1000 #km/s, from Ramella et al.'''

integral1, err1 = quad(LuminosityFunction,-300,M_lim)   #another calculated constant
######################################################################################################################################################################################################################
# Diagnositc Plotting #
#######################

'''v_avgs=np.linspace(1,13000,10000)
M12_s=M12(v_avgs)
M12_s_sort=np.sort(M12_s)
arg=M12_s.argsort()
rev=arg.argsort()

plt.plot(v_avgs,M12_s,lw=3)
plt.xlabel('v_avg [km/s]',fontsize=30)
plt.ylabel('M12 [abs mag]',fontsize=30)
plt.axvline(vf,color='k',ls='dashed',lw=3)
plt.text(vf-1,1,'vf=1000',fontsize=12)
plt.gca().invert_yaxis()
plt.show()


LMs=LuminosityFunction(M12_s)
plt.plot(M12_s,LMs,lw=3,c='m')
plt.xlabel('M [abs mag]',fontsize=30)
plt.ylabel('$\Phi$ [Mpc$^{-3}$ abs$^{-1}$]',fontsize=30)
plt.axvline(M12(vf),color='k',lw=3,ls='dashed',c='k')
plt.text(M12(vf)-0.2,0.0003,'Mlim = '+str(round(M12(vf),2)),fontsize=12)
plt.gca().invert_xaxis()
plt.show()


M12_si=np.append(np.array([-32]),M12_s_sort)
intLMs=cumtrapz(LuminosityFunction(M12_si),M12_si)
intLMs=intLMs[rev]


plt.plot(M12_s,intLMs,lw=3,c='r')
plt.xlabel('M12 [abs mag]',fontsize=30)
plt.ylabel('$\int_{-\infty}^{M12} \Phi (M) dM$ [Mpc$^{-1/3}$]',fontsize=30)
plt.gca().invert_xaxis()
plt.axvline(M12(vf),color='k',lw=3,ls='dashed',c='k')
plt.text(M12(vf)-0.2,-0.017,'Mlim = '+str(round(M12(vf),2)),fontsize=12)
plt.show()

dls=D_0*((intLMs/integral1)**(-1./3))


plt.plot(M12_s,dls,lw=3,c='g')
plt.axvline(M12(vf),color='k',lw=3,ls='dashed',c='k')
plt.xlabel('M12 [abs mag]',fontsize=30)
plt.ylabel('Dl [Mpc$^{-1}$]',fontsize=30)
plt.text(M12(vf)-0.2,-1,'Mlim = '+str(round(M12(vf),2)),fontsize=12)
plt.gca().invert_xaxis()
plt.show()

def dlv(v):
	M12_s=M12(v_avgs)
	M12_s_sort=np.sort(M12_s)
	arg=M12_s.argsort()
	rev=arg.argsort()
	M12_si=np.sort(np.append(np.array([-32]),M12_s_sort))
	intLMs=cumtrapz(LuminosityFunction(M12_si),M12_si)
	intLMs=intLMs[rev]
	dls=D_0*((intLMs/integral1)**(-1./3))
	return dls

plt.plot(v_avgs,dlv(v_avgs),lw=3,c='c')
plt.xlabel('v_avg [km/s]',fontsize=30)
plt.ylabel('Dl [Mpc$^{-1}$]',fontsize=30)
plt.axvline(vf,c='k',lw=3,ls='dashed')
plt.show()

'''
######################################################################################################################################################################################################################

def angsep(ra1,ra2,dec1,dec2):
	try:
		faq=(np.pi/180)
		ra_1,ra_2,dec_1,dec_2=faq*ra1,faq*ra2,faq*dec1,faq*dec2
		cosl=np.cos((np.pi/2)-dec_2)*np.cos((np.pi/2)-dec_1)+np.sin((np.pi/2)-dec_2)*np.sin((np.pi/2)-dec_1)*np.cos(ra_2-ra_1)
		val=(1./faq)*np.arccos(cosl)
	except RuntimeWarning:
		c1=SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
		c2=SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
		sep=c1.separation(c2)
		val=sep.value

	return val


def FindFriends(galaxy_index,v,ra,dec,v0,d0): #add checked to the end if need be
	
	v_cut=np.where((v>=v[galaxy_index]-v0) & (v<=v[galaxy_index]+v0))[0]

	ra_friends=ra[v_cut]
	dec_friends=dec[v_cut]
	v_friends=v[v_cut]

	separations=angsep(ra_friends,ra[galaxy_index],dec_friends,dec[galaxy_index])
	theta=(np.pi/180)*(separations/2)
	v_average=(v[galaxy_index]+v_friends)/2
	D12=(np.sin(theta)*(v_average/H0))

	M12s=M12(v_average)
	M12sort=np.sort(M12s)
	arg=M12s.argsort()
	rev=arg.argsort()

	M12sort=np.append(np.array([-32]),M12sort) #add first step
	yts=LuminosityFunction(M12sort)
	integrals=cumtrapz(yts,M12sort)
	integrals=integrals[rev] #return to correct order

	D_L=d0*((integrals/integral1)**(-1./3))
	
	checklimit=np.where(D_L>projected_limit)[0]				#making sure that the linking length can never be larger than the projected limits (at high cz)
	if len(checklimit)<>0:
		D_L[checklimit]=projected_limit-0.1

	pos_distances=D_L-D12
	pos=np.where(pos_distances>=0)[0]   ## MADE AN EDIT HERE DON't KNOW IF THIS SHOULD ACTUALLY BE GREATER THAN AND EQUAL TOO###
	return v_cut[pos]



def FindGroupQuick(galaxy_index,v,ra,dec,v0,d0,group_vlim,lim_rad):
	friends_after=FindFriends(galaxy_index,v,ra,dec,v0,d0)
	friends_before=[]
	new=np.setdiff1d(friends_after,np.array([galaxy_index]))
	iterations=0	
	while list(friends_after) <> list(friends_before) and iterations<20:
		iterations+=1
		friends_before=friends_after

		#update the group center
		group_ra=WrapMean(ra[friends_before])
		group_dec=np.mean(dec[friends_before])
		group_v=np.mean(v[friends_before])
		#print 'groupv= ', group_v
		#print 'vs: ', v[friends_before]
		#print 'vlim', group_vlim
		V_cut=np.where((np.abs(v[friends_before]-group_v)<=group_vlim))[0]  #returns which array indicies of the current group satisy this condition
		RA_friends=ra[friends_before][V_cut]      #excludes all the galaxies which fall out of the velocity limits 
		DEC_friends=dec[friends_before][V_cut]
		V_friends=v[friends_before][V_cut]

		Separations=angsep(RA_friends,group_ra,DEC_friends,group_dec)
		Theta=(np.pi/180)*(Separations)
		V_average=(group_v+V_friends)/2      #keep the average here since the center is bound to change. Using the center only would be less accurate. 
		DD12=(np.sin(Theta)*(V_average/H0))  #projected distances for each galaxy to the center of the group

		Pos=np.where(DD12<=lim_rad)[0]
		friends_before=friends_before[V_cut][Pos]
		new=np.intersect1d(friends_before,new)
		variable=[]
		for i in range(len(new)):
			variable+=list(FindFriends(new[i],v,ra,dec,v0,d0))
		friends_after=np.unique(np.append(friends_before,variable)).astype(int)
		new=np.setdiff1d(friends_after,friends_before)

		print '\t', len(friends_before), len(new), len(friends_after)
	return friends_after

#v0=350.
#d0=0.56
#vlim=1000.
#dlim=1.


def run_fof(Ra,Dec,v,v0,d0,vlim,dlim):

	checked=np.zeros(len(Ra))
	Groups=[]
	tic=datetime.datetime.now() 
	local_cut=checked
	local_cut=np.where(checked==0)[0]    #new added
	while len(local_cut)>1:    #<>1 or >1?
		local_cut=np.where(checked==0)[0]
		local_p=np.random.randint(0,len(local_cut))
		print len(local_cut)
		local_ra=Ra[local_cut]
		local_dec=Dec[local_cut]
		local_v=v[local_cut]
		local_group=FindGroupQuick(int(local_p),local_v,local_ra,local_dec,v0,d0,vlim,dlim)
		group=local_cut[local_group]
		Groups.append(group)       
		checked[group]=1
		#local_cut=np.where(checked==0)[0]
	toc=datetime.datetime.now()
	print toc-tic

	groups=[]
	notgroups=[]
	for i in range(len(Groups)):
		if len(Groups[i])>2:
			groups.append(Groups[i])
		else:
			notgroups.append(Groups[i])

	if len(notgroups)>1:    #need at least one array to concatenate or else it wont work
		notgroups=np.concatenate(notgroups)

	return groups,notgroups



runs=100
velocity_limit=1500.#2000.#1500.
global projected_limit
projected_limit=2.01
d0_i,d0_f=0.46,0.66#,0.56
dd0=d0_f-d0_i
it=dd0/runs
f=open('run_results.txt','w')
f.close()
start_time=datetime.datetime.now()
for i in range(runs):
	f=open('run_results.txt','a')
	print '\t\t RUN NUMBER ', i 

	v0=350.#250.+i*20.
	d0=d0_i+it*i
	print '\t\t d0 = ', d0
	t=run_fof(Ra,Dec,v,v0,d0,velocity_limit,projected_limit)
	#run_results+=t[0]
	for trial in range(len(t[0])):
		for incident in range(len(t[0][trial])):
			f.write(str(t[0][trial][incident])+' ')
		f.write('\n')
	f.close()

#Reading in the results of the numerous runs and putting them back into a list of arrays to be used for the graph theory 
f=open('run_results.txt') #open the text file 
lines=f.readlines()   #make a list of each row  (which is a string)
run_results=[]
for i in range(len(lines)):
	run_results.append(np.array([int(x) for x in lines[i].split()]))  #make the array of our list which splits each string into a list of integers (galaxy labels)


#done reading in 
####

groups,edges,weighting,weighting_normed,sub_groupings=GM.Stabalize(run_results,0.5,runs)
notgroups=np.setdiff1d(np.arange(0,len(Ra)),np.concatenate(groups))
end_time=datetime.datetime.now()

print 'Program Time:', end_time-start_time
print
print 'Writing to File'
prefix_string='2MRS'


#Write edge metadata to a file quickly (May want to move this into a proper place )
edge_1,edge_2,weight=np.zeros(len(edges)),np.zeros(len(edges)),np.zeros(len(edges)) #split the edges tuples triplates into 3 different arrays
for i in range(len(edges)):
	edge_1[i]=edges[i][0]
	edge_2[i]=edges[i][1]
	weight[i]=edges[i][2]






'''non_Ra=Ra[t[1]]
non_Dec=Dec[t[1]]

g_Ra=Ra[np.concatenate(t[0])]
g_Dec=Dec[np.concatenate(t[0])]

c=SkyCoord(ra=non_Ra*u.degree,dec=non_Dec*u.degree,frame='icrs')
c1=SkyCoord(ra=g_Ra*u.degree,dec=g_Dec*u.degree,frame='icrs')

nl=c.galactic.l.value
nb=c.galactic.b.value
gl=c1.galactic.l.value
gb=c1.galactic.b.value

rl=np.random.random(len(nl))*360
a=np.random.random(len(nl))*2-1
rb=np.arccos(a)-np.pi/2
rb=rb*(180./np.pi)
print min(rb),max(rb)

plt.plot(rl,rb,'bo')
plt.show()

fig=plt.figure()
ax=fig.add_subplot(111,projection="aitoff")
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.scatter(-(np.pi/180)*nl,(np.pi/180)*nb,s=1)
ax.scatter(-(np.pi/180)*nl+2*np.pi,(np.pi/180)*nb,s=1)
plt.show()

fig=plt.figure()
ax=fig.add_subplot(111,projection="aitoff")
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.scatter(-(np.pi/180)*rl,(np.pi/180)*rb,s=1)
ax.scatter(-(np.pi/180)*rl+2*np.pi,(np.pi/180)*rb,s=1)
plt.show()

fig=plt.figure()
ax=fig.add_subplot(111,projection="aitoff")
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.scatter(-(np.pi/180)*gl,(np.pi/180)*gb,s=1)
ax.scatter(-(np.pi/180)*gl+2*np.pi,(np.pi/180)*gb,s=1)
plt.show()

fig=plt.figure()
ax=fig.add_subplot(111,projection="aitoff")
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.scatter(-(np.pi/180)*nl,(np.pi/180)*nb,s=1)
ax.scatter(-(np.pi/180)*nl+2*np.pi,(np.pi/180)*nb,s=1)
ax.scatter(-(np.pi/180)*gl,(np.pi/180)*gb,s=1)
ax.scatter(-(np.pi/180)*gl+2*np.pi,(np.pi/180)*gb,s=1)
plt.show()'''


'''singles_ra=Ra[t[1]]
singles_dec=Dec[t[1]]

thing=np.concatenate(t[0])
groups_ra=Ra[thing]
groups_dec=Dec[thing]

uni_ra=np.random.random(len(Ra))*360
a=np.random.random(len(nl))*2-1
uni_dec=np.arccos(a)-np.pi/2
uni_dec=uni_dec*(180./np.pi)

singles_b=get_lb(singles_ra,singles_dec)[1]
groups_b=get_lb(groups_ra,groups_dec)[1]
uni_b=get_lb(uni_ra,uni_dec)[1]

north=np.where(b>0)[0]
south=np.where(b<0)[0]
Ra_north=Ra[north]
Ra_south=Ra[south]
Dec_north=Dec[north]
Dec_south=Dec[south]
DD_north=Generate_2p(Ra_north,Dec_north,100,30)
DD_south=Generate_2p(Ra_south,Dec_south,100,30)
DD=Generate_2p(Ra,Dec,100,30)

north=np.where(singles_b>0)[0]
south=np.where(singles_b<0)[0]
singles_ra_north=singles_ra[north]
singles_ra_south=singles_ra[south]
singles_dec_north=singles_dec[north]
singles_dec_south=singles_dec[south]
DDs_south=Generate_2p(singles_ra_south,singles_dec_south,100,30)
DDs_north=Generate_2p(singles_ra_north,singles_dec_north,100,30)
DDs=Generate_2p(singles_ra,singles_dec,100,30)

north=np.where(groups_b>0)[0]
south=np.where(groups_b<0)[0]
groups_ra_north=groups_ra[north]
groups_ra_south=groups_ra[south]
groups_dec_north=groups_dec[north]
groups_dec_south=groups_dec[south]
DDg_south=Generate_2p(groups_ra_south,groups_dec_south,100,30)
DDg_north=Generate_2p(groups_ra_north,groups_dec_north,100,30)
DDg=Generate_2p(groups_ra,groups_dec,100,30)

north=np.where(uni_b>0)[0]
south=np.where(uni_b<0)[0]
uni_ra_north=uni_ra[north]
uni_ra_south=uni_ra[south]
uni_dec_north=uni_dec[north]
uni_dec_south=uni_dec[south]
DDu_south=Generate_2p(uni_ra_south,uni_dec_south,100,30)
DDu_north=Generate_2p(uni_ra_north,uni_dec_north,100,30)
DDu=Generate_2p(uni_ra,uni_dec,100,30)


x=np.linspace(1,100,1000)
s=10
plt.plot(DDs[0],DDs[1],'-bo',lw=3,ms=s,label='Left Overs')
plt.plot(DDg[0],DDg[1],'-rs',lw=3,ms=s,label='Structure')
plt.plot(DDu[0],DDu[1],'-g^',lw=3,ms=s,label='Uniform')
plt.plot(DD[0],DD[1],'-mv',lw=3,ms=s,label='2MRS')
plt.plot(x,x**(-0.8),'--k',lw=2,label='Toms Line')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta$ [deg]',fontsize=30)
plt.ylabel(r'$\omega (\theta)$',fontsize=30)
plt.show()


x=np.linspace(1,100,1000)
s=10
plt.plot(DDs_south[0],DDs_south[1],'-bo',lw=3,ms=s,label='Left Overs')
plt.plot(DDg_south[0],DDg_south[1],'-rs',lw=3,ms=s,label='Structure')
plt.plot(DDu_south[0],DDu_south[1],'-g^',lw=3,ms=s,label='Uniform')
plt.plot(DD_south[0],DD_south[1],'-mv',lw=3,ms=s,label='2MRS')
plt.plot(x,x**(-0.8),'--k',lw=2,label='Toms Line')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta$ [deg]',fontsize=30)
plt.ylabel(r'$\omega (\theta)$',fontsize=30)
plt.show()

x=np.linspace(1,100,1000)
s=10
plt.plot(DDs_north[0],DDs_north[1],'-bo',lw=3,ms=s,label='Left Overs')
plt.plot(DDg_north[0],DDg_north[1],'-rs',lw=3,ms=s,label='Structure')
plt.plot(DDu_north[0],DDu_north[1],'-g^',lw=3,ms=s,label='Uniform')
plt.plot(DD_north[0],DD_north[1],'-mv',lw=3,ms=s,label='2MRS')
plt.plot(x,x**(-0.8),'--k',lw=2,label='Toms Line')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta$ [deg]',fontsize=30)
plt.ylabel(r'$\omega (\theta)$',fontsize=30)
plt.show()'''
####################################################################################################################
## At this point the algorithm is done running and we have an array of arrays with the groups and not groups
##  This is formating into the various different formats that we want. 
##  !NEEDS AN ARRAY of ARRAYS WHICH CONTAIN THE INDICIES OF THE GROUPS. 
####################################################################################################################
#cosmology section. Defining these functions to get heliocentric velocities into CMB velocities. 
def VCMB(glon,glat,v):
	Lapex=264.14
	Bapex=48.26
	Vapex=371.0

	dL=glon-Lapex
	T1=np.sin(glat*(np.pi/180))*np.sin(Bapex*(np.pi/180))
	T2=np.cos(glat*(np.pi/180))*np.cos(Bapex*(np.pi/180))*np.cos(dL*(np.pi/180))

	return v+(Vapex*(T1+T2))

h=0.73
Om_m=0.3 #Omegas 
Om_e=0.7
Om_k=0.0
Dh=3000*(1./h)
def E(z):
	term1=Om_m*((1+z)**3)
	term2=Om_k*((1+z)**2)
	term3=Om_e
	return np.sqrt(term1+term2+term3)

def integrand(z):
	k=E(z)
	return 1./k

def CoMoving(cz):
	cm=[]
	for i in range(len(cz)):
		integral,err=quad(integrand,0,float(cz[i]/300000.))
		cm.append(Dh*integral)
	return np.array(cm)

def Get_Dist(glon,glat,v_helio):
	cmb_velocity=VCMB(glon,glat,v_helio)
	cmb_redshift=cmb_velocity/300000. 
	cmb_comoving_distance=CoMoving(cmb_velocity)
	return cmb_velocity,cmb_redshift,cmb_comoving_distance

def Convert_to_XYZ(glon,glat,v_helio):
	cmb_v,cmb_z,cmb_d,=Get_Dist(glon,glat,v_helio)
	c=SkyCoord(l=glon*u.degree,b=glat*u.degree,distance=cmb_d*u.Mpc,frame='galactic')
	val_x=c.cartesian.x.value
	val_y=c.cartesian.y.value
	val_z=c.cartesian.z.value
	return cmb_v,cmb_z,cmb_d,val_x,val_y,val_z
#######################################################################################3
# Writing the edges meta data, hopedully using the same cosmology as everything else
########################################################################################

weight=weight.astype(str)

l1=l[edge_1.astype(int)]
l2=l[edge_2.astype(int)]
b1=b[edge_1.astype(int)]
b2=b[edge_2.astype(int)]
vh1=v[edge_1.astype(int)]
vh2=v[edge_2.astype(int)]

vcmb1=VCMB(l1,b1,vh1)
vcmb2=VCMB(l2,b2,vh2)

dist1=CoMoving(vcmb1)
dist2=CoMoving(vcmb2)

c=SkyCoord(l=l1*u.degree,b=b1*u.degree,distance=dist1*u.Mpc,frame='galactic')
x1=c.cartesian.x.value.astype(str) 
y1=c.cartesian.y.value.astype(str)
z1=c.cartesian.z.value.astype(str)

c=SkyCoord(l=l2*u.degree,b=b2*u.degree,distance=dist2*u.Mpc,frame='galactic')
x2=c.cartesian.x.value.astype(str) 
y2=c.cartesian.y.value.astype(str)
z2=c.cartesian.z.value.astype(str)


#write as an ascii format
labels=['X1','Y1','Z1','X2','Y2','Z2','Weight']
types=['r','r','r','r','r','r','r']

data=[x1,y1,z1,x2,y2,z2,weight]
widths=np.zeros(len(data))
for i in range(len(data)):
	max_val_data=len(max(data[i],key=len))+1
	max_val_label=len(labels[i])+1
	if max_val_label>=max_val_data:
		widths[i]=max_val_label
	else:
		widths[i]=max_val_data
widths=widths.astype(int)		

f=open('edges.txt','w')
file_string='|'.join('%*s' % p for p in zip(widths,labels))
f.write('|'+file_string+'| \n')
file_string='|'.join('%*s' % p for p in zip(widths,types))
f.write('|'+file_string+'| \n')
for i in range(len(x1)):
	headers=[x1[i],y1[i],z1[i],x2[i],y2[i],z2[i],weight[i]]
	file_string=' '.join('%*s' % p for p in zip(widths,headers))
	f.write(' '+file_string+' \n')
f.close()

#write for partiview.

levels=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]
color_vals=[2,5,8,11,13,16,18,21,23,25]
for i in range(len(levels)-1):
	weight=weight.astype(float)
	cut_lower=np.where(weight>levels[i])
	cut_upper=np.where(weight<=levels[i+1])
	bins=np.intersect1d(cut_lower,cut_upper)
	f=open('Edges_'+str(levels[i+1])+'.speck','w')
	f.write('#Edges from ('+str(levels[i])+','+str(levels[i+1])+'] \n')
	f.write(' \n ')
	print 'Writing Partiview File: '+str(i+1)+'/10'
	for k in tqdm(range(len(x1[bins]))):
		f.write('mesh -c '+ str(color_vals[i])+' -s wire { \n')
		f.write('1 2 \n')
		f.write(x1[bins][k]+' '+y1[bins][k]+' '+z1[bins][k]+' \n')
		f.write(x2[bins][k]+' '+y2[bins][k]+' '+z2[bins][k]+' \n')
		f.write('} \n')
		f.write(' \n')
	f.close()

#######################################################
# writing the not in groups to file since easiest.
#######################################################
ng_vcmb=VCMB(l[notgroups],b[notgroups],v[notgroups])
ng_zcmb=ng_vcmb/300000
ng_dist=CoMoving(ng_vcmb)
c=SkyCoord(l=l[notgroups]*u.degree,b=b[notgroups]*u.degree,distance=ng_dist*u.Mpc,frame='galactic')
nx=(c.cartesian.x.value)
ny=(c.cartesian.y.value)
nz=(c.cartesian.z.value)

n_K=K[notgroups]-5.*np.log10(ng_dist*1e6)+5
ng_lum_scale=-1.*(K[notgroups]-5.*np.log10(ng_dist*1e6)+5+20.) #Adding 20 to the current magnitudes and then inverting them to make (most of them) positive
ng_ls_floor=-3. #setting a floor level for the lum_scale (important for partiview to keep some dynamic range)
floor=np.where(ng_lum_scale<ng_ls_floor)[0]   #finding where all the lumscale values are less than the floor 
ng_lum_scale[floor]=ng_ls_floor     		  #assigning those values less than the floor value to the floor value

data=[MASXJIDs[notgroups],Ra[notgroups],Dec[notgroups],l[notgroups],b[notgroups],v[notgroups],ng_vcmb,ng_zcmb,K[notgroups],n_K,ng_lum_scale,ng_dist,nx,ny,nz,other_names[notgroups]]
for i in range(len(data)):  #convert all the arrays into string arrays
	data[i]=data[i].astype(str)

labels=['2MASXJID','Ra','Dec','l','b','vh','vcmb','Zcmb','k_app','K_abs','lum_scale','dist','X','Y','Z','other_names']
units=['s','r','r','r','r','r','r','r','r','r','r','r','r','r','r','s']


widths=np.zeros(len(data))
for i in range(len(data)):
	max_val_data=len(max(data[i],key=len))+1
	max_val_label=len(labels[i])+1
	if max_val_label>=max_val_data:
		widths[i]=max_val_label
	else:
		widths[i]=max_val_data
widths=widths.astype(int)

f=open(prefix_string+'_notGroups.txt','w')
file_string='|'.join('%*s' % p for p in zip(widths,labels))
f.write('|'+file_string+'| \n')
file_string='|'.join('%*s' % p for p in zip(widths,units))
f.write('|'+file_string+'| \n')



for i in range(len(notgroups)):
	headers=[MASXJIDs[notgroups[i]],Ra[notgroups[i]],Dec[notgroups[i]],l[notgroups[i]],b[notgroups[i]],v[notgroups[i]],ng_vcmb[i],ng_zcmb[i],K[notgroups[i]],n_K[i],ng_lum_scale[i],ng_dist[i],nx[i],ny[i],nz[i],other_names[notgroups[i]]]
	file_string=' '.join('%*s' % p for p in zip(widths,headers))
	f.write(' '+file_string+' \n')
f.close()

####################################################################################
# Writing the file for subgroups. This will only be the ones that have splintered  #
# The whole data set will be needed however, to do the redshift distortion correction #
########################################################################################

parent_labels=[]
sub_groups=[]
normal_groups=[]
for i in range(len(sub_groupings)):
	if len(sub_groupings[i])>1:
		for k in range(len(sub_groupings[i])):
			sub_groups.append(sub_groupings[i][k])
			parent_labels.append(i+1)
	else:
		normal_groups.append(sub_groupings[i])

subs_l,subs_b,subs_vh=[],[],[]
for i in range(len(sub_groups)):
	subs_l.append(WrapMean(l[sub_groups[i]]))
	subs_b.append(np.mean(b[sub_groups[i]]))
	subs_vh.append(np.mean(v[sub_groups[i]]))

subs_l,subs_b,subs_vh=np.array(subs_l),np.array(subs_b),np.array(subs_vh)
subs_vcmb,subs_zcmb,subs_d,subs_x,subs_y,subs_z=Convert_to_XYZ(subs_l,subs_b,subs_vh)

##################################################################
# Will do the correction right away, need it for the sub_raddii
##################################################################

subs_radii=[]
subs_corrected_d=[]
subs_galaxies=[]
for i in range(len(sub_groups)):
	group_ra=WrapMean(Ra[sub_groups[i]])
	group_dec=np.mean(Dec[sub_groups[i]])
	group_dist=subs_d[i]

	group_members_ra=Ra[sub_groups[i]]
	group_members_dec=Dec[sub_groups[i]]


	seps=angsep(group_members_ra,group_ra,group_members_dec,group_dec)		
	Theta=(np.pi/180)*(seps)												
	projected_seps=np.sin(Theta)*(group_dist)							   

	disp=np.std(projected_seps)								

	mu,sigma = 0.,disp 										
	corrected_comoving=np.random.normal(mu,sigma,len(seps))			
	corrected_comoving=group_dist+corrected_comoving

	subs_corrected_d+=list(corrected_comoving)
	subs_radii.append(np.median(projected_seps))   ################ Change over here ####################
	subs_galaxies+=list(sub_groups[i])
subs_galaxies=np.array(subs_galaxies)
subs_galaxies_labels=MASXJIDs[subs_galaxies]


sub_labels=[]
counter=0
parent_labels=np.array(parent_labels)
while counter<len(parent_labels):
	val=np.where(parent_labels==parent_labels[counter])[0]
	subval=np.arange(len(val))+1
	for i in range(len(val)):
		for k in range(len(sub_groups[val[i]])):
			sub_labels.append(str(parent_labels[counter])+'-'+str(subval[i])) 
	counter+=len(val)  


f=open('helloworld.txt','w')
f.write('# 2MRS_label   SubGroupID  \n')
for i in range(len(subs_galaxies_labels)):
	f.write(subs_galaxies_labels[i]+' '+sub_labels[i]+' \n')
f.close()




#################################################################
# Calculating required fields for uncorrected galaxies in groups#
#################################################################

#Calculating the v_cmb
gal_groups=np.concatenate(groups) #gives a list (in order) of all the galaxies which are in groups
gal_groups_vcmb=VCMB(l[gal_groups],b[gal_groups],v[gal_groups]) #calculating the vcmb
gal_groups_Zcmb=gal_groups_vcmb/300000   #calculating the redshift in the CMB reference frame (assuming c=500000 km/s)

#Calculate the comoving distance
gal_groups_dist=CoMoving(gal_groups_vcmb) #using the CMB velocities not the heliocentric!

#generate XYZ coordinates ussing the GALACTIC! Reference frame!
c=SkyCoord(l=l[gal_groups]*u.degree,b=b[gal_groups]*u.degree,distance=gal_groups_dist*u.Mpc,frame='galactic')
gal_groups_X=c.cartesian.x.value
gal_groups_Y=c.cartesian.y.value
gal_groups_Z=c.cartesian.z.value


#Calculate the Absolite magnitude and the Luminosity scale. 
#since these are extinction corrected magnitudes we don't have an extinction correction
gal_groups_K=K[gal_groups]-5.*np.log10(gal_groups_dist*1e6)+5
gal_groups_lum_scale=-1.*(K[gal_groups]-5.*np.log10(gal_groups_dist*1e6)+5+20.) #Adding 20 to the current magnitudes and then inverting them to make (most of them) positive
gal_groups_ls_floor=-3. #setting a floor level for the lum_scale (important for partiview to keep some dynamic range)
floor=np.where(gal_groups_lum_scale<gal_groups_ls_floor)[0]   #finding where all the lumscale values are less than the floor 
gal_groups_lum_scale[floor]=gal_groups_ls_floor     		  #assigning those values less than the floor value to the floor value



################################################################################################################################################################################
#Generating the first instance of the groupsgal catalog. This can then be read back in to get all the arrays we need as strings and to work out the group data on the spot. 
################################################################################################################################################################################
counter=0 #counter to make sure that entries found in one array eg. gal_groups, dont get repeated in the nested for loop. 
file=open(prefix_string+'_GroupGalaxies.txt','w')
#file.write(\n')
#file.write('[string] \t [deg] \t [deg] \t [deg] \t [deg] \t [km/s] \t [km/s] \t [ap_mag] \t [abs_mag] \t [Mpc] \t [Mpc] \t [Mpc] \t [Mpc] \t [string] \t [deg] \t [deg] \t [deg] \t [deg] \t [km/s] \t [number] \t [number] \n')
for i in range(len(groups)):    #run through all of the groups
	for k in range(len(groups[i])): #go through every galaxy in every group
		file_2masxj=str(MASXJIDs[gal_groups[counter]])+' '
		file_ra=str(Ra[groups[i][k]])+' '
		file_dec=str(Dec[groups[i][k]])+' '
		file_l=str(l[groups[i][k]])+' '
		file_b=str(b[groups[i][k]])+' '
		file_v=str(v[groups[i][k]])+' '
		file_vcmb=str(gal_groups_vcmb[counter])+' '
		file_zcmb=str(gal_groups_Zcmb[counter])+' '
		file_mag=str(K[groups[i][k]])+' '
		file_abs=str(gal_groups_K[counter])+' '
		file_dist=str(gal_groups_dist[counter])+' '
		file_X=str(gal_groups_X[counter])+' '
		file_Y=str(gal_groups_Y[counter])+' '
		file_Z=str(gal_groups_Z[counter])+' '
		file_othernames=str(other_names[gal_groups[counter]])+' '
		file_group_ra=str(WrapMean(Ra[groups[i]]))+' '
		file_group_dec=str(np.mean(Dec[groups[i]]))+' '
		file_group_l=str(WrapMean(l[groups[i]]))+' '
		file_group_b=str(np.mean(b[groups[i]]))+' '
		file_group_v=str(np.mean(v[groups[i]]))+' '
		file_group_no=str(len(groups[i]))+' '
		file_group_ID=str(i+1)+ ' \n'
		file_string=file_2masxj+file_ra+file_dec+file_l+file_b+file_v+file_vcmb+file_zcmb+file_mag+file_abs+file_dist+file_X+file_Y+file_Z+file_othernames+file_group_ra+file_group_dec+file_group_l+file_group_b+file_group_v+file_group_no+file_group_ID
		file.write(file_string)
		counter+=1  #increasing counter so that the next value in the 1d calculated arrays, eg. gal_groups will move on to the correct next incriment
file.close()

###################################
# Dealing with the Gropup deatails
###################################


#read back in the file and get all the arrays as strings (which we can use to format everything in the correct manner with string widths)
infile=prefix_string+'_GroupGalaxies.txt'
g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gd,gX,gY,gZ,gother,Gra,Gdec,Gl,Gb,Gvh,Gno,GID=np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21),unpack=True,dtype=str)

#this method does have redundancy but thats irrelavant based on how fast the numpy arrays can be manipulated
#calculating the group properties that we need. 
Gvcmb=VCMB(Gl.astype(float),Gb.astype(float),Gvh.astype(float))  #getting the VCMB of the galaxy groups
GZcmb=Gvcmb/300000 #calculating the cmb redshift of the group
Gd=CoMoving(Gvcmb) #calculating the comving distance of the group
c=SkyCoord(l=Gl.astype(float)*u.degree,b=Gb.astype(float)*u.degree,distance=Gd*u.Mpc,frame='galactic') #getting the X,Y,Z's of the groups. 
Gx=c.cartesian.x.value
Gy=c.cartesian.y.value
Gz=c.cartesian.z.value

Gvcmb=Gvcmb.astype(str)
GZcmb=GZcmb.astype(str)
Gd=Gd.astype(str)
Gx=Gx.astype(str)
Gy=Gy.astype(str)
Gz=Gz.astype(str)

gal_groups_lum_scale=gal_groups_lum_scale.astype(str)  #make a string to be read in

gal_groups_wc=np.concatenate(weighting).astype(str)
gal_groups_wcn=np.concatenate(weighting_normed).astype(str)

data=[g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gal_groups_lum_scale,gd,gX,gY,gZ,gal_groups_wc,gal_groups_wcn,gother,Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID] #a list of arrays which will make writing data easier. 
labels=['2MASXJID','Ra','Dec','l','b','vh','vcmb','Zcmb','k_app','K_abs','lum_scale','dist','X','Y','Z','degree','degree_normed','other_names','Group_Ra','Group_Dec','Group_l','Group_b','Group_vh','Group_vcmb','Group_Zcmb','Group_dist','Group_X','Group_Y','Group_Z','no_Group_Members','Group_ID']
units=['s','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','s','r','r','r','r','r','r','r','r','r','r','r','r','r']


widths=np.zeros(len(data))
for i in range(len(data)):
	max_val_data=len(max(data[i],key=len))+1
	max_val_label=len(labels[i])+1
	if max_val_label>=max_val_data:
		widths[i]=max_val_label
	else:
		widths[i]=max_val_data
widths=widths.astype(int)

f=open(prefix_string+'_GroupGalaxies.txt','w')
file_string='|'.join('%*s' % p for p in zip(widths,labels))
f.write('|'+file_string+'| \n')
file_string='|'.join('%*s' % p for p in zip(widths,units))
f.write('|'+file_string+'| \n')



for i in range(len(g2masid)):
	headers=[g2masid[i],gra[i],gdec[i],gl[i],gb[i],gvh[i],gvcmb[i],gzcmb[i],gk[i],gK[i],gal_groups_lum_scale[i],gd[i],gX[i],gY[i],gZ[i],gal_groups_wc[i],gal_groups_wcn[i],gother[i],Gra[i],Gdec[i],Gl[i],Gb[i],Gvh[i],Gvcmb[i],GZcmb[i],Gd[i],Gx[i],Gy[i],Gz[i],Gno[i],GID[i]]
	file_string=' '.join('%*s' % p for p in zip(widths,headers))
	f.write(' '+file_string+' \n')
f.close()

#############################################################################
# Uncorrected complete, will now create the group catalog (just the groups)
#############################################################################

#just need to find the first instance of each new group.
  #making an int so that the ordering can be ignored with string which would go 1, 10, 100, 1000, 11 etc..
idx=np.unique(GID.astype(int),return_index=True)[1]    #Find where all the new integers are in the labels and return their indecies. 
Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID=Gra[idx],Gdec[idx],Gl[idx],Gb[idx],Gvh[idx],Gvcmb[idx],GZcmb[idx],Gd[idx],Gx[idx],Gy[idx],Gz[idx],Gno[idx],GID[idx] #apply to all the fields
Gdata=data[18:]      #Only take the last few fields which belong to the Groups
Glabels=labels[18:]
Gwidths=widths[18:]

Glabels[8]='X'
Glabels[9]='Y'
Glabels[10]='Z'

units=['r','r','r','r','r','r','r','r','r','r','r','r','r']

f=open(prefix_string+'_Groups.txt','w')
file_string='|'.join('%*s' % p for p in zip(Gwidths,Glabels))
f.write('|'+file_string+'| \n')
file_string='|'.join('%*s' % p for p in zip(Gwidths,units))
f.write('|'+file_string+'| \n')


for i in range(len(Gra)):
	headers=[Gra[i],Gdec[i],Gl[i],Gb[i],Gvh[i],Gvcmb[i],GZcmb[i],Gd[i],Gx[i],Gy[i],Gz[i],Gno[i],GID[i]]
	file_string=' '.join('%*s' % p for p in zip(Gwidths,headers))
	f.write(' '+file_string+' \n')
f.close()

##########################################################################################
# Groups and GalaxyGroups done, using both to generate the redshift distortion correction
##########################################################################################

infile=prefix_string+'_GroupGalaxies.txt'  #reading in the group galaxies file
g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gls,gd,gX,gY,gZ,gw,gwn,gother,Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID=np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),unpack=True,dtype=str,skiprows=2)

corrected_d=[]
radii=[]
gX=[]
gY=[]
gZ=[]
for i in range(len(idx)):   #idx still represents the places in the main group galaxies catalog where the new groups are
	group_members=np.where(GID.astype(int)==i+1)[0]   #taking advantage of the labels being i+1. This gives us the places where the current group members are
	group_members_ra=gra[group_members].astype(float)		#get the ra's of all the group members
	group_members_dec=gdec[group_members].astype(float)		#get the dec's of all the group members 

	group_ra=float(Gra[idx[i]])									#get the systematic ra of the group
	group_dec=float(Gdec[idx[i]])									#get the systematic dec of the group
	group_dist=float(Gd[idx[i]]) 								    #get the systematic !CMB! velocity of the group

	seps=angsep(group_members_ra,group_ra,group_members_dec,group_dec)		#calculate the angular separations (on sky) of all the members from the systematic
	Theta=(np.pi/180)*(seps)												#convert angular separations into radians (to use in trig function) and divide by 2 (instead of dealing with it later)
	projected_seps=np.sin(Theta)*(group_dist)							    #prjected separations at the systemic comoving distance

	disp=np.std(projected_seps)								#get the projected dispersion (assuming a gaussian profile) (using this standard deviation as the characteristic when we draw random samples from a gaussian PDF)

	mu,sigma = 0.,disp 										#set up parameters of the PDF
	corrected_comoving=np.random.normal(mu,sigma,len(group_members))			#generate random comoving
	corrected_comoving=group_dist+corrected_comoving

	corrected_d.append(corrected_comoving)
	radii.append(np.median(projected_seps)) ################ Change over here ######################
	gd[group_members]=corrected_comoving

'''for i in range(len(subs_galaxies_labels)):
	idx=np.where(g2masid==str(subs_galaxies_labels[i]))   ############################## This will convert the distances to the sub groups I THINK ######################################s
	gd[idx]=str(subs_corrected_d[i])'''

radii=np.array(radii).astype(str)
c=SkyCoord(l=gl.astype(float)*u.degree,b=gb.astype(float)*u.degree,distance=gd.astype(float)*u.Mpc,frame='galactic')
gX=c.cartesian.x.value
gY=c.cartesian.y.value
gZ=c.cartesian.z.value


gd=gd.astype(str)  #convert arrays back to strings
gX=gX.astype(str)
gY=gY.astype(str)
gZ=gZ.astype(str)

data=[g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gls,gd,gX,gY,gZ,gw,gwn,gother,Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID] #a list of arrays which will make writing data easier. 
labels=['2MASXJID','Ra','Dec','l','b','vh','vcmb','Zcmb','k_app','K_abs','lum_scale','dist','X','Y','Z','degree','degree_normed','other_names','Group_Ra','Group_Dec','Group_l','Group_b','Group_vh','Group_vcmb','Group_Zcmb','Group_dist','Group_X','Group_Y','Group_Z','no_Group_Members','Group_ID']
units=['s','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','s','r','r','r','r','r','r','r','r','r','r','r','r','r']


widths=np.zeros(len(data))
for i in range(len(data)):
	max_val_data=len(max(data[i],key=len))+1
	max_val_label=len(labels[i])+1
	if max_val_label>=max_val_data:
		widths[i]=max_val_label
	else:
		widths[i]=max_val_data
widths=widths.astype(int)

f=open(prefix_string+'_GroupGalaxies_Corrected.txt','w')
file_string='|'.join('%*s' % p for p in zip(widths,labels))
f.write('|'+file_string+'| \n')
file_string='|'.join('%*s' % p for p in zip(widths,units))
f.write('|'+file_string+'| \n')

for i in range(len(g2masid)):
	headers=[g2masid[i],gra[i],gdec[i],gl[i],gb[i],gvh[i],gvcmb[i],gzcmb[i],gk[i],gK[i],gls[i],gd[i],gX[i],gY[i],gZ[i],gw[i],gwn[i],gother[i],Gra[i],Gdec[i],Gl[i],Gb[i],Gvh[i],Gvcmb[i],GZcmb[i],Gd[i],Gx[i],Gy[i],Gz[i],Gno[i],GID[i]]
	file_string=' '.join('%*s' % p for p in zip(widths,headers))
	f.write(' '+file_string+' \n')
f.close()

#####################################################################
#####################################################################
# Writing the Ascii files to different formats
#####################################################################
#####################################################################


################
# Parti View
###############

#Group Galaxies

infile=prefix_string+'_GroupGalaxies.txt'

g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gls,gd,gX,gY,gZ,gw,gwc,gother,Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID=np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),unpack=True,dtype=str,skiprows=2)
labels=['Ra','Dec','l','b','vh','vcmb','Zcmb','k_app','K_abs','lum_scale','dist','degree','degree_normed','Group_ID']

f=open(prefix_string+'_GroupGalaxies.speck','w')
f.write('datavar 0 texnum \n')
for i in range(len(labels)):
	f.write('datavar '+str(i+1)+' '+labels[i]+' \n')
f.write('texturevar texnum \n')
f.write('texture -M 1 halo.rgb \n')

for i in range(len(gX)):
	file_string=gX[i]+' '+gY[i]+' '+gZ[i]+' 1 '+gra[i]+' '+gdec[i]+' '+gl[i]+' '+gb[i]+' '+gvh[i]+' '+gvcmb[i]+' '+gzcmb[i]+' '+gk[i]+' '+gK[i]+' '+gls[i]+' '+gd[i]+' '+gw[i]+' '+gwn[i]+' '+GID[i]+' # '+gother[i]+' '+g2masid[i]+' \n'
	f.write(file_string)
f.close()

f=open(prefix_string+'_GroupGalaxies_labels.speck','w')
for i in range(len(gX)):
	file_string=' '+gX[i]+' '+gY[i]+' '+gZ[i]+' text '+gother[i]+' \n'
	f.write(file_string)
f.close()

##################################################
#############################
#Corrected group galaxies
##########################

infile=prefix_string+'_GroupGalaxies_Corrected.txt'

g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gls,gd,gX,gY,gZ,gw,gwn,gother,Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID=np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),unpack=True,dtype=str,skiprows=2)
labels=['Ra','Dec','l','b','vh','vcmb','Zcmb','k_app','K_abs','lum_scale','dist','degree','degree_normed','Group_ID']

f=open(prefix_string+'_GroupGalaxies_Corrected.speck','w')
f.write('datavar 0 texnum \n')
for i in range(len(labels)):
	f.write('datavar '+str(i+1)+' '+labels[i]+' \n')
f.write('texturevar texnum \n')
f.write('texture -M 1 halo.rgb \n')

for i in range(len(gX)):
	file_string=gX[i]+' '+gY[i]+' '+gZ[i]+' 1 '+gra[i]+' '+gdec[i]+' '+gl[i]+' '+gb[i]+' '+gvh[i]+' '+gvcmb[i]+' '+gzcmb[i]+' '+gk[i]+' '+gK[i]+' '+gls[i]+' '+gd[i]+' '+gw[i]+' '+gwn[i]+' '+GID[i]+' # '+gother[i]+' '+g2masid[i]+' \n'
	f.write(file_string)
f.close()

f=open(prefix_string+'_GroupGalaxies_Corrected_labels.speck','w')
for i in range(len(gX)):
	file_string=' '+gX[i]+' '+gY[i]+' '+gZ[i]+' text '+gother[i]+' \n'
	f.write(file_string)
f.close()
#####################################################
######################
# Not Groups
######################
#### Note that if there is a single entry in not groups.txt, we will have a problem reading it in. 
infile=prefix_string+'_notGroups.txt'
g2masid,gra,gdec,gl,gb,gvh,gvcmb,gzcmb,gk,gK,gls,gd,gX,gY,gZ,gother=np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),unpack=True,dtype=str,skiprows=2)
labels=['Ra','Dec','l','b','vh','vcmb','Zcmb','k_app','K_abs','lum_scale','dist']


f=open(prefix_string+'_notGroups.speck','w')
f.write('datavar 0 texnum \n')
for i in range(len(labels)):
	f.write('datavar '+str(i+1)+' '+labels[i]+' \n')
f.write('texturevar texnum \n')
f.write('texture -M 1 halo.rgb \n')

for i in range(len(gX)):
	file_string=gX[i]+' '+gY[i]+' '+gZ[i]+' 1 '+gra[i]+' '+gdec[i]+' '+gl[i]+' '+gb[i]+' '+gvh[i]+' '+gvcmb[i]+' '+gzcmb[i]+' '+gk[i]+' '+gK[i]+' '+gls[i]+' '+gd[i]+' # '+gother[i]+' '+g2masid[i]+' \n'
	f.write(file_string)
f.close()

f=open(prefix_string+'_notGroups_labels.speck','w')
for i in range(len(gX)):
	file_string=' '+gX[i]+' '+gY[i]+' '+gZ[i]+' text '+ gother[i]+' \n'
	f.write(file_string)
f.close()
##############################################################
#################
# Groups
#################

infile=prefix_string+'_Groups.txt'
Gra,Gdec,Gl,Gb,Gvh,Gvcmb,GZcmb,Gd,Gx,Gy,Gz,Gno,GID=np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12),unpack=True,dtype=str,skiprows=2)

#have to make the fancy mesh pv file. 
tom_radii=np.sqrt(0.96+(Gno.astype(float)/50.))  #scale radius to some value set by Tom
tom_radii=tom_radii.astype(str)

f=open(prefix_string+'_Groups.speck','w')
for i in range(len(Gx)):
	file_string=Gx[i]+' '+Gy[i]+' '+Gz[i]+' ellipsoid -r '+radii[i]+' -c 10 -s wire -n 24 \n' 
	f.write(file_string)
f.close()

f=open(prefix_string+'_Groups_labels.speck','w')
for i in range(len(Gx)):
	file_string=' '+Gx[i]+' '+Gy[i]+' '+Gz[i]+' text '+GID[i]+' \n'
	f.write(file_string)
f.close()
#############################################################

#################
# Sub Groups
#################
#have to make the fancy mesh pv file. 
tom_radii=np.sqrt(0.96+(Gno.astype(float)/50.))  #scale radius to some value set by Tom
tom_radii=tom_radii.astype(str)

f=open(prefix_string+'_SubGroups.speck','w')
for i in range(len(subs_x)):
	file_string=str(subs_x[i])+' '+str(subs_y[i])+' '+str(subs_z[i])+' ellipsoid -r '+str(subs_radii[i])+' -c 20 -s wire -n 24 \n' 
	f.write(file_string)
f.close()


#############################################################


data=[np.array(subs_x).astype(str),np.array(subs_y).astype(str),np.array(subs_z).astype(str),np.array(subs_radii).astype(str)] #a list of arrays which will make writing data easier. 
labels=['X','Y','Z','Radius']
units=['r','r','r','r']

widths=np.zeros(len(data))
for i in range(len(data)):
	max_val_data=len(max(data[i],key=len))+1
	max_val_label=len(labels[i])+1
	if max_val_label>=max_val_data:
		widths[i]=max_val_label
	else:
		widths[i]=max_val_data
widths=widths.astype(int)

f=open(prefix_string+'_SubGroups.txt','w')
file_string='|'.join('%*s' % p for p in zip(widths,labels))
f.write('|'+file_string+'| \n')
file_string='|'.join('%*s' % p for p in zip(widths,units))
f.write('|'+file_string+'| \n')

for i in range(len(subs_x)):
	headers=[subs_x[i],subs_y[i],subs_z[i],subs_radii[i]]
	file_string=' '.join('%*s' % p for p in zip(widths,headers))
	f.write(' '+file_string+' \n')
f.close()
##########################################################################################################

