#################################################
#
# Script to Calculate Various group properties
# Mainly holding functions 
#
# Trystan Lambert 
# 09/07/2019
# 
#################################################

import numpy as np 
import pylab as plt
import IAPUC 
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.integrate import quad
from astropy.constants import G

h,MagLim,vf,red_start,redlim,alpha,M_star,Phi_star=np.loadtxt('Params.txt',usecols=(1))
Phi_star=Phi_star*(h**3)
H0 = 100*h
M_lim = MagLim-25-5*np.log10(vf/H0)
lum_const=0.4*np.log(10)*Phi_star


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

def WrapMean(array):
	if np.max(array)-np.min(array)>=180  and len(np.where((array>90) & (array<270))[0])==0:
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
	cmb_redshift=cmb_velocity.astype(float)/300000. 
	cmb_comoving_distance=CoMoving(cmb_velocity)
	return cmb_velocity,cmb_redshift,cmb_comoving_distance

def Convert_to_XYZ(glon,glat,v_helio):
	cmb_v,cmb_z,cmb_d,=Get_Dist(glon,glat,v_helio)
	c=SkyCoord(l=glon*u.degree,b=glat*u.degree,distance=cmb_d*u.Mpc,frame='galactic')
	val_x=c.cartesian.x.value
	val_y=c.cartesian.y.value
	val_z=c.cartesian.z.value
	return cmb_v,cmb_z,cmb_d,val_x,val_y,val_z




infile='2MRS_GroupGalaxies.txt'
GroupRa,GroupDEC,Group_ID,velocities,k_app,Group_no = np.loadtxt(infile,usecols=(1,2,-1,6,8,-2),skiprows=2,unpack=True)
tmrs_id=np.loadtxt(infile,usecols=(0),unpack=True,skiprows=2,dtype=str)

infile='helloworld.txt'
check_tmrs_id,subgroup_id=np.loadtxt(infile,skiprows=2,unpack=True,dtype=str)
parent_ids,children_ids=[],[]
for i in range(len(subgroup_id)):
	val=subgroup_id[i].split('-')
	parent_ids.append(val[0])
	children_ids.append(val[1])
parent_ids=np.array(parent_ids).astype(float)
children_ids=np.array(children_ids).astype(float)
#function for calculating the velocity dispersion using the "gapper" method

'''disp=250
vel=np.random.normal(1000,disp,10)
error=0.1*vel'''


#Function for determining the gap dispersion of a group making use of the velocity error of the mesurement
def group_disp(velocity_array,velocity_error_array):
	velocities=np.sort(velocity_array)
	gaps=velocities[1:]-velocities[:-1]  #calculate the distances between the distances. 
	N=len(velocities)
	i_val=np.arange(N)[1:]
	weights=i_val*(N-i_val)    #calculate gaussian weightings

	gw=weights*gaps
	gw_sum=np.sum(gw)  				#sum of the product of weights and gaps

	sigma_gap=(np.sqrt(np.pi)/(N*(N-1)))*gw_sum 		# calculate sigma gap
	error_quad=np.sqrt(np.sum(velocity_error_array**2)) 			#add errors in quadrature
	sigma=np.sqrt((N/(N-1))*(sigma_gap)**2-(error_quad)**2) #final sigma value
	return sigma
'''
x=group_disp(vel,error)
print 'True: ',disp
print 'gap: ', x 
print 'std: ', np.std(vel) '''

#using tempel et. al 2016 and 2018 to try calculate the dispersion to see if we do a better job.
def tempel_disp(velocity_array):
	mean_z=np.mean(velocity_array/300000)
	mean_v=np.mean(velocity_array)

	term1 = np.sum((velocity_array-mean_v)**2)
	term2=1./((1+mean_z)*(len(velocity_array)-1))

	return np.sqrt(term1*term2)





h=0.73
OL=0.7
O0=0.3

#standard function for calculating the disperson of a galaxy group using Toms paper from 2017
# https://arxiv.org/pdf/1607.01190.pdf
def tom_disp(velocity_array):
	mean_v=np.median(velocity_array)     				############ Change here #####################
	denominator = 1 + velocity_array/300000
	distribution= (velocity_array-mean_v)/denominator

	return np.round(np.sqrt(np.mean(distribution**2)),3)

#function for calulating the R200 radius for the galaxu groups from Poggianti et. al 2010
def R200(velocity_array):
	disp=tom_disp(velocity_array)
	z=np.mean(velocity_array)/300000
	R2=1.73*(float(disp)/1000)*(1./np.sqrt(OL+O0*(1+z)**3))*(h**(-1))
	return np.round(R2,3) 

#function for calulating the Mass for the galaxu groups from Poggianti et. al 2010
def M200(velocity_array):
	disp=tom_disp(velocity_array)
	z=np.mean(velocity_array)/300000
	M2=1.2e15*((disp/1000)**3)*(1./np.sqrt(OL+O0*(1+z)**3))*(h**(-1))  
	return np.round(M2/10e11,3)

#function to calculate our sudo-random definiton of richnes. Finding all the galaxies with R<-23.5 Absolute magnitudes. 
def Richness(k_app_array,v_cmb_array):
	v_cmb_value=np.mean(v_cmb_array)
	co_moving_value=CoMoving([v_cmb_value])
	lum_dist_value=(1+v_cmb_value/300000)*co_moving_value
	
	KABS=k_app_array-5*np.log10(lum_dist_value*1e6)+5
	if lum_dist_value<110:
		return str(len(np.where(KABS<-23.5)[0]))
	else:
		return '>'+str(len(np.where(KABS<-23.5)[0]))#-999


def LuminosityFunction(M):  #Used in all the integrals 
	t2=10**(0.4*(alpha+1)*(M_star-M))
	t3=np.exp(-10**(0.4*(M_star-M)))
	return lum_const*t2*t3

mlim=11.75
Mmin=-24.13
def Ltot(k_app_array,velocity_array):
	mean_vel_array=[np.mean(velocity_array)]
	LD=(1+(np.mean(velocity_array)/300000))*CoMoving(mean_vel_array)

	K=k_app_array-5*np.log10(LD*1e6)+5

	Lis=10**(-0.4*(K-3.28))
	Lobs=np.sum(Lis)

	integral1, err1 = quad(LuminosityFunction,-300,Mmin)
	x=np.linspace(-300,50)
	y=LuminosityFunction(x)

	integral2, err2 = quad(LuminosityFunction,-300,-5*(np.log10(LD))+mlim+5)

	return (Lobs*(integral1/integral2))/1e11


### Tempel 2016 inluded this in his paper. Maybe worth citing
def max_radius(ra_array,dec_array,velocity_array):
	mean_vel_array=[np.mean(velocity_array)]
	LD=(1+(np.mean(velocity_array)/300000))*CoMoving(mean_vel_array)  #get the luminosity distance

	group_ra=WrapMean(ra_array)
	group_dec=np.mean(dec_array)

	thetas=angsep(group_ra,ra_array,group_dec,dec_array)
	thetas=(np.pi/180)*thetas

	thetas=np.sin(thetas)*LD
	return max(thetas)


################ Crook Virial Mass function (requires several minor functions) #####################
#function to calculate the projected distance between two galaxies
def proj_radii(ra_array,dec_array,velocity_array,a,b):
	thetas=angsep(ra_array[a],ra_array[b],dec_array[a],dec_array[b])
	thetas=(np.pi/180)*thetas
	return 2*(np.mean(velocity_array)/H0)*np.tan(thetas/2)

def Virial_Radius_func(ra_array,dec_array,velocity_array):
	N=len(ra_array)
	term1=N*(N-1)
	radiis=[]
	for i in range(len(ra_array)):
		for j in range(len(ra_array)):
			if j>i:
				radiis.append(proj_radii(ra_array,dec_array,velocity_array,i,j))
	radiis=np.array(radiis)
	term2=np.sum(radiis)
	return term1/term2

def projected_dispersion(velocity_array):
	Vg=np.mean(velocity_array)
	N=len(ra_array)
	term1=np.sum((velocity_array-Vg)**2)
	term2=N-1
	return term1/term2

def Virial_properties(ra_array,dec_array,velocity_array):
	Rp=Virial_Radius_func(ra_array,dec_array,velocity_array)
	sigmasquared=projected_dispersion(velocity_array)
	Mv=((3*np.pi)/2)*((Rp*sigmasquared)/(G.value))
	return Mv/1e11, np.sqrt(sigmasquared), Rp

###################################################################################################



all_groups=np.unique(Group_ID)
#f=open('Group_Properties.txt')
labels=['Group_ID','Other Names' ,'Ra','Dec','l','b','Nmem','Nmem_Corrected','Velocity_h','Co_Moving_Distance','Luminosity_Distance','Tom_Dispersion','Standard_Deviation','Virial_Dispersion','Tempel_Dispersion','R200','max_radius','Virial_Radius','M200','Virial_Mass','L_tot','Number_of_Subs','Richness']
#units=['i','r','r','r','r','r','r','r','i']
units=['string','string','deg','deg','deg','deg','number','number','km/s','Mpc','Mpc','km/s','km/s','km/s','km/s','Mpc','Mpc','Mpc','10^11 Msol','Msol?','10^11 L_sol','number','number']
dispersions=np.zeros(len(all_groups))
velocities1=np.zeros(len(all_groups))
masses=np.zeros(len(all_groups))
radii=np.zeros(len(all_groups))
mean_vels=np.zeros(len(all_groups))
numbers=np.zeros(len(all_groups))
RA=np.zeros(len(all_groups))
DEC=np.zeros(len(all_groups))
normal=np.zeros(len(all_groups))
rich=[]#np.zeros(len(all_groups))

total_luminosity = np.zeros(len(all_groups))
max_radii=np.zeros(len(all_groups))
Virial_Mass=np.zeros(len(all_groups))
Virial_dispersion=np.zeros(len(all_groups))
Virial_Radius=np.zeros(len(all_groups))
tempel_dispersion=np.zeros(len(all_groups))

for i in range(len(all_groups)):
	cut=np.where(Group_ID==all_groups[i])[0]
	vel_array=velocities[cut]
	k_app_array=k_app[cut]
	ra_array=GroupRa[cut]
	dec_array=GroupDEC[cut]
	

	l_tot=Ltot(k_app_array,vel_array)
	dispersion_tempel=tempel_disp(vel_array)
	proj_radius=max_radius(ra_array,dec_array,vel_array)
	Mv,sv,Rv = Virial_properties(ra_array,dec_array,vel_array)
	

	disp=tom_disp(vel_array)
	radius=R200(vel_array)
	mass=M200(vel_array)
	number=len(cut)#int(np.unique(Group_no[cut]))
	mean_vel=np.mean(vel_array)

	dispersions[i]=disp
	velocities1[i]=round(mean_vel,2)

	masses[i]=mass
	numbers[i]=len(cut)#int(number)
	radii[i]=radius

	RA[i]=WrapMean(GroupRa[cut])
	DEC[i]=np.mean(GroupDEC[cut])

	normal[i]=np.std(vel_array)
	rich.append(Richness(k_app_array,vel_array))

	total_luminosity[i]=l_tot
	max_radii[i]=proj_radius
	Virial_Mass[i]=Mv
	Virial_dispersion[i]=sv
	Virial_Radius[i]=Rv
	tempel_dispersion[i]=dispersion_tempel

rich=np.array(rich)
order=numbers.argsort()[::-1]   #sort by maxium members, descending
order = RA.argsort()

#Galactic coords plus Luminosity and co_moving distance
c=SkyCoord(ra=RA*u.deg,dec=DEC*u.deg,frame='icrs')
L=np.around(c.galactic.l.value,5)
B=np.around(c.galactic.b.value,5)

Co_Moving=CoMoving(velocities1)
Lum_dist=(1+velocities1/300000)*Co_Moving

Co_Moving=np.around(Co_Moving,1)
Lum_dist=np.around(Lum_dist,1)


#getting the number of subgroups
#######################################################################
number_sub_groups=[]
for i in range(len(all_groups[order])):
	val=np.where(parent_ids==all_groups[order][i])[0]
	if len(val)>0:
		number_sub_groups.append(np.max(children_ids[val]))
	else:
		number_sub_groups.append(0)
number_sub_groups=np.array(number_sub_groups).astype(int)
#######################################################################

########################################################################################################################
#### Gettting the cluster names that we can ############################################################################
########################################################################################################################

infile='big_clusters.tbl'
big_ra,big_dec,big_l,big_b,big_z = np.loadtxt(infile,usecols=(1,2,3,4,5),unpack=True)
big_name=np.loadtxt(infile,usecols=(0),unpack=True,dtype=str)

big_cz=big_z*300000 

big_d=CoMoving(big_cz)

c=SkyCoord(l=big_l*u.degree,b=big_b*u.degree,distance=big_d*u.Mpc,frame='galactic')
big_x=c.cartesian.x.value
big_y=c.cartesian.y.value
big_z=c.cartesian.z.value



c1=SkyCoord(l=L[order]*u.degree,b=B[order]*u.degree,distance=Co_Moving[order]*u.Mpc,frame='galactic')
groups_x=c1.cartesian.x.value 
groups_y=c1.cartesian.y.value 
groups_z=c1.cartesian.z.value 

big_cluster_labels=np.zeros(len(RA))
big_cluster_labels=big_cluster_labels.astype(str)

for i in range(len(big_z)):
	looking_distance = np.sqrt((big_x[i]-groups_x)**2 + (big_y[i]-groups_y)**2 + (big_z[i]-groups_z)**2)
	minimum=np.min(looking_distance)
	if minimum < 10.:
		sorted_pos = np.argsort(looking_distance)
		lowest_five=sorted_pos[:5]
		correct_group = np.where(numbers[order][lowest_five]==np.max(numbers[order][lowest_five]))[0]
		if len(correct_group)>1:
			correct_group=correct_group[0]
		position=lowest_five[correct_group]
		if looking_distance[position] < 10.:
			print big_name[i]
			big_cluster_labels[position] = big_name[i]

neaten = np.where(big_cluster_labels=='0.0')
big_cluster_labels[neaten]= '-'
#########################################################################################################################

#########################################################################################################################
#################################### Adding the corrected Nmem ##########################################################
#########################################################################################################################
infile='groups_Nmem-correction.tbl'
matching_id,n_corrected = np.loadtxt(infile,usecols=(0,2),unpack=True,skiprows=1,dtype=int)

N_corrected=[]
test_matching_id=[]
match_to_order=all_groups[order]
for i in range(len(order)):
	idx = np.where(matching_id==match_to_order[i])[0]
	N_corrected.append(n_corrected[idx][0])
	test_matching_id.append(matching_id[idx][0])
N_corrected = np.array(N_corrected)
#########################################################################################################################


data=[all_groups[order].astype(int),big_cluster_labels,np.around(RA[order],2),np.around(DEC[order],2),np.around(L[order],2),np.around(B[order],2),numbers[order].astype(int),N_corrected,np.around(velocities1[order]).astype(int),Co_Moving[order],Lum_dist[order],np.around(dispersions[order]).astype(int),np.around(normal[order]).astype(int),np.around(Virial_dispersion[order],2),np.around(tempel_dispersion[order],2),np.around(radii[order],2),np.around(max_radii[order],2),np.around(Virial_Radius[order],2),np.around(masses[order],1),np.around(Virial_Mass[order],2),np.around(total_luminosity[order],2),number_sub_groups,rich[order]]
#data=[all_groups.astype(int),np.around(RA,5),np.around(DEC,5),np.around(L,5),np.around(B,5),numbers.astype(int),np.around(velocities1).astype(int),np.around(dispersions).astype(int),np.around(normal).astype(int),np.around(radii,2),np.around(masses,1)]


############ Only using to match what is in the paper for MNRAS online table. Un hash the first data = instance above. ##########################################
#################################################################################################################################################################
data=[all_groups[order].astype(int),big_cluster_labels,np.around(RA[order],2),np.around(DEC[order],2),np.around(L[order],2),np.around(B[order],2),numbers[order].astype(int),N_corrected,rich[order],np.around(velocities1[order]).astype(int),Co_Moving[order],Lum_dist[order],np.around(dispersions[order]).astype(int),np.around(normal[order]).astype(int),np.around(radii[order],2),np.around(masses[order],1),number_sub_groups]
labels=['Group_ID','Other Names' ,'Ra','Dec','l','b','Nmem','Nmem_Corrected','Richness','Velocity_cmb','Co_Moving_Distance','Luminosity_Distance','RMS','Standard_Deviation','R200','M200','Number_of_Subs']
units=['string','string','deg','deg','deg','deg','number','number','number','km/s','Mpc','Mpc','km/s','km/s','Mpc','10^11 Msol','number']
IAPUC.Generate_IAPUC('MNRAS_Table_1.txt',data,labels,units)
#################################################################################################################################################################
#################################################################################################################################################################

IAPUC.Generate_IAPUC('Group_Properties.txt',data,labels,units)


#### Writing A csv file so that the latex table can be generated ####
row_max=20  #number of rows to include use len(data[0]) for the entire set
row_max=len(data[0])

f=open('Group_Properties.csv','w')
for i in range(len(labels)):
	if i < len(labels)-1:
		f.write(str(labels[i])+',')
	else:
		f.write(str(labels[i]))
f.write(' \n')

for i in range(row_max):
	for k in range(len(data)):
		if k < len(data)-1:
			f.write(str(data[k][i])+',')
		else:
			f.write(str(data[k][i]))
	f.write(' \n')
f.close()

###################################################################


plt.style.use('seaborn-paper')

params = {'legend.fontsize': 30,
          'legend.handlelength': 2}
plt.rcParams.update(params)



plot_rich=[]
for i in range(len(rich)):
	if rich[i][0] <> '>' and rich[i][0] <> '0':
		plot_rich.append(i)
plot_rich=np.array(plot_rich)

plt.scatter(radii[plot_rich],np.log10(rich[plot_rich].astype(float)),color='k',s=20)
plt.xlabel(r'$R_{200}$ [Mpc]',fontsize=40)
plt.ylabel(r'$\log({\rm Richness})$',fontsize=40)
plt.tick_params(axis='both',which='major',labelsize=40)
plt.tick_params(axis='both',which='minor',labelsize=40)
#plt.yscale('log')
plt.show()


plt.scatter(radii[order],np.log10(N_corrected.astype(float)),color='k',s=20)
plt.xlabel(r'$R_{200}$ [Mpc]',fontsize=40)
plt.ylabel(r'$\log(N_{\rm corr})$',fontsize=40)
plt.tick_params(axis='both',which='major',labelsize=40)
plt.tick_params(axis='both',which='minor',labelsize=40)
#plt.yscale('log')
plt.show()


plt.scatter(velocities1,numbers.astype(float),color='k',s=20)
plt.xlabel('cz [km s$^{-1}$]',fontsize=40)
plt.ylabel('Number of Members',fontsize=40)
plt.tick_params(axis='both',which='major',labelsize=40)
plt.tick_params(axis='both',which='minor',labelsize=40)
plt.show()


########################## Making Contour Plots for comparison ######################################
'''
from scipy.stats import kde

data=np.array(zip(radii[order],np.log10(N_corrected.astype(float))))
x,y=data.T

k=kde.gaussian_kde(data.T)
nbins=10
xi,yi=np.mgrid[x.min():x.max():nbins*1j,y.min():y.max():nbins*1j]
zi=k(np.vstack([xi.flatten(),yi.flatten()]))

plt.pcolormesh(xi,yi,zi.reshape(xi.shape),shading='gouraud')
plt.show()'''