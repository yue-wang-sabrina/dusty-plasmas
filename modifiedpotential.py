import numpy
import scipy
import matplotlib.pyplot as plt
import math
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd
import mpl_toolkits.mplot3d.axes3d as p3
import itertools
import time
from scipy import special
from tqdm import tqdm
from scipy.optimize import fsolve
from sklearn import preprocessing
from numba import jit
import pickle
from numpy import isclose
from scipy.interpolate import interp2d
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev


radd=1.5*10**(-6) #radius of dust particle
e=1.60217662*10**(-19) #electron charge
e0=8.85418782*10**(-12) #perimittivity free space
Te=46000 #Electrons usually a few EV #Non-thermal plasma, E field couples with electron more efficiently, higher temperature
Ti=310 #Ion thermal temperature, room temperature
kb=1.38064852*10**(-23) #boltzmans constant
ne0=1.*10**(15) #initial electron density
ni0=ne0 #Quasi-neutral initial conditions
md=1000.*(4/3)*math.pi*radd**3 #Mass of dust particle
mi=39.948*1.66053904*10**(-27) #argon mass ion
me=9.10938356*10**(-31) #electron mass
cs=numpy.sqrt(kb*Te/mi) #Bohm speed = sound speed
lambdade=(kb*Te*e0/(ne0*e**2))**0.5 #Electron debye rad
lambdadi=(kb*Ti*e0/(ni0*e**2))**0.5 #Ion debye rad
lambdaD= 1./(1./(lambdadi**2)+1./(lambdade**2))**0.5  #dusty debye radius
boxr=523*lambdaD #cylinder radius
boxz=0.001 #cylinder height
g=-9.8 #gravitational acceleration
dt=0.0001 #Much longer than charging time which is order nanoseconds
sheathd=10*lambdaD
electrodeV=abs((kb*Te/(2*e))*(numpy.log(2*math.pi*me/mi))) #potential at electrode
wallV=electrodeV #cylindrical sides of wall same potential
radinfluence=5*lambdaD
dipolea=boxr/100.
mu0=4*math.pi*10**(-7) #Permeaility free space
Bmom=((2*math.pi*(0.003)**3)*0.014/mu0)*numpy.array([0,0,1]) #Nm/T #At 1cm away I want the B to be 0.014T

##OML to find phia
def OLMsol(): #Solve for dust grain surface potential and dust grain charge
	x0=3.3*Te/e
	xnew=x0

	def f(x):
		return math.sqrt(me*Ti/(mi*Te))*(1.-(Te/Ti)*x)-numpy.exp(x) #x is the normalised potential
	def fprime(x):
		return math.sqrt(me*Ti/(mi*Te))*(-Te/Ti)-numpy.exp(x)

	x0=0. #Value approximated as number given in M.Lampe paper
	xnew=x0

	for i in numpy.arange(1000): #number of iterations
		xnew=xnew-f(xnew)/fprime(xnew)
	
	surfpot=xnew #normalised potential
	chargedust=(surfpot*radd*4*math.pi*e0)*(Te*kb/e)/e
	return surfpot*(Te*kb/e), chargedust #This is in Volts and Coloumbs, i.e. SI units

phia,Zd = OLMsol() #phi at surface of dust and the corresponding dust charge


#Moving the magnet up and down will change the strength of the dipole at the bottom
dipolepos=[0,0,-0.00] #If [0,0,0] then dipole is at original position where finite size of approx 0.003m with field strength 0.014 at edge
def pointdipoleB(r):
	r=numpy.array(r)+numpy.array([0,0,-0.003])
	magr=numpy.sqrt(r[0]**2+r[1]**2+r[2]**2)
	rhat=numpy.array(r)/magr
	magBmom=numpy.sqrt(Bmom[0]**2+Bmom[1]**2+Bmom[2]**2)
	Bmomhat=numpy.array(Bmom)/magBmom
	return (mu0*(3*numpy.dot(rhat,numpy.array(Bmom))*rhat-Bmom)/(4*math.pi*magr**3))
B=pointdipoleB(dipolepos)
Bmom=((2*math.pi*(0.003)**3)*B/mu0)*numpy.array([0,0,1]) #Change magnitude of magnetic moment assuming magnet just got weaker
magBmom=numpy.sqrt(Bmom[0]**2+Bmom[1]**2+Bmom[2]**2)
Bmomhat=numpy.array(Bmom)/magBmom 

#Constants following from Joe's thesis
LAMBDA=numpy.sqrt(magBmom*mu0*e/(4*math.pi*me*numpy.sqrt(kb*Te*3/me))) #Distance at which electron is influenced from Joe's thesis
curlyM=magBmom*mu0/(4*math.pi)
VeT=numpy.sqrt(kb*Te*3/me) #Thermal velocity of electrons
critp0=2*numpy.sqrt(e*me*VeT*curlyM) #Critical p0 value determining regimes

##Variables that can change when calculating the region of electron accessibility
##NB: e drift velocity is used for p0 calculation whilst VeT+edrift is used for Z calculation 
electrondriftpos=[boxr,0,0] #Position to calculate the drift in electron velocity
rvalp0 = numpy.sqrt(electrondriftpos[0]**2+electrondriftpos[1]**2+electrondriftpos[2]**2)#The r value at which we are calculating p0


#Calculate E/B drift for electrons rough estimate Assuming B field is 0.014 vertical at all positions
def electrondrift(r):
	omega=abs(e*0.014/me) #Ion cyclotron frequency
	tau=2./omega #electron-neutral collision time #Ratio Taken from konopkas paper on experimentalbasis crystals
	V=-wallV*(r[0]**2+r[1]**2)/(boxr**2)
	mag=(wallV/(boxr**2))*2*math.sqrt(r[0]**2+r[1]**2)
	unitr=[-r[0],-r[1],0]
	if (unitr[0]**2+unitr[1]**2)==0.:
		return numpy.array([0,0,0])
	else:
		Eapplied=numpy.array([-r[0],-r[1],0])*mag/numpy.sqrt(r[0]**2+r[1]**2)
		return numpy.array(Eapplied)*(omega*tau)**2/((1+(omega*tau)**2)*0.014)

edriftvel=electrondrift(r=electrondriftpos)
electrondriftvel=numpy.sqrt(edriftvel[0]**2+edriftvel[1]**2+edriftvel[2]**2)

#Define z+/- function for determining accessible/inaccessible electron regions
def Zpm(r,p,pm):
	if math.isnan(e*curlyM*r**2/(p+pm*me*r*(VeT+electrondriftvel))**(2./3.))==True:
		return 0
	elif (e*curlyM*r**2/(p+pm*me*r*(VeT+electrondriftvel)))**(2./3.)-r**2<=0:
		return 0
	else:
		return numpy.sqrt((e*curlyM*r**2/(p+pm*me*r*(VeT+electrondriftvel)))**(2./3.)-r**2)

# ##Test different values of r effect on p0

# term1=me*rvalp0*electrondriftvel #As boxr is the maximum r that the electron can come from
# term2=e*mu0*magBmom/(4*math.pi*rvalp0) #when electron is at boxr this term is at its minimum
# p0=term1+term2 #As p0 is always very tiny, P0=0 essentially 
# print("Difference p0-critp0 = ",p0-critp0)

# ##Get inaccessble zone at current p0
# r2=numpy.arange(lambdaD,boxr,boxr/10000)
# zinaccesspos=[]
# zinaccessneg=[]
# for i in r2:
# 	zinaccessneg.append(Zpm(i,p0,-1))
# 	zinaccesspos.append(Zpm(i,p0,+1))

# rnorm=numpy.array(r2)/LAMBDA
# znormpos=numpy.array(zinaccesspos)/LAMBDA
# znormneg=numpy.array(zinaccessneg)/LAMBDA

# #Get inaccessible zone at any p0

r=numpy.arange(0,boxr,boxr/10000)
zcrit=[]
for i in r:
	zcrit.append(Zpm(i,critp0,1))
rnormalize=numpy.array(r)/LAMBDA
znormalize=numpy.array(zcrit)/LAMBDA

# ##Plot inaccessibility regions

# fig, ax = plt.subplots(1,1)
# ax.plot(rnorm,znormpos,'b.-',label='Inaccessble below this for this p0, zinaccesspos')
# ax.plot(rnorm,znormneg,'k.-',label='Inaccessble above this for this p0, zinaccessneg')
# ax.fill_between(rnorm,znormpos,znormneg,where=znormneg>znormpos,interpolate=True, color='pink')
# ax.plot(rnormalize,znormalize,'r.-',label="inaccessible for all p0 below curve")
# ax.plot(rnormalize,numpy.ones(len(rnormalize))*sheathd/LAMBDA,'m-',label="Sheath top")
# ax.fill_between(rnormalize,0,znormalize,interpolate=True,color='red')
# plt.xlabel("r/LAMBDA")
# plt.ylabel("z/LAMBDA")
# plt.xlim([0,max(rnormalize)])
# plt.ylim([0,1])
# plt.title("Regions of electron acessibility with B field presense")
# plt.legend()
# plt.show()

############################################################################################
##Find the total charge inside the void assuming only protons present
##Add radial and (vertical?) fields in the sheath
def radialfield(r):#Apply radial E field at the sheath, r = position of dust
	if r[1]<0:
		raise ValueError("particle fell out of cylinder!!!!")
	elif r[1]<sheathd:
		V=-wallV*(r[0]**2+r[1]**2)/(boxr**2)
		mag=(wallV/(boxr**2))*2*math.sqrt(r[0]**2+r[1]**2)
		unitr=[-r[0],0]
		if (unitr[0]**2+unitr[1]**2)==0.:
			return numpy.array([0,0])
		else:
			Eapplied=numpy.array([-1,0])*mag
			return Eapplied
	else:
		return numpy.array([0,0])

def sheathfield(r):
	if r[1]<sheathd and r[1]>=0:
		V=electrodeV*(1.-r[1]/sheathd)**2
		field=abs(2*(1.-r[1]/sheathd)*(electrodeV/sheathd))
		return numpy.array([0,field])
	else:
		return numpy.array([0,0])


def voidvol(rmax):
	intfunc = lambda r : 2*math.pi*r*numpy.sqrt((e*curlyM*r**2/(2*numpy.sqrt(e*curlyM*VeT*me)+me*r*VeT))**(2./3.)-r**2)
	voidvol=scipy.integrate.quad(intfunc,0,rmax)
	return voidvol[0]

def voidQ(rmax):
	Qfunc = lambda r: e*ne0*numpy.exp((e/(kb*Te))*electrodeV*((1.-(1./sheathd)*numpy.sqrt((e*curlyM*r**2/(2*numpy.sqrt(e*curlyM*VeT*me)+me*VeT*r))**(2./3.)-r**2))**2+(r/boxr)**2))*2*math.pi*r*numpy.sqrt((e*curlyM*r**2/(2*numpy.sqrt(e*curlyM*VeT*me)+me*VeT*r))**(2./3.)-r**2)
	voidcharge = scipy.integrate.quad(Qfunc,0,rmax)
	return voidcharge

##Create function for modified potential and hence electric field
def pospos(): #positive charge positions
	separation=lambdaD/2
	chargepos=[] #positive charge positions
	intcharge=int(voidvol(rmax)) #integer number of positive particles
	rows=int(rmax/separation) #number of positive particles separated by lambdaD along the r axis
	for i in numpy.linspace(0,rmax,rows):
		columns=int(Zpm(i,critp0,1)/separation)
		if columns==0:
			chargepos.append([i,Zpm(i,critp0,1)])
		else:
			for j in numpy.arange(columns):
				chargepos.append([i,j*separation]) #number of partices along the z axis separated by lambdaD for specific r value
	return numpy.array(chargepos).T

def magEcharge(r): #magnitude of E field due to charge
	normr=numpy.linalg.norm(r)
	rhat=r/normr
	return abs((voidchargeguess/numbcharge)/(4*math.pi*e0*normr**2))*rhat

@jit
def gridcheck(chargepos):
	Evalsr=[]
	Evalsz=[]
	Evalsradial=[]
	Evalsheath=[]
	Evalsmodifyall=[]
	for i in tqdm(numpy.arange(len(gridr))):
		Evalstempr=[]
		Evalstempz=[]
		Evalsradialtemp=[]
		Evalsheathtemp=[]
		Evalsmodifyalltemp=[]
		for j in tqdm(numpy.arange(len(gridz[i]))):
			totalE=numpy.array([0,0])
			radialE=radialfield([gridr[i][j],gridz[i][j]])
			Evalsradialtemp.append(radialE[0])
			Evalsheathtemp.append(sheathfield([gridr[i][j],gridz[i][j]])[1])
			for k in tqdm(numpy.arange(len(chargepos[0]))):
				r=[chargepos[0][k]-gridr[i][j],chargepos[1][k]-gridz[i][j]]
				if numpy.linalg.norm(r)==0.:
					pass #Ignore the points that fall directly on a positive charge to avoid infinities
				else:
					totalE=totalE+numpy.array(magEcharge(numpy.array(r)))
			Evalstempr.append(totalE[0])
			Evalstempz.append(totalE[1])
			Evalsmodifyalltemp.append(totalE)
		Evalsr.append(Evalstempr)
		Evalsz.append(Evalstempz)
		Evalsradial.append(Evalsradialtemp)
		Evalsheath.append(Evalsheathtemp)
		Evalsmodifyall.append(Evalsmodifyalltemp)
	return Evalsr,Evalsz,Evalsradial,Evalsheath, Evalsmodifyall

def interpolate(r):
	r0=r[0]
	z0=r[1]
	if r0<rmax:
		Eresult=numpy.array([0,0])
		indleft=(r0-firstpoint)/separationhor1
		indlow=(z0-firstpoint)/separationsheath
		s1=isclose((indleft**3) ** (1.0/3), int(indleft))
		s2=isclose((indlow**3) ** (1.0/3), int(indlow))
		if s1==True and s2==True:
			print("case 1")
			gridEs=[[Evalsr[int(indlow)],Evalsz[int(indleft)]]]
			gridpoints=[[gridr[int(indlow)],gridz[int(indleft)]]]
		elif s1 ==True and s2 == False:
			print("case 2")
			indlow=int(indlow)
			gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow+1][indleft],Evalsz[indlow+1][indleft]]]
			gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow+1][indleft],gridz[indlow+1][indleft]]]

		elif s1 ==False and s2 == True:
			print("case 3")
			indleft=int(indleft)
			gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow][indleft+1],Evalsz[indlow][indleft+1]]]
			gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow][indleft+1],gridz[indlow][indleft+1]]]

		else:
			indlow=int(indlow)
			indleft=int(indleft)
			gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow+1][indleft],Evalsz[indlow+1][indleft]],[Evalsr[indlow+1][indleft+1],Evalsz[indlow+1][indleft+1]],[Evalsr[indlow][indleft+1],Evalsz[indlow][indleft+1]]]
			gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow+1][indleft],gridz[indlow+1][indleft]],[gridr[indlow+1][indleft+1],gridz[indlow+1][indleft+1]],[gridr[indlow][indleft+1],gridz[indlow][indleft+1]]]
	elif r0>rmax:
		Eresult=numpy.array([0,0])
		indleft=(boxr-firstpoint-rmax)/separationhor2
		indlow=(z0-firstpoint)/separationsheath
		s1=isclose((indleft**3) ** (1.0/3), int(indleft))
		s2=isclose((indlow**3) ** (1.0/3), int(indlow))
		indleft=int(indleft)
		indlow=int(indlow)
		if s1!=False or s2!=False:
			print("Case 4")
		gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow+1][indleft],Evalsz[indlow+1][indleft]],[Evalsr[indlow+1][indleft+1],Evalsz[indlow+1][indleft+1]],[Evalsr[indlow][indleft+1],Evalsz[indlow][indleft+1]]]
		gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow+1][indleft],gridz[indlow+1][indleft]],[gridr[indlow+1][indleft+1],gridz[indlow+1][indleft+1]],[gridr[indlow][indleft+1],gridz[indlow][indleft+1]]]

	d1=(abs(r0-gridEs[0][0]))
	d4=(abs(r0-gridEs[0][1]))
	d2=(abs(r0-gridEs[2][0]))
	d3=(abs(r0-gridEs[2][1]))
	dhorsq=d1**2+d2**2
	dvertsq=d3**2+d4**2
	Etop=(d2**2/dhorsq)*numpy.array(gridEs[1])+(d1**2/dhorsq)*numpy.array(gridEs[2])
	Ebottom=(d2**2/dhorsq)*numpy.array(gridEs[0])+(d1**2/dhorsq)*numpy.array(gridEs[3])
	Efinal=(d3**2/dvertsq)*Ebottom+(d4**2/dvertsq)*Etop
	return Efinal, gridpoints
		




func = lambda rmax : 1. - e*curlyM/(rmax*(2*numpy.sqrt(e*curlyM*me*VeT)+me*rmax*VeT))
rmax=fsolve(func,0.4*LAMBDA)
voidchargeguess=voidvol(rmax)*ne0*e
voidcharge = voidQ(rmax)[0]
chargepos=pospos()
numbcharge=len(chargepos[0])

separationsheath=lambdaD#separation distance between grid points
separationhor1=lambdaD
separationhor2=lambdaD*10
firstpoint=separationsheath/4
gridlinesheath=numpy.arange(firstpoint,sheathd*1.5,separationsheath)
gridlinehor=numpy.arange(firstpoint,rmax,separationhor1)
gridlinehor2=numpy.arange(rmax+firstpoint,boxr*1.5,separationhor2)
gridhor=list(gridlinehor)+list(gridlinehor2)
gridr,gridz=numpy.meshgrid(gridhor,gridlinesheath)

Evalsr, Evalsz,Evalsradial,Evalsheath,Evalsmodifyall=gridcheck(chargepos)
Evalsradialz=[numpy.zeros(len(Evalsradial[0]))]*len(Evalsradial)
Evalsheathr=[numpy.zeros(len(Evalsheath[0]))]*len(Evalsheath)
def checkedenofm(): #Run this function to calculate ratio of charge in the void that is displaced compared to outside void in rest of sheath
	print("Ratio of charge inside to outside region of inaccessibility=", voidcharge/(ne0*e*(2*math.pi*boxr*sheathd-voidvol(rmax))))
################################Plotting################################################

# plt.figure()
# plt.plot(numpy.array(rnormalize)*LAMBDA,numpy.array(znormalize)*LAMBDA,'r-',label='Inaccessible region for all p0')
# plt.quiver(gridr,gridz,Evalsr,Evalsz,color='b',label='modified field')
# plt.quiver(gridr,gridz,Evalsradial,Evalsradialz,color='k',label='radial field')
# plt.quiver(gridr,gridz,Evalsheathr,Evalsheath,color='g',label='sheath field')
# plt.plot(chargepos.T[:,0],chargepos.T[:,1],'r.')
# plt.plot(numpy.arange(len(gridr[0])),numpy.ones(len(gridr[0]))*0.00038257340070632558,'m-',label='crystal plane')
# plt.legend()
# plt.xlim([0,rmax])
# plt.ylim([0,sheathd*1.5])
# plt.show()

########################################################################################
##Generate particles in their equilibrium position
filehandler = open("crystalpositions2,5K.obj",'rb') ##2k particles 5.5hrs to run 2500 iterations just for settling down
xinit=pickle.load(filehandler)
yinit=pickle.load(filehandler)
zinit=pickle.load(filehandler)
filehandler.close()
initpositions=[[i,j,k] for i,j,k in zip(xinit,yinit,zinit)]

## Test scipy's inbuilt interpolator
tck = interp2d(gridhor, gridlinesheath, Evalsr,kind='cubic') #Using interp2d
tck2 = interp2d(gridhor, gridlinesheath, Evalsz,kind='cubic')
# tck = bisplrep(gridr,gridz,Evalsr) #Using bisplrev
# tck2 = bisplrep(gridr,gridz,Evalsz)
plt.figure()
for p in numpy.arange(len(xinit)):
	if numpy.sqrt(xinit[p]**2+yinit[p]**2)<rmax:
		Etestr=tck(xinit[p], yinit[p]) #Using interp2d
		Etestz=tck2(xinit[p], yinit[p])
		# Etestr = bisplev(xinit[p], yinit[p],tck) #Using bisplrev
		# Etestz = bisplev(xinit[p], yinit[p],tck2)
		#plt.scatter(numpy.array(numpy.mat(gridpointstest).T[0]),numpy.array(numpy.mat(gridpointstest).T[1]),color='k',s=5)
		plt.quiver(numpy.sqrt(xinit[p]**2+yinit[p]**2),zinit[p],Etestr,Etestz,width=0.002,color='m')
	else:
		pass
plt.plot(numpy.array(rnormalize)*LAMBDA,numpy.array(znormalize)*LAMBDA,'r-',label='Inaccessible region for all p0')
plt.quiver(gridr,gridz,Evalsr,Evalsz,color='b',label='modified field')
# plt.quiver(gridr,gridz,Evalsradial,Evalsradialz,color='k',label='radial field')
# plt.quiver(gridr,gridz,Evalsheathr,Evalsheath,color='g',label='sheath field')
# plt.plot(chargepos.T[:,0],chargepos.T[:,1],'r.')
# plt.plot(numpy.arange(len(gridr[0])),numpy.ones(len(gridr[0]))*0.00038257340070632558,'m-',label='crystal plane')
plt.xlabel("r")
plt.ylabel("z")
plt.title("Purple arrows = Interpolated E fields at crystal equilibrium positions")
plt.legend()
plt.xlim([0,rmax])
plt.ylim([0,sheathd*1.5])
plt.show()


##Test my interpolator
plt.figure()
for p in numpy.arange(len(xinit)):
	if numpy.sqrt(xinit[p]**2+yinit[p]**2)<rmax:
		Etest,gridpointstest=interpolate([numpy.sqrt(xinit[p]**2+yinit[p]**2),zinit[p]])
		#plt.scatter(numpy.array(numpy.mat(gridpointstest).T[0]),numpy.array(numpy.mat(gridpointstest).T[1]),color='k',s=5)
		plt.quiver(numpy.sqrt(xinit[p]**2+yinit[p]**2),zinit[p],Etest[0],Etest[1],width=0.002,color='m')
	else:
		pass
plt.plot(numpy.array(rnormalize)*LAMBDA,numpy.array(znormalize)*LAMBDA,'r-',label='Inaccessible region for all p0')
plt.quiver(gridr,gridz,Evalsr,Evalsz,color='b',label='modified field')
# plt.quiver(gridr,gridz,Evalsradial,Evalsradialz,color='k',label='radial field')
# plt.quiver(gridr,gridz,Evalsheathr,Evalsheath,color='g',label='sheath field')
# plt.plot(chargepos.T[:,0],chargepos.T[:,1],'r.')
# plt.plot(numpy.arange(len(gridr[0])),numpy.ones(len(gridr[0]))*0.00038257340070632558,'m-',label='crystal plane')
plt.xlabel("r")
plt.ylabel("z")
plt.title("Purple arrows = Interpolated E fields at crystal equilibrium positions")
plt.legend()
plt.xlim([0,rmax])
plt.ylim([0,sheathd*1.5])
plt.show()

##Save the modified potential in pickle

# filehandler = open(b"modifiedfield.obj",'wb') ##2k particles 5.5hrs to run 2500 iterations just for settling down
# pickle.dump(gridr,filehandler)
# pickle.dump(gridz,filehandler)
# pickle.dump(Evalsr,filehandler)
# pickle.dump(Evalsz,filehandler)
# pickle.dump(rmax,filehandler)
# pickle.dump(separationsheath,filehandler)
# pickle.dump(separationhor1,filehandler)
# pickle.dump(separationhor2,filehandler)
# pickle.dump(firstpoint,filehandler)
# filehandler.close()














