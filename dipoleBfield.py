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


#Some constants
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
radinfluence=10*lambdaD
dipolea=boxr/100.
mu0=4*math.pi*10**(-7) #Permeaility free space
Bmom=((2*math.pi*(0.003)**3)*0.014/mu0)*numpy.array([0,0,1]) #Nm/T #At 1cm away I want the B to be 0.014T
magBmom=numpy.sqrt(Bmom[0]**2+Bmom[1]**2+Bmom[2]**2)
Bmomhat=numpy.array(Bmom)/magBmom
dipolepos=[0,0,-0.0005] 

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



def pointdipoleB(r):
	r=numpy.array(r)-numpy.array(dipolepos)
	magr=numpy.sqrt(r[0]**2+r[1]**2+r[2]**2)
	rhat=numpy.array(r)/magr
	magBmom=numpy.sqrt(Bmom[0]**2+Bmom[1]**2+Bmom[2]**2)
	Bmomhat=numpy.array(Bmom)/magBmom
	return (mu0*(3*numpy.dot(rhat,numpy.array(Bmom))*rhat-Bmom)/(4*math.pi*magr**3))

def finitesolenoid(r,N=50,L=0.5,rad=0.1,I=1000): #N coils over length of solenoid L, rad is radius of solenoid cross section, current I
	magr = numpy.sqrt(r[0]**2+r[1]**2)
	sigmaplus = r[2]+L/2.
	sigmaminus = r[2]-L/2.
	n = N/L
	ksquarep = 4*radd*magr/(sigmaplus**2+(rad+magr)**2)
	kp=numpy.sqrt(ksquarep)
	ksquarem = 4*radd*magr/(sigmaminus**2+(rad+magr)**2)
	km=numpy.sqrt(ksquarem)
	phip=numpy.arctan(abs(sigmaplus/(rad-magr)))
	phim=numpy.arctan(abs(sigmaminus/(rad-magr)))
	def zeta(phi,ksquare):
		return special.ellipeinc(phi,ksquare)-special.ellipe(ksquare)*special.ellipkinc(phi,ksquare)/special.ellipk(ksquare)
	def heuman(phi,ksquare):
		return (special.ellipkinc(phi,(1.-ksquare))/special.ellipk((1-ksquare)))+(2./math.pi)*special.ellipk(ksquare)*zeta(phi,(1-ksquare))

	Br=(mu0*n*I/math.pi)*numpy.sqrt(radd/magr)*\
		((((2-ksquarep)/(2*kp))*special.ellipk(kp)-special.ellipe(kp)/kp)\
		-(((2-ksquarem)/(2*km))*special.ellipk(km)-special.ellipe(km)/km))
	x,y=((numpy.array(r)/magr)*Br)[0:2]
	Bz=(mu0*n*I/4.)*((sigmaplus*kp*special.ellipk(kp)/(math.pi*numpy.sqrt(magr*rad))+(rad-magr)*sigmaplus*heuman(phip,kp)/abs((rad-magr)*sigmaplus))\
		-(sigmaminus*km*special.ellipk(km)/(math.pi*numpy.sqrt(magr*rad))+(rad-magr)*sigmaminus*heuman(phim,km)/abs((rad-magr)*sigmaminus)))
	print(Br,Bz)
	return numpy.array([x,y,0])


##Plot 3D quiver B field
heightbox=sheathd
widthbox=0.002
dx=0.
stepsize=widthbox/3
stepsizeheight=heightbox/3

#X, Y, Z = numpy.meshgrid(numpy.arange(-0.5, 0.5 , stepsize),numpy.arange(-0.5, 0.5 , stepsize), numpy.arange(-0.5, 0.5, stepsize))
X, Y, Z = numpy.meshgrid(numpy.arange(-widthbox, widthbox , stepsize),numpy.arange(-widthbox, widthbox, stepsize), numpy.arange(-heightbox+dx, heightbox+dx, stepsizeheight))

U=[]
V=[]
W=[]
for i in numpy.arange(len(Y)):
	rowU=[]
	rowV=[]
	rowW=[]
	for j in numpy.arange(len(Y[i])):
		rowU2=[]
		rowV2=[]
		rowW2=[]
		for k in numpy.arange(len(Y[i][j])):
			if abs(numpy.sqrt(X[i][j][k]**2+Y[i][j][k]**2+Z[i][j][k]**2))<0.05*boxz:
				rowU2.append(0)
				rowV2.append(0)
				rowW2.append(0)
			else:
				newpoint=pointdipoleB([X[i][j][k],Y[i][j][k],Z[i][j][k]])
				rowU2.append(newpoint[0])
				rowV2.append(newpoint[1])
				rowW2.append(newpoint[2])
		rowU.append(rowU2)
		rowV.append(rowV2)
		rowW.append(rowW2)
	U.append(numpy.array(rowU))
	V.append(numpy.array(rowV))
	W.append(numpy.array(rowW))
U=numpy.array(U)
V=numpy.array(V)
W=numpy.array(W)

# fig=plt.figure()
# ax=fig.gca(projection='3d')
# #ax.view_init(elev=90., azim=90)
# Ulist=abs(numpy.hstack(numpy.hstack(U)))
# Vlist=abs(numpy.hstack(numpy.hstack(V)))
# Wlist=abs(numpy.hstack(numpy.hstack(W)))
# maximumarrow=max(max(Ulist),max(Vlist),max(Wlist))
# Q = ax.quiver3D(X,Y, Z, U, V, W,length=0.0001,arrow_length_ratio=0.5,normalize=True) #length is just a multiplication factor
# #P = ax.scatter(X,Y,Z)
# # maxx=max(Ulist)
# # maxy=max(Vlist)
# # maxz=max(Wlist)
# # ax.set_xlim([-widthbox,widthbox])
# # ax.set_ylim([-widthbox,widthbox])
# # ax.set_zlim([-heightbox,heightbox+dx])
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.show()


# ################2D quiver plot?
# plt.figure()
# X=np.mean(X,axis=2)
# plt.quiver(X,Y,U,V)
# set_xlabel("x")
# set_ylabel("y")
# plt.show()













##Plot for pointdipole
'''
heightbox=0.00038
widthbox=boxr
dx=1.
stepsize=widthbox/3

stepsizeheight=heightbox/3
X, Y, Z = numpy.meshgrid(numpy.arange(-0.5, 0.5 , stepsize),numpy.arange(-0.5, 0.5 , stepsize), numpy.arange(-0.5, 0.5, stepsize))
#X, Y, Z = numpy.meshgrid(numpy.arange(-widthbox, widthbox , stepsize),numpy.arange(-widthbox, widthbox, stepsize), numpy.arange(-heightbox+dx, heightbox+dx, stepsizeheight))

U=[]
V=[]
W=[]
for i in numpy.arange(len(Y)):
	rowU=[]
	rowV=[]
	rowW=[]
	for j in numpy.arange(len(Y[i])):
		rowU2=[]
		rowV2=[]
		rowW2=[]
		for k in numpy.arange(len(Y[i][j])):
			if abs(numpy.sqrt(X[i][j][k]**2+Y[i][j][k]**2+Z[i][j][k]**2))<100*dipolea:
				rowU2.append(0)
				rowV2.append(0)
				rowW2.append(0)
			else:
				newpoint=pointdipoleB([X[i][j][k],Y[i][j][k],Z[i][j][k]])
				rowU2.append(newpoint[0])
				rowV2.append(newpoint[1])
				rowW2.append(newpoint[2])
		rowU.append(rowU2)
		rowV.append(rowV2)
		rowW.append(rowW2)
	U.append(numpy.array(rowU))
	V.append(numpy.array(rowV))
	W.append(numpy.array(rowW))
U=numpy.array(U)
V=numpy.array(V)
W=numpy.array(W)

fig=plt.figure()
ax=fig.gca(projection='3d')
Ulist=abs(numpy.hstack(numpy.hstack(U)))
Vlist=abs(numpy.hstack(numpy.hstack(V)))
Wlist=abs(numpy.hstack(numpy.hstack(W)))
maximumarrow=max(max(Ulist),max(Vlist),max(Wlist))
Q = ax.quiver3D(X,Y, Z, U, V, W,length=0.0015,linewidth=1.5,arrow_length_ratio=0.5,normalize=False) #length is just a multiplication factor
#P = ax.scatter(X,Y,Z)
# maxx=max(Ulist)
# maxy=max(Vlist)
# maxz=max(Wlist)
# ax.set_xlim([-widthbox,widthbox])
# ax.set_ylim([-widthbox,widthbox])
# ax.set_zlim([-heightbox,heightbox+dx])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
'''


