##This file is mainly for looking at B field effects of crystals that are already formed.
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
import pickle
from tqdm import tqdm
from numpy import isclose
from scipy import interpolate

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
dt=0.000001 #Much longer than charging time which is order nanoseconds
sheathd=10*lambdaD
electrodeV=abs((kb*Te/(2*e))*(numpy.log(2*math.pi*me/mi))) #potential at electrode
wallV=electrodeV #cylindrical sides of wall same potential
radinfluence=0.002#10*lambdaD
dipolea=boxr/100.
mu0=4*math.pi*10**(-7) #Permeaility free space
Bmom=((2*math.pi*(0.003)**3)*0.014/mu0)*numpy.array([0,0,1]) #Nm/T #At 1cm away I want the B to be 0.014T
magBmom=numpy.sqrt(Bmom[0]**2+Bmom[1]**2+Bmom[2]**2)
Bmomhat=numpy.array(Bmom)/magBmom
dipolepos=[0,0,-0.0005] 

omegatau=0.01 #Referenced by konopka in his experiments for ion collisions
class ion:
	def __init__(self,pos,vel,acc):
		self.pos=numpy.array(pos)
		self.vel=numpy.array(vel)
		self.acc=numpy.array(acc)
		self.pos1=numpy.array(pos)
		self.vel1=numpy.array(vel)
		self.charge=e
		self.dr=0
	def constB(self):
		return numpy.array([0,0.014,0])
	def constE(self,E=[0,0,1]):
		return numpy.array(E)
	def updateeulerforward(self,B,E=[0,0,1]):
		##Euler forward - should be unstable!
		self.pos = self.pos + dt*self.vel
		self.vel = ((self.charge/mi)*numpy.cross(self.vel,B) + numpy.array(E)*(self.charge/mi))*dt+self.vel
	def updateeulerback(self,B,E=[0,0,1]): #Also unstable?
		self.pos = self.pos - dt*self.vel
		self.vel = self.vel - ((self.charge/mi)*numpy.cross(self.vel,B)*dt + numpy.array(E)*(self.charge/mi))
	def updateRK4(self,B,E=[0,0,1]):
		##RK4 integration
		fv1=(self.charge/mi)*numpy.cross(self.vel1,B) + numpy.array(E)*(self.charge/mi)
		fy1=self.vel1
		v1=self.vel1+0.5*dt*fv1
		y1=self.pos1+0.5*dt*fy1
		fv2=(self.charge/mi)*numpy.cross(v1,B) + numpy.array(E)*(self.charge/mi)
		fy2=v1
		v2=self.vel1+0.5*dt*fv2
		y2=self.pos1+0.5*dt*fy2
		fv3=(self.charge/mi)*numpy.cross(v2,B) + numpy.array(E)*(self.charge/mi)
		fy3=v2
		v3=self.vel1+dt*fv3
		y3=self.pos1+dt*fy3
		fv4=(self.charge/mi)*numpy.cross(v3,B) + numpy.array(E)*(self.charge/mi)
		fy4=v3
		self.vel=self.vel1+(1./6.)*dt*(fv1+2.*fv2+2.*fv3+fv4)
		self.pos=self.pos1+(1./6.)*dt*(fy1+2.*fy2+2.*fy3+fy4)
		self.vel1=self.vel
		self.dr = numpy.array(self.pos)-numpy.array(self.pos1)
		self.pos1=self.pos
	def updateboris(self,B):
		pass

	def getselfvel(self):
		return self.vel
	def getselfpos(self):
		return self.pos

def runsim(iterations=200):
	iterations = iterations
	vinit=numpy.sqrt(kb*Ti/mi)
	ion1=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])
	# ion2=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])
	# ion3=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])

	position=[]
	velocity=[]
	# velocity2=[]
	# position2=[]
	# position3=[]
	for i in numpy.arange(iterations):
		ion1.updateRK4(B=ion1.constB())
		# ion2.updateeulerforward(B=ion2.constB())
		# ion3.updateeulerback(B=ion3.constB())
		position.append(ion1.getselfpos())
		velocity.append(numpy.linalg.norm(ion1.getselfvel()))
		# velocity2.append(numpy.linalg.norm(ion2.getselfvel()))
		# position2.append(ion2.getselfpos())
		# position3.append(ion3.getselfpos())

	newx=[i[0] for i in position]
	newy=[i[1] for i in position]
	newz=[i[2] for i in position]
	# newx2=[i[0] for i in position2]
	# newy2=[i[1] for i in position2]
	# newz2=[i[2] for i in position2]
	# newx3=[i[0] for i in position3]
	# newy3=[i[1] for i in position3]
	# newz3=[i[2] for i in position3]
	direction=numpy.cross(ion1.constE(),ion1.constB())/(numpy.linalg.norm(numpy.cross(ion1.constE(),ion1.constB())))
	return newx,newy,newz,numpy.array(direction),ion1

def plotrunsim(newx,newy,newz):
	newx=newx
	newy=newy
	newz=newz
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot(newx,newy,newz,'r-',label="RK4")
	#ax.plot(newx2,newy2,newz2,'b-',label="Euler Forwards")
	#ax.plot(newx3,newy3,newz3,'m-',label="Euler backwards")

	plt.legend()
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")
	plt.show()

def checkdrift():	
	iterations = numpy.arange(200,20000,500)
	particledrift=[]
	for i in tqdm(iterations):
		x,y,z,direction,ion1=runsim(i)
		EXBdrift = numpy.linalg.norm(ion1.constE())/numpy.linalg.norm(ion1.constB())
		iterationtime = dt*i
		totald=EXBdrift*iterationtime
		theorydrift=numpy.linalg.norm(direction*totald)
		particlepos=numpy.linalg.norm(numpy.array([x[-1],y[-1],z[-1]])*direction)
		particledrift.append([theorydrift,particlepos])
	particledrift = numpy.array(numpy.mat(particledrift).T)
	plt.figure()
	plt.scatter(iterations,particledrift[1]/particledrift[0])
	plt.ylabel("r_simulation/r_theory")
	plt.xlabel("Iteration number")
	plt.title("Ratio of simulated drift position vs theoretical drift position")
	plt.show()



#checkdrift()

def compareFEM(iterations=500):
	iterations = iterations
	vinit=numpy.sqrt(kb*Ti/mi)
	ion1=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])
	ion2=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])
	ion3=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])

	KERK4=[]
	KEEF=[]
	KEEB=[]
	for i in numpy.arange(iterations):
		ion1.updateRK4(B=ion1.constB(),E=[0,0,0])
		ion2.updateeulerforward(B=ion2.constB(),E=[0,0,0])
		ion3.updateeulerback(B=ion3.constB(),E=[0,0,0])
		KERK4.append(0.5*mi*numpy.linalg.norm(ion1.getselfvel())**2)
		KEEF.append(0.5*mi*numpy.linalg.norm(ion2.getselfvel())**2)
		KEEB.append(0.5*mi*numpy.linalg.norm(ion3.getselfvel())**2)

	#Compare Kinetic energy NB: when E field is [0,0,0] between RK4 an Euler forward
	plt.figure()
	plt.plot(numpy.arange(len(KERK4)),KERK4,label="RK4")
	plt.plot(numpy.arange(len(KEEB)),KEEB,'r*',label="Euler backward")
	plt.plot(numpy.arange(len(KEEF)),KEEF,label="Euler forward")
	plt.legend()
	plt.ylim([0,max([max(KEEF)*1.1,max(KEEB)*1.1])])
	plt.title("No E field present, iterations = %s"%iterations)
	plt.xlabel("Iteration number")
	plt.ylabel("Kinetic energy of particle")
	plt.show()

#compareFEM()

def thermalkickexb():
	pass


