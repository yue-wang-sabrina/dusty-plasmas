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
dipolepos=[0,0,-0.01] 

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

class Dust: 
	def __init__(self,mass=0,radius=0,lambdaD=0,phia=phia, initcharge=0,initpos=[0,0,0],initvel=[0,0,0],initacc=[0,0,0]):
		self.m=mass
		self.rad=radius
		self.lambdaD=lambdaD
		self.charge=initcharge
		self.pos=initpos #[rcostheta,rsintheta,z]
		self.vel=initvel
		self.acc=initacc
		self.phia=phia
		self.sheathd=10*lambdaD
		self.multifields=numpy.array([0,0,0])
		self.Bswitch=False
		#Extra testing stuff for diff numerical schemes and checks: can be deleted if not used later
		self.pos1=initpos
		self.pos2=initpos
		self.acc1=initacc
		self.vel1=initvel
		self.vel2=initvel
		self.check=0
		self.mach=numpy.array([0,0,0])

	def getselfpos(self):
		return self.pos
	def getselfvel(self):
		return self.vel
	def getselfacc(self):
		return self.acc
	def scatterwalls(self):
		if abs(numpy.array(self.pos[0])+numpy.array(self.vel[0])*dt+0.5*numpy.array(self.acc[0])*dt**2)+abs(numpy.array(self.pos[1])+numpy.array(self.vel[1])*dt+0.5*numpy.array(self.acc[1])*dt**2)>boxr:
			print("Particle is scattered on cylinder walls at position", self.getselfpos())
			self.vel[0]*=-1
			self.vel[1]*=-1
		if numpy.array(self.pos[2])+numpy.array(self.vel[2])*dt+0.5*numpy.array(self.acc[2])*dt**2<0:
			self.vel[2]=abs(self.vel[2])
			print("Particle is vertically scattered")

	def step(self): #Equation of motion involving forces, to update position and velocitiy and acceleration in time
		##Add gravitational acceleration
		self.acc=numpy.array([0,0,g])
		##Add linearly varying E field caused acceleration
		accsheath=self.sheathfield()*abs(self.charge/self.m)
		accradial=self.radialfield()*abs(self.charge/self.m)
		self.acc=numpy.array(self.acc)+accsheath+accradial+(E*numpy.array(unitr)*self.charge/self.m)+self.verticalion()
		self.vel=numpy.array(self.vel)+dt*numpy.array(self.acc)
		self.pos=numpy.array(self.pos)+numpy.array(self.vel)*dt+0.5*numpy.array(self.acc)*dt**2
		
	
	def steptest(self,Emultifield=numpy.array([0,0,0])): #Test function for interactions, alter incooporate into step function
		#self.scatterwalls() #Check if hit walls
		##Add gravitational acceleration
		self.acc=numpy.array([0,0,g])
		##Add linearly varying E field caused acceleration
		accsheath=self.sheathfield()*abs(self.charge/self.m)
		accradial=self.radialfield()*abs(self.charge/self.m)
		accB=self.EXBacchybrid(B=self.dipoleB(r=dipolepos))
		self.acc=numpy.array(self.acc)+accsheath+(Emultifield*abs(self.charge/self.m))+self.damping()+accradial+accB
		self.vel=numpy.array(self.vel)+dt*numpy.array(self.acc)
		self.pos1=self.pos #Save previous position
		self.pos=numpy.array(self.pos)+numpy.array(self.vel)*dt+0.5*numpy.array(self.acc)*dt**2
		self.multifields=numpy.array([0,0,0])

		##Leapfrog - unstable
		# self.vel=numpy.array(self.vel2)+2*dt*numpy.array(self.acc1)
		# self.pos=numpy.array(self.pos2)+2*dt*numpy.array(self.vel1)
		# self.acc1=self.acc
		# self.vel2=self.vel1
		# self.vel1=self.vel
		# self.pos2=self.pos1
		# self.pos1=self.pos
		##RK4 - also unstable?
		# fv1=numpy.array(self.acc1)
		# fy1=numpy.array(self.vel1)
		# v1=numpy.array(self.vel1)+0.5*dt*fv1
		# y1=numpy.array(self.pos1)+0.5*dt*fy1
		# fv2=numpy.array(self.acc1)
		# fy2=v1
		# v2=self.vel1+0.5*dt*fv2
		# y2=self.pos1+0.5*dt*fy2
		# fv3=fv2
		# fy3=v2
		# v3=self.vel1+dt*fv3
		# y3=self.pos1+dt*fy3
		# fv4=fv3
		# fy4=v3
		# self.vel=self.vel1+(1./6.)*dt*(fv1+2.*fv2+2.*fv3+fv4)
		# self.pos=self.pos1+(1./6.)*dt*(fy1+2.*fy2+2.*fy3+fy4)
		# self.acc1=self.acc
		# self.vel1=self.vel
		# self.pos1=self.pos


	def damping(self):
		gamma=5000.
		magvel=numpy.sqrt(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
		if magvel==0:
			return numpy.array([0,0,0])
		else:
			acc=-1*gamma*(numpy.array(self.vel)/magvel)*magvel**2
			return acc
	def verticalion(self):
		if self.pos[2]<sheathd:
			mach=5.
			beta=abs(self.charge*e/(Ti*lambdade))
			return numpy.array([0,0,(((Ti/e)**2)*numpy.log(lambdade*mach**2/(beta*lambdadi))*beta**2/mach**2)/self.m])
		else:
			return numpy.array([0,0,0])

	def selffield(self,g2): #E field due to dust g2. NB: Assuming all particles have same size!
		##Energy between dressed grains (see Shukla book)
		#energy=(self.charge*grain2.charge/r)*math.exp(-r/self.lambdaD)*(1.-r/(2*self.lambdaD))
		#return energy
		##Electric field at some distance r from the dust particle
		d=numpy.sqrt((g2.pos[0]-self.pos[0])**2+(g2.pos[1]-self.pos[1])**2+(g2.pos[2]-self.pos[2])**2)
		di=numpy.sqrt((g2.pos1[0]-self.pos1[0])**2+(g2.pos1[1]-self.pos1[1])**2+(g2.pos1[2]-self.pos1[2])**2)
		if di==0.:
			pass#print("di=0 PROBLEM")
		connecti=numpy.array([self.pos1[0]-g2.pos1[0],self.pos1[1]-g2.pos1[1],self.pos1[2]-g2.pos1[2]])/di
		connectf=numpy.array([self.pos[0]-g2.pos[0],self.pos[1]-g2.pos[1],self.pos[2]-g2.pos[2]])/d
		if connecti[0]!=connectf[0] or connecti[1]!=connectf[1] or connecti[2]!=connectf[2]:
			pass#raise ValueError("Particles passed each other")
		if d-2*self.rad<0.:
			pass#raise ValueError("Particles are inside each other!?")
		r=[self.pos[0]-g2.pos[0],self.pos[1]-g2.pos[1],self.pos[2]-g2.pos[2]]
		unitr=numpy.array(r)/d
		E=abs((self.phia/self.lambdaD)*numpy.exp(-(d-2*self.rad)/self.lambdaD))
		return E*unitr

	def selffieldmany(self,E):
		self.multifields=self.multifields+numpy.array(E)

	def sheathfield(self):
		if self.pos[2]>=self.sheathd and self.pos[2]>0:
			return numpy.array([0,0,0])
		elif self.pos[2]<self.sheathd and self.pos[2]>0:
			V=electrodeV*(1.-self.pos[2]/self.sheathd)**2
			field=abs(2*(1.-self.pos[2]/sheathd)*(electrodeV/sheathd))
			return numpy.array([0,0,field])
		else:
			raise ValueError("dust fell out of cylinder!!With position",self.pos)

	def radialfield(self):#Apply radial E field at the sheath
		if self.pos[2]<0:
			raise ValueError("particle fell out of cylinder!!!!")
		elif self.pos[2]<sheathd:
			V=-wallV*(self.pos[0]**2+self.pos[1]**2)/(boxr**2)
			mag=(wallV/(boxr**2))*2*math.sqrt(self.pos[0]**2+self.pos[1]**2)
			unitr=[-self.pos[0],-self.pos[1],0]
			if (unitr[0]**2+unitr[1]**2)==0.:
				return numpy.array([0,0,0])
			else:
				Eapplied=numpy.array([-self.pos[0],-self.pos[1],0])*mag/numpy.sqrt(self.pos[0]**2+self.pos[1]**2)
				return numpy.array(Eapplied)
		else:
			return numpy.array([0,0,0])

	def intergraind(self,g2):
		if numpy.sqrt((self.pos[0]-g2.pos[0])**2+(self.pos[1]-g2.pos[1])**2+(self.pos[2]-g2.pos[2])**2)<=radinfluence:
			return True
		else:
			return False
	def EXBacc(self,B=[0,0,0.014],machmag=0.): ##Ion drag force due to Constant vertical B field causing ExB drift of ions
		if self.pos[2]<sheathd and self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vT=numpy.sqrt(kb*Ti/mi)#Thermal speed ions
			E=self.radialfield()#+self.sheathfield()
			vdrift=numpy.cross(E,numpy.array(B))/(magB**2) #vdrift for constant vertical B field
			
			# ##Ion-neutral collision change to vdrift
			# omega=abs(e*magB/mi) #Ion cyclotron frequency
			# tau=0.005/omega #ion-neutral collision time
			# vdrift[0]*=(omega*tau)**2/(1+(omega*tau)**2) #neglect v_r component
			# vdrift[1]*=(omega*tau)**2/(1+(omega*tau)**2)
			# vdrift[2]*=(omega*tau)**2/(1+(omega*tau)**2)

			
			##Ion-neutral collision drift velocity new D.D. Millar 1976
			omega=abs(e*magB/mi)
			tau=0.05/omega
			Br=numpy.sqrt(B[0]**2+B[1]**2)
			Er=numpy.sqrt(E[0]**2+E[1]**2)
			k=((e/mi)**2)*(B[2]*Er-Br*E[2])
			p=Er-B[2]*k/omega**2
			s=E[2]+Br*k/omega**2
			drift0=(-k/(omega**2+1./tau**2))
			drift1=(e/(tau**2*mi))*(tau**3*p-(B[2]*k/omega**4)*(-tau**3*omega**2/(1+tau**2*omega**2)))
			drift2=(e/(tau**2*mi))*(tau**3*s+(Br*k/omega**4)*(-tau**3*omega**3/(1+tau**2*omega**2)))

			r=numpy.sqrt(self.pos[0]**2+self.pos[1]**2)
			theta=numpy.arctan(self.pos[2]/self.pos[0])
			vdrift[0]=abs((drift1)*r*numpy.sin(theta))*numpy.sign(vdrift[0])
			vdrift[1]=abs((drift1)*r*numpy.cos(theta))*numpy.sign(vdrift[1])
			#print(vdrift[0],vdrift[1])

			##Calculate total acceleration on dust particle due to ion and neutral drag
			mach=numpy.array([vdrift[0]-self.vel[0],vdrift[1]-self.vel[1],vdrift[2]-self.vel[2]])/vT
			machmag=numpy.sqrt(mach[0]**2+mach[1]**2+mach[2]**2)
			self.mach=mach
			if machmag<=2:
				LAMBDA=numpy.sqrt(1./(numpy.exp(-machmag**2/2.)*lambdadi**(-2)+lambdade**(-2)))
				beta=abs(self.charge*e/(4*math.pi*e0*kb*Ti*LAMBDA))
				coloumblog=5.
				forceS=(1./3.)*numpy.sqrt(32*math.pi)*((Ti*kb/e)**2)*coloumblog*beta**2*e0*mach
				phi=abs(self.charge*e/(kb*Te*4*math.pi*e0*self.rad))
				forceC=mach*vT*4*math.pi*self.rad**2*ni0*mi*numpy.sqrt(kb*Te/(2*math.pi*me))*numpy.exp(-phi)
				forceN=-numpy.array(self.vel)*mi*numpy.sqrt(4*math.pi)*ni0*numpy.sqrt(kb*Ti/(2*math.pi*mi))*math.pi*radd**2
				nn= 7./(kb*Ti)#Number density neutrals, 7Pa used by Saitou
				forceionneutraldrag=mach*mi*numpy.sqrt(4*math.pi)*nn*numpy.sqrt(kb*Ti/(2*math.pi*mi))*math.pi*radd**2
				forcetotal=forceS+forceC+forceN
				#print(forceionneutraldrag,forceS,forceN,forceC)
				return forcetotal/self.m
			else: 
				return mach*machmag*math.pi*Ti*kb*self.rad**2*ni0/self.m
		else:
			return numpy.array([0,0,0])#raise ValueError("Dust particle is not in sheath")

	def EXBaccfortov(self,B=[0,0,0.014],machmag=0.): ##Test other literature theories
		if self.pos[2]<sheathd and self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vT=numpy.sqrt(kb*Ti/mi)#Thermal speed ions
			vdrift=numpy.cross(self.radialfield()+self.sheathfield(),numpy.array(B))/(magB**2) #vdrift for constant vertical B field
			mach=numpy.array([vdrift[0]-self.vel[0],vdrift[1]-self.vel[1],vdrift[2]-self.vel[2]])/vT
			machmag=numpy.sqrt(mach[0]**2+mach[1]**2+mach[2]**2)

			LAMBDA=numpy.sqrt(1./(numpy.exp(-machmag**2/2.)*lambdadi**(-2)+lambdade**(-2)))
			beta=abs(Zd*e**2/(LAMBDA*Ti*kb))
			if machmag<=1 and beta<=13:
				coloumblog=5.
				force=abs((1./3.)*numpy.sqrt(2./math.pi)*(beta**2)*coloumblog)*mach

				return force/self.m
			elif machmag<=1 and beta>13:
				force=abs((2./3.)*numpy.sqrt(2/math.pi)*((numpy.log(beta))**2+2*numpy.log(beta)+2))*mach

				return force/self.m
			elif machmag>1:
				force=abs((beta**2*numpy.log(lambdade*machmag**2/(lambdaD*beta))*(1./machmag)**2))*mach/(machmag)

				return force/self.m
			else:
				return numpy.array([0,0,0])
		else:
			return numpy.array([0,0,0])

	def EXBacckinetic(self,B=[0,0,0.014],machmag=0.):
		if self.pos[2]<sheathd and self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vT=numpy.sqrt(kb*Ti/mi)#Thermal speed ions
			vdrift=numpy.cross(self.radialfield()+self.sheathfield(),numpy.array(B))/(magB**2) #vdrift for constant vertical B field
			mach=numpy.array([vdrift[0]-self.vel[0],vdrift[1]-self.vel[1],vdrift[2]-self.vel[2]])/vT
			machmag=numpy.sqrt(mach[0]**2+mach[1]**2+mach[2]**2)

			LAMBDA=numpy.sqrt(1./(numpy.exp(-machmag**2/2.)*lambdadi**(-2)+lambdade**(-2)))
			beta=abs(Zd*e**2/(LAMBDA*Ti*kb))
			E=self.radialfield()
			incollfreq=e*numpy.sqrt(E[0]**2+E[1]**2+E[2]**2)/(mi*numpy.sqrt(vdrift[0]**2+vdrift[1]**2+vdrift[2]**2))
			l=vT/incollfreq
			LAMBDA=1./numpy.sqrt(lambdade**(-2)+lambdadi**(-2))
			if machmag<=1:
				force=abs((1./3.)*numpy.sqrt(2/math.pi)*beta**2*(numpy.log(1./beta)+(1./numpy.sqrt(2*math.pi))*K(LAMBDA/l)))*mach

				return force/self.m
			else:
				force=abs(numpy.sqrt(2./math.pi)*beta**2*numpy.log(4*l*machmag/(LAMBDA*beta))*(1./machmag))*mach/machmag
				return force/self.m
		else:
			return numpy.array([0,0,0])#raise ValueError("Dust particle is not in sheath")

	def EXBacchybrid(self,B=[0,0,0.014],machmag=0.):
		if self.pos[2]<sheathd and self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vT=numpy.sqrt(kb*Ti/mi)#Thermal speed ions
			vdrift=numpy.cross(self.radialfield()+self.sheathfield(),numpy.array(B))/(magB**2) #vdrift for constant vertical B field
			E=self.radialfield()+self.sheathfield()

			##Just blindly multiplying by factor
			omega=abs(e*magB/mi)
			tau=0.05/omega			
			vdrift[0]*=(omega*tau)**2/(1+(omega*tau)**2)
			vdrift[1]*=(omega*tau)**2/(1+(omega*tau)**2)
			vdrift[2]*=(omega*tau)**2/(1+(omega*tau)**2)

			##Ion-neutral collision drift velocity new D.D. Millar 1976
			# omega=abs(e*magB/mi)
			# tau=0.05/omega
			# Br=numpy.sqrt(B[0]**2+B[1]**2)
			# Er=numpy.sqrt(E[0]**2+E[1]**2)
			# k=((e/mi)**2)*(B[2]*Er-Br*E[2])
			# p=Er-B[2]*k/omega**2
			# s=E[2]+Br*k/omega**2
			# drift0=(-k/(omega**2+1./tau**2))
			# drift1=(e/(tau**2*mi))*(tau**3*p-(B[2]*k/omega**4)*(-tau**3*omega**2/(1+tau**2*omega**2)))
			# drift2=(e/(tau**2*mi))*(tau**3*s+(Br*k/omega**4)*(-tau**3*omega**3/(1+tau**2*omega**2)))
			# r=numpy.sqrt(self.pos[0]**2+self.pos[1]**2)
			# theta=numpy.arctan(self.pos[2]/self.pos[0])
			# vdrift[0]=abs((drift1)*r*numpy.sin(theta))*numpy.sign(vdrift[0])
			# vdrift[1]=abs((drift1)*r*numpy.cos(theta))*numpy.sign(vdrift[1])
			# vdrift[2]=abs(drift2)*numpy.sign(vdrift[2])

			mach=numpy.array([vdrift[0]-self.vel[0],vdrift[1]-self.vel[1],vdrift[2]-self.vel[2]])/vT
			machmag=numpy.sqrt(mach[0]**2+mach[1]**2+mach[2]**2)
			LAMBDA=numpy.sqrt(1./(numpy.exp(-machmag**2/2.)*lambdadi**(-2)+lambdade**(-2)))
			beta=abs(Zd*e**2/(LAMBDA*Ti*kb))
			u=numpy.sqrt(vdrift[0]**2+vdrift[1]**2+vdrift[2]**2)/vT
			z=abs(Zd)*e**2/(radd*Te*kb)
			tau=Te/Ti
			coloumblog=5.
			force=numpy.sqrt(2*math.pi)*radd**2*ni0*mi*vT**2*\
				(numpy.sqrt(math.pi/2)*special.erf(u/numpy.sqrt(2))*\
					(1+u**2+(1-u**2)*(1+2*z*tau)+4*z**2*tau**2*u**(-2)*numpy.log(coloumblog))+\
						(u**(-1)*(1+2*z*tau+u**2-4*z**2*tau**2*numpy.log(coloumblog))*numpy.exp(-u**2/2.)))*mach/machmag
			return force/self.m

		else:
			return numpy.array([0,0,0])#raise ValueError("Dust particle is not in sheath")	

	def EXBaccDUSTT(self,B=[0.,0.,0.014],machmag=0.):
		if self.pos[2]<sheathd and self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vT=numpy.sqrt(kb*Ti/mi)#Thermal speed ions
			vdrift=numpy.cross(self.radialfield()+self.sheathfield(),numpy.array(B))/(magB**2) #vdrift for constant vertical B field
			mach=numpy.array([vdrift[0]-self.vel[0],vdrift[1]-self.vel[1],vdrift[2]-self.vel[2]])/vT
			machmag=numpy.sqrt(mach[0]**2+mach[1]**2+mach[2]**2)

			LAMBDA=numpy.sqrt(1./(numpy.exp(-machmag**2/2.)*lambdadi**(-2)+lambdade**(-2)))
			beta=abs(Zd*e**2/(LAMBDA*Ti*kb))
			u=numpy.sqrt(vdrift[0]**2+vdrift[1]**2+vdrift[2]**2)/vT
			z=abs(Zd)*e**2/(radd*Te*kb)
			tau=Te/Ti
			coloumblog=5.
			force=2*math.pi*radd**2*mi*ni0*vT*numpy.sqrt(2)*(vdrift-numpy.array(self.vel))*numpy.log(coloumblog)*(-phia)*chandra(machmag/numpy.sqrt(2))/(Ti*(machmag/numpy.sqrt(2))/Te)
			return force/self.m
		else:
			return numpy.array([0,0,0])#raise ValueError("Dust particle is not in sheath")	

	def EXBaccShukla(self,B=[0.,0.,0.014],machmag=0.):
		if self.pos[2]<sheathd and self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vT=numpy.sqrt(kb*Ti/mi)#Thermal speed ions
			vdrift=numpy.cross(self.radialfield()+self.sheathfield(),numpy.array(B))/(magB**2) #vdrift for constant vertical B field
			mach=numpy.array([vdrift[0]-self.vel[0],vdrift[1]-self.vel[1],vdrift[2]-self.vel[2]])/vT
			machmag=numpy.sqrt(mach[0]**2+mach[1]**2+mach[2]**2)
			vdriftmagsquare=vdrift[0]**2+vdrift[1]**2+vdrift[2]**2
			V_it=(vdriftmagsquare+8*kb*Ti/(mi*math.pi))**0.5
			sigmacoll=math.pi*radd**2*(1-2*e*phia/(mi*vdriftmagsquare))
			b0=radd*e*phia/(mi*vdriftmagsquare)
			bc=radd*numpy.sqrt(1.-2*e*phia/(mi*vdriftmagsquare))
			sigmacoul=2*math.pi*b0**2*numpy.log((b0**2+lambdade**2)/(b0**2+bc**2))
			force=ni0*mi*V_it*vdrift*(sigmacoul+sigmacoll)
			return force/self.m
		else:
			return numpy.array([0,0,0]) #raise ValueError("Dust particle is not in sheath")	


	def Btest(self):
		if self.pos[2]<sheathd:
			B=0.014
			vdrift=numpy.cross(self.radialfield()+self.multifields+self.sheathfield(),numpy.array([0,0,B]))/(B**2)
			return vdrift/numpy.sqrt(vdrift[0]**2+vdrift[1]**2+vdrift[2]**2)
		else:
			return numpy.array([0,0,0])
	def combinedrift(self,B=[0,0,0.014]):
		if self.Bswitch==True:
			magB=numpy.sqrt(B[0]**2+B[1]**2+B[2]**2)
			vdrift=numpy.cross(self.radialfield()+self.sheathfield(),numpy.array(B))/(magB**2)
			vperp=numpy.sqrt(2*kb*Ti/mi)
			vpa=vperp*numpy.sqrt(1./2.)
			rL=mi*vperp/(e*magB)
			r=numpy.sqrt(self.pos[0]**2+self.pos[1]**2+self.pos[2]**2)
			theta=numpy.arctan(self.pos[1]/self.pos[0])
			gradB=(-3.*mu0*numpy.sqrt(9*self.pos[2]**2+magBmom**2)/(4*math.pi*r**4))*numpy.array([numpy.cos(theta),numpy.sin(theta),0.])\
			+(mu0*18*self.pos[2]/(8*math.pi*r**3*numpy.sqrt(9*self.pos[2]**2+magBmom**2)))*numpy.array([0,0,1])
			totaldrift=(numpy.cross(B,gradB)/(magB**2))*(mi/(e*magB))*(vpa**2+0.5*vperp**2)
			return totaldrift
		else:
			return numpy.array([0,0,0])

	def dipoleB(self,r):
		if self.pos[2]<sheathd and self.Bswitch==True:
			magr=numpy.sqrt((self.pos[0]-r[0])**2+(self.pos[1]-r[1])**2+(self.pos[2]-r[2])**2)
			rhat=(numpy.array(self.pos)-numpy.array(r))/magr
			return (mu0*(3*numpy.dot(rhat,numpy.array(Bmom))*rhat-Bmom)/(4*math.pi*magr**3))
		else:
			return numpy.array([0,0,0])
	def thermalmotion(self):
		return 0.0001*numpy.array(self.acc)*numpy.array([numpy.random.random(),numpy.random.random(),numpy.random.random()]*numpy.array([(-1)**numpy.random.randint(2),(-1)**numpy.random.randint(2),(-1)**numpy.random.randint(2)]))


##Generate particles in their equilibrium position
filehandler = open("crystalpositions2K.obj",'rb') ##2k particles 5.5hrs to run 2500 iterations just for settling down
xinit=pickle.load(filehandler)
yinit=pickle.load(filehandler)
zinit=pickle.load(filehandler)
filehandler.close()
initpositions=[[i,j,k] for i,j,k in zip(xinit,yinit,zinit)]

## Quick plot to view initial configuration
# plt.figure()
# plt.scatter(xinit,yinit)
# plt.xlim((min(xinit),max(xinit)))
# plt.ylim((min(yinit),max(yinit)))
# plt.title("Initial positions of 2K particles")
# plt.show()

# #Tile the hexagonal shape
# lentile=max(xinit)-min(xinit)
# heightile=max(yinit)-min(yinit)
# xnumtiles=boxr/lentile
# ynumtiles=boxr/heightile
# print(lentile,heightile)

# lentile=max(xinit)-min(xinit) #Length of each tile
# heightile=max(yinit)-min(yinit) #height of each tile
# hexlayers=int(boxr/lentile) #Number of hexagonal layers

##Prepare modified B field
filehandler = open("modifiedfield.obj",'rb') ##2k particles 5.5hrs to run 2500 iterations just for settling down
gridr=pickle.load(filehandler)
gridz=pickle.load(filehandler)
Evalsr=pickle.load(filehandler)
Evalsz=pickle.load(filehandler)
rmax=pickle.load(filehandler)
separationsheath=pickle.load(filehandler)
separationhor1=pickle.load(filehandler)
separationhor2=pickle.load(filehandler)
firstpoint=pickle.load(filehandler)
filehandler.close()

def interpolate(r):
	x=r[0]
	y=r[1]
	r0=numpy.sqrt(r[0]**2+r[1]**2)
	z0=r[2]
	if z0>sheathd:
		return numpy.array([0,0,0])
	else:
		if r0<rmax:
			Eresult=numpy.array([0,0])
			indleft=(abs(r0)-abs(firstpoint))/separationhor1
			indlow=(abs(z0)-abs(firstpoint))/separationsheath
			s1=isclose((indleft**10) ** (1.0/10), int(indleft))
			s2=isclose((indlow**10) ** (1.0/10), int(indlow))
			# if s1==True and s2==True:
			# 	indlow=int(indlow)
			# 	indleft=int(indleft)
			# 	print("case 1")
			# 	gridEs=[[Evalsr[int(indlow)],Evalsz[int(indleft)]]]
			# 	gridpoints=[[gridr[int(indlow)],gridz[int(indleft)]]]
			# elif s1 ==True and s2 == False:
			# 	print("case 2")
			# 	indlow=int(indlow)
			# 	indleft=int(indleft)
			# 	gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow+1][indleft],Evalsz[indlow+1][indleft]]]
			# 	gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow+1][indleft],gridz[indlow+1][indleft]]]

			# elif s1 ==False and s2 == True:
			# 	print("case 3")
			# 	indleft=int(indleft)
			# 	indlow=int(indlow)
			# 	gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow][indleft+1],Evalsz[indlow][indleft+1]]]
			# 	gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow][indleft+1],gridz[indlow][indleft+1]]]

				# else:
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
		theta=numpy.arctan(abs(r[1]/r[0]))
		if r[0]>0 and r[1]>0: #Quadrant 1
			pass
		elif r[0]<0 and r[1]>0: #Quadrant 2
			theta=math.pi-theta
		elif r[0]<0 and r[1]<0: #Quadrant 3
			theta+=math.pi
		elif r[0]>0 and r[1]<0: #Quadrant 4
			theta=2*math.pi-theta
		else:
			print("Problem: dust grain at x and y positions", [r[0],r[1]])
		Efinal=[Efinal[0]*numpy.cos(theta),Efinal[0]*numpy.sin(theta),Efinal[1]]
		return numpy.array([Efinal[0],Efinal[1],0])*diminishfactor



##Create dictionary of particles from pickle object
diminishfactor=10**(0)
position=[]
numparticles=100
names=[]
for i in numpy.arange(numparticles):
	names.append('g%s'%i)
dustdict = {name: Dust(md,radd,lambdaD,phia,Zd*e,[0,0,0],[0,0,0],[0,0,0]) for name in names}
for i,j in zip(dustdict,numpy.arange(len(dustdict))):
	dustdict[i].pos=initpositions[j]

#Create pairs
##Check every pair of particle interaction
list1=list(dustdict.keys())
list2=list1
pairlist=list(itertools.product(list1,list2))
pairs=set()
for x,y in pairlist:
	if x!=y:
		pairs.add(frozenset((x,y)))
pairs=list(pairs)
for i in numpy.arange(len(pairs)):
	pairs[i]=list(pairs[i])

removelist=[]
for i in pairs:
	if i[0]==i[1]:
		removelist.append(i)
pairs=[i for i in pairs if i not in removelist]


##Interact and iterate 
iterationsB=20000
inititerations=1000
g9velcheck=[]
g9poscheck=[]
g9acccheck=[]


for i in tqdm(numpy.arange(inititerations)):
	pairsfinal=[]
	for b in pairs:
		if dustdict[b[0]].intergraind(dustdict[b[1]])==True:
			pairsfinal.append(b)
		else:
			pass #pairsfinal.append(b)
	for j in pairsfinal:
		interactfield=dustdict[j[0]].selffield(dustdict[j[1]])	
		dustdict[j[0]].selffieldmany(interactfield)
		dustdict[j[1]].selffieldmany(-interactfield)
	for k in dustdict:	
		dustdict[k].steptest(numpy.array(dustdict[k].multifields))
		# acc.append(dustdict[k].getselfacc())
		# vel.append(dustdict[k].getselfvel())
		position.append(dustdict[k].getselfpos())	

for l in dustdict:
	dustdict[l].Bswitch=True

for i in tqdm(numpy.arange(iterationsB)):
	pairsfinal=[]
	for b in pairs:
		if dustdict[b[0]].intergraind(dustdict[b[1]])==True:
			pairsfinal.append(b)
		else:
			pass #pairsfinal.append(b)
	# g9velcheck.append(dustdict['g29'].getselfvel())
	# g9poscheck.append(dustdict['g29'].getselfpos())
	# g9acccheck.append(dustdict['g29'].getselfacc())
	##Pick out the pairs that are less than 5 lambdaDs apart
	for j in pairsfinal:
		interactfield=dustdict[j[0]].selffield(dustdict[j[1]])	
		dustdict[j[0]].selffieldmany(interactfield)
		dustdict[j[1]].selffieldmany(-interactfield)
	for k in dustdict:	
		dustdict[k].steptest(numpy.array(dustdict[k].multifields)+interpolate(dustdict[k].getselfpos()))
		# acc.append(dustdict[k].getselfacc())
		# vel.append(dustdict[k].getselfvel())
		position.append(dustdict[k].getselfpos())



#Sort the positions of the particles
newx=[i[0] for i in position]
newy=[i[1] for i in position]
newz=[i[2] for i in position]


t = numpy.array(numpy.sort([list(numpy.arange(int(len(newx)/numparticles)))*numparticles]))[0]
df = pd.DataFrame({"time": t ,"x" : newx, "y" : newy, "z" : newz})

def update_graph(num):
	data=df[df['time']==num]
	point.set_data(data.x, data.y)
	point.set_3d_properties(data.z)		
	return point, 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=90., azim=90)
if min(newx)==max(newx):
	ax.set_xlim([-10*lambdaD,10*lambdaD])
else:
	ax.set_xlim([min(newx),max(newx)])
if min(newy)==max(newy):
	ax.set_ylim([-10*lambdaD,10*lambdaD])
else:
	ax.set_ylim([min(newy),max(newy)])
if min(newz)==max(newz):
	ax.set_zlim([-10*lambdaD,10*lambdaD])
else:
	ax.set_zlim([min(newz)-0.5*min(newz),max(newz)])

ax.set_xlim([-rmax*1.5,rmax*1.5])
ax.set_ylim([-rmax*1.5,rmax*1.5])


data=df[df['time']==0]
point, = ax.plot(data.x, data.y, data.z, linestyle="", marker=".")
#plt.cla()
# for i in [[0,0],[0,1],[1,0],[1,1]]:
# 	circx=((-1)**(i[0]))*numpy.arange(100)*boxr*0.9
# 	circy=((-1)**(i[1]))*numpy.sqrt(boxr**2-circx**2)
# 	ax.plot(circx,circy,'r-')
ani = matplotlib.animation.FuncAnimation(fig, update_graph, frames=iterationsB+inititerations,interval=1, blit=True)
ani.save('voidplusrotationdimfactor1.mp4', fps=30,extra_args=['-vcodec', 'libx264'])
plt.show()



