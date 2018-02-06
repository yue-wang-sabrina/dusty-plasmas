import crystalformation
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
##Test generating list of particles
##Preparing initialization of dust particles
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
boxr=10000*lambdaD #cylinder radius
boxz=0.001 #cylinder height
g=-9.8 #gravitational acceleration
dt=0.0001 #Much longer than charging time which is order nanoseconds
sheathd=10*lambdaD
electrodeV=abs((kb*Te/(2*e))*(numpy.log(2*math.pi*me/mi))) #potential at electrode
wallV=electrodeV #cylindrical sides of wall same potential
radinfluence=5*lambdaD
phia=crystalformation.phia
Zd=crystalformation.Zd


numparticles=100
names=[]
initpositions=[]
for i in numpy.arange(numparticles):
	names.append('g%s'%(i+1))
dustdict = {'g0':crystalformation.Dust(md,radd,lambdaD,phia,Zd*e,[lambdaD*numpy.random.randint(-6,6),lambdaD*numpy.random.randint(-6,6),boxz],[0,0,0],[0,0,0])}#{name: Dust(md,radd,lambdaD,phia,Zd*e,[lambdaD*numpy.random.randint(-6,6),lambdaD*numpy.random.randint(-6,6),0.00038257340070632547],[0,0,0],[0,0,0]) for name in names}
simtime=7 #simulation time in seconds
iterations=simtime/dt
timedrop=[] #Time at which the ith particle is dropped
for i in numpy.arange(numparticles):
	timedrop.append(i*0.005)
position=[]
vel=[]
acc=[]

##Drop 1 by 1
currentime=0.
def findpairs(dustdict):
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
	pairsfinal=[]
	for b in pairs:
		if dustdict[b[0]].intergraind(dustdict[b[1]])==True:
			pairsfinal.append(b)
		else:
			pass #pairsfinal.append(b)
	return pairsfinal
starttime=time.time()
while currentime<simtime:
	acctemp=[]
	veltemp=[]
	positiontemp=[]
	#pairsfinal=findpairs(dustdict)
	pairs=findpairs(dustdict)
	pairsfinal=[]
	for b in pairs:
		if dustdict[b[0]].intergraind(dustdict[b[1]])==True:
			pairsfinal.append(b)
		else:
			pass #pairsfinal.append(b)
	if len(timedrop)!=0:
		if abs(currentime-timedrop[0])<0.0000001:
			del timedrop[0]
			dustdict.update({names[0]:crystalformation.Dust(md,radd,lambdaD,phia,Zd*e,[lambdaD*numpy.random.randint(-6,6),lambdaD*numpy.random.randint(-6,6),boxz],[0,0,0],[0,0,0])})
			del names[0]
			#pairsfinal=findpairs(dustdict)
			pairs=findpairs(dustdict)
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
				dustdict[k].steptest(dustdict[k].multifields)
				acctemp.append(dustdict[k].getselfacc())
				veltemp.append(dustdict[k].getselfvel())
				positiontemp.append(dustdict[k].getselfpos())	
		else:
			for j in pairsfinal:
				interactfield=dustdict[j[0]].selffield(dustdict[j[1]])	
				dustdict[j[0]].selffieldmany(interactfield)
				dustdict[j[1]].selffieldmany(-interactfield)
			for k in dustdict:
				dustdict[k].steptest(dustdict[k].multifields)
				acctemp.append(dustdict[k].getselfacc())
				veltemp.append(dustdict[k].getselfvel())
				positiontemp.append(dustdict[k].getselfpos())	
	else:
		for j in pairsfinal:
			interactfield=dustdict[j[0]].selffield(dustdict[j[1]])	
			dustdict[j[0]].selffieldmany(interactfield)
			dustdict[j[1]].selffieldmany(-interactfield)
		for k in dustdict:
			dustdict[k].steptest(dustdict[k].multifields)
			acctemp.append(dustdict[k].getselfacc())
			veltemp.append(dustdict[k].getselfvel())
			positiontemp.append(dustdict[k].getselfpos())	
	currentime+=dt
	acc.append(acctemp)
	vel.append(veltemp)
	position.append(positiontemp)
stoptime=time.time()

print("Time taken for computation=",stoptime-starttime," seconds")

newx=[]
newy=[]
newz=[]
for i in numpy.arange(len(position)):
	tempx=[]
	tempy=[]
	tempz=[]
	for j in numpy.arange(len(position[i])):
		tempx.append(position[i][j][0])
		tempy.append(position[i][j][1])
		tempz.append(position[i][j][2])
	newx.append(tempx)
	newy.append(tempy)
	newz.append(tempz)

lenperstep=[len(l) for l in newx]
t=[numpy.ones(k)*ind for k,ind in zip(lenperstep,numpy.arange(iterations+1))]
t=numpy.array([float(j) for i in t for j in i])
newx=numpy.array([j for i in newx for j in i])
newy=numpy.array([j for i in newy for j in i])
newz=numpy.array([j for i in newz for j in i])

df=pd.DataFrame({'time':t,"x":newx,"y":newy,"z":newz})

def update_graph(num):
	data=df[df['time']==num]
	point.set_data(data.x, data.y)
	point.set_3d_properties(data.z)
	return point, 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

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
	ax.set_zlim([min(newz),max(newz)])

data=df[df['time']==0]
point, = ax.plot(data.x, data.y, data.z, linestyle="", marker="o")
ani = matplotlib.animation.FuncAnimation(fig, update_graph,frames=int(iterations), interval=40, blit=True,repeat=False)
#ani.save('video.mp4', fps=30,extra_args=['-vcodec', 'libx264'])
plt.show()



