import numpy
import math
import matplotlib.pyplot as plt

from tqdm import tqdm

import constants as const
from particles.ions import Ion

dt = 0.0000001

class IonAnalysis:
    def __init__(self):
        self.position = None
        self.iterations = None
        self.direction = None
        self.method = None
        self.ion = None
        self.particledrift = None

    def runsim(self,iterations,method): # Run a single ion motion through E and B field and make a plot
        self.iterations = iterations
        self.method = method
        vinit = numpy.sqrt(const.kb * const.Ti / const.mi)
        self.ion = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])

        position = []
        velocity = []

        if method == 'RK4':
            for i in numpy.arange(self.iterations):
                self.ion.updateRK4(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        elif method == 'EF':
            for i in numpy.arange(self.iterations):
                self.ion.updateeulerforward(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        elif method == 'EB':
            for i in numpy.arange(self.iterations):
                self.ion.updateeulerback(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        elif method == 'BORIS':
            for i in numpy.arange(self.iterations):
                self.ion.updateboris(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        else:
            raise ValueError("{} does not exist".format(method))

        self.position = numpy.array(position)
        self.direction = numpy.cross(self.ion.constE(), self.ion.constB()) / (
            numpy.linalg.norm(numpy.cross(self.ion.constE(), self.ion.constB())))

    def plotrunsim(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.position[:,0], self.position[:,1], self.position, 'r-', label="{}".format(self.method))
        plt.legend()
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_zlabel("z (m)")
        fig.show()

    def checkdrift(self,method,iterations=numpy.arange(200, 20000, 500)):
        drift = []
        for i in tqdm(iterations):
            self.runsim(i,method=method)
            EXBdrift = numpy.linalg.norm(self.ion.constE()) / numpy.linalg.norm(self.ion.constB())
            iterationtime = dt * i
            totald = EXBdrift * iterationtime
            theorydrift = numpy.linalg.norm(self.direction * totald)
            particlepos = numpy.linalg.norm(
                numpy.array([self.position[-1, 0], self.position[-1,1], self.position[-1,2]]) * self.direction
            )
            drift.append([theorydrift, particlepos])
        self.particledrift = numpy.array(numpy.mat(drift).T)

    def plotcheckdrift(self):
        plt.figure()
        plt.scatter(self.iterations, self.particledrift[1] / self.particledrift[0])
        plt.ylabel("r_simulation/r_theory")
        plt.xlabel("Iteration number")
        plt.title("Ratio of simulated drift position vs theoretical drift position")
        plt.show()

    def

