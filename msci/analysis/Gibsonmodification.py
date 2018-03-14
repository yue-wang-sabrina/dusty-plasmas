import numpy
import math
import constants as const
import scipy
import pickle
from numpy import isclose

from tqdm import tqdm
from msci.utils.utils import pointdipoleB
from scipy.optimize import fsolve


# Variables that can change when calculating the region of electron accessibility are electrondriftpos and rvalp0
# Note pointdipoleB function has the original dipole position at -0.003m.

class ModifyBFieldAnalysis:
    def __init__(self, constants=const, electrondriftpos=[const.boxr, 0, 0]):  # Position to calculate the drift in electron velocity
        self.const = constants
        self.electrondriftpos = electrondriftpos
        self.rvalp0 = numpy.linalg.norm(electrondriftpos)  # the r value at which we are calculating p0
        self.B = None
        self.Bmom = const.Bmom
        self.magBmom = const.magBmom
        self.Bmomhat = const.Bmomhat
        self.LAMBDA = None
        self.curlyM = None
        self.VeT = None
        self.critp0 = None
        self.edriftvel = None
        self.electrondriftvel = None
        self.rnorm = None
        self.znormneg = None
        self.znormpos = None
        self.rnormalize = None
        self.znormalize = None
        self.sheathfield = None
        self.voidcharge = None
        self.voidchargeguess = None
        self.chargepos = None
        self.rmax = None
        self.numbcharge = None
        self.gridr = None
        self.gridz = None
        self.Evalsheath = None
        self.Evalsheathr = None
        self.Evalsz = None
        self.Evalsr = None
        self.Evalsradial = None
        self.Evalsradialz = None
        self.voidvol = None
        self.Evalsmodifyall = None
        self.magEcharge = None
        self.separationsheath = None
        self.separationhor1 = None
        self.separationhor2 = None
        self.firstpoint = None

    def Gibsonconstants(self):  # Constants following from Joe's thesis, Run this when initialise
        self.LAMBDA = numpy.sqrt(self.magBmom * const.mu0 * const.e / (
            4 * math.pi * const.me * numpy.sqrt(
                const.kb * const.Te * 3 / const.me)))  # Distance at which electron is influenced from Joe's thesis
        self.curlyM = self.magBmom * const.mu0 / (4 * math.pi)
        self.VeT = numpy.sqrt(const.kb * const.Te * 3 / const.me)  # Thermal velocity of electrons
        self.critp0 = 2 * numpy.sqrt(
            const.e * const.me * self.VeT * self.curlyM)  # Critical p0 value determining regimes

    def electrondrift(
            self):  # Calculate E/B drift for electrons rough estimate Assuming B field is 0.014 vertical at all positions. Run this when initialise
        r = self.electrondriftpos
        omega = abs(const.e * 0.014 / const.me)  # Ion cyclotron frequency
        tau = 2. / omega  # electron-neutral collision time #Ratio Taken from konopkas paper on experimentalbasis crystals
        V = -const.wallV * (r[0] ** 2 + r[1] ** 2) / (const.boxr ** 2)
        mag = (const.wallV / (const.boxr ** 2)) * 2 * math.sqrt(r[0] ** 2 + r[1] ** 2)
        unitr = [-r[0], -r[1], 0]
        if (unitr[0] ** 2 + unitr[1] ** 2) == 0.:
            self.edriftvel = numpy.array([0, 0, 0])
        else:
            Eapplied = numpy.array([-r[0], -r[1], 0]) * mag / numpy.sqrt(r[0] ** 2 + r[1] ** 2)
            self.edriftvel = numpy.array(Eapplied) * (omega * tau) ** 2 / ((1 + (omega * tau) ** 2) * 0.014)
            self.electrondriftvel = numpy.linalg.norm(self.edriftvel)

    def getrmax(self):
        func = lambda rmax: 1. - const.e * self.curlyM / (
            rmax * (2 * numpy.sqrt(const.e * self.curlyM * const.me * self.VeT) + const.me * rmax * self.VeT))
        self.rmax = fsolve(func, 0.4 * self.LAMBDA)

    def VoidVol(self):
        intfunc = lambda r: 2 * math.pi * r * numpy.sqrt(
            (const.e * self.curlyM * r ** 2 / (
                2 * numpy.sqrt(const.e * self.curlyM * self.VeT * const.me) + const.me * r * self.VeT)) ** (
                2. / 3.) - r ** 2)
        voidvol = scipy.integrate.quad(intfunc, 0, self.rmax)
        self.voidvol = voidvol[0]

    def voidQ(self):
        Qfunc = lambda r: const.e * const.ne0 * numpy.exp(
            (const.e / (const.kb * const.Te)) * const.electrodeV * ((1. - (1. / const.sheathd) * numpy.sqrt(
                (const.e * self.curlyM * r ** 2 / (
                    2 * numpy.sqrt(const.e * self.curlyM * self.VeT * const.me) + const.me * self.VeT * r)) ** (
                    2. / 3.) - r ** 2)) ** 2 + (
                                                                        r / const.boxr) ** 2)) * 2 * math.pi * r * numpy.sqrt(
            (const.e * self.curlyM * r ** 2 / (
                2 * numpy.sqrt(const.e * self.curlyM * self.VeT * const.me) + const.me * self.VeT * r)) ** (
                2. / 3.) - r ** 2)
        self.voidcharge = scipy.integrate.quad(Qfunc, 0, self.rmax)[0]

    def Zpm(self, r, p, pm):  # Define z+/- function for determining accessible/inaccessible electron regions
        if math.isnan(const.e * self.curlyM * r ** 2 / (p + pm * const.me * r * (self.VeT + self.electrondriftvel)) ** (
                    2. / 3.)):
            return 10
        elif (const.e * self.curlyM * r ** 2 / (p + pm * const.me * r * (self.VeT + self.electrondriftvel))) ** (
                    2. / 3.) - r ** 2 <= 0:
            return 0
        else:
            return numpy.sqrt(
                (const.e * self.curlyM * r ** 2 / (p + pm * const.me * r * (self.VeT + self.electrondriftvel))) ** (
                    2. / 3.) - r ** 2)

    def inaccessiblecurrentp0(self):  # Get inaccessible zone at current p0
        r2 = numpy.arange(const.lambdaD, const.boxr, const.boxr / 10000)
        zinaccesspos = []
        zinaccessneg = []
        term1 = const.me * self.rvalp0 * self.electrondriftvel  # As boxr is the maximum r that the electron can come from
        term2 = const.e * const.mu0 * const.magBmom / (
            4 * math.pi * self.rvalp0)  # when electron is at boxr this term is at its minimum
        p0 = term1 + term2  # As p0 is always very tiny, P0=0 essentially
        for i in r2:
            zinaccessneg.append(self.Zpm(i, p0, -1))
            zinaccesspos.append(self.Zpm(i, p0, +1))

        self.rnorm = numpy.array(r2) / self.LAMBDA
        self.znormpos = numpy.array(zinaccesspos) / self.LAMBDA
        self.znormneg = numpy.array(zinaccessneg) / self.LAMBDA

    def inaccessibleanyp0(self):
        r = numpy.arange(0, const.boxr, const.boxr / 10000)
        zcrit = []
        for i in r:
            zcrit.append(self.Zpm(i, self.critp0, 1))
        self.rnormalize = numpy.array(r) / self.LAMBDA
        self.znormalize = numpy.array(zcrit) / self.LAMBDA

    def radialfield(self, r):  # Apply radial E field at the sheath, r = position of dust
        if r[1] < 0:
            raise ValueError("particle fell out of cylinder!!!!")
        elif r[1] < const.sheathd:
            V = -const.wallV * (r[0] ** 2 + r[1] ** 2) / (const.boxr ** 2)
            mag = (const.wallV / (const.boxr ** 2)) * 2 * math.sqrt(r[0] ** 2 + r[1] ** 2)
            unitr = [-r[0], 0]
            if (unitr[0] ** 2 + unitr[1] ** 2) == 0.:
                self.Eapplied = numpy.array([0, 0])
            else:
                self.Eapplied = numpy.array([-1, 0]) * mag
        else:
            self.Eapplied = numpy.array([0, 0])

    def SheathField(self, r):
        if r[1] < const.sheathd and r[1] >= 0:
            V = const.electrodeV * (1. - r[1] / const.sheathd) ** 2
            field = abs(2 * (1. - r[1] / const.sheathd) * (const.electrodeV / const.sheathd))
            self.sheathfield = numpy.array([0, field])
        else:
            self.sheathfield = numpy.array([0, 0])

    def pospos(self):  # positive charge positions
        separation = const.lambdaD / 2
        chargepos = []  # positive charge positions
        intcharge = int(self.voidvol)  # integer number of positive particles
        rows = int(self.rmax / separation)  # number of positive particles separated by lambdaD along the r axis
        for i in numpy.linspace(0, self.rmax, rows):
            columns = int(self.Zpm(i, self.critp0, 1) / separation)
            if columns == 0:
                chargepos.append([i, self.Zpm(i, self.critp0, 1)])
            else:
                for j in numpy.arange(columns):
                    chargepos.append([i,
                                      j * separation])  # number of partices along the z axis separated by lambdaD for specific r value
        self.chargepos = numpy.array(chargepos).T

    def MagEcharge(self, r):  # magnitude of E field due to charge
        normr = numpy.linalg.norm(r)
        rhat = r / normr
        self.magEcharge = abs((self.voidcharge / self.numbcharge) / (4 * math.pi * const.e0 * normr ** 2)) * rhat

    def gridcheck(self, chargepos):
        Evalsr = []
        Evalsz = []
        Evalsradial = []
        Evalsheath = []
        Evalsmodifyall = []
        for i in tqdm(numpy.arange(len(self.gridr))):
            Evalstempr = []
            Evalstempz = []
            Evalsradialtemp = []
            Evalsheathtemp = []
            Evalsmodifyalltemp = []
            for j in tqdm(numpy.arange(len(self.gridz[i]))):
                totalE = numpy.array([0, 0])
                self.radialfield([self.gridr[i][j], self.gridz[i][j]])
                Evalsradialtemp.append(self.Eapplied[0])
                self.SheathField([self.gridr[i][j], self.gridz[i][j]])
                Evalsheathtemp.append(self.sheathfield[1])
                for k in tqdm(numpy.arange(len(chargepos[0]))):
                    r = [chargepos[0][k] - self.gridr[i][j], chargepos[1][k] - self.gridz[i][j]]
                    if numpy.linalg.norm(r) == 0.:
                        pass  # Ignore the points that fall directly on a positive charge to avoid infinities
                    else:
                        self.MagEcharge(numpy.array(r))
                        totalE = totalE + numpy.array(self.magEcharge)
                Evalstempr.append(totalE[0])
                Evalstempz.append(totalE[1])
                Evalsmodifyalltemp.append(totalE)
            Evalsr.append(Evalstempr)
            Evalsz.append(Evalstempz)
            Evalsradial.append(Evalsradialtemp)
            Evalsheath.append(Evalsheathtemp)
            Evalsmodifyall.append(Evalsmodifyalltemp)
        self.Evalsr = Evalsr
        self.Evalsz = Evalsz
        self.Evalsradial = Evalsradial
        self.Evalsheath = Evalsheath
        self.Evalsmodifyall = Evalsmodifyall

    def checkedenofm(
            self):  # Run this function to calculate ratio of charge in the void that is displaced compared to outside void in rest of sheath
        print("Ratio of charge inside to outside region of inaccessibility=",
              self.voidcharge / (const.ne0 * const.e * (2 * math.pi * const.boxr * const.sheathd - self.voidvol)))
        print(
            "Approximate volume using half doughnut calculation is",
            2 * math.pi * (self.rmax / 2) * math.pi * self.rmax ** 2 / 2,
            "should be bigger than integral of my solid, which is", self.voidvol)

    def getgrids(self):
        self.voidchargeguess = self.voidvol * const.ne0 * const.e
        self.pospos()
        self.numbcharge = len(self.chargepos[
                                  0]) * 2 * math.pi * const.boxr / const.lambdaD  # The factor multiplied by is the number of 2D slices in the r-z plane there are in the 3D volume

        separationsheath = const.lambdaD  # separation distance between grid points
        self.separationsheath = separationsheath
        separationhor1 = const.lambdaD
        separationhor2 = const.lambdaD * 10
        self.separationhor1 = separationhor1
        self.separationhor2 = separationhor2
        firstpoint = separationsheath / 4
        self.firstpoint = firstpoint
        gridlinesheath = numpy.arange(firstpoint, const.sheathd * 1.5, separationsheath)
        gridlinehor = numpy.arange(firstpoint, self.rmax, separationhor1)
        gridlinehor2 = numpy.arange(self.rmax + firstpoint, const.boxr * 1.5, separationhor2)
        gridhor = list(gridlinehor) + list(gridlinehor2)
        self.gridr, self.gridz = numpy.meshgrid(gridhor, gridlinesheath)

        self.gridcheck(self.chargepos)
        self.Evalsradialz = [numpy.zeros(len(self.Evalsradial[0]))] * len(self.Evalsradial)
        self.Evalsheathr = [numpy.zeros(len(self.Evalsheath[0]))] * len(self.Evalsheath)

    def interpolate(self, r):
        r0 = r[0]
        z0 = r[1]
        if r0 < self.rmax:
            Eresult = numpy.array([0, 0])
            indleft = (r0 - self.firstpoint) / self.separationhor1
            indlow = (z0 - self.firstpoint) / self.separationsheath
            s1 = isclose((indleft ** 3) ** (1.0 / 3), int(indleft))
            s2 = isclose((indlow ** 3) ** (1.0 / 3), int(indlow))
            if s1 == True and s2 == True:
                print("case 1")
                gridEs = [[self.Evalsr[int(indlow)], self.Evalsz[int(indleft)]]]
                gridpoints = [[self.gridr[int(indlow)], self.gridz[int(indleft)]]]
            elif s1 == True and s2 == False:
                print("case 2")
                indlow = int(indlow)
                gridEs = [[self.Evalsr[indlow][indleft], self.Evalsz[indlow][indleft]],
                          [self.Evalsr[indlow + 1][indleft], self.Evalsz[indlow + 1][indleft]]]
                gridpoints = [[self.gridr[indlow][indleft], self.gridz[indlow][indleft]],
                              [self.gridr[indlow + 1][indleft], self.gridz[indlow + 1][indleft]]]

            elif s1 == False and s2 == True:
                print("case 3")
                indleft = int(indleft)
                gridEs = [[self.Evalsr[indlow][indleft], self.Evalsz[indlow][indleft]],
                          [self.Evalsr[indlow][indleft + 1], self.Evalsz[indlow][indleft + 1]]]
                gridpoints = [[self.gridr[indlow][indleft], self.gridz[indlow][indleft]],
                              [self.gridr[indlow][indleft + 1], self.gridz[indlow][indleft + 1]]]

            else:
                indlow = int(indlow)
                indleft = int(indleft)
                gridEs = [[self.Evalsr[indlow][indleft], self.Evalsz[indlow][indleft]],
                          [self.Evalsr[indlow + 1][indleft], self.Evalsz[indlow + 1][indleft]],
                          [self.Evalsr[indlow + 1][indleft + 1], self.Evalsz[indlow + 1][indleft + 1]],
                          [self.Evalsr[indlow][indleft + 1], self.Evalsz[indlow][indleft + 1]]]
                gridpoints = [[self.gridr[indlow][indleft], self.gridz[indlow][indleft]],
                              [self.gridr[indlow + 1][indleft], self.gridz[indlow + 1][indleft]],
                              [self.gridr[indlow + 1][indleft + 1], self.gridz[indlow + 1][indleft + 1]],
                              [self.gridr[indlow][indleft + 1], self.gridz[indlow][indleft + 1]]]
        elif r0 > self.rmax:
            Eresult = numpy.array([0, 0])
            indleft = (const.boxr - self.firstpoint - self.rmax) / self.separationhor2
            indlow = (z0 - self.firstpoint) / self.separationsheath
            s1 = isclose((indleft ** 3) ** (1.0 / 3), int(indleft))
            s2 = isclose((indlow ** 3) ** (1.0 / 3), int(indlow))
            indleft = int(indleft)
            indlow = int(indlow)
            if s1 != False or s2 != False:
                print("Case 4")
            gridEs = [[self.Evalsr[indlow][indleft], self.Evalsz[indlow][indleft]],
                      [self.Evalsr[indlow + 1][indleft], self.Evalsz[indlow + 1][indleft]],
                      [self.Evalsr[indlow + 1][indleft + 1], self.Evalsz[indlow + 1][indleft + 1]],
                      [self.Evalsr[indlow][indleft + 1], self.Evalsz[indlow][indleft + 1]]]
            gridpoints = [[self.gridr[indlow][indleft], self.gridz[indlow][indleft]],
                          [self.gridr[indlow + 1][indleft], self.gridz[indlow + 1][indleft]],
                          [self.gridr[indlow + 1][indleft + 1], self.gridz[indlow + 1][indleft + 1]],
                          [self.gridr[indlow][indleft + 1], self.gridz[indlow][indleft + 1]]]

        d1 = (abs(r0 - gridEs[0][0]))
        d4 = (abs(r0 - gridEs[0][1]))
        d2 = (abs(r0 - gridEs[2][0]))
        d3 = (abs(r0 - gridEs[2][1]))
        dhorsq = d1 ** 2 + d2 ** 2
        dvertsq = d3 ** 2 + d4 ** 2
        Etop = (d2 ** 2 / dhorsq) * numpy.array(gridEs[1]) + (d1 ** 2 / dhorsq) * numpy.array(gridEs[2])
        Ebottom = (d2 ** 2 / dhorsq) * numpy.array(gridEs[0]) + (d1 ** 2 / dhorsq) * numpy.array(gridEs[3])
        Efinal = (d3 ** 2 / dvertsq) * Ebottom + (d4 ** 2 / dvertsq) * Etop
        return Efinal, gridpoints

    def savetopickle(self, name, security=False):
        if security:
            filehandler = open(b"modifiedfield{}.obj".format(name),
                               'wb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
            pickle.dump(self.gridr, filehandler)
            pickle.dump(self.gridz, filehandler)
            pickle.dump(self.Evalsr, filehandler)
            pickle.dump(self.Evalsz, filehandler)
            pickle.dump(self.rmax, filehandler)
            pickle.dump(self.separationsheath, filehandler)
            pickle.dump(self.separationhor1, filehandler)
            pickle.dump(self.separationhor2, filehandler)
            pickle.dump(self.firstpoint, filehandler)
            filehandler.close()
            print("Saved modified fields to {}".format("modifiedbfield{}.obj".format(name)))
