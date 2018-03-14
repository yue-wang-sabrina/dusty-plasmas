import math

import numpy
from scipy import special


class Dust:
    def __init__(self, const, mass=0, radius=0, lambdaD=0, phia=0, initcharge=0, initpos=[0, 0, 0], initvel=[0, 0, 0],
                 initacc=[0, 0, 0]):
        """
        Initiate the class

        :param mass: (float) mass of the dust particle
        :param radius:
        :param lambdaD:
        :param phia:
        :param initcharge:
        :param initpos:
        :param initvel:
        :param initacc:
        """
        self.m = mass
        self.rad = radius
        self.lambdaD = lambdaD
        self.charge = initcharge
        self.pos = initpos
        self.vel = initvel
        self.acc = initacc
        self.phia = phia
        self.sheathd = 10 * lambdaD
        self.multifields = numpy.array([0, 0, 0])
        self.Bswitch = False
        self.const = const

        # Extra testing stuff for diff numerical schemes and checks: can be deleted if not used later
        self.pos1 = initpos
        self.pos2 = initpos
        self.acc1 = initacc
        self.vel1 = initvel
        self.vel2 = initvel
        self.check = 0
        self.mach = numpy.array([0, 0, 0])

    def getselfpos(self):
        return self.pos

    def getselfvel(self):
        return self.vel

    def getselfacc(self):
        return self.acc

    def scatterwalls(self):
        if abs(numpy.array(self.pos[0]) + numpy.array(self.vel[0]) * self.const.dt + 0.5 * numpy.array(
                self.acc[0]) * self.const.dt ** 2) + abs(
                            numpy.array(self.pos[1]) + numpy.array(self.vel[1]) * self.const.dt + 0.5 * numpy.array(
                self.acc[1]) * self.const.dt ** 2) > self.const.boxr:
            print("Particle is scattered on cylinder walls at position", self.getselfpos())
            self.vel[0] *= -1
            self.vel[1] *= -1

        if numpy.array(self.pos[2]) + numpy.array(self.vel[2]) * dt + 0.5 * numpy.array(self.acc[2]) * dt ** 2 < 0:
            self.vel[2] = abs(self.vel[2])
            print("Particle is vertically scattered")

    def steptest(self, Emultifield=numpy.array([0, 0, 0]), method='integral'):
        # Test function for interactions, alter incooporate into step function
        # self.scatterwalls() #Check if hit walls

        # Add gravitational acceleration
        self.acc = numpy.array([0, 0, self.const.g])

        # Add linearly varying E field caused acceleration
        if method == 'integral':
            accsheath = self.sheathfield() * abs(self.charge / self.m)
            accradial = self.radialfield() * abs(self.charge / self.m)
            accB = self.EXBacchybrid(B=self.dipoleB(r=self.const.dipolepos))
            self.acc = numpy.array(self.acc) + accsheath + (
                Emultifield * abs(self.charge / self.m)) + self.damping() + accradial + accB
            self.vel = numpy.array(self.vel) + self.const.dt * numpy.array(self.acc)
            self.pos1 = self.pos  # Save previous position
            self.pos = numpy.array(self.pos) + numpy.array(self.vel) * self.const.dt + 0.5 * numpy.array(
                self.acc) * self.const.dt ** 2
            self.multifields = numpy.array([0, 0, 0])

        elif method == 'leap frog':
            # Leapfrog - unstable
            self.vel = numpy.array(self.vel2) + 2 * self.const.dt * numpy.array(self.acc1)
            self.pos = numpy.array(self.pos2) + 2 * self.const.dt * numpy.array(self.vel1)
            self.acc1 = self.acc
            self.vel2 = self.vel1
            self.vel1 = self.vel
            self.pos2 = self.pos1
            self.pos1 = self.pos

        elif method == 'rk4':
            # RK4 - also unstable?
            fv1 = numpy.array(self.acc1)
            fy1 = numpy.array(self.vel1)
            v1 = numpy.array(self.vel1) + 0.5 * self.const.dt * fv1
            y1 = numpy.array(self.pos1) + 0.5 * self.const.dt * fy1
            fv2 = numpy.array(self.acc1)
            fy2 = v1
            v2 = self.vel1 + 0.5 * self.const.dt * fv2
            y2 = self.pos1 + 0.5 * self.const.dt * fy2
            fv3 = fv2
            fy3 = v2
            v3 = self.vel1 + self.const.dt * fv3
            y3 = self.pos1 + self.const.dt * fy3
            fv4 = fv3
            fy4 = v3
            self.vel = self.vel1 + (1. / 6.) * self.const.dt * (fv1 + 2. * fv2 + 2. * fv3 + fv4)
            self.pos = self.pos1 + (1. / 6.) * self.const.dt * (fy1 + 2. * fy2 + 2. * fy3 + fy4)
            self.acc1 = self.acc
            self.vel1 = self.vel
            self.pos1 = self.pos

        else:
            raise Exception('Method does not exist!')

    def damping(self):
        magvel = numpy.sqrt(self.vel[0] ** 2 + self.vel[1] ** 2 + self.vel[2] ** 2)
        if magvel == 0:
            return numpy.array([0, 0, 0])
        else:
            acc = -1 * self.const.gamma * (numpy.array(self.vel) / magvel) * magvel ** 2
            return acc

    def verticalion(self):
        if self.pos[2] < self.const.sheathd:
            mach = 5.
            beta = abs(self.charge * self.const.e / (self.const.Ti * self.const.lambdade))
            return numpy.array([0, 0, (
                ((self.const.Ti / self.const.e) ** 2) * numpy.log(
                    self.const.lambdade * mach ** 2 / (beta * self.const.lambdadi)) * beta ** 2 / mach ** 2) / self.m])
        else:
            return numpy.array([0, 0, 0])

    def selffield(self, g2):
        # E field due to dust g2. NB: Assuming all particles have same size!
        # Energy between dressed grains (see Shukla book)
        # energy=(self.charge*grain2.charge/r)*math.exp(-r/self.lambdaD)*(1.-r/(2*self.lambdaD))
        # return energy
        # Electric field at some distance r from the dust particle
        d = numpy.sqrt((g2.pos[0] - self.pos[0]) ** 2 + (g2.pos[1] - self.pos[1]) ** 2 + (g2.pos[2] - self.pos[2]) ** 2)
        di = numpy.sqrt(
            (g2.pos1[0] - self.pos1[0]) ** 2 + (g2.pos1[1] - self.pos1[1]) ** 2 + (g2.pos1[2] - self.pos1[2]) ** 2)
        if di == 0.:
            pass  # print("di=0 PROBLEM")
        connecti = numpy.array([self.pos1[0] - g2.pos1[0], self.pos1[1] - g2.pos1[1], self.pos1[2] - g2.pos1[2]]) / di
        connectf = numpy.array([self.pos[0] - g2.pos[0], self.pos[1] - g2.pos[1], self.pos[2] - g2.pos[2]]) / d
        if connecti[0] != connectf[0] or connecti[1] != connectf[1] or connecti[2] != connectf[2]:
            pass  # raise ValueError("Particles passed each other")
        if d - 2 * self.rad < 0.:
            pass  # raise ValueError("Particles are inside each other!?")
        r = [self.pos[0] - g2.pos[0], self.pos[1] - g2.pos[1], self.pos[2] - g2.pos[2]]
        unitr = numpy.array(r) / d
        E = abs((self.phia / self.lambdaD) * numpy.exp(-(d - 2 * self.rad) / self.lambdaD))
        return E * unitr

    def selffieldmany(self, E):
        self.multifields = self.multifields + numpy.array(E)

    def sheathfield(self):
        if self.pos[2] >= self.sheathd and self.pos[2] > 0:
            return numpy.array([0, 0, 0])
        elif 0 < self.pos[2] < self.sheathd:
            self.V = self.const.electrodeV * (1. - self.pos[2] / self.sheathd) ** 2
            field = abs(2 * (1. - self.pos[2] / self.const.sheathd) * (self.const.electrodeV / self.const.sheathd))
            return numpy.array([0, 0, field])
        else:
            raise ValueError("dust fell out of cylinder!!With position", self.pos)

    def radialfield(self):  # Apply radial E field at the sheath
        if self.pos[2] < 0:
            raise ValueError("particle fell out of cylinder!!!!")
        elif self.pos[2] < self.const.sheathd:
            V = -self.const.wallV * (self.pos[0] ** 2 + self.pos[1] ** 2) / (self.const.boxr ** 2)
            mag = (self.const.wallV / (self.const.boxr ** 2)) * 2 * math.sqrt(self.pos[0] ** 2 + self.pos[1] ** 2)
            unitr = [-self.pos[0], -self.pos[1], 0]
            if (unitr[0] ** 2 + unitr[1] ** 2) == 0.:
                return numpy.array([0, 0, 0])
            else:
                Eapplied = numpy.array([-self.pos[0], -self.pos[1], 0]) * mag / numpy.sqrt(
                    self.pos[0] ** 2 + self.pos[1] ** 2)
                return numpy.array(Eapplied)
        else:
            return numpy.array([0, 0, 0])

    def intergraind(self, g2):
        if numpy.sqrt((self.pos[0] - g2.pos[0]) ** 2 + (self.pos[1] - g2.pos[1]) ** 2 + (
                    self.pos[2] - g2.pos[2]) ** 2) <= self.const.radinfluence:
            return True
        else:
            return False

    def EXBacc(self, B=[0, 0, 0.014], machmag=0., method='factor'):
        # Ion drag force due to Constant vertical B field causing ExB drift of ions

        if self.pos[2] < self.const.sheathd and self.Bswitch == True:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vT = numpy.sqrt(self.const.kb * self.const.Ti / self.const.mi)  # Thermal speed ions
            E = self.radialfield()  # +self.sheathfield()
            vdrift = numpy.cross(E, numpy.array(B)) / (magB ** 2)  # vdrift for self.constant vertical B field
            omega = abs(self.const.e * magB / self.const.mi)  # Ion cyclotron frequency
            tau = 0.01 / omega  # ion-neutral collision time
            theta = 0;
            if abs(self.pos[0]) == 0.:
                if self.pos[1]>0.:
                    theta = math.pi/2
                elif self.pos[1] < 0.:
                    theta = 3*math.pi/2
                else:
                    print("1Problem with theta calculation EXB drift")
            elif abs(self.pos[1]) == 0.:
                if self.pos[0] > 0:
                    theta = 0
                elif self.pos[0] < 0:
                    theta = math.pi
                else:
                    print("2Problem with theta calculation EXB drift")
            else:
                if self.pos[0] > 0. and self.pos[1] > 0.:
                    theta = numpy.arctan(self.pos[1]/self.pos[0])
                elif self.pos[0]>0. and self.pos[1] < 0.:
                    theta = 2*math.pi - numpy.arctan(abs(self.pos[1] /self.pos[0]))
                elif self.pos[0] < 0. and self.pos[1] > 0.:
                    theta = math.pi - numpy.arctan(abs(self.pos[1]/self.pos[0]))
                elif self.pos[0] < 0. and self.pos[1] < 0:
                    theta = math.pi + numpy.arctan(abs(self.pos[1]/self.pos[0]))
                else:
                    print("3Problem with theta calculation EXB drift")

            # Ion-neutral collision change to vdrift
            if method == 'factor':
                vdrift[0] *= (omega * tau) ** 2 / (1 + (omega * tau) ** 2)  # neglect v_r component
                vdrift[1] *= (omega * tau) ** 2 / (1 + (omega * tau) ** 2)
                vdrift[2] *= (omega * tau) ** 2 / (1 + (omega * tau) ** 2)

            # Ion-neutral collision drift velocity new D.D. Millar 1976
            elif method == 'derivation':
                # omega = abs(self.const.e * magB / self.const.mi)
                # tau = 0.05 / omega
                Br = numpy.sqrt(B[0] ** 2 + B[1] ** 2)
                Er = numpy.sqrt(E[0] ** 2 + E[1] ** 2)
                # k = ((self.const.e / self.const.mi) ** 2) * (B[2] * Er - Br * E[2])
                # p = Er - B[2] * k / omega ** 2
                # s = E[2] + Br * k / omega ** 2
                # drift0 = (-k / (omega ** 2 + 1. / tau ** 2))
                # drift1 = (self.const.e / (tau ** 2 * self.const.mi)) * (
                #     tau ** 3 * p - (B[2] * k / omega ** 4) * (-tau ** 3 * omega ** 2 / (1 + tau ** 2 * omega ** 2)))
                # drift2 = (self.const.e / (tau ** 2 * self.const.mi)) * (
                #     tau ** 3 * s + (Br * k / omega ** 4) * (-tau ** 3 * omega ** 3 / (1 + tau ** 2 * omega ** 2)))
                #
                # r = numpy.sqrt(self.pos[0] ** 2 + self.pos[1] ** 2)
                # theta = numpy.arctan(self.pos[2] / self.pos[0])
                # vdrift[0] = abs(drift1 * r * numpy.sin(theta)) * numpy.sign(vdrift[0])
                # vdrift[1] = abs(drift1 * r * numpy.cos(theta)) * numpy.sign(vdrift[1])

                Factor = (self.const.e * tau**2 / self.const.mi**2 * (1 + (omega*tau)**2)) * (Er * B[2] - E[2] * Br)
                vdrift[0] = Factor * -1 * numpy.sin(theta)
                vdrift[1] = Factor * numpy.cos(theta)
            else:
                raise Exception('Method does not exist!')

            # Calculate total acceleration on dust particle due to ion and neutral drag
            mach = numpy.array([vdrift[0] - self.vel[0], vdrift[1] - self.vel[1], vdrift[2] - self.vel[2]]) / vT
            machmag = numpy.sqrt(mach[0] ** 2 + mach[1] ** 2 + mach[2] ** 2)
            self.mach = mach
            if machmag <= 2:
                LAMBDA = numpy.sqrt(
                    1. / (numpy.exp(-machmag ** 2 / 2.) * self.const.lambdadi ** (-2) + self.const.lambdade ** (-2)))
                beta = abs(self.charge * self.const.e / (4 * math.pi * self.const.e0 * self.const.kb * self.const.Ti * LAMBDA))
                coloumblog = 5.
                forceS = (1. / 3.) * numpy.sqrt(32 * math.pi) * (
                    (self.const.Ti * self.const.kb / self.const.e) ** 2) * coloumblog * beta ** 2 * self.const.e0 * mach
                phi = abs(self.charge * self.const.e / (self.const.kb * self.const.Te * 4 * math.pi * self.const.e0 * self.rad))
                forceC = mach * vT * 4 * math.pi * self.rad ** 2 * self.const.ni0 * self.const.mi * numpy.sqrt(
                    self.const.kb * self.const.Te / (2 * math.pi * self.const.me)) * numpy.exp(-phi)
                forceN = -numpy.array(self.vel) * self.const.mi * numpy.sqrt(4 * math.pi) * self.const.ni0 * numpy.sqrt(
                    self.const.kb * self.const.Ti / (2 * math.pi * self.const.mi)) * math.pi * self.const.radd ** 2
                nn = 7. / (self.const.kb * self.const.Ti)  # Number density neutrals, 7Pa used by Saitou
                forceionneutraldrag = mach * self.const.mi * numpy.sqrt(4 * math.pi) * nn * numpy.sqrt(
                    self.const.kb * self.const.Ti / (2 * math.pi * self.const.mi)) * math.pi * self.const.radd ** 2
                forcetotal = forceS + forceC + forceN

                return forcetotal / self.m
            else:
                return mach * machmag * math.pi * self.const.Ti * self.const.kb * self.rad ** 2 * self.const.ni0 / self.m
        else:
            return numpy.array([0, 0, 0])  # raise ValueError("Dust particle is not in sheath")

    def EXBaccfortov(self, B=[0, 0, 0.014], machmag=0.):  ##Test other literature theories
        if self.pos[2] < self.const.sheathd and self.Bswitch == True:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vT = numpy.sqrt(self.const.kb * self.const.Ti / self.const.mi)  # Thermal speed ions
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (
                magB ** 2)  # vdrift for self.constant vertical B field
            mach = numpy.array([vdrift[0] - self.vel[0], vdrift[1] - self.vel[1], vdrift[2] - self.vel[2]]) / vT
            machmag = numpy.sqrt(mach[0] ** 2 + mach[1] ** 2 + mach[2] ** 2)

            LAMBDA = numpy.sqrt(1. / (numpy.exp(-machmag ** 2 / 2.) * self.const.lambdadi ** (-2) + self.const.lambdade ** (-2)))
            beta = abs(self.const.Zd * self.const.e ** 2 / (LAMBDA * self.const.Ti * self.const.kb))
            if machmag <= 1 and beta <= 13:
                coloumblog = 5.
                force = abs((1. / 3.) * numpy.sqrt(2. / math.pi) * (beta ** 2) * coloumblog) * mach

                return force / self.m
            elif machmag <= 1 and beta > 13:
                force = abs(
                    (2. / 3.) * numpy.sqrt(2 / math.pi) * ((numpy.log(beta)) ** 2 + 2 * numpy.log(beta) + 2)) * mach

                return force / self.m
            elif machmag > 1:
                force = abs((beta ** 2 * numpy.log(self.const.lambdade * machmag ** 2 / (self.const.lambdaD * beta)) * (
                    1. / machmag) ** 2)) * mach / (machmag)

                return force / self.m
            else:
                return numpy.array([0, 0, 0])
        else:
            return numpy.array([0, 0, 0])

    def EXBacckinetic(self, B=[0, 0, 0.014], machmag=0.):
        if self.pos[2] < self.const.sheathd and self.Bswitch == True:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vT = numpy.sqrt(self.const.kb * self.const.Ti / self.const.mi)  # Thermal speed ions
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (
                magB ** 2)  # vdrift for self.constant vertical B field
            mach = numpy.array([vdrift[0] - self.vel[0], vdrift[1] - self.vel[1], vdrift[2] - self.vel[2]]) / vT
            machmag = numpy.sqrt(mach[0] ** 2 + mach[1] ** 2 + mach[2] ** 2)

            LAMBDA = numpy.sqrt(1. / (numpy.exp(-machmag ** 2 / 2.) * self.const.lambdadi ** (-2) + self.const.lambdade ** (-2)))
            beta = abs(self.const.Zd * self.const.e ** 2 / (LAMBDA * self.const.Ti * self.const.kb))
            E = self.radialfield()
            incollfreq = self.const.e * numpy.sqrt(E[0] ** 2 + E[1] ** 2 + E[2] ** 2) / (
                self.const.mi * numpy.sqrt(vdrift[0] ** 2 + vdrift[1] ** 2 + vdrift[2] ** 2))
            l = vT / incollfreq
            LAMBDA = 1. / numpy.sqrt(self.const.lambdade ** (-2) + self.const.lambdadi ** (-2))
            if machmag <= 1:
                force = abs((1. / 3.) * numpy.sqrt(2 / math.pi) * beta ** 2 * (
                    numpy.log(1. / beta) + (1. / numpy.sqrt(2 * math.pi)) * K(LAMBDA / l))) * mach

                return force / self.m
            else:
                force = abs(numpy.sqrt(2. / math.pi) * beta ** 2 * numpy.log(4 * l * machmag / (LAMBDA * beta)) * (
                    1. / machmag)) * mach / machmag
                return force / self.m
        else:
            return numpy.array([0, 0, 0])  # raise ValueError("Dust particle is not in sheath")

    def EXBacchybrid(self, B=[0, 0, 0.014], machmag=0., method='factor',combinedrifts=False):
        if self.pos[2] < self.const.sheathd and self.Bswitch:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vT = numpy.sqrt(self.const.kb * self.const.Ti / self.const.mi)  # Thermal speed ions
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (
                magB ** 2)  # vdrift for self.constant vertical B field
            if combinedrifts:
                vdrift = vdrift + numpy.array(self.combinedrift(B))
            E = self.radialfield() + self.sheathfield()
            r = numpy.sqrt(self.pos[0] ** 2 + self.pos[1] ** 2)
            if abs(self.pos[0]) == 0.:
                if self.pos[1] > 0:
                    theta = math.pi/2
                elif self.pos[1]<0:
                    theta = 3*math.pi/2
                else:
                    print("Error in finding theta")
            elif abs(self.pos[1]) == 0:
                if self.pos[0] >0:
                    theta = 0;
                elif self.pos[0]<0:
                    theta =math.pi
            else:
                if self.pos[0] >0 and self.pos[1] >0:
                    theta = numpy.arctan(self.pos[1] / self.pos[0])
                elif self.pos[0]>0 and self.pos[1]<0:
                    theta = 2*math.pi - numpy.arctan(abs(self.pos[1] / self.pos[0]))
                elif self.pos[0]<0 and self.pos[1]>0:
                    theta = math.pi - numpy.arctan(abs(self.pos[1] / self.pos[0]))
                elif self.pos[0]<0 and self.pos[1] <0:
                    theta = math.pi + numpy.arctan(abs(self.pos[1] / self.pos[0]))

            # Just blindly multiplying by factor
            if method == 'factor':
                omega = abs(self.const.e * magB / self.const.mi)
                tau = 0.01 / abs(self.const.e * 0.014/self.const.mi)
                vdrift[0] = vdrift[0] * (omega * tau) ** 2 / (
                    1 + (omega * tau) ** 2) # Multiplication by factor just to make simulation faster
                vdrift[1] = vdrift[1]*(omega * tau) ** 2 / (1 + (omega * tau) ** 2)
                vdrift[2] = 0#abs(vdrift[2])*((omega * tau) ** 2 / (1 + (omega * tau) ** 2))
            # Ion-neutral collision drift velocity new D.D. Millar 1976
            elif method == 'derivation':
                omega = abs(self.const.e * magB / self.const.mi)
                tau = 0.01 / omega
                Br = numpy.sqrt(B[0] ** 2 + B[1] ** 2)
                Er = numpy.sqrt(E[0] ** 2 + E[1] ** 2)
                # k = ((self.const.e / self.const.mi) ** 2) * (B[2] * Er - Br * E[2])
                # p = Er - B[2] * k / omega ** 2
                # s = E[2] + Br * k / omega ** 2
                # drift0 = (-k / (omega ** 2 + 1. / tau ** 2))
                # drift1 = (self.const.e / (tau ** 2 * self.const.mi)) * (
                #     tau ** 3 * p - (B[2] * k / omega ** 4) * (-tau ** 3 * omega ** 2 / (1 + tau ** 2 * omega ** 2)))
                # drift2 = (self.const.e / (tau ** 2 * self.const.mi)) * (
                #     tau ** 3 * s + (Br * k / omega ** 4) * (-tau ** 3 * omega ** 3 / (1 + tau ** 2 * omega ** 2)))
                # vdrift[0] = abs((drift1) * r * numpy.sin(theta)) * numpy.sign(vdrift[0])
                # vdrift[1] = abs((drift1) * r * numpy.cos(theta)) * numpy.sign(vdrift[1])
                # vdrift[2] = abs(drift2) * numpy.sign(vdrift[2])
                factor = (self.const.e**2 * tau**2 /(self.const.mi**2 * (1+(omega*tau)**2))) * (Er * B[2] - E[2]*Br)
                vdrift[0] = factor * -numpy.sin(theta)
                vdrift[1] = factor * numpy.cos(theta)
                vdrift[2]=0
            else:
                raise Exception('Method does not exist!')
            mach = numpy.array([vdrift[0] - self.vel[0], vdrift[1] - self.vel[1], vdrift[2] - self.vel[2]]) / vT
            machmag = numpy.sqrt(mach[0] ** 2 + mach[1] ** 2 + mach[2] ** 2)
            LAMBDA = numpy.sqrt(1. / (numpy.exp(-machmag ** 2 / 2.) * self.const.lambdadi ** (-2) + self.const.lambdade ** (-2)))
            beta = abs(self.const.Zd * self.const.e ** 2 / (LAMBDA * self.const.Ti * self.const.kb))
            u = numpy.sqrt(vdrift[0] ** 2 + vdrift[1] ** 2 + vdrift[2] ** 2) / vT
            z = abs(self.const.Zd) * self.const.e ** 2 / (4 * math.pi * self.const.e0 * self.const.radd * self.const.Te * self.const.kb)
            tau = self.const.Te / self.const.Ti
            coloumblog = 5.
            force = numpy.sqrt(2 * math.pi) * self.const.radd ** 2 * self.const.ni0 * self.const.mi * vT ** 2 * \
                    (numpy.sqrt(math.pi / 2) * special.erf(u / numpy.sqrt(2)) *
                     (1 + u ** 2 + (1 - u ** (-2)) * (1 + 2 * z * tau) + 4 * z ** 2 * tau ** 2 * u ** (-2) * numpy.log(
                         coloumblog)) +
                     (
                         u ** (-1) * (
                             1 + 2 * z * tau + u ** 2 - 4 * z ** 2 * tau ** 2 * numpy.log(coloumblog)) * numpy.exp(
                             -u ** 2 / 2.))) * mach / machmag
            return force / self.m

        else:
            return numpy.array([0, 0, 0])  # raise ValueError("Dust particle is not in sheath")

    def EXBaccDUSTT(self, B=[0., 0., 0.014], machmag=0.):
        if self.pos[2] < self.const.sheathd and self.Bswitch:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vT = numpy.sqrt(self.const.kb * self.const.Ti / self.const.mi)  # Thermal speed ions
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (
                magB ** 2)  # vdrift for self.constant vertical B field
            mach = numpy.array([vdrift[0] - self.vel[0], vdrift[1] - self.vel[1], vdrift[2] - self.vel[2]]) / vT
            machmag = numpy.sqrt(mach[0] ** 2 + mach[1] ** 2 + mach[2] ** 2)

            LAMBDA = numpy.sqrt(1. / (numpy.exp(-machmag ** 2 / 2.) * self.const.lambdadi ** (-2) + self.const.lambdade ** (-2)))
            beta = abs(self.const.Zd * self.const.e ** 2 / (LAMBDA * self.const.Ti * self.const.kb))
            u = numpy.sqrt(vdrift[0] ** 2 + vdrift[1] ** 2 + vdrift[2] ** 2) / vT
            z = abs(self.const.Zd) * self.const.e ** 2 / (self.const.radd * self.const.Te * self.const.kb)
            tau = self.const.Te / self.const.Ti
            coloumblog = 5.
            force = 2 * math.pi * self.const.radd ** 2 * self.const.mi * self.const.ni0 * vT * numpy.sqrt(2) * (
                vdrift - numpy.array(self.vel)) * numpy.log(coloumblog) * (-self.const.phia) * self.const.chandra(
                machmag / numpy.sqrt(2)) / (
                        self.const.Ti * (machmag / numpy.sqrt(2)) / self.const.Te)
            return force / self.m
        else:
            return numpy.array([0, 0, 0])  # raise ValueError("Dust particle is not in sheath")

    def EXBaccShukla(self, B=[0., 0., 0.014], machmag=0.):
        if self.pos[2] < self.const.sheathd and self.Bswitch:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vT = numpy.sqrt(self.const.kb * self.const.Ti / self.const.mi)  # Thermal speed ions
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (
                magB ** 2)  # vdrift for self.constant vertical B field
            mach = numpy.array([vdrift[0] - self.vel[0], vdrift[1] - self.vel[1], vdrift[2] - self.vel[2]]) / vT
            machmag = numpy.sqrt(mach[0] ** 2 + mach[1] ** 2 + mach[2] ** 2)
            vdriftmagsquare = vdrift[0] ** 2 + vdrift[1] ** 2 + vdrift[2] ** 2
            V_it = (vdriftmagsquare + 8 * self.const.kb * self.const.Ti / (self.const.mi * math.pi)) ** 0.5
            sigmacoll = math.pi * self.const.radd ** 2 * (1 - 2 * self.const.e * self.const.phia / (self.const.mi * vdriftmagsquare))
            b0 = self.const.radd * self.const.e * self.const.phia / (self.const.mi * vdriftmagsquare)
            bc = self.const.radd * numpy.sqrt(1. - 2 * self.const.e * self.const.phia / (self.const.mi * vdriftmagsquare))
            sigmacoul = 2 * math.pi * self.const.b0 ** 2 * numpy.log(
                (self.const.b0 ** 2 + self.const.lambdade ** 2) / (self.const.b0 ** 2 + self.const.bc ** 2))
            force = self.const.ni0 * self.const.mi * V_it * vdrift * (sigmacoul + sigmacoll)
            return force / self.m
        else:
            return numpy.array([0, 0, 0])  # raise ValueError("Dust particle is not in sheath")

    def Btest(self):
        if self.pos[2] < self.const.sheathd:
            B = 0.014
            vdrift = numpy.cross(self.radialfield() + self.multifields + self.sheathfield(), numpy.array([0, 0, B])) / (
                B ** 2)
            return vdrift / numpy.sqrt(vdrift[0] ** 2 + vdrift[1] ** 2 + vdrift[2] ** 2)
        else:
            return numpy.array([0, 0, 0])

    def combinedrift(self, B=[0, 0, 0.014]):
        if self.Bswitch:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (magB ** 2)
            vperp = numpy.sqrt(2 * self.const.kb * self.const.Ti / self.const.mi)
            vpa = vperp * numpy.sqrt(1. / 2.)
            rL = self.const.mi * vperp / (self.const.e * magB)
            r = numpy.sqrt(self.pos[0] ** 2 + self.pos[1] ** 2 + self.pos[2] ** 2)
            rplane = numpy.sqrt(self.pos[0] ** 2 + self.pos[1] ** 2)
            if abs(self.pos[0]) == 0.:
                if self.pos[1] > 0:
                    theta = math.pi/2
                elif self.pos[1]<0:
                    theta = 3*math.pi/2
                else:
                    print("Error in finding theta")
            elif abs(self.pos[1]) == 0:
                if self.pos[0] >0:
                    theta = 0
                elif self.pos[0]<0:
                    theta =math.pi
            else:
                if self.pos[0] >0 and self.pos[1] >0:
                    theta = numpy.arctan(self.pos[1] / self.pos[0])
                elif self.pos[0]>0 and self.pos[1]<0:
                    theta = 2*math.pi - numpy.arctan(abs(self.pos[1] / self.pos[0]))
                elif self.pos[0]<0 and self.pos[1]>0:
                    theta = math.pi - numpy.arctan(abs(self.pos[1] / self.pos[0]))
                elif self.pos[0]<0 and self.pos[1] <0:
                    theta = math.pi + numpy.arctan(abs(self.pos[1] / self.pos[0]))

            tempz = self.const.dipolepos[2] + self.pos[2]
            sqrtterm = numpy.sqrt(9 * self.const.magBmom ** 2 * rplane ** 2 * tempz ** 2 + (2 * tempz ** 2 * self.const.magBmom - self.const.magBmom * rplane ** 2) ** 2)
            partial_r1 = (-5. / 2.) * (2 * rplane) * sqrtterm / r ** 7
            partial_r2 = 0.5 * (18 * self.const.magBmom ** 2 * rplane * tempz ** 2 - 4 * self.const.magBmom * rplane * (2 * self.const.magBmom * tempz ** 2 - self.const.magBmom * rplane ** 2)) / (r ** 5 * sqrtterm)
            partial_r = (self.const.mu0 / (4 * math.pi)) * (partial_r1 + partial_r2)
            partial_z1 = (-5. / 2.) * (2 * tempz) * sqrtterm / r ** 7
            partial_z2 = 0.5 * (18 * self.const.magBmom ** 2 * tempz * rplane ** 2 + 8 * self.const.magBmom * tempz * (2 * self.const.magBmom * tempz ** 2 - self.const.magBmom * rplane ** 2)) / (r ** 5 * sqrtterm)
            partial_z = (self.const.mu0 / (4 * math.pi)) * (partial_z1 + partial_z2)
            gradB = partial_r * numpy.array([numpy.cos(theta), numpy.sin(theta), 0.]) + partial_z * numpy.array([0, 0, 1])
            totaldrift = (numpy.cross(B, gradB) / (magB ** 2)) * (self.const.mi / (self.const.e * magB)) * (vpa ** 2 + 0.5 * vperp ** 2)
            return totaldrift
        else:
            return numpy.array([0, 0, 0])

    def EXBDRIFT(self,B):
        if self.pos[2] < self.const.sheathd and self.Bswitch:
            magB = numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
            vdrift = numpy.cross(self.radialfield() + self.sheathfield(), numpy.array(B)) / (
                magB ** 2)  # vdrift for self.constant vertical B field
            return vdrift
        else:
            return numpy.array([0,0,0])

    def vxBforce(self):
        force = self.charge*numpy.cross(self.dipoleB(self.const.dipolepos),self.vel)
        return force

    def dipoleB(self, r):
        if self.pos[2] < self.const.sheathd and self.Bswitch:
            magr = numpy.sqrt((self.pos[0] - r[0]) ** 2 + (self.pos[1] - r[1]) ** 2 + (self.pos[2] - r[2]) ** 2)
            rhat = (numpy.array(self.pos) - numpy.array(r)) / magr
            return (
                self.const.mu0 * (3 * numpy.dot(rhat, numpy.array(self.const.Bmom)) * rhat - self.const.Bmom) / (
                    4 * math.pi * magr ** 3))
        else:
            return numpy.array([0, 0, 0])

    def thermalmotion(self):
        return 0.0001 * numpy.array(self.acc) * numpy.array(
            [numpy.random.random(), numpy.random.random(), numpy.random.random()] * numpy.array(
                [(-1) ** numpy.random.randint(2), (-1) ** numpy.random.randint(2), (-1) ** numpy.random.randint(2)]))
