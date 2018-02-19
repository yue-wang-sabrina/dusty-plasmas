import numpy

import constants as const


class Ion:
    def __init__(self, pos, vel, acc):
        self.dt = 0.0000001
        self.pos = numpy.array(pos)
        self.vel = numpy.array(vel)
        self.acc = numpy.array(acc)
        self.pos1 = numpy.array(pos)
        self.vel1 = numpy.array(vel)
        self.charge = const.e
        self.dr = 0

    def constB(self):
        return numpy.array([0, 0.014, 0])

    def constE(self, E=[0, 0, 1]):
        return numpy.array(E)

    def updateeulerforward(self, B, E=[0, 0, 1]):
        # Euler forward - should be unstable!
        self.pos += self.dt * self.vel
        self.vel = ((self.charge / const.mi) * numpy.cross(self.vel, B) + numpy.array(E) * (
            self.charge / const.mi
        )) * const.dt + self.vel

    def updateeulerback(self, B, E=[0, 0, 1]):  # Also unstable?
        self.pos -= self.dt * self.vel
        self.vel -= (
            (self.charge / const.mi) * numpy.cross(self.vel, B) * self.dt + numpy.array(E) * (self.charge / const.mi)
        )

    def updateRK4(self, B, E=[0, 0, 1]):
        # RK4 integration
        fv1 = (self.charge / const.mi) * numpy.cross(self.vel1, B) + numpy.array(E) * (self.charge / const.mi)
        fy1 = self.vel1
        v1 = self.vel1 + 0.5 * self.dt * fv1
        y1 = self.pos1 + 0.5 * self.dt * fy1
        fv2 = (self.charge / const.mi) * numpy.cross(v1, B) + numpy.array(E) * (self.charge / const.mi)
        fy2 = v1
        v2 = self.vel1 + 0.5 * self.dt * fv2
        y2 = self.pos1 + 0.5 * self.dt * fy2
        fv3 = (self.charge / const.mi) * numpy.cross(v2, B) + numpy.array(E) * (self.charge / const.mi)
        fy3 = v2
        v3 = self.vel1 + self.dt * fv3
        y3 = self.pos1 + self.dt * fy3
        fv4 = (self.charge / const.mi) * numpy.cross(v3, B) + numpy.array(E) * (self.charge / const.mi)
        fy4 = v3
        self.vel = self.vel1 + (1. / 6.) * self.dt * (fv1 + 2. * fv2 + 2. * fv3 + fv4)
        self.pos = self.pos1 + (1. / 6.) * self.dt * (fy1 + 2. * fy2 + 2. * fy3 + fy4)
        self.vel1 = self.vel
        self.dr = numpy.array(self.pos) - numpy.array(self.pos1)
        self.pos1 = self.pos

    def updateboris(self, B):
        pass

    def getselfvel(self):
        return self.vel

    def getselfpos(self):
        return self.pos
