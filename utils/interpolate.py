import numpy
import math

from numpy import isclose


import constants as const


def interpolate(r, gridr, gridz, Evalsr, Evalsz, rmax, separationsheath, separationhor1, separationhor2, firstpoint):
    x = r[0]
    y = r[1]
    r0 = numpy.sqrt(r[0] ** 2 + r[1] ** 2)
    z0 = r[2]
    if z0 > const.sheathd:
        return numpy.array([0, 0, 0])
    else:
        if r0 < rmax:
            Eresult = numpy.array([0, 0])
            indleft = (abs(r0) - abs(firstpoint)) / separationhor1
            indlow = (abs(z0) - abs(firstpoint)) / separationsheath
            indlow = int(indlow)
            indleft = int(indleft)
            gridEs = [[Evalsr[indlow][indleft], Evalsz[indlow][indleft]],
                      [Evalsr[indlow + 1][indleft], Evalsz[indlow + 1][indleft]],
                      [Evalsr[indlow + 1][indleft + 1], Evalsz[indlow + 1][indleft + 1]],
                      [Evalsr[indlow][indleft + 1], Evalsz[indlow][indleft + 1]]]
            gridpoints = [[gridr[indlow][indleft], gridz[indlow][indleft]],
                          [gridr[indlow + 1][indleft], gridz[indlow + 1][indleft]],
                          [gridr[indlow + 1][indleft + 1], gridz[indlow + 1][indleft + 1]],
                          [gridr[indlow][indleft + 1], gridz[indlow][indleft + 1]]]

        elif r0 > rmax:
            Eresult = numpy.array([0, 0])
            indleft = (const.boxr - firstpoint - rmax) / separationhor2
            indlow = (z0 - firstpoint) / separationsheath
            s1 = isclose((indleft ** 3) ** (1.0 / 3), int(indleft))
            s2 = isclose((indlow ** 3) ** (1.0 / 3), int(indlow))
            indleft = int(indleft)
            indlow = int(indlow)
            if s1 != False or s2 != False:
                print("Case 4")
            gridEs = [[Evalsr[indlow][indleft], Evalsz[indlow][indleft]],
                      [Evalsr[indlow + 1][indleft], Evalsz[indlow + 1][indleft]],
                      [Evalsr[indlow + 1][indleft + 1], Evalsz[indlow + 1][indleft + 1]],
                      [Evalsr[indlow][indleft + 1], Evalsz[indlow][indleft + 1]]]
            gridpoints = [[gridr[indlow][indleft], gridz[indlow][indleft]],
                          [gridr[indlow + 1][indleft], gridz[indlow + 1][indleft]],
                          [gridr[indlow + 1][indleft + 1], gridz[indlow + 1][indleft + 1]],
                          [gridr[indlow][indleft + 1], gridz[indlow][indleft + 1]]]

        d1 = (abs(r0 - gridEs[0][0]))
        d4 = (abs(r0 - gridEs[0][1]))
        d2 = (abs(r0 - gridEs[2][0]))
        d3 = (abs(r0 - gridEs[2][1]))
        dhorsq = d1 ** 2 + d2 ** 2
        dvertsq = d3 ** 2 + d4 ** 2
        Etop = (d2 ** 2 / dhorsq) * numpy.array(gridEs[1]) + (d1 ** 2 / dhorsq) * numpy.array(gridEs[2])
        Ebottom = (d2 ** 2 / dhorsq) * numpy.array(gridEs[0]) + (d1 ** 2 / dhorsq) * numpy.array(gridEs[3])
        Efinal = (d3 ** 2 / dvertsq) * Ebottom + (d4 ** 2 / dvertsq) * Etop
        theta = numpy.arctan(abs(r[1] / r[0]))
        if r[0] > 0 and r[1] > 0:  # Quadrant 1
            pass
        elif r[0] < 0 and r[1] > 0:  # Quadrant 2
            theta = math.pi - theta
        elif r[0] < 0 and r[1] < 0:  # Quadrant 3
            theta += math.pi
        elif r[0] > 0 and r[1] < 0:  # Quadrant 4
            theta = 2 * math.pi - theta
        else:
            print("Problem: dust grain at x and y positions", [r[0], r[1]])
        Efinal = [Efinal[0] * numpy.cos(theta), Efinal[0] * numpy.sin(theta), Efinal[1]]
        return numpy.array([Efinal[0], Efinal[1], 0])