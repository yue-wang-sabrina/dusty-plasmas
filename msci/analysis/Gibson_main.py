from msci.analysis.Gibsonmodification import ModifyBFieldAnalysis
import numpy
import math

Bmomstrength = numpy.arange(0.0121,0.018,0.00045)

def norm(x):
    return numpy.sqrt(x[0]**2+x[1]**2+x[2]**2)

for i in numpy.arange(len(Bmomstrength)):
    gibtemp = ModifyBFieldAnalysis()
    gibtemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * Bmomstrength[i] / gibtemp.const.mu0) * numpy.array([0, 0, 1])
    gibtemp.const.magBmom = norm(gibtemp.const.Bmom)
    gibtemp.const.Bmomhat = numpy.array(gibtemp.const.Bmom) / gibtemp.const.magBmom
    gibtemp.Gibsonconstants()
    gibtemp.electrondrift()
    gibtemp.getrmax()
    gibtemp.VoidVol()
    gibtemp.voidQ()
    gibtemp.inaccessiblecurrentp0()
    gibtemp.inaccessibleanyp0()
    gibtemp.pospos()
    gibtemp.getgrids()
    gibtemp.gridcheck(gibtemp.chargepos)
    gibtemp.savetopickle(Bmomstrength[i], decimals=5, security=False)
