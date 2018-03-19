from msci.analysis.Gibsonmodification import ModifyBFieldAnalysis
from msci.utils.interpolate import interpolate
import msci.plots.Gibsonmodificatonplots as Gplot

from IPython import get_ipython
import msci.analysis.constants as const

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

g = ModifyBFieldAnalysis(electrondriftpos=[const.boxr, 0, 0])

g.Gibsonconstants()
g.electrondrift()
g.getrmax()
g.VoidVol()
g.voidQ()
g.inaccessibleanyp0()
g.inaccessiblecurrentp0()
g.getgrids()

# Plot inaccessibility region
# Gplot.plotinaccessibility(g.rnorm, g.znormpos, g.znormneg, g.rnormalize, g.znormalize, g.LAMBDA, save=False)

# Plot modified electric field inside region inaccessible to electrons
Gplot.plotmodifiedefield(g.rnormalize, g.LAMBDA, g.znormalize, g.gridr, g.gridz, g.Evalsr, g.Evalsz, g.Evalsradial,
                         g.Evalsradialz,
                         g.Evalsheathr, g.Evalsheath, g.chargepos, g.rmax)

# Plot interpolated electric field at dust equilibrium positions
# Gplot.testinterpolator(g, g.rmax, g.rnormalize, g.znormalize, g.LAMBDA, g.gridr, g.gridz, g.Evalsr, g.Evalsz,
#                        save=False)
