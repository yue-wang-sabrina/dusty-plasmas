# This file is mainly for looking at B field effects of crystals that are already formed.

from msci.analysis.analysis_ion import IonAnalysis
from msci.plots.ionplots import *
from msci.utils.utils import bootstrap
from msci.plots import ionplots

from IPython import get_ipython
import msci.analysis.constants as const

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

# Define constants wanted for analysis
tau = 10 ** (-5)
omega = const.e * 0.014 / const.mi

# ioneffect1 = IonAnalysis(tau=tau, omega=omega)

#
# # Analysis different tau values
#
# taulist = [0.01 / ioneffect1.omega, 0.05 / ioneffect1.omega]
# ioneffect1.changetau(taulist, 1000, 1, 'testchangetau.obj', True)
# plotratiodriftmanytau('testchangetau.obj', len(taulist), ioneffect1.omega, ioneffect1.tau)
#
# # Analysis different time values
#
# iterationlist = [0.001 / ioneffect1.dt]
# ioneffect1.changetime(tau=0.275 / (const.e * 0.014 / const.mi), iterationlist=iterationlist, runs=1,
#                       filename='testchangetime.obj', security=True)
# plotratiomanytime('testchangetime.obj', len(iterationlist), ioneffect1.omega, ioneffect1.tau)
#
# # Analysis old drift from long runs of different omega*tau for time = 3s

# plotolddrifts(fontsize = 20,figsize0=50, figsize1=50)

