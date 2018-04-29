from msci.analysis.Gibsonmodification import ModifyBFieldAnalysis
from msci.utils.interpolate import interpolate
import msci.plots.Gibsonmodificatonplots as Gplot
from IPython import get_ipython
import msci.analysis.constants as const
import matplotlib.pyplot as plt
ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')
import numpy
import math


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



#Plot inaccessibility region for report
# f, axarr = plt.subplots(2, 2)
# rnorm=[]
# znormpos=[]
# znormneg=[]
# thermalvel=numpy.sqrt(8*const.kb*const.Te/(math.pi*const.me))
# initialvels=[thermalvel*0.5,thermalvel*6,thermalvel*12,thermalvel*18]
# rnormalize =[]
# znormalize = []
# LAMBDA = []
# p0list=[]
# for i in numpy.arange(4):
#     g = ModifyBFieldAnalysis(electrondriftpos=[const.boxr, 0, 0],initialvel=initialvels[i])
#     g.Gibsonconstants()
#     g.electrondrift()
#     g.getrmax()
#     g.VoidVol()
#     g.voidQ()
#     g.inaccessibleanyp0()
#     g.inaccessiblecurrentp0()
#     rnorm.append(g.rnorm)
#     znormpos.append(g.znormpos)
#     znormneg.append(g.znormneg)
#     rnormalize.append(g.rnormalize)
#     znormalize.append(g.znormalize)
#     LAMBDA.append(g.LAMBDA)
#     p0list.append(g.p0)
# j=0
# # axarr[0, 0].plot(rnorm[j], znormpos[j], 'b.-', label='Region below inaccessble to electrons')  # for this p0, zinaccesspos')
# # axarr[0, 0].plot(rnorm[j], znormneg[j], 'k.-', label='Region above inaccessible to electrons')  # for this p0, zinaccessneg')
# axarr[0, 0].fill_between(rnorm[j], znormpos[j], znormneg[j], where=znormneg[j] > znormpos[j], interpolate=True, color='pink')
# axarr[0, 0].plot(rnormalize[j], znormalize[j], 'r-')#,label="Region below is inaccessible for electrons of all velocities")  # p0 below curve")
# axarr[0, 0].plot(rnormalize[j], numpy.ones(len(rnormalize[j])) * const.sheathd / LAMBDA[j], 'm-', label="Region of interest is below this line (Sheath top)")
# axarr[0, 0].fill_between(rnormalize[j], 0, znormalize[j], interpolate=True, color='red')
# # axarr[0, 0].set_xlabel(r"$\frac{r}{\Lambda}$", fontsize=35)
# # axarr[0, 0].set_ylabel(r"$\frac{z}{\Lambda}$", fontsize=35)
# axarr[0, 0].set_xlim([0, 4])
# axarr[0, 0].set_ylim([0, 1])
# axarr[0,0].set_title(r"$p_0$=%s"%"%0.1f"%p0list[j])
#
# j=1
# # axarr[0, 1].plot(rnorm[j], znormpos[j], 'b.-', label='Region below inaccessble to electrons')  # for this p0, zinaccesspos')
# # axarr[0, 1].plot(rnorm[j], znormneg[j], 'k.-', label='Region above inaccessible to electrons')  # for this p0, zinaccessneg')
# axarr[0, 1].fill_between(rnorm[j], znormpos[j], znormneg[j], where=znormneg[j] > znormpos[j], interpolate=True, color='pink')
# axarr[0, 1].plot(rnormalize[j], znormalize[j], 'r-')#,label="Region below is inaccessible for electrons of all velocities")  # p0 below curve")
# axarr[0, 1].plot(rnormalize[j], numpy.ones(len(rnormalize[j])) * const.sheathd / LAMBDA[j], 'm.-', label="Top of sheath")
# axarr[0, 1].fill_between(rnormalize[j], 0, znormalize[j], interpolate=True, color='red')
# # axarr[0, 1].set_xlabel(r"$\frac{r}{\Lambda}$", fontsize=35)
# # axarr[0, 1].set_ylabel(r"$\frac{z}{\Lambda}$", fontsize=35)
# axarr[0, 1].set_xlim([0, 4])
# axarr[0, 1].set_ylim([0, 1])
# axarr[0,1].set_title(r"$p_0$=%s"%"%0.1f"%p0list[j])
#
# j=2
# # axarr[1, 0].plot(rnorm[j], znormpos[j], 'b.-', label='Region below inaccessble to electrons')  # for this p0, zinaccesspos')
# # axarr[1, 0].plot(rnorm[j], znormneg[j], 'k.-', label='Region above inaccessible to electrons')  # for this p0, zinaccessneg')
# axarr[1, 0].fill_between(rnorm[j], znormpos[j], znormneg[j], where=znormneg[j] > znormpos[j], interpolate=True, color='pink')
# axarr[1, 0].plot(rnormalize[j], znormalize[j], 'r-')#,label="Region below is inaccessible for electrons of all velocities")  # p0 below curve")
# axarr[1, 0].plot(rnormalize[j], numpy.ones(len(rnormalize[j])) * const.sheathd / LAMBDA[j], 'm-', label="Top of sheath")
# axarr[1, 0].fill_between(rnormalize[j], 0, znormalize[j], interpolate=True, color='red')
# # axarr[1, 0].set_xlabel(r"$\frac{r}{\Lambda}$", fontsize=35)
# # axarr[1, 0].set_ylabel(r"$\frac{z}{\Lambda}$", fontsize=35)
# axarr[1, 0].set_xlim([0, 4])
# axarr[1, 0].set_ylim([0, 1])
# axarr[1,0].set_title(r"$p_0$=%s"%"%0.1f"%p0list[j])
#
# j=3
# # axarr[1, 1].plot(rnorm[j], znormpos[j], 'b.-', label='Region below inaccessble to electrons')  # for this p0, zinaccesspos')
# # axarr[1, 1].plot(rnorm[j], znormneg[j], 'k.-', label='Region above inaccessible to electrons')  # for this p0, zinaccessneg')
# axarr[1, 1].fill_between(rnorm[j], znormpos[j], znormneg[j], where=znormneg[j] > znormpos[j], interpolate=True, color='pink')
# axarr[1, 1].plot(rnormalize[j], znormalize[j], 'r-')#,label="Region below is inaccessible for electrons of all velocities")  # p0 below curve")
# axarr[1, 1].fill_between(rnormalize[j], 0, znormalize[j], interpolate=True, color='red')
# axarr[1, 1].plot(rnormalize[j], numpy.ones(len(rnormalize[j])) * const.sheathd / LAMBDA[j], 'm-', label="Top of sheath")
# # axarr[1, 1].set_xlabel(r"$\frac{r}{\Lambda}$", fontsize=35)
# # axarr[1, 1].set_ylabel(r"$\frac{z}{\Lambda}$", fontsize=35)
# axarr[1, 1].set_xlim([0, 4])
# axarr[1, 1].set_ylim([0, 1])
# axarr[1,1].set_title(r"$p_0$=%s"%"%0.1f"%p0list[j])
#
# # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
# plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
# plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
# from pylab import rcParams
# rcParams['figure.figsize'] = 13, 12
# f.text(0.5, 0.04, r"$\frac{r}{\Lambda}$", ha='center', va='center',fontsize=35)
# f.text(0.06, 0.5, r"$\frac{z}{\Lambda}$", ha='center', va='center', rotation='vertical',fontsize=35)
# f.savefig('forbidden.pdf')
# f.show()