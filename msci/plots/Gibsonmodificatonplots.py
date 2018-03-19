# # ##Plot inaccessibility regions
import matplotlib
import matplotlib.pyplot as plt
import numpy
import msci.analysis.constants as const
import pickle

font = {'family': 'normal',
        'weight': 'bold',
        'size': 20}

matplotlib.rc('font', **font)


def plotinaccessibility(rnorm, znormpos, znormneg, rnormalize, znormalize, LAMBDA, save=False):
    fig, ax = plt.subplots(1, 1)
    ax.plot(rnorm, znormpos, 'b.-', label='Region below inaccessble to electrons')  # for this p0, zinaccesspos')
    ax.plot(rnorm, znormneg, 'k.-', label='Region above inaccessible to electrons')  # for this p0, zinaccessneg')
    ax.fill_between(rnorm, znormpos, znormneg, where=znormneg > znormpos, interpolate=True, color='pink')
    ax.plot(rnormalize, znormalize, 'r.-',
            label="Region below is inaccessible for electrons of all velocities")  # p0 below curve")
    ax.plot(rnormalize, numpy.ones(len(rnormalize)) * const.sheathd / LAMBDA, 'm-',
            label="Region of interest is below this line (Sheath top)")
    ax.fill_between(rnormalize, 0, znormalize, interpolate=True, color='red')
    plt.xlabel(r"$\frac{r}{\Lambda}$", fontsize=35)
    plt.ylabel(r"$\frac{z}{\Lambda}$", fontsize=35)
    plt.xlim([0, max(rnormalize)])
    plt.ylim([0, 1])
    plt.title("Regions accessible to electrons \n with magnetic field presence")
    plt.legend()
    from pylab import rcParams
    rcParams['figure.figsize'] = 13, 12
    if save:
        plt.savefig('inaccessible', format='eps')
    fig.show()


def plotmodifiedefield(rnormalize, LAMBDA, znormalize, gridr, gridz, Evalsr, Evalsz, Evalsradial, Evalsradialz,
                       Evalsheathr, Evalsheath, chargepos, rmax):
    plt.figure()
    plt.plot(numpy.array(rnormalize) * LAMBDA, numpy.array(znormalize) * LAMBDA, 'r-',
             label='Inaccessible region for all p0')
    plt.quiver(gridr, gridz, Evalsr, Evalsz, color='b', label='modified field')
    # plt.quiver(gridr, gridz, Evalsradial, Evalsradialz, color='k', label='radial field')
    # plt.quiver(gridr, gridz, Evalsheathr, Evalsheath, color='g', label='sheath field')
    # plt.plot(chargepos.T[:, 0], chargepos.T[:, 1], 'r.')
    # plt.plot(numpy.arange(len(gridr[0])), numpy.ones(len(gridr[0])) * 0.00038257340070632558, 'm-',
    #          label='crystal plane')
    plt.legend()
    plt.xlim([0, rmax])
    plt.ylim([0, const.sheathd * 1.5])
    plt.show()


def testinterpolator(object, rmax, rnormalize, znormalize, LAMBDA, gridr, gridz, Evalsr, Evalsz, save=False):
    # Generate particles in their equilibrium position
    filehandler = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/crystalpositions2,5K.obj",
                       'rb')
    xinit = pickle.load(filehandler)
    yinit = pickle.load(filehandler)
    zinit = pickle.load(filehandler)
    filehandler.close()
    initpositions = [[i, j, k] for i, j, k in zip(xinit, yinit, zinit)]
    plt.figure()
    for p in numpy.arange(len(xinit)):
        if numpy.sqrt(xinit[p] ** 2 + yinit[p] ** 2) < rmax:
            Etest, gridpointstest = object.interpolate([numpy.sqrt(xinit[p] ** 2 + yinit[p] ** 2), zinit[p]])
            plt.scatter(numpy.array(numpy.mat(gridpointstest).T[0]), numpy.array(numpy.mat(gridpointstest).T[1]),
                        color='k', s=5)
            plt.quiver(numpy.sqrt(xinit[p] ** 2 + yinit[p] ** 2), zinit[p], Etest[0], Etest[1], width=0.002, color='m')
        else:
            pass
    plt.plot(numpy.array(rnormalize) * LAMBDA, numpy.array(znormalize) * LAMBDA, 'r-',
             label='Inaccessible region for all p0')
    plt.quiver(gridr, gridz, Evalsr, Evalsz, color='b', label='modified field')
    # plt.quiver(gridr,gridz,Evalsradial,Evalsradialz,color='k',label='radial field')
    # plt.quiver(gridr,gridz,Evalsheathr,Evalsheath,color='g',label='sheath field')
    # plt.plot(chargepos.T[:,0],chargepos.T[:,1],'r.')
    # plt.plot(numpy.arange(len(gridr[0])),numpy.ones(len(gridr[0]))*0.00038257340070632558,'m-',label='crystal plane')
    plt.xlabel("r (m) ")
    plt.ylabel("z (m) ")
    plt.title("Modified E fields inside \n electron inaccessible region")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.legend(loc=2)
    plt.xlim([0, rmax])
    plt.ylim([0, const.sheathd * 1.5])
    if save:
        from pylab import rcParams
        rcParams['figure.figsize'] = 16, 10
        plt.savefig('fig', format='eps')
    plt.show()
