import matplotlib
import numpy
import pickle
import scipy
import matplotlib.pyplot as plt
import matplotlib.animation
import mpl_toolkits.mplot3d.axes3d as p3

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

from particles.ions import Ion
from tqdm import tqdm
import constants as const


def plotthermalkick(position, tau):
    position = numpy.array(position)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(position[:, 0], position[:, 1], position[:, 2], 'r--', label="RK4")
    ax.scatter(position[:, 0][-1], position[:, 1][-1], position[:, 2][-1], 'ro', label="End", s=50)
    ax.scatter(position[:, 0][0], position[:, 1][0], position[:, 2][0], 'bo', label="Start", s=50)
    plt.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.title("Thermal kick included at every dt=%s seconds" % tau)
    plt.show()


def plotratiodriftsingletau(drifts): # correspond to function averagekickeffect in analysis_ion
    fig = plt.figure()
    plt.scatter(numpy.arange(len(drifts)), drifts)
    plt.xlabel("Iteration number")
    plt.ylabel("driftlength_collision/driftlength_nocollision")
    plt.title("Plot of ratio of drift with thermal kicks to drift using pure EXB drift")
    fig.show()


def plotratiodriftmanytau(filename, lentaulist, omega, tau): # correspond to function changetau in analysis_ion
    ion = Ion(pos=[0, 0, 0], vel=[0, 0, 0], acc=[0, 0, 0], omega=omega, tau=tau)
    filehandler = open(filename, 'rb')
    drifts = []
    bs = []

    for i in numpy.arange(lentaulist):
        drifts.append(pickle.load(filehandler))
        bs.append(pickle.load(filehandler))
    omegataulist = numpy.array(pickle.load(filehandler))
    time = pickle.load(filehandler)
    filehandler.close()
    driftsav = [i[0] for i in bs]
    driftsav = numpy.array(driftsav) / (-ion.dt * (ion.Eabs / ion.Babs) * time)
    poserr = [i[1][0] for i in bs]
    negerr = [i[1][1] for i in bs]
    poserr = numpy.array(poserr) / (ion.dt * (ion.Eabs / ion.Babs) * time)
    negerr = numpy.array(negerr) / (ion.dt * (ion.Eabs / ion.Babs) * time)

    func = lambda n: (omega * tau) ** n / (1 + (omega * tau) ** n) - driftsav[-1]
    n = scipy.optimize.fsolve(func, 2)

    fig = plt.figure()
    plt.errorbar(omegataulist, driftsav, yerr=[poserr, negerr], fmt='o', ecolor='g', label='Average drifts')
    plt.xlabel("omega*tau (s)")
    plt.ylabel("Simulated drift/theoretical drift (m)")
    plt.legend()
    plt.title(
        r"Different runtime, (omega*tau)^n for n=%s" % n)
    fig.show()

def plotratiomanytime(filename, lentimelist, omega, tau):
    ion = Ion(pos=[0, 0, 0], vel=[0, 0, 0], acc=[0, 0, 0], omega=omega, tau=tau)
    filehandler = open(filename, 'rb')
    drifts = []
    bs = []

    for i in numpy.arange(lentimelist):
        drifts.append(pickle.load(filehandler))
        bs.append(pickle.load(filehandler))
    timelist = pickle.load(filehandler)
    filehandler.close()
    timelist = numpy.array(timelist)
    driftsav = [i[0] for i in bs]
    driftsav = numpy.array(driftsav) / ((ion.Eabs / ion.Babs) * numpy.array(timelist))
    poserr = [i[1][0] for i in bs]
    negerr = [i[1][1] for i in bs]
    poserr = numpy.array(poserr) / ((ion.Eabs / ion.Babs) * numpy.array(timelist))
    negerr = numpy.array(negerr) / ((ion.Eabs / ion.Babs) * numpy.array(timelist))

    func = lambda n: (omega * tau) ** n / (1 + (omega * tau) ** n) - driftsav[-1]
    n = scipy.optimize.fsolve(func, 2)

    fig = plt.figure()
    plt.errorbar(timelist, driftsav, yerr=[poserr, negerr], fmt='o', ecolor='g', label='Average drifts')
    plt.xlabel("Run time (s)")
    plt.ylabel("Simulated drift/theoretical drift (m)")
    plt.legend()
    plt.title(
        r"Different runtime, (omega*tau)^n for n=%s" % n)
    fig.show()


def plotolddrifts():
    ion = Ion(pos=[0, 0, 0], vel=[0, 0, 0], acc=[0, 0, 0], omega=None, tau=None)
    # filehandler1 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt1.obj",
    #     'rb')
    # filehandler2 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEtauE6.obj", 'rb')
    # filehandler3 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithE.obj",
    #                     'rb')
    # filehandler4 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt15.obj",
    #     'rb')
    # filehandler5 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt25.obj",
    #     'rb')
    # filehandler6 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt3.obj",
    #     'rb')
    # filehandler7 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt1real.obj",
    #     'rb')
    # filehandler8 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt2.obj",
    #     'rb')
    # filehandler9 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt07.obj",
    #     'rb')
    # filehandler10 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt325.obj",
    #     'rb')
    # filehandler11 = open(
    #     "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/driftsconsttauwithEomegataupt275.obj",
    #     'rb')

    driftobjects = [
        {'name': 'driftobject0', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt1.obj'},
        {'name': 'driftobject1', 'drifts': [], 'bs': [], 'num_of_runs': 4, 'filename': 'driftsconsttauwithEtauE6.obj'},
        {'name': 'driftobject2', 'drifts': [], 'bs': [], 'num_of_runs': 2, 'filename': 'driftsconsttauwithE.obj'},
        {'name': 'driftobject3', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt15.obj'},
        {'name': 'driftobject4', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt25.obj'},
        {'name': 'driftobject5', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt3.obj'},
        {'name': 'driftobject6', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt1real.obj'},
        {'name': 'driftobject7', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt2.obj'},
        {'name': 'driftobject8', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt07.obj'},
        {'name': 'driftobject9', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt325.obj'},
        {'name': 'driftobject10', 'drifts': [], 'bs': [], 'num_of_runs': 1,
         'filename': 'driftsconsttauwithEomegataupt275.obj'},

    ]

    for driftobject in tqdm(driftobjects):
        filehandler = open(
            "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/objects_oldruns/{}".format(driftobject['filename']),
            'rb'
        )

        for i in numpy.arange(driftobject['num_of_runs']):
            driftobject['drifts'].append(pickle.load(filehandler))
            driftobject['bs'].append(pickle.load(filehandler))

        driftobject['iterations'] = pickle.load(filehandler)
        filehandler.close()
        if driftobject['name']=='driftobject0':
            driftobject['timelist'] = [0.5, 1]
        elif driftobject['name']=='driftobject1' or driftobject['name']=='driftobject2':
            driftobject['timelist'] = numpy.array(driftobject['iterations']) * ion.dt/(const.e*0.014/const.mi)
        else:
            driftobject['timelist'] = numpy.array(driftobject['iterations']) * ion.dt
        print(driftobject['timelist'])
        driftobject['avx'] = [i[0] for i in driftobject['bs']]
        driftobject['av'] = numpy.array(driftobject['avx']) / (
        -ion.dt * (1 / 0.014) * numpy.array(driftobject['timelist']) / ion.dt)
        driftobject['poserr'] = [i[1][0] for i in driftobject['bs']]
        driftobject['negerr'] = [i[1][1] for i in driftobject['bs']]
        driftobject['poserr'] = numpy.array(driftobject['poserr']) / (
        ion.dt * (1 / 0.014) * numpy.array(driftobject['timelist']) / ion.dt)
        driftobject['negerr'] = numpy.array(driftobject['negerr']) / (
        ion.dt * (1 / 0.014) * numpy.array(driftobject['timelist']) / ion.dt)

    font = {'family': 'normal',
            'weight': 'bold',
            'size': 25}

    matplotlib.rc('font', **font)
    fig = plt.figure()
    finalratios = []
    for i in numpy.arange(len(driftobjects)):
        finalratios.append(driftobjects[i]['av'][-1])
    omegatau = [0.01, 0.033813824472513944, 0.3381382447251395, 0.15, 0.25, 0.3, 0.1, 0.2, 0.07, 0.325, 0.275]
    plt.errorbar(omegatau, finalratios, yerr=11 * [driftobjects[1]['poserr'][-1]],
                 xerr=11 * [driftobjects[1]['negerr'][-1]], fmt='o', capthick=2)

    theoryx = numpy.arange(0, max(omegatau), 10 ** (-2))
    plt.plot(theoryx, theoryx ** 3 / (1 + theoryx) ** 3, 'r-', label=r'$y = \frac{(\omega*\tau)^3}{1+(\omega*\tau)^3}$')
    plt.xlabel(r"$\omega*\tau$", fontsize=25)
    plt.ylabel("simulated drift velocity / theoretical drift velocity", fontsize=25)
    plt.title("Reduction factor for ion-drift velocity \n due to ion-neutral collisions")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.ylim([0,max(finalratios)])
    plt.legend(loc=2)
    from pylab import rcParams
    rcParams['figure.figsize'] = 5, 2
    fig.show()
    # plt.savefig('fig', format='eps')
