
import numpy
from msci.analysis.analysis_dust import BEffectsAnalysis
from msci.utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field, generate_particle_equilibrium_positions2
from msci.plots import dustplots
import msci.dustyplasma_cpp.dustcpp_wrapper as dcpp
import matplotlib.pyplot as plt
from IPython import get_ipython
import msci.analysis.constants as const
import math
import decimal
import pickle
from msci.utils.interpolate import interpolate
import matplotlib
import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(family='Times New Roman',style='normal', size=30)
csfont = {'fontname':'Times New Roman'}
import dill
ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

beffect1 = BEffectsAnalysis(const)

METHOD = "velocitiesforreport"

if METHOD == "dustvelocities":
    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    # beffect1.const.Bmom = ((2 * math.pi * (0.001) ** 3) * 0.014 / beffect1.const.mu0) * numpy.array(
    #     [0, 0, 1])
    # beffect1.const.dipolepos = [0,0,.0003] #This is for the plot with differential rotation in report

    beffect1.const.dipolepos = [0,0,-0.001] #This is for plot with no differential rotation in report

    beffect1.create_particles(
        numparticles=250,
        initpositions=generate_particle_equilibrium_positions2()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=2000,
        init_iterations=10,
        method='Gibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj')
    )
    beffect1.sort_positions_of_particles()
    dustplots.pplot(beffect1)
    vellist=[]
    rlist=[]
    for i in beffect1.dustdict:
        vellist.append(numpy.sqrt(beffect1.dustdict[i].vel[0]**2+beffect1.dustdict[i].vel[1]**2))
        rlist.append(numpy.sqrt(beffect1.dustdict[i].pos[0]**2+beffect1.dustdict[i].pos[1]**2))
    angvel = numpy.array(vellist)/numpy.array(rlist)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(rlist,vellist)
    ax.set_xlim([0,max(rlist)*1.1])
    ax.set_ylim([0.1*min(vellist),max(vellist)*1.1])
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.set_xlabel("r (m)",fontsize = 25)
    ax.set_ylabel(r"Dust rotational speed (m/s)", fontsize = 25)
    fontsize = 25
    font = {'family': 'normal',
            'weight': 'bold',
            'size': fontsize}
    matplotlib.rc('font', **font)
    fig.show()


elif METHOD == "velocitiesforreport":
    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    beffect1.const.Bmom = ((2 * math.pi * (0.002) ** 3) * 0.014 / beffect1.const.mu0) * numpy.array(
        [0, 0, 1])
    beffect1.const.dipolepos = [0,0,-.0000] #This is for the plot with differential rotation in report

    # beffect1.const.dipolepos = [0, 0, -0.0001]  # This is for plot with no differential rotation in report

    numparticles = 250
    beffect1.create_particles(
        numparticles=numparticles,
        initpositions=generate_particle_equilibrium_positions2()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=1500,
        init_iterations=10,
        method='Gibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj')
    )


    beffect1.sort_positions_of_particles()
    dustplots.pplot(beffect1)

    vels = []
    rs=[]

    for i in beffect1.dustdict:
        vels.append(beffect1.dustdict[i].vel)
        rs.append(numpy.sqrt(beffect1.dustdict[i].pos[0]**2+beffect1.dustdict[i].pos[1]**2))

    positions = beffect1.position[-numparticles:]
    positionsx=[i[0] for i in positions]
    positionsy=[i[1] for i in positions]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 16))
    ax.scatter(positionsx,positionsy)
    ax.set_xlim([min(positionsx),max(positionsx)])
    ax.set_ylim([min(positionsy),max(positionsy)])
    ax.quiver(positionsx,positionsy,[i[0] for i in vels],[i[1] for i in vels],norm = True)
    ax.set_xlabel("x (m)",fontsize=25)
    ax.set_ylabel("y (m)",fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=25)
    fontsize = 25
    font = {'family': 'normal',
            'weight': 'bold',
            'size': fontsize}
    matplotlib.rc('font', **font)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    fig.show()

    vellist=[]
    rlist=[]
    for i in beffect1.dustdict:
        if numpy.sqrt(beffect1.dustdict[i].pos[0]**2+beffect1.dustdict[i].pos[1]**2)>0.7*10**(-3):
            vellist.append(-numpy.sqrt(beffect1.dustdict[i].vel[0] ** 2 + beffect1.dustdict[i].vel[1] ** 2))
        else:
            vellist.append(numpy.sqrt(beffect1.dustdict[i].vel[0] ** 2 + beffect1.dustdict[i].vel[1] ** 2))
        rlist.append(numpy.sqrt(beffect1.dustdict[i].pos[0] ** 2 + beffect1.dustdict[i].pos[1] ** 2))
    angvel = numpy.array(vellist)/numpy.array(rlist)
    fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16))
    ax2.scatter(rlist,vellist)
    ax2.set_xlim([0,max(rlist)*1.1])
    ax2.set_ylim([min(vellist),max(vellist)*1.1])
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax2.set_xlabel("r (m)",fontsize = 25)
    ax2.set_ylabel(r"Dust rotational velocity (m/s)", fontsize = 25)
    fontsize = 20
    font = {'family': 'normal',
            'weight': 'bold',
            'size': fontsize}
    matplotlib.rc('font', **font)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    fig2.show()

elif METHOD == "velvsB": #The pickled object is every velocity vector of every particle in every beffectlist (i.e. every numparticles block is a single beffect run)
    Bmomstrength = numpy.arange(0.014, 2, 0.05)

    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    filehandler = open(b'velsvsB.obj','wb')
    beffectlist=dill.load(beffectlist)
    filehandler.close()
    speeds = []
    distaway = 0.0006
    positions = []
    for i in beffectlist:
        minlist = []
        for j in i.dustdict:
            POS = i.dustdict[j].pos
            absdist = abs(numpy.sqrt(POS[0] ** 2 + POS[1] ** 2) - distaway)
            minlist.append(absdist)
        index = minlist.index(min(minlist))
        veltemp = i.dustdict[list(i.dustdict.keys())[index]].vel
        postemp = i.dustdict[list(i.dustdict.keys())[index]].pos
        speeds.append(norm(veltemp))
        positions.append(numpy.sqrt(postemp[0] ** 2 + postemp[1] ** 2))

    figBmom, axBmom = plt.subplots(nrows=1, ncols=1, figsize=(16, 8))

    axBmom.plot(Bmomstrength, speeds, 'o')#, label=r'Dust particle $%sm$ away' % distaway)
    axBmom.set_xlabel("Magnetic strength (T) ", fontsize=30)
    axBmom.set_ylabel("Speeds (m/s)", fontsize=30)
    axBmom.legend(fontsize=25);
    axBmom.tick_params(labelsize=25)
    axBmom.set_xticklabels(axBmom.get_xticks().astype(int), csfont)
    axBmom.set_yticklabels(axBmom.get_yticks().astype(int), csfont);
    figBmom.show()
