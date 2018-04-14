# This file is mainly for looking at B field effects of crystals that are already formed.

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

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

beffect1 = BEffectsAnalysis(const)

METHOD = "velvsB"


def norm(x):
    return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


if METHOD == "ALLATONCE":  # Drop all particles at once
    beffect1.create_particles(
        numparticles=50,
        initpositions=generate_particle_equilibrium_positions2()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=600,
        init_iterations=100,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj')
    )
    beffect1.sort_positions_of_particles()
    dustplots.pplot(beffect1)


elif METHOD == "DROP":  # Drop 1 by 1
    beffect1.numparticles = 250
    beffect1.interact_and_iterate_drop_method(
        iterationsB=0,
        init_iterations=5000,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj'))
    beffect1.sort_posititions_drop_method()
    dustplots.pplot(beffect1,name="drop250",save=True,jn=True)


elif METHOD == "CPP":  # Use C++ to iterate
    particlenum = 50
    itB = 2000
    initit = 0
    beffect2 = dcpp.DustAnalysisCpp(initit, itB, particlenum, 1)

    beffect2.get_modified_field()
    beffect2.get_equilibrium_positions()
    beffect2.run()
    beffect1.numparticles = particlenum
    beffect1.iterationsB = itB
    beffect1.init_iterations = initit
    beffect1.method = "Gibs"
    beffect1.position = [[i, j, k] for i, j, k in zip(beffect2.positions_x, beffect2.positions_y, beffect2.positions_z)]
    beffect1.position_array = numpy.array(beffect1.position)
    beffect1.sort_positions_of_particles()
    beffect1.modified_b_field = prepare_modified_b_field('modifiedfield.obj')
    dustplots.pplot(beffect1)

elif METHOD == "Finddriftvelocities":
    from msci.particles.dust import Dust

    positions = numpy.arange(const.lambdaD, 200 * const.lambdaD, 3 * const.lambdaD)
    magcombinevel = []
    magexbvel = []


    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    for i in positions:
        g0 = Dust(const, const.md, const.radd, const.lambdaD, const.phia, const.Zd * const.e, [0, 0, 0], [0, 0, 0], [0, 0, 0])
        g0.Bswitch = True
        g0.pos = [i, 0, 0.0003825734]
        magcombinevel.append(norm(g0.combinedrift(B=g0.dipoleB(r=const.dipolepos))))
        magexbvel.append(norm(g0.EXBDRIFT(B=g0.dipoleB(r=const.dipolepos))))

    plt.plot(positions, magcombinevel, 'r-', label='combine drifts')
    plt.plot(positions, magexbvel, 'b-', label='ExB drifts')
    plt.xlabel("r (m)")
    plt.ylabel(r'$v_{drift}$ (m/s)')
    plt.legend()
    plt.show()

elif METHOD == "comparevxbwithotherdrifts":
    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    beffect1.create_particles(
        numparticles=10,
        initpositions=generate_particle_equilibrium_positions()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=1000,
        init_iterations=100,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj')
    )
    beffect1.sort_positions_of_particles()
    dustplots.pplot(beffect1)
    vxbforcelist = []
    positions = []
    velocities = []
    radialfield = []
    cross = []
    Bfield = []
    driftforcelist = []

    for i in beffect1.dustdict:
        positions.append(numpy.sqrt(beffect1.dustdict[i].pos[0] ** 2 + beffect1.dustdict[i].pos[1] ** 2))
        velocities.append(norm(beffect1.dustdict[i].vel))
        radialfield.append(norm(beffect1.dustdict[i].radialfield()))
        normB = norm(beffect1.dustdict[i].dipoleB(const.dipolepos))
        Bfield.append(beffect1.dustdict[i].dipoleB(const.dipolepos)[2])
        cross.append(norm(numpy.cross(beffect1.dustdict[i].radialfield(),
                                      beffect1.dustdict[i].dipoleB(const.dipolepos)) / normB ** 2))
        vxbforcelist.append(norm(beffect1.dustdict[i].vxBforce()))
        driftforcelist.append(norm(
            beffect1.dustdict[i].EXBacchybrid(B=beffect1.dustdict[i].dipoleB(const.dipolepos), combinedrifts=True)))

    plt.plot(positions,driftforcelist,'bo',label='driftforce')
    plt.plot(positions,vxbforcelist,'ro',label='vxB force')
    plt.xlabel("r (m)")
    plt.ylabel("Force (N)")
    plt.legend()
    plt.show()

elif METHOD == "dustvelocities":
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
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.set_xlabel("r (m)",fontsize = 25)
    ax.set_ylabel(r"Dust rotational speed (m/s)", fontsize = 25)
    fontsize = 20
    font = {'family': 'normal',
            'weight': 'bold',
            'size': fontsize}
    matplotlib.rc('font', **font)
    fig.show()


elif METHOD == "velocitiesforreport":
    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    beffect1.const.Bmom = ((2 * math.pi * (0.0009) ** 3) * 0.014 / beffect1.const.mu0) * numpy.array(
        [0, 0, 1])
    beffect1.const.dipolepos = [0,0,-.0001] #This is for the plot with differential rotation in report

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
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(positionsx,positionsy)
    ax.set_xlim([min(positionsx),max(positionsx)])
    ax.set_ylim([min(positionsy),max(positionsy)])
    ax.quiver(positionsx,positionsy,[i[0] for i in vels],[i[1] for i in vels],norm = True)
    ax.set_xlabel("x (m)",fontsize=25)
    ax.set_ylabel("y (m)",fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20)
    fontsize = 20
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
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
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


elif METHOD == "compareforces": #Plot in report
    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    beffect1.create_particles(
        numparticles=120,
        initpositions=generate_particle_equilibrium_positions2()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=5000,
        init_iterations=1000,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj')
    )
    # beffect1.sort_positions_of_particles()
    # dustplots.pplot(beffect1)

    neutraldragforcelist = [norm(i) for i in beffect1.neutraldragforcelist]
    vxbforces = [norm(i) for i in beffect1.vxbforces]
    exbforcelist = [norm(i) for i in beffect1.exbforcelist]
    otherdriftforcelist = [norm(i) for i in beffect1.otherdriftforcelist]
    time = beffect1.const.dt*numpy.arange(len(beffect1.neutraldragforcelist))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, neutraldragforcelist, 'b-', label='Neutral Drag Force')
    ax.plot(time, vxbforces, 'r-', label='Lorentz Force')
    ax.plot(time, exbforcelist, 'g-',label='EXB Force')
    ax.plot(time, otherdriftforcelist,'m-',label='Curvature + GradB Forces')
    ax.set_xlabel("Time (s)",fontsize=25)
    ax.set_ylabel("Force (N)",fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20)
    fontsize = 20
    font = {'family': 'normal',
            'weight': 'bold',
            'size': fontsize}
    matplotlib.rc('font', **font)
    ax.legend(loc=5)
    fig.show()


elif METHOD == "velvsB": #The pickled object is every velocity vector of every particle in every beffectlist (i.e. every numparticles block is a single beffect run)
    beffectlist = []
    Bmomstrength = numpy.arange(0.014, 2, 0.05)

    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    for i in numpy.arange(len(Bmomstrength)):
        beffecttemp = BEffectsAnalysis(const)
        beffecttemp.const.Bmom = ((2 * math.pi * (0.00038) ** 3) * Bmomstrength[i] / beffecttemp.const.mu0) * numpy.array(
            [0, 0, 1])
        beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
        beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
        beffecttemp.create_particles(
            numparticles=250,
            initpositions=generate_particle_equilibrium_positions2()
        )
        filename = 'modifiedfield0.014.obj'
        beffecttemp.create_pairs()
        beffecttemp.interact_and_iterate(
            iterationsB=2000,
            init_iterations=10,
            method='NoGibs',
            modified_b_field=prepare_modified_b_field(filename),
            combinedrifts=True
        )
        beffectlist.append(beffecttemp)

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
    axBmom.set_xlabel("Magnetic strength (T) ", fontsize=15)
    axBmom.set_ylabel("Speeds (m/s)", fontsize=15)
    axBmom.legend(fontsize=15);
    figBmom.show()


elif METHOD == "voidsizewrtN":
    lenNlist = 25
    beffectlist = []
    Nlist = 20 * numpy.arange(1,lenNlist)
    voidsize = []

    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    for i in numpy.arange(len(Nlist)):
        beffecttemp = BEffectsAnalysis(const)
        beffecttemp.create_particles(
            numparticles=Nlist[i],
            initpositions=generate_particle_equilibrium_positions()
        )
        beffecttemp.create_pairs()
        beffecttemp.interact_and_iterate(
            iterationsB=500,
            init_iterations=100,
            method='Gibs',
            modified_b_field=prepare_modified_b_field('modifiedfield.obj'),
            combinedrifts=True
        )
        pos=[]
        for j in beffecttemp.dustdict.keys():
            dust=beffecttemp.dustdict[j]
            pos.append(numpy.sqrt(dust.pos[0]**2+dust.pos[1]**2))
        voidsize.append(min(pos))

elif METHOD == "voidsizewrtB":
    beffectlist = []
    Bmomstrength1 = numpy.arange(0.005, 0.019, 0.001)
    Bmomstrength2 = numpy.arange(0.0121, 0.01660, 0.00045)
    numberpart = 10
    initit= 100
    Bit=1000
    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    for i in numpy.arange(len(Bmomstrength1)):
        threeplaces = decimal.Decimal(10) ** (-3)
        namenumber = decimal.Decimal(Bmomstrength1[i]).quantize(threeplaces)
        filename = 'modifiedfield{}.obj'.format(namenumber)
        beffecttemp = BEffectsAnalysis(const)
        beffecttemp.const.Bmom = (
                                 (2 * math.pi * (0.003) ** 3) * Bmomstrength1[i] / beffecttemp.const.mu0) * numpy.array(
            [0, 0, 1])
        beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
        beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
        beffecttemp.create_particles(
            numparticles=numberpart,
            initpositions=generate_particle_equilibrium_positions()
        )
        beffecttemp.create_pairs()
        beffecttemp.interact_and_iterate(
            iterationsB=Bit,
            init_iterations=initit,
            method='Gibs',
            modified_b_field=prepare_modified_b_field(filename),
            combinedrifts=True
        )
        beffectlist.append(beffecttemp)



    for i in numpy.arange(len(Bmomstrength2)):
        namenumber = format(Bmomstrength2[i], '.{}f'.format(5))
        filename = 'modifiedfield{}.obj'.format(namenumber)
        beffecttemp = BEffectsAnalysis(const)
        beffecttemp.const.Bmom = (
                                 (2 * math.pi * (0.003) ** 3) * Bmomstrength2[i] / beffecttemp.const.mu0) * numpy.array(
            [0, 0, 1])
        beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
        beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
        beffecttemp.create_particles(
            numparticles=numberpart,
            initpositions=generate_particle_equilibrium_positions()
        )
        beffecttemp.create_pairs()
        beffecttemp.interact_and_iterate(
            iterationsB=Bit,
            init_iterations=initit,
            method='Gibs',
            modified_b_field=prepare_modified_b_field(filename),
            combinedrifts=True
        )
        beffectlist.append(beffecttemp)

    voidsize = []

    for i in beffectlist:
        pos = []
        for j in i.dustdict.keys():
            dust = i.dustdict[j]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        voidsize.append(min(pos))

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 8))

    ax.plot(Bmomstrength1, voidsize[0:len(Bmomstrength1)], 'bo', label='%s particles in crystal' % numberpart)
    ax.plot(Bmomstrength2, voidsize[len(Bmomstrength1):], 'ro', label='%s particles in crystal' % numberpart)
    ax.set_xlabel("Magnetic moment (Nm)", fontsize=15)
    ax.set_ylabel("Void size (m)", fontsize=15)
    ax.legend(fontsize=15)
    fig.show()

elif METHOD == "testspecificmodifiedEfield":
    filename = 'modifiedfield0.014.obj'
    beffecttemp = BEffectsAnalysis(const)
    beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * 0.014 / beffecttemp.const.mu0) * numpy.array(
        [0, 0, 1])
    beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
    beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
    beffecttemp.const.dipolepos = numpy.array([0,0,-0.0005])
    beffecttemp.const.boxr = 0.04
    beffecttemp.create_particles(
        numparticles=250,
        initpositions=generate_particle_equilibrium_positions()
    )
    beffecttemp.create_pairs()
    beffecttemp.interact_and_iterate(
        iterationsB=1000,
        init_iterations=2000,
        method='Gibs',
        modified_b_field=prepare_modified_b_field(filename),
        combinedrifts=True
    )
    beffecttemp.sort_positions_of_particles()
    dustplots.pplot(beffecttemp)# ,name="Dust_rotation_NoGibs",save=True,jn=True)

    # positions = []
    # velocities = []
    # angvels = []
    #
    # for i in beffecttemp.dustdict:
    #     rad = numpy.sqrt(beffecttemp.dustdict[i].pos[0] ** 2 + beffecttemp.dustdict[i].pos[1] ** 2)
    #     positions.append(rad)
    #     speed=norm(beffecttemp.dustdict[i].vel)
    #     velocities.append(speed)
    #     angvels.append(speed/rad)
    #
    # plt.figure()
    # plt.plot(positions,velocities,'o')
    # plt.xlabel("r(m)")
    # plt.ylabel("speed (m/s)")
    # plt.show()
    # plt.figure()
    # plt.plot(positions,angvels,'o')
    # plt.xlabel("r(m)")
    # plt.ylabel(r"Angular velocity $\omega$ (rad/s)")
    # plt.show()

elif METHOD == "ChecktotalradialEfield":
    numpart=200
    filename = 'modifiedfield0.005.obj'
    beffecttemp = BEffectsAnalysis(const)
    beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * 0.005 / beffecttemp.const.mu0) * numpy.array(
        [0, 0, 1])
    beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
    beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
    beffecttemp.create_particles(
        numparticles=numpart,
        initpositions=generate_particle_equilibrium_positions()
    )
    positions = numpy.arange(5*beffecttemp.const.lambdaD,0.0035,(0.0035-5*beffecttemp.const.lambdaD)/numpart)
    ind=0
    radialfieldlist=[]
    for i in beffecttemp.dustdict:
        beffecttemp.dustdict[i].pos=[positions[ind],0.000001,0.0003825734]
        ind+=1
        modified_b_field = prepare_modified_b_field(filename)
        Efield=interpolate(
            beffecttemp.dustdict[i].getselfpos(),
            modified_b_field['gridr'], modified_b_field['gridz'], modified_b_field['Evalsr'],
            modified_b_field['Evalsz'], modified_b_field['rmax'], modified_b_field['separationsheath'],
            modified_b_field['separationhor1'], modified_b_field['separationhor2'],
            modified_b_field['firstpoint']
        )
        radialfield = beffecttemp.dustdict[i].radialfield()+numpy.array(Efield)
        radialfieldmag = numpy.sqrt(radialfield[0]**2+radialfield[2]**2)
        if radialfield[0]<0:
            radialfieldmag *= -1
        radialfieldlist.append(radialfieldmag)
    plt.plot(positions,radialfieldlist)
    plt.plot(positions,numpy.zeros(len(positions)))
    plt.xlabel("r (m)")
    plt.ylabel("Radial E field strength (+ x-axis)")
    plt.show()

else:
    pass

def generate_particle_equilibrium_positions():
    # Generate particles in their equilibrium position
    filehandler = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/crystalpositions2,5K.obj",
                       'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
    xinit = pickle.load(filehandler)
    yinit = pickle.load(filehandler)
    zinit = pickle.load(filehandler)
    filehandler.close()
    initpositions = numpy.array([[i, j,k] for i, j, k in zip(xinit, yinit,zinit)])
    return initpositions

# initpos = generate_particle_equilibrium_positions()
# plt.plot(initpos[:,0],initpos[:,1],'ro',markersize=2)
# plt.show()



