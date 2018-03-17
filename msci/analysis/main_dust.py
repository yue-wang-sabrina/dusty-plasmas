# This file is mainly for looking at B field effects of crystals that are already formed.

import numpy
from msci.analysis.analysis_dust import BEffectsAnalysis
from msci.utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
from msci.plots import dustplots
import msci.dustyplasma_cpp.dustcpp_wrapper as dcpp
import matplotlib.pyplot as plt
from IPython import get_ipython
import msci.analysis.constants as const
import math
import decimal

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

beffect1 = BEffectsAnalysis(const)

METHOD = "testspecificmodifiedEfield"


def norm(x):
    return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


if METHOD == "ALLATONCE":  # Drop all particles at once
    beffect1.create_particles(
        numparticles=50,
        initpositions=generate_particle_equilibrium_positions()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=2000,
        init_iterations=0,
        method='Gibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj')
    )
    beffect1.sort_positions_of_particles()
    dustplots.pplot(beffect1)


elif METHOD == "DROP":  # Drop 1 by 1
    beffect1.numparticles = 10
    beffect1.interact_and_iterate_drop_method(
        iterationsB=500,
        init_iterations=500,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field('modifiedfield.obj'))
    beffect1.sort_posititions_drop_method()
    dustplots.pplot(beffect1)


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
        g0 = Dust(const.md, const.radd, const.lambdaD, const.phia, const.Zd * const.e, [0, 0, 0], [0, 0, 0], [0, 0, 0])
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
        # plt.plot(positions,driftforcelist,'bo',label='driftforce')
        # plt.plot(positions,vxbforcelist,'ro',label='vxB force')
        # plt.xlabel("r (m)")
        # plt.ylabel("Force (N)")
        # plt.legend()
        # plt.show()
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
    filename = 'modifiedfield0.01255.obj'
    beffecttemp = BEffectsAnalysis(const)
    beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * 0.01255 / beffecttemp.const.mu0) * numpy.array(
        [0, 0, 1])
    beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
    beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
    beffecttemp.create_particles(
        numparticles=250,
        initpositions=generate_particle_equilibrium_positions()
    )
    beffecttemp.create_pairs()
    beffecttemp.interact_and_iterate(
        iterationsB=500,
        init_iterations=10,
        method='Gibs',
        modified_b_field=prepare_modified_b_field(filename),
        combinedrifts=True
    )
    beffecttemp.sort_positions_of_particles()
    dustplots.pplot(beffecttemp)




