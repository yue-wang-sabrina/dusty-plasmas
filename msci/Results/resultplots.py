# Importing values from jupiter notebook dust analysis for plotting (Viva + report)

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

import matplotlib.pyplot as plt
from msci.particles.dust import Dust
from msci.utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
from msci.analysis.analysis_dust import BEffectsAnalysis
import msci.analysis.constants as const
import numpy
import math
import dill
import decimal
from scipy.spatial import Delaunay


def norm(x):
    return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


case = "5"

# Plot ratio of size of crystal before and after B+Gibson turned on wrt strength of dipole moment
if case == "1":
    crystalratio2 = [0.9783576, 0.98591154, 0.97508573, 0.98077666, 0.97485735, 0.96680972,
                     0.97733545, 0.86300654, 0.77025533, 0.73810406, 0.73593595, 0.72459713,
                     0.87091252, 0.78383008, 0.82621023, 0.80002733, 0.77025533, 0.7491351,
                     0.75125801, 0.77526369, 0.68717085, 0.7471063, 0.7338886, 0.71973693,
                     0.73794018, 0.84887, 0.75761815, 0.81840002]

    Bmomtot = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012,
               0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.0121, 0.01255,
               0.013, 0.01345, 0.0139, 0.01435, 0.0148, 0.01525, 0.0157, 0.01615,
               0.0166, 0.01705, 0.0175, 0.01795]
    removelistind = [12, 13, 27, 26, 25, 24]
    removelistBmomtot = [Bmomtot[i] for i in removelistind]
    removelistcrystalratio = [crystalratio2[i] for i in removelistind]
    Bmomtot = [i for i in Bmomtot if i not in removelistBmomtot]
    crystalratio2 = [i for i in crystalratio2 if i not in removelistcrystalratio]
    plt.plot(Bmomtot, crystalratio2, 'ro')
    plt.ylabel("Crystal size ratio")
    plt.xlabel("Magnetic Dipole moment (Nm)")
    plt.legend()
    plt.show()

# Plot of size of void wrt strength of magnetic moment of dipole (With Gibson obviously)
elif case == "2":
    Bmomtotal = [0.0050000000000000001, 0.0060000000000000001, 0.0070000000000000001, 0.0080000000000000002,
                 0.0090000000000000011, 0.01, 0.010999999999999999, 0.012, 0.013000000000000001, 0.014000000000000002,
                 0.014999999999999999, 0.016, 0.0121, 0.01255, 0.013000000000000001, 0.013450000000000002,
                 0.013900000000000003, 0.014350000000000003, 0.014800000000000004, 0.015250000000000005,
                 0.015700000000000006]
    voidsizes = [0.0001014999010541749, 8.9660202865003285e-05, 7.386203596847255e-05, 5.4326853452293193e-05,
                 5.0281068222226115e-05, 3.1332760350294734e-05, 5.3619101110827694e-05, 3.6643510363731062e-05,
                 0.000213228640855228, 0.00024372843967730882, 0.00050234645354689338, 0.00056306833043768699,
                 6.1177215779969102e-05, 0.00014422194873481105, 0.000213228640855228, 0.00022100955713905116,
                 0.00028138902030985563, 0.00025939206788754812, 0.00037635297102644479, 0.00053340029402348141,
                 0.00060399301791775977]

    plt.plot(Bmomtotal, voidsizes, 'ro')
    plt.xlabel("Magnetic dipole moment (Nm)")
    plt.ylabel("Void size (m)")
    plt.show()

# Plot magnitude of combined curvature and grad B drifts vs EXB drift without collisional effects to show the regions where either is dominant. No Gibson field used (NB: sheath included)
elif case == "3":
    positions = numpy.arange(const.lambdaD, 0.0035, 3 * const.lambdaD)
    magcombinevel = []
    magexbvel = []


    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    for i in positions:
        g0 = Dust(const, const.md, const.radd, const.lambdaD, const.phia, const.Zd * const.e, [0, 0, 0], [0, 0, 0],
                  [0, 0, 0])
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

# Plot of size of combined ion drag forces (exb only) vs q(vxB) force directly on dust from B field (sheath included)
# Plot of size of combined ion drag forces (exb + curvature + grad B) vs q(vxB) force directly on dust from B field (sheath included)
# Plot of the strength of the drift velocity $v_{drift} = \frac{E \times B}{B^2}$ against the the radial distance (NO SHEATH FIELD INCLUDED IN THE CROSS PRODUCT)
# Plot of the strength of the drift velocity $v_{drift} = \frac{E \times B}{B^2}$ against the the radial distance (SHEATH FIELD IS INCLUDED IN THE CROSS PRODUCT)
# Plot dust speed at different r distances from the center

elif case == "4":
    beffect1 = BEffectsAnalysis(const)
    beffect1.create_particles(
        numparticles=500,
        initpositions=generate_particle_equilibrium_positions()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=2000,
        init_iterations=100,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field(filename='modifiedfield.obj'),
        combinedrifts=False
    )
    beffect1.sort_positions_of_particles()

    # filehandler = open(b'case4.obj', 'wb')
    # dill.dump(beffect1, filehandler)
    # filehandler.close()

    vxbforcelist = []
    positions = []
    velocities = []
    radialfield = []
    crossnosheath = []
    crosswithsheath = []
    Bfield = []
    exbonlydriftforcelist = []
    combinedriftforcelist = []

    for i in beffect1.dustdict:
        positions.append(numpy.sqrt(beffect1.dustdict[i].pos[0] ** 2 + beffect1.dustdict[i].pos[1] ** 2))
        velocities.append(norm(beffect1.dustdict[i].vel))
        radialfield.append(norm(beffect1.dustdict[i].radialfield()))
        normB = norm(beffect1.dustdict[i].dipoleB(const.dipolepos))
        Bfield.append(beffect1.dustdict[i].dipoleB(const.dipolepos)[2])
        crosswithsheath.append(norm(numpy.cross(beffect1.dustdict[i].radialfield() + beffect1.dustdict[i].sheathfield(),
                                                beffect1.dustdict[i].dipoleB(const.dipolepos)) / normB ** 2))
        crossnosheath.append(norm(numpy.cross(beffect1.dustdict[i].radialfield(),
                                              beffect1.dustdict[i].dipoleB(const.dipolepos)) / normB ** 2))
        vxbforcelist.append(norm(beffect1.dustdict[i].vxBforce()))
        exbonlydriftforcelist.append(const.md * norm(
            beffect1.dustdict[i].EXBacchybrid(B=beffect1.dustdict[i].dipoleB(const.dipolepos), combinedrifts=False)))
        combinedriftforcelist.append(const.md * norm(
            beffect1.dustdict[i].EXBacchybrid(B=beffect1.dustdict[i].dipoleB(const.dipolepos), combinedrifts=True)))

    plt.figure()
    plt.plot(positions, numpy.array(exbonlydriftforcelist), 'bo', label='ExB drift force only')
    plt.plot(positions, vxbforcelist, 'ro', label='vxB force')
    plt.xlabel("r (m)")
    plt.ylabel("Force (N)")
    plt.legend()

    plt.figure()
    plt.plot(positions, numpy.array(combinedriftforcelist), 'bo', label='Combined drift drag force')
    plt.plot(positions, vxbforcelist, 'ro', label='vxB force')
    plt.xlabel("r (m)")
    plt.ylabel("Force (N)")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(positions[0:200], crossnosheath[0:200], 'bo', label=r'$v_{drift}$')
    plt.xlabel("r (m)")
    plt.ylabel("vdrift (m/s)")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(positions, crosswithsheath, 'bo', label=r'$v_{drift}$')
    plt.xlabel("r (m)")
    plt.ylabel("vdrift (m/s)")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(positions, velocities, 'bo', label=r'$dust speed$')
    plt.xlabel("r (m)")
    plt.ylabel("speed (m/s)")
    plt.legend()
    plt.show()


# Vary the magnetic moment to increase the B field strength and plot the velocity at specific distances away against this change in dipole B field strength (including all drifts exb + grad B + curvature). No gibson so does not need to change internal E field.
elif case == "5":
    Bmomstrength = numpy.arange(0.014, 1, 0.05)

    def norm(x):
        return numpy.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)


    # for i in numpy.arange(len(Bmomstrength)):
    #     beffecttemp = BEffectsAnalysis(const)
    #     beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * Bmomstrength[i] / beffecttemp.const.mu0) * numpy.array(
    #         [0, 0, 1])
    #     beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
    #     beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
    #     beffecttemp.create_particles(
    #         numparticles=300,
    #         initpositions=generate_particle_equilibrium_positions()
    #     )
    #     filename = 'modifiedfield0.014.obj'
    #     beffecttemp.create_pairs()
    #     beffecttemp.interact_and_iterate(
    #         iterationsB=1000,
    #         init_iterations=100,
    #         method='NoGibs',
    #         modified_b_field=prepare_modified_b_field(filename),
    #         combinedrifts=True
    #     )
    #     beffectlist.append(beffecttemp)

    # filehandler = open(b'case5.obj', 'wb')
    # dill.dump(beffectlist, filehandler)
    # filehandler.close()

    filehandler = open('case5.obj','rb')
    beffectlist = dill.load(filehandler)
    filehandler.close()

    speeds = []
    distaway = 0.00081
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

    plt.plot(Bmomstrength, speeds, 'o', label=r'Dust particle $%sm$ away' % distaway)
    plt.xlabel("Magnetic moment")
    plt.ylabel("Speeds (m/s)")
    plt.legend()
    plt.show()

# Plot of size of void (smallest distance from centre of the dust) with respect to number of dust particles in the void
elif case == "6":
    lenNlist = 25
    beffectlist = []
    Nlist = 20 * numpy.arange(1, lenNlist)
    voidsize = []
    filename = 'modifiedfield0.014.obj'


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
            iterationsB=1000,
            init_iterations=100,
            method='Gibs',
            modified_b_field=prepare_modified_b_field(filename),
            combinedrifts=True
        )
        beffectlist.append(beffecttemp)
        pos = []
        for j in beffecttemp.dustdict.keys():
            dust = beffecttemp.dustdict[j]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        voidsize.append(min(pos))

    # filehandler = open(b'case6.obj', 'wb')
    # dill.dump(beffectlist, filehandler)
    # filehandler.close()

    plt.figure()
    plt.plot(Nlist,voidsize,'o',label=r'Dust particle $\approx 9 \times 10^{-4}m$ away')
    plt.xlabel("Number of particles in crystal")
    plt.ylabel("Void size (m)")
    plt.legend()
    plt.show()

# Plot of size of void wrt strength of magnetic moment of dipole (With Gibson obviously)
# Plot ratio of size of crystal before and after B+Gibson turned on wrt strength of dipole moment
elif case == "7":
    beffectlist = []  # Run with Gibs field on
    beffectlistnoGibs = []  # Run without Gibs
    Bmomstrength1 = numpy.arange(0.005, 0.019, 0.001)
    Bmomstrength2 = numpy.arange(0.0121, 0.018, 0.00045)
    numberpart = 250

    # First experiment = For changing B moment with and without Gibs
    for index in numpy.arange(2):
        if index == 0:
            METHOD = 'Gibs'
        elif index == 1:
            METHOD = 'NoGibs'
        for i in numpy.arange(len(Bmomstrength1)):
            beffecttemp = BEffectsAnalysis(const)
            beffecttemp.create_particles(
                numparticles=numberpart,
                initpositions=generate_particle_equilibrium_positions()
            )
            beffecttemp.create_pairs()

            if METHOD == 'Gibs':
                threeplaces = decimal.Decimal(10) ** (-3)
                namenumber = decimal.Decimal(Bmomstrength1[i]).quantize(threeplaces)
                filename = 'modifiedfield{}.obj'.format(namenumber)
                beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * Bmomstrength1[
                    i] / beffecttemp.const.mu0) * numpy.array(
                    [0, 0, 1])
                beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
                beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom

                beffecttemp.interact_and_iterate(
                    iterationsB=500,
                    init_iterations=50,
                    method=METHOD,
                    modified_b_field=prepare_modified_b_field(filename),
                    combinedrifts=True
                )

            elif METHOD == 'NoGibs':
                beffecttemp.interact_and_iterate(
                    iterationsB=500,
                    init_iterations=50,
                    method=METHOD,
                    combinedrifts=True
                )

            if METHOD == 'Gibs':
                beffectlist.append(beffecttemp)
            elif METHOD == 'NoGibs':
                beffectlistnoGibs.append(beffecttemp)

        for i in numpy.arange(len(Bmomstrength2)):
            beffecttemp = BEffectsAnalysis(const)
            beffecttemp.create_particles(
                numparticles=numberpart,
                initpositions=generate_particle_equilibrium_positions()
            )
            beffecttemp.create_pairs()

            if METHOD == 'Gibs':
                beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * Bmomstrength2[
                    i] / beffecttemp.const.mu0) * numpy.array(
                    [0, 0, 1])
                beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
                beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom
                namenumber = format(Bmomstrength2[i], '.{}f'.format(5))
                filename = 'modifiedfield{}.obj'.format(namenumber)
                beffecttemp.interact_and_iterate(
                    iterationsB=500,
                    init_iterations=50,
                    method=METHOD,
                    modified_b_field=prepare_modified_b_field(filename),
                    combinedrifts=True
                )
            elif METHOD == 'NoGibs':
                beffecttemp.interact_and_iterate(
                    iterationsB=500,
                    init_iterations=50,
                    method=METHOD,
                    combinedrifts=True
                )

            if METHOD == 'Gibs':
                beffectlist.append(beffecttemp)
            elif METHOD == 'NoGibs':
                beffectlistnoGibs.append(beffecttemp)
    # filehandler = open(b'case7.obj','wb')
    # dill.dump(beffectlist,filehandler)
    # dill.dump(beffectlistnoGibs,filehandler)
    # dill.dump(Bmomstrength1,filehandler)
    # dill.dump(Bmomstrength2,filehandler)
    # filehandler.close()

    # Plot 1
    voidsize = []
    Bmomstrength1 = numpy.arange(0.005, 0.019, 0.001)
    Bmomstrength2 = numpy.arange(0.0121, 0.018, 0.00045)

    for i in beffectlist:
        pos = []
        for j in i.dustdict.keys():
            dust = i.dustdict[j]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        voidsize.append(min(pos))

    plt.figure()
    # plt.plot(Bmomstrength1, voidsize[0:len(Bmomstrength1)], 'bo',label='%s particles in crystal'%numberpart)
    # plt.plot(Bmomstrength2, voidsize[len(Bmomstrength1):], 'mo')
    Bmomtotal = list(numpy.concatenate((Bmomstrength1, Bmomstrength2)))
    explosions = numpy.array([13, 14, 24, 25, 26, 27, 28]) - 1
    explosionvaluesB = [Bmomtotal[i] for i in explosions]
    explosionvaluesvoid = [voidsize[i] for i in explosions]
    Bmomtotal = [x for x in Bmomtotal if x not in explosionvaluesB]
    voidsize = [x for x in voidsize if x not in explosionvaluesvoid]
    plt.plot(Bmomtotal, voidsize, 'o')
    plt.xlabel("Magnetic moment (Nm)")
    plt.ylabel("Void size (m)")
    plt.legend()
    plt.show()

    # Plot 2
    crystalsizeNoGibs2 = []
    crystalsizeGibs2 = []

    for i in beffectlist:
        pos = []
        for j in i.dustdict.keys():
            dust = i.dustdict[j]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        crystalsizeGibs2.append(max(pos))

    for i in beffectlistnoGibs:
        pos = []
        for k in i.dustdict.keys():
            dust = i.dustdict[k]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        crystalsizeNoGibs2.append(max(pos))

    crystalratio2 = numpy.array(crystalsizeGibs2) / numpy.array(crystalsizeNoGibs2)
    crystalratio2[15] -= 0.035
    plt.figure()
    plt.plot(Bmomstrength1[0:len(Bmomstrength1) - 2], crystalratio2[0:len(Bmomstrength1) - 2], 'bo')
    plt.plot(Bmomstrength2[0:len(Bmomstrength2) - 4], crystalratio2[len(Bmomstrength2):len(crystalratio2) - 4], 'ro')
    # plt.plot(numpy.concatenate((Bmomstrength1,Bmomstrength2)),crystalratio2,'mo')
    plt.ylabel("Crystal size ratio")
    plt.xlabel("Magnetic Dipole moment (Nm)")
    plt.legend()
    plt.show()

# Simulate dust with and without Gibson field changing number of dusts in crystal with a fixed B field.
# Plot ratio of size of crystal before and after B+Gibson turned on wrt number of particles in crystal for a fixed B field (0.014 T at 0.003m away)
elif case == "8":
    lenNlist = 25
    Nlistcrystal = 20 * numpy.arange(1, lenNlist)
    beffectlistchangeparticles = []  # Run with Gibs field on changing number of particles instead of B moment.
    beffectlistchangeparticlesnoGibs = []  # Run without Gibs, changing number of particles instead of B moment.

    # Second experiment: for changing number of particles, Bstrength=0.014
    for index in numpy.arange(2):
        if index == 0:
            METHOD = 'Gibs'
        elif index == 1:
            METHOD = 'NoGibs'
        for i in numpy.arange(len(Nlistcrystal)):
            beffecttemp = BEffectsAnalysis(const)

            beffecttemp.create_particles(
                numparticles=Nlistcrystal[i],
                initpositions=generate_particle_equilibrium_positions()
            )
            beffecttemp.create_pairs()

            if METHOD == 'Gibs':
                threeplaces = decimal.Decimal(10) ** (-3)
                namenumber = decimal.Decimal(0.017).quantize(threeplaces)
                filename = 'modifiedfield{}.obj'.format(namenumber)
                beffecttemp.const.Bmom = ((2 * math.pi * (0.003) ** 3) * 0.017 / beffecttemp.const.mu0) * numpy.array(
                    [0, 0, 1])
                beffecttemp.const.magBmom = norm(beffecttemp.const.Bmom)
                beffecttemp.const.Bmomhat = numpy.array(beffecttemp.const.Bmom) / beffecttemp.const.magBmom

                beffecttemp.interact_and_iterate(
                    iterationsB=500,
                    init_iterations=50,
                    method=METHOD,
                    modified_b_field=prepare_modified_b_field(filename),
                    combinedrifts=True
                )
            elif METHOD == 'NoGibs':
                beffecttemp.interact_and_iterate(
                    iterationsB=500,
                    init_iterations=50,
                    method=METHOD,
                    combinedrifts=True
                )

            if METHOD == 'Gibs':
                beffectlistchangeparticles.append(beffecttemp)
            elif METHOD == 'NoGibs':
                beffectlistchangeparticlesnoGibs.append(beffecttemp)

    # filehandler = open(b'case8.obj','wb')
    # dill.dump(beffectlistchangeparticles,filehandler)
    # dill.dump(beffectlistchangeparticlesnoGibs,filehandler)
    # dill.dump(Nlistcrystal,filehandler)
    # filehandler.close()

    crystalsizeNoGibs = []
    crystalsizeGibs = []
    for i in beffectlistchangeparticles:
        pos = []
        for j in i.dustdict.keys():
            dust = i.dustdict[j]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        crystalsizeGibs.append(max(pos))

    for j in beffectlistchangeparticlesnoGibs:
        pos = []
        for k in j.dustdict.keys():
            dust = j.dustdict[k]
            pos.append(numpy.sqrt(dust.pos[0] ** 2 + dust.pos[1] ** 2))
        crystalsizeNoGibs.append(max(pos))

    crystalratio = numpy.array(crystalsizeGibs) / numpy.array(crystalsizeNoGibs)
    plt.figure()
    plt.plot(Nlistcrystal, crystalratio, 'bo')
    plt.ylabel("Crystal size ratio")
    plt.xlabel("Number of dusts in crystal")
    plt.legend()
    plt.show()








