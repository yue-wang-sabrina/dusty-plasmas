import numpy as np
import matplotlib.pyplot as plt
from paircorrelationutilities import *
from paircorrelation import pairCorrelationFunction_2D
import msci.analysis.constants as const
import pickle
import numpy
import matplotlib

method = "250"


if method == "1Kparticles":
    filehandler = open("/Users/yuewang/Documents/Repos/dusty-plasmas/msci/Results/40mmcylinder1Kdust.obj",
                           'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
    initpositions = pickle.load(filehandler)
    filehandler.close()
    # Particle setup
    domain_size = const.boxr
    num_particles = 1000

    # Calculation setup
    dr = 1 * const.lambdaD

    ### Random arrangement of particles ###
    particle_radius = const.radd
    rMax = domain_size / 20
    # x = np.random.uniform(low=0, high=domain_size, size=num_particles)
    # y = np.random.uniform(low=0, high=domain_size, size=num_particles)
    x = numpy.array([i[0] for i in initpositions])
    y = numpy.array([i[1] for i in initpositions])


elif method == "250":
    filehandler = open("/Users/yuewang/Documents/Repos/dusty-plasmas/msci/Results/rigidrotation.obj",
                       'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
    initpositions = pickle.load(filehandler)
    filehandler.close()
    # Particle setup
    domain_size = const.boxr
    num_particles = 250

    # Calculation setup
    dr = 0.4 * const.lambdaD

    ### Random arrangement of particles ###
    particle_radius = const.radd
    rMax = domain_size / 50
    # x = np.random.uniform(low=0, high=domain_size, size=num_particles)
    # y = np.random.uniform(low=0, high=domain_size, size=num_particles)
    x = numpy.array([i[0] for i in initpositions])
    y = numpy.array([i[1] for i in initpositions])


    filehandler2 = open("/Users/yuewang/Documents/Repos/dusty-plasmas/msci/Results/differentialrotation.obj",
                       'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
    initpositions2 = pickle.load(filehandler2)
    filehandler2.close()

    x2 = numpy.array([i[0] for i in initpositions2])
    y2 = numpy.array([i[1] for i in initpositions2])


# Compute pair correlation
g_r, r, reference_indices = pairCorrelationFunction_2D(x, y, domain_size, rMax, dr)
g_r2, r2, reference_indices2 = pairCorrelationFunction_2D(x2, y2, domain_size, rMax, dr)

# Visualize
plt.figure()
plt.plot(r2,g_r2, '-',color='red',label="Differential rotation")

plt.plot(r, g_r, '-',color='blue',label="Rigid rotation")
plt.xlabel('r',fontsize=25)
plt.ylabel('g(r)',fontsize=25)
plt.xlim( (0, rMax) )
plt.ylim( (0, 1.05 * g_r.max()) )
plt.tick_params(axis='both', which='major', labelsize=20)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
fontsize = 20
font = {'family': 'normal',
        'weight': 'bold',
        'size': fontsize}
matplotlib.rc('font', **font)
plt.legend()
# plot_adsorbed_circles(x, y, particle_radius, domain_size, reference_indices=reference_indices)

# ### Hexagonal circle packing ###
# particle_radius = 1.0
# domain_size = 50.0
# rMax = domain_size / 3
#
# x, y, domain_width, domain_height = generate_hex_circle_packing(particle_radius, domain_size)
#
# # Compute pair correlation
# g_r, r, reference_indices = pairCorrelationFunction_2D(x, y, domain_size,
#         rMax, dr)
#
# # Visualize
# plt.figure()
# plt.plot(r, g_r, color='black')
# plt.xlabel('r')
# plt.ylabel('g(r)')
# plt.xlim( (0, rMax) )
# plt.ylim( (0, 1.05 * g_r.max()) )
#
# plot_adsorbed_circles(x, y, particle_radius, domain_size, reference_indices=reference_indices)
#
plt.show()