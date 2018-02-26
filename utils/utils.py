import pickle
import numpy
import math
import matplotlib.pyplot as plt


def generate_particle_equilibrium_positions():
    # Generate particles in their equilibrium position
    filehandler = open("objects/crystalpositions2K.obj",
                       'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
    xinit = pickle.load(filehandler)
    yinit = pickle.load(filehandler)
    zinit = pickle.load(filehandler)
    filehandler.close()
    initpositions = [[i, j, k] for i, j, k in zip(xinit, yinit, zinit)]
    return initpositions


def plot_particle_equilibrium_positions(initpositions):
    # Quick plot to view initial configuration
    plt.figure()
    plt.scatter(initpositions[:, 0], initpositions[:, 1])
    plt.xlim((min(initpositions[:, 0]), max(initpositions[:, 0])))
    plt.ylim((min(initpositions[:, 1]), max(initpositions[:, 1])))
    plt.title("Initial positions of 2K particles")
    plt.show()


def prepare_modified_b_field():
    # Prepare modified B field
    filehandler = open("objects/modifiedfield.obj",
                       'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
    gridr = pickle.load(filehandler)
    gridz = pickle.load(filehandler)
    Evalsr = pickle.load(filehandler)
    Evalsz = pickle.load(filehandler)
    rmax = pickle.load(filehandler)
    separationsheath = pickle.load(filehandler)
    separationhor1 = pickle.load(filehandler)
    separationhor2 = pickle.load(filehandler)
    firstpoint = pickle.load(filehandler)
    filehandler.close()
    modified_b_field = {
        'gridr': gridr,
        'gridz': gridz,
        'Evalsr': Evalsr,
        'Evalsz': Evalsz,
        'rmax': rmax,
        'separationsheath': separationsheath,
        'separationhor1': separationhor1,
        'separationhor2': separationhor2,
        'firstpoint': firstpoint,
    }
    return modified_b_field


def bootstrap(drifts, bsit=2000):  # Bootstrap iteration is to take bsit resamplings.
    drifts.sort()
    mean = numpy.mean(drifts)
    samplesind = numpy.random.choice(len(drifts), (bsit, len(drifts)))
    samples = []
    sampleavs = []
    for i in numpy.arange(bsit):
        sample = [drifts[k] for k in samplesind[i]]
        samples.append(sample)
        sampleavs.append(numpy.mean(sample))
    delta = [j - mean for j in sampleavs]
    delta.sort()
    tenpt = math.ceil(0.1 * len(drifts))
    nintypt = math.ceil(0.9 * (len(drifts)))
    error = [delta[tenpt - 1], delta[nintypt - 1]]
    return mean, error


def pointdipoleB(mu0, Bmom, r):
    r = numpy.array(r) + numpy.array([0, 0, -0.003])
    magr = numpy.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    rhat = numpy.array(r) / magr
    magBmom = numpy.sqrt(Bmom[0] ** 2 + Bmom[1] ** 2 + Bmom[2] ** 2)
    Bmomhat = numpy.array(Bmom) / magBmom
    return (mu0 * (3 * numpy.dot(rhat, numpy.array(Bmom)) * rhat - Bmom) / (4 * math.pi * magr ** 3))
