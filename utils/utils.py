import pickle

import matplotlib.pyplot as plt


def generate_particle_equilibrium_positions():
    # Generate particles in their equilibrium position
    filehandler = open("objects/crystalpositions2K.obj", 'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
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
    plt.xlim((min(initpositions[:, 0]),max(initpositions[:, 0])))
    plt.ylim((min(initpositions[:, 1]),max(initpositions[:, 1])))
    plt.title("Initial positions of 2K particles")
    plt.show()


def prepare_modified_b_field():
    # Prepare modified B field
    filehandler = open("objects/modifiedfield.obj", 'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
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
