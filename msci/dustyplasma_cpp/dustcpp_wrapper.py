import ctypes
import os
import numpy as np
import pickle

dir_path = os.path.dirname(os.path.realpath(__file__))


class DustAnalysisCpp:

    def __init__(self, init_iterations, iterationsB, n_particles):
        self.dustlib = ctypes.cdll.LoadLibrary(dir_path + "/dust_wrapper.so")

        self.init_iterations = init_iterations
        self.iterationsB = iterationsB
        self.n_particles = n_particles
        self.xinit = None
        self.yinit = None
        self.zinit = None

        self.positions_x = np.zeros((self.init_iterations + self.iterationsB) * self.n_particles, np.float)
        self.positions_y = np.zeros((self.init_iterations + self.iterationsB) * self.n_particles, np.float)
        self.positions_z = np.zeros((self.init_iterations + self.iterationsB) * self.n_particles, np.float)

    def run(self):
        self.dustlib.simulate.restype = None
        self.dustlib.simulate(
            ctypes.c_int(self.init_iterations),
            ctypes.c_int(self.iterationsB),
            ctypes.c_int(self.n_particles),
            np.ctypeslib.as_ctypes(self.xinit[0:self.n_particles]),
            np.ctypeslib.as_ctypes(self.yinit[0:self.n_particles]),
            np.ctypeslib.as_ctypes(self.zinit[0:self.n_particles]),
            np.ctypeslib.as_ctypes(self.positions_x),
            np.ctypeslib.as_ctypes(self.positions_y),
            np.ctypeslib.as_ctypes(self.positions_z)

        )

    def get_equilibrium_positions(self):
        filehandler = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/crystalpositions2K.obj",
                       'rb')
        self.xinit = np.array(pickle.load(filehandler))
        self.yinit = np.array(pickle.load(filehandler))
        self.zinit = np.array(pickle.load(filehandler))
        filehandler.close()


