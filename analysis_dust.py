import itertools

import numpy
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import mpl_toolkits.mplot3d.axes3d as p3

from tqdm import tqdm

import constants as const
from particles.dust import Dust
from utils.interpolate import interpolate


class BEffectsAnalysis:
    def __init__(self):
        self.dustdict = None
        self.pairs = None
        self.position = []
        self.positions_df = None
        self.t = None
        self.particles = None
        self.numparticles = None
        self.iterationsB = None
        self.init_iterations = None
        self.position_array = None
        self.modified_b_field = None

    def create_particles(self, numparticles, initpositions):
        # Create dictionary of particles from pickle object

        self.numparticles = numparticles

        names = []
        for int_name in numpy.arange(numparticles):
            names.append('g%s' % int_name)

        self.dustdict = {
            name: Dust(const.md, const.radd, const.lambdaD, const.phia, const.Zd * const.e, [0, 0, 0], [0, 0, 0], [0, 0, 0])
            for name in names
            }

        for i, j in zip(self.dustdict, numpy.arange(len(self.dustdict))):
            self.dustdict[i].pos = initpositions[j]

    def create_pairs(self):
        # Create pairs
        # Check every pair of particle interaction
        list1 = list(self.dustdict.keys())
        list2 = list1
        pairlist = list(itertools.product(list1, list2))
        pairs = set()

        for x, y in pairlist:
            if x != y:
                pairs.add(frozenset((x, y)))

        pairs = list(pairs)

        for pair_index in numpy.arange(len(pairs)):
            pairs[pair_index] = list(pairs[pair_index])

        removelist = []
        for i in pairs:
            if i[0] == i[1]:
                removelist.append(i)

        self.pairs = [i for i in pairs if i not in removelist]

    def interact_and_iterate(self, iterationsB, init_iterations, method='Gibs', modified_b_field=None):
        # Interact and iterate
        self.iterationsB = iterationsB
        self.init_iterations = init_iterations
        self.modified_b_field = modified_b_field

        for i in tqdm(numpy.arange(init_iterations)):
            pairsfinal = []
            for b in self.pairs:
                if self.dustdict[b[0]].intergraind(self.dustdict[b[1]]):
                    pairsfinal.append(b)
                else:
                    pass  # pairsfinal.append(b)
            for j in pairsfinal:
                interactfield = self.dustdict[j[0]].selffield(self.dustdict[j[1]])
                self.dustdict[j[0]].selffieldmany(interactfield)
                self.dustdict[j[1]].selffieldmany(-interactfield)
            for k in self.dustdict:
                self.dustdict[k].steptest(numpy.array(self.dustdict[k].multifields))
                self.position.append(self.dustdict[k].getselfpos())

        for l in self.dustdict:
            self.dustdict[l].Bswitch = True

        for i in tqdm(numpy.arange(iterationsB)):
            pairsfinal = []
            for b in self.pairs:
                if self.dustdict[b[0]].intergraind(self.dustdict[b[1]]):
                    pairsfinal.append(b)
                else:
                    pass
            # Pick out the pairs that are less than 5 lambdaDs apart
            for j in pairsfinal:
                interactfield = self.dustdict[j[0]].selffield(self.dustdict[j[1]])
                self.dustdict[j[0]].selffieldmany(interactfield)
                self.dustdict[j[1]].selffieldmany(-interactfield)
            for k in self.dustdict:
                if method == 'NoGibs':
                    self.dustdict[k].steptest(numpy.array(self.dustdict[k].multifields))
                elif method == 'Gibs':
                    Efield = numpy.array(self.dustdict[k].multifields) + interpolate(
                        self.dustdict[k].getselfpos(),
                        modified_b_field['gridr'], modified_b_field['gridz'], modified_b_field['Evalsr'],
                        modified_b_field['Evalsz'], modified_b_field['rmax'], modified_b_field['separationsheath'],
                        modified_b_field['separationhor1'], modified_b_field['separationhor2'],
                        modified_b_field['firstpoint']
                    )
                    self.dustdict[k].steptest(Efield)
                    self.position.append(self.dustdict[k].getselfpos())

        self.position_array = numpy.array(self.position)

    def sort_positions_of_particles(self):
        # Sort the positions of the particles

        self.t = numpy.array(
            numpy.sort([list(numpy.arange(int(len(self.position_array[:, 0]) / self.numparticles))) * self.numparticles])
        )[0]

        self.positions_df = pd.DataFrame({
            "time": self.t,
            "x": self.position_array[:, 0],
            "y": self.position_array[:, 1],
            "z": self.position_array[:, 2]
        })

    def plot(self):

        def _update_graph(n_iter):
            data = self.positions_df[self.positions_df['time'] == n_iter]
            point.set_data(data.x, data.y)
            point.set_3d_properties(data.z)
            return point,

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(elev=90., azim=90)
        if min(self.position_array[:, 0]) == max(self.position_array[:, 0]):
            ax.set_xlim([-10 * const.lambdaD, 10 * const.lambdaD])
        else:
            ax.set_xlim([min(self.position_array[:, 0]), max(self.position_array[:, 0])])
        if min(self.position_array[:, 1]) == max(self.position_array[:, 1]):
            ax.set_ylim([-10 * const.lambdaD, 10 * const.lambdaD])
        else:
            ax.set_ylim([min(self.position_array[:, 1]), max(self.position_array[:, 1])])
        if min(self.position_array[:, 2]) == max(self.position_array[:, 2]):
            ax.set_zlim([-10 * const.lambdaD, 10 * const.lambdaD])
        else:
            ax.set_zlim([min(self.position_array[:, 2]) - 0.5 * min(self.position_array[:, 2]), max(self.position_array[:, 2])])

        ax.set_xlim([-self.modified_b_field['rmax'] * 1.5, self.modified_b_field['rmax'] * 1.5])
        ax.set_ylim([-self.modified_b_field['rmax'] * 1.5, self.modified_b_field['rmax'] * 1.5])

        data = self.positions_df[self.positions_df['time'] == 0]
        point, = ax.plot(data.x, data.y, data.z, linestyle="", marker=".")
        plt.cla()

        ani = matplotlib.animation.FuncAnimation(
            fig,
            _update_graph,
            frames=self.iterationsB + self.init_iterations,
            interval=1,
            blit=True
        )

        # ani.save('RotationnoGibsvdrifttimes1.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        fig.show()