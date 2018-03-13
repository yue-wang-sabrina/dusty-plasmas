import itertools

import numpy
import pandas as pd

from tqdm import tqdm

import msci.analysis.constants as const
from msci.particles.dust import Dust
from msci.utils.interpolate import interpolate


class BEffectsAnalysis:
    def __init__(self):
        self.dustdict = None
        self.pairs = None
        self.position = []
        self.positions_df = None
        self.t = []
        self.particles = None
        self.numparticles = None
        self.iterationsB = None
        self.init_iterations = None
        self.position_array = None
        self.modified_b_field = {}
        self.const = const
        self.tindex = 0;
        self.vxbforces=[]
        self.hybridforces=[]
        self.radialEforces=[]
        self.intergrainforces=[]
        self.combinedriftvel=[]
        self.exbdriftvel=[]

    def create_particles(self, numparticles, initpositions):
        # Create dictionary of particles from pickle object

        self.numparticles = numparticles

        names = []
        for int_name in numpy.arange(numparticles):
            names.append('g%s' % int_name)

        self.dustdict = {
            name: Dust(const.md, const.radd, const.lambdaD, const.phia, const.Zd * const.e, [0, 0, 0], [0, 0, 0],
                       [0, 0, 0])
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

    def interact_only(self, method='Gibs', modified_b_field=None):
        positionstemp = []
        pairsfinal = []
        for b in self.pairs:
            if self.dustdict[b[0]].intergraind(self.dustdict[b[1]]):
                pairsfinal.append(b)
            else:
                pass
        for j in pairsfinal:
            interactfield = self.dustdict[j[0]].selffield(self.dustdict[j[1]])
            self.dustdict[j[0]].selffieldmany(interactfield)
            self.dustdict[j[1]].selffieldmany(-interactfield)
        for k in self.dustdict:
            if self.dustdict[k].Bswitch:
                if method == "NoGibs":
                    self.dustdict[k].steptest(numpy.array(self.dustdict[k].multifields))
                    positionstemp.append(self.dustdict[k].getselfpos())
                elif method == "Gibs":
                    Efield = numpy.array(self.dustdict[k].multifields) + interpolate(
                        self.dustdict[k].getselfpos(),
                        modified_b_field['gridr'], modified_b_field['gridz'], modified_b_field['Evalsr'],
                        modified_b_field['Evalsz'], modified_b_field['rmax'], modified_b_field['separationsheath'],
                        modified_b_field['separationhor1'], modified_b_field['separationhor2'],
                        modified_b_field['firstpoint']
                    )
                    self.dustdict[k].steptest(Efield)
                    positionstemp.append(self.dustdict[k].getselfpos())
            else:
                self.dustdict[k].steptest(numpy.array(self.dustdict[k].multifields))
                positionstemp.append(self.dustdict[k].getselfpos())

        self.position.extend(positionstemp)
        self.t.extend(len(positionstemp) * [self.tindex])
        self.tindex += 1

    def interact_and_iterate_drop_method(self, iterationsB, init_iterations, method='Gibs', modified_b_field=None):
        self.iterationsB = iterationsB
        self.init_iterations = init_iterations
        self.modified_b_field = modified_b_field

        names = []
        for int_name in numpy.arange(self.numparticles):
            names.append('g%s' % int_name)

        self.dustdict = {}

        nameindex = 0;

        self.dustdict.update({
            names[nameindex]: Dust(
                const.md, const.radd, const.lambdaD, const.phia,
                const.Zd * const.e,
                [
                    const.lambdaD * numpy.random.randint(-50, 50),
                    const.lambdaD * numpy.random.randint(-50, 50),
                    const.boxz
                ], [0, 0, 0], [0, 0, 0])
        })

        nameindex += 1

        self.dustdict.update({
            names[nameindex]: Dust(
                const.md, const.radd, const.lambdaD, const.phia,
                const.Zd * const.e,
                [
                    const.lambdaD * numpy.random.randint(-50, 50),
                    const.lambdaD * numpy.random.randint(-50, 50),
                    const.boxz
                ], [0, 0, 0], [0, 0, 0])
        })

        nameindex += 1

        for i in tqdm(numpy.arange(init_iterations)):
            if i % 10 == 0 and nameindex <= self.numparticles - 1:
                self.dustdict.update({
                    names[nameindex]: Dust(
                        const.md, const.radd, const.lambdaD, const.phia,
                        const.Zd * const.e,
                        [
                            const.lambdaD * numpy.random.randint(-50, 50),
                            const.lambdaD * numpy.random.randint(-50, 50),
                            const.boxz
                        ],[0, 0, 0], [0, 0, 0])
                })
                nameindex += 1
            self.create_pairs()
            self.interact_only(method, modified_b_field)

        for l in self.dustdict:
            self.dustdict[l].Bswitch = True

        self.create_pairs()

        for k in tqdm(numpy.arange(iterationsB)):
            self.interact_only(method, modified_b_field)

        self.position_array = numpy.array(self.position)

    def sort_posititions_drop_method(self):
        # Sort the positions of the particles for drop method
        self.positions_df = pd.DataFrame({
            "time": self.t,
            "x": self.position_array[:, 0],
            "y": self.position_array[:, 1],
            "z": self.position_array[:, 2]
        })

    def interact_and_iterate(self, iterationsB, init_iterations, method='Gibs', modified_b_field=None):
        # Interact and iterate non-drop method
        self.iterationsB = iterationsB
        self.init_iterations = init_iterations
        self.modified_b_field = modified_b_field

        for i in tqdm(numpy.arange(init_iterations)):
            pairsfinal = []
            for b in self.pairs:
                if self.dustdict[b[0]].intergraind(self.dustdict[b[1]]):
                    pairsfinal.append(b)
                else:
                    pass
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
            MASSDUST = self.dustdict['g0'].m
            dustwanted = self.dustdict[list(self.dustdict.keys())[-1]]
            self.vxbforces.append(MASSDUST*dustwanted.vxBforce())
            self.hybridforces.append(MASSDUST*dustwanted.EXBacchybrid(B=dustwanted.dipoleB(const.dipolepos),method='factor',combinedrifts=True))
            self.radialEforces.append(MASSDUST*dustwanted.radialfield())
            self.combinedriftvel.append(dustwanted.combinedrift(B=dustwanted.dipoleB(r=const.dipolepos)))
            self.exbdriftvel.append(dustwanted.EXBDRIFT(B=dustwanted.dipoleB(r=const.dipolepos)))

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

            self.intergrainforces.append(MASSDUST*numpy.array(dustwanted.multifields))

            for k in self.dustdict:
                if method == 'NoGibs':
                    self.dustdict[k].steptest(numpy.array(self.dustdict[k].multifields))
                    self.position.append(self.dustdict[k].getselfpos())
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

    def sort_positions_of_particles(self):  # works with interact_and_iterate function
        # Sort the positions of the particles of non-drop method

        self.t = numpy.array(
            numpy.sort(
                [list(numpy.arange(int(len(self.position_array[:, 0]) / self.numparticles))) * self.numparticles])
        )[0]

        self.positions_df = pd.DataFrame({
            "time": self.t,
            "x": self.position_array[:, 0],
            "y": self.position_array[:, 1],
            "z": self.position_array[:, 2]
        }
        )
