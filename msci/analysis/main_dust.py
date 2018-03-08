# This file is mainly for looking at B field effects of crystals that are already formed.

import numpy
from msci.analysis.analysis_dust import BEffectsAnalysis
from msci.utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
from msci.plots import dustplots
import msci.dustyplasma_cpp.dustcpp_wrapper as dcpp

from IPython import get_ipython

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')

beffect1 = BEffectsAnalysis()

METHOD = "CPP"  # "ALLATONCE" OR "DROP"

if METHOD == "ALLATONCE":
    beffect1.create_particles(
        numparticles=10,
        initpositions=generate_particle_equilibrium_positions()
    )
    beffect1.create_pairs()
    beffect1.interact_and_iterate(
        iterationsB=1000,
        init_iterations=500,
        method='Gibs',
        modified_b_field=prepare_modified_b_field()
    )
    beffect1.sort_positions_of_particles()
    dustplots.pplot(beffect1)

elif METHOD == "DROP":
    beffect1.numparticles = 10
    beffect1.interact_and_iterate_drop_method(
        iterationsB=500,
        init_iterations=500,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field())
    beffect1.sort_posititions_drop_method()
    dustplots.pplot(beffect1)


elif METHOD == "CPP":
    particlenum = 10;
    itB=1000
    initit = 200
    beffect2 = dcpp.DustAnalysisCpp(initit, itB, particlenum, 1)

    beffect2.get_modified_field()
    beffect2.get_equilibrium_positions()
    beffect2.run()
    beffect1.numparticles = particlenum
    beffect1.iterationsB = itB
    beffect1.init_iterations = initit
    beffect1.method = "Gibs"
    beffect1.position=[[i,j,k] for i,j,k in zip(beffect2.positions_x,beffect2.positions_y,beffect2.positions_z)]
    beffect1.position_array = numpy.array(beffect1.position)
    beffect1.sort_positions_of_particles()
    beffect1.modified_b_field=prepare_modified_b_field()
    dustplots.pplot(beffect1)