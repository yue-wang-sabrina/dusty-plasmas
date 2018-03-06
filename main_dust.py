# This file is mainly for looking at B field effects of crystals that are already formed.

from analysis_dust import BEffectsAnalysis
from utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
from plots import dustplots

from IPython import get_ipython


ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')


beffect1 = BEffectsAnalysis()
beffect1.create_particles(
    numparticles=10,
    initpositions=generate_particle_equilibrium_positions()
)
beffect1.create_pairs()
beffect1.interact_and_iterate(
    iterationsB=5000,
    init_iterations=500,
    method='Gibs',
    modified_b_field=prepare_modified_b_field()
)
beffect1.sort_positions_of_particles()
dustplots.plot(beffect1, False)

