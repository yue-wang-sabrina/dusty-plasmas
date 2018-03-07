# This file is mainly for looking at B field effects of crystals that are already formed.

from msci.analysis.analysis_dust import BEffectsAnalysis
from msci.utils.utils import generate_particle_equilibrium_positions, prepare_modified_b_field
from msci.plots import dustplots

from IPython import get_ipython


ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')


beffect1 = BEffectsAnalysis()

METHOD = "DROP" # "ALLATONCE" OR "DROP"

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
        iterationsB = 500,
        init_iterations = 500,
        method='NoGibs',
        modified_b_field=prepare_modified_b_field())
    beffect1.sort_posititions_drop_method()
    dustplots.pplot(beffect1)


