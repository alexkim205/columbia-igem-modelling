from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import tetR_decay as m
import os

t = pl.linspace(0, 60000)
network_f = "tetR_network"
graph_f = "tetR_graph"
filename = "tetR_decay.py"

simres = ScipyOdeSimulator(m.model, tspan=t).run()
yout = simres.all

pl.ion()
pl.figure()
pl.plot(t, yout['obs_dox'], label="dox")
pl.plot(t, yout['obs_translate_ready'], label="induced transgenes")
pl.plot(t, yout['obs_rtTA_i'], label="inactive rtTA")
pl.plot(t, yout['obs_rtTA_a'], label="active rtTA")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")
pl.show()

os.system("python -m pysb.tools.render_reactions " + filename + " | dot -T pdf -o " + network_f + "_dot.pdf")
os.system("python -m pysb.tools.render_reactions " + filename + " | circo -T pdf -o " + network_f + "_circo.pdf")
os.system("python -m pysb.tools.render_species " + filename + " | ccomps -x | dot | gvpack -m0 | neato -n2 -T pdf -o " + network_f + "_species.pdf")

pl.savefig(graph_f + ".png", bbox_inches='tight')
