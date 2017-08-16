from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import shRNA_knockdown as m
import os

t = pl.linspace(0, 20000)
network_f = "shRNA_network"
graph_f = "shRNA_graph"
filename = "shRNA_knockdown.py"

simres = ScipyOdeSimulator(m.model, tspan=t).run()
yout = simres.all

pl.ion()
pl.figure()
pl.plot(t, yout['obs_shRNA'], label="shRNA")
pl.plot(t, yout['obs_siRNA'], label="siRNA")
pl.plot(t, yout['obs_amRNA'], label="active mRNA")
pl.plot(t, yout['obs_imRNA'], label="inactive mRNA")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")
pl.show()

os.system("python -m pysb.tools.render_reactions " + filename + " | dot -T pdf -o " + network_f + "_dot.pdf")
os.system("python -m pysb.tools.render_reactions " + filename + " | circo -T pdf -o " + network_f + "_circo.pdf")
os.system("python -m pysb.tools.render_species " + filename + " | ccomps -x | dot | gvpack -m0 | neato -n2 -T pdf -o " + network_f + "_species.pdf")

pl.savefig(graph_f + ".png", bbox_inches='tight')
