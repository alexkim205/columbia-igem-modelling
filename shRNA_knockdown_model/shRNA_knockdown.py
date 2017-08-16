from pysb import *
from pysb.integrate import Solver
from pysb.macros import *

Model()

# RNA
# maybe i should separate sRNA into separate shRNA and siRNA
Monomer('sRNA', ['b', 'S'], {'S':['sh', 'si']}) # States: shRNA, siRNA
Monomer('mRNA', ['b', 'S'], {'S':['a', 'i']}) # States: active, inactive

# Enzymes
Monomer('Dicer', ['b'])
Monomer('RISC', ['b', 'b_', 'S'], {'S':['a', 'i']}) # States: active, inactive

# Reaction Rates
Parameter('kf_mRNA_deg', 1.0e-03)
Parameter('kf_siRNA_deg', 1.0e-03)
Parameter('kf_dicer', 1.0e-03)
Parameter('kr_dicer', 1.0e-03)
Parameter('kc_dicer', 1.0e-07)
Parameter('kc_siRNA_RISC_bind', 1.0e-02)
Parameter('kc_siRNA_RISC_mRNA_bind', 1.0e-02)
Parameter('kc_cleave', 1)

## Rules

# mRNA degradation
degrade(mRNA(b=None), kf_mRNA_deg)

# siRNA degradation
degrade(sRNA(b=None, S='si'), kf_siRNA_deg)

# Dicer cuts shRNA --> siRNA
# shRNA + Dicer <>(kf, kr) shRNA - Dicer
# shRNA - Dicer -> Dicer + siRNA
catalyze_state(Dicer(), 'b', sRNA(), 'b', 'S', 'sh', 'si', (kf_dicer, kr_dicer, kc_dicer))

# bind siRNA and RISC
# siRNA + RISC -> siRNA - RISC
Rule('siRNA_RISC_bind', RISC(b=None, b_=None) + sRNA(b=None, S='si') >> RISC(b=1, b_=None) % sRNA(b=1, S='si'), kc_siRNA_RISC_bind)

# bind siRNA-RISC complex with mRNA
# siRNA - RISC + mRNA -> mRNA - siRNA - RISC
Rule('siRNA_RISC_mRNA_bind', RISC(b=1, b_=None) % sRNA(b=1, S='si') + mRNA(b=None, S='a') >> RISC(b=1, b_=2) % sRNA(b=1, S='si') % mRNA(b=2, S='a'), kc_siRNA_RISC_mRNA_bind)

# cleave mRNA with siRNA-guided active RISC complex
Rule('cleave_mRNA', RISC(b=1, b_=2) % sRNA(b=1, S='si') % mRNA(b=2, S='a') >> mRNA(b=None, S='i') + RISC(b=1, b_=None) % sRNA(b=1, S='si'), kc_cleave)

# initial conditions
Parameter('sRNA_0', 1000)
Parameter('mRNA_0', 1000)
Parameter('Dicer_0', 10000)
Parameter('RISC_0', 10000)
Initial(sRNA(b=None, S='sh'), sRNA_0)
Initial(mRNA(b=None, S='a'), mRNA_0)
Initial(Dicer(b=None), Dicer_0)
Initial(RISC(b=None, b_=None, S='i'), RISC_0)

# Observables
Observable('obs_shRNA', sRNA(b=None, S='sh'))
Observable('obs_siRNA', sRNA(b=None, S='si'))
Observable('obs_amRNA', mRNA(b=None, S='a'))
Observable('obs_imRNA', mRNA(b=None, S='i'))
#Observable('obs_Dicer', Bid(b=None, S='t'))
#Observable('obs_RISC', Bid(b=None, S='t'))
