from pysb import *
from pysb.integrate import Solver
from pysb.macros import *

Model()

## Monomers
# Tc (b) <- can bond to tetR to make it detach from operon
# tetO (b) <- tet operon can be bonded to tetR and be repressed, or not
# tetR (b, S) S: a, i <- tetR % tetO --> repress OR tetR % Tc + tetC --> allow induce
# gene <- synthesize

## doxycycline analogue of Tc
Monomer('dox', ['t'])                        

# tetO + promoter + transgene
# bound + active = allows gene expression
# unbound + inactive = represses gene expression
Monomer('tetO_pro_gene', ['t', 'S'], {'S':['a', 'i']}) 

# reverse tetR + transactivator domain (rtTA)
# bonding and transactivator site
# active = rtetR isn't bound to dox, so it'll be attached to tetO --> inducing gene expression
# inactive = rtetR is bound to dox, so it'll be unbound from tetO --> repressing gene expression
Monomer('rtTA', ['b', 't', 'S'], {'S':['a', 'i']}) 

# gene to be regulated
# Monomer('gene')

## Rates
# reaction rates
Parameter('kf_rtTA_translation', .001)
Parameter('kf_rtTA_dox_bind', .01)
Parameter('kf_rtTA_tetO_bind', .01)

# degradation/synthesis rates
Parameter('kf_rtTA_deg', .0001)
Parameter('kf_dox_deg', .0001)
Parameter('kf_dox_syn', .7)

## Rules

# dox antibiotic degradation
degrade(dox(t=None), kf_dox_deg)

# rtTA degradation all unbound
degrade(rtTA(b=None, t=None, S='i'), kf_rtTA_deg)
degrade(dox(t=1) % rtTA(b=None, t=1, S='a'), kf_rtTA_deg)

# injection/synthesis of dox
synthesize(dox(t=None), kf_dox_syn) 

# tetR production/translation
synthesize(rtTA(b=None, t=None, S='i'), kf_rtTA_translation)

# Bind dox to rtTA
Rule('dox_rtTA_bind', dox(t=None) + rtTA(b=None, t=None, S='i') >> 
                      dox(t=1) % rtTA(b=None, t=1, S='a'), kf_rtTA_dox_bind)

# Bind dox+rtTA complex to tetO
Rule('rtTA_tetO_bind', rtTA(b=None, t=1, S='a') % dox(t=1) + tetO_pro_gene(t=None, S='i') >> 
                       rtTA(b=2, t=1, S='a') % dox(t=1) % tetO_pro_gene(t=2, S='a'), kf_rtTA_tetO_bind)

# Unbind rtTA if inactivated by absence of dox??

## initial conditions
# Monomer('dox', ['t'])                        
# Monomer('tetO_pro_gene', ['b', 'S'], {'S':['a', 'i']}) 
# Monomer('rtTA', ['b', 't', 'S'], {'S':['a', 'i']}) 

Parameter('dox_0', 0)
Parameter('tetO_pro_gene_0', 500)
Parameter('rtTA_0', 50)

Initial(dox(t=None), dox_0)
Initial(tetO_pro_gene(t=None, S='i'), tetO_pro_gene_0)
Initial(rtTA(b=None, t=None, S='i'), rtTA_0)

## Observables
Observable('obs_dox', dox(t=None))
Observable('obs_translate_ready', rtTA(b=2, t=1, S='a') % dox(t=1) % tetO_pro_gene(t=2, S='a'))
Observable('obs_rtTA_i', rtTA(b=None, t=None, S='i'))
Observable('obs_rtTA_a', dox(t=1) % rtTA(b=None, t=1, S='a'))