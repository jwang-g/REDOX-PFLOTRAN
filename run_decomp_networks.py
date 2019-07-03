import decomp_network
import os
import plot_pf_output

pflotran_exe='../pflotran-interface/src/pflotran/pflotran'

# Simple decomposition network 1: 1 POM pool, one DOM pool, first order
POM={'name':'POM','CN':15,'initval':1e3}
DOM={'name':'DOM1','CN':15,'initval':1e-3,'immobile':False}
simple_network_DOM1 = decomp_network.decomp_network(pools=[POM,DOM])
simple_network_DOM1.add_reaction('POM','DOM1',CUE=1.0,turnover_val=0.1)

DOM1_result,DOM1_units=simple_network_DOM1.run_simulation('SOMdecomp_template.txt','DOM1',pflotran_exe)

# Simple decomposition network 2: 1 POM pool, one DOM pool, inhibited by DOM
POM={'name':'POM','CN':15,'initval':1e3}
DOM={'name':'DOM1','CN':15,'initval':1e-3,'immobile':False}
simple_network_DOM2 = decomp_network.decomp_network(pools=[POM,DOM])
simple_network_DOM2.add_reaction('POM','DOM1',CUE=1.0,turnover_val=0.1,
            inhibition_terms=[{'species':'DOM1','type':'MONOD','k':1e-2}])

DOM2_result,DOM2_units=simple_network_DOM2.run_simulation('SOMdecomp_template.txt','DOM2',pflotran_exe)

# Microbial decomposition network with priming effects
microbes={'name':'MICROBES','CN':10,'initval':1e-3}
SOM1={'name':'SOM1','CN':10,'initval':1.1e3}
CO2={'name':'CO2','immobile':False}

priming_network = decomp_network.decomp_network(pools=[POM,microbes,SOM1,CO2])
priming_network.add_reaction('POM','MICROBES',CUE=0.6,turnover_val=0.2,monod_terms=[{'species':'MICROBES','k':1e0}])
priming_network.add_reaction('SOM1','MICROBES',CUE=1e-1,turnover_val=10,monod_terms=[{'species':'MICROBES','k':1e0}])
priming_network.add_reaction('MICROBES','CO2',CUE=1.0,turnover_val=0.05,inhibition_terms=[{'species':'MICROBES','k':1e-5,'type':'INVERSE_MONOD'}])
priming_result,priming_units=priming_network.run_simulation('SOMdecomp_template.txt','priming',pflotran_exe)

from pylab import *
import pandas
figure('Simulation results',figsize=(10,8));clf()
ax=subplot(221)
DOM1_result,DOM1_units=plot_pf_output.convert_units(DOM1_result,DOM1_units,'M')
DOM2_result,DOM2_units=plot_pf_output.convert_units(DOM2_result,DOM2_units,'M')
plot(DOM1_result['Total DOM1'],c='C0',label='DOM')
plot(DOM1_result['POM'],c='C1',label='POM')
plot(DOM2_result['Total DOM1'],c='C0',label='DOM (with inhibition k=%1.1g)'%simple_network_DOM2.edges[('POM','DOM1')]['inhibition_terms'][0]['k'],ls='--')
plot(DOM2_result['POM'],c='C1',label='POM (with inhibition k=%1.1g)'%simple_network_DOM2.edges[('POM','DOM1')]['inhibition_terms'][0]['k'],ls='--')
title('DOM simulations')
ylabel('Concentration (M)')
legend()

ax=subplot(222)
plot(priming_result['POM'],c='C1')
plot(priming_result['SOM1'],c='C2')
plot(priming_result['MICROBES'],c='C3')
axins=ax.inset_axes([0.5,0.2,0.4,0.3])
axins.semilogy(priming_result['MICROBES'],c='C3')
title('Priming simulations')
ylabel('Concentration (mol/m3)')
legend()

ax=subplot(223)
pos=decomp_network.nx.layout.circular_layout(simple_network_DOM1)
decomp_network.nx.draw_networkx(simple_network_DOM1,pos=pos,with_labels=True,ax=ax,nodes=simple_network_DOM1.nodes,node_color=decomp_network.make_nodecolors(simple_network_DOM1.nodes),arrowsize=20)

ax=subplot(224)
pos=decomp_network.nx.layout.spectral_layout(priming_network)
decomp_network.nx.draw_networkx(priming_network,pos=pos,with_labels=True,ax=ax,nodes=priming_network.nodes,node_color=decomp_network.make_nodecolors(priming_network.nodes),arrowsize=20)

tight_layout()
show()

# A more complex decomposition network approximating the structure of the CORPSE model

decomp_network_CORPSE=decomp_network.decomp_network()
decomp_network_CORPSE.add_pool('LABILE_POM',CN=    15.0)
decomp_network_CORPSE.add_pool('RESISTANT_POM',CN= 15.0)
decomp_network_CORPSE.add_pool('DOM1',CN=    15.0,immobile=False)
decomp_network_CORPSE.add_pool('DOM2',CN= 15.0,immobile=False)
decomp_network_CORPSE.add_pool('SOIL_MICROBES',CN= 10.0)
decomp_network_CORPSE.add_pool('NECROMASS',CN=     10.0)
decomp_network_CORPSE.add_pool('DOM3',CN= 10.0,immobile=False)
decomp_network_CORPSE.add_pool('MAOM',CN=          10.0)

# Microbial decomposition of unprotected SOM
decomp_network_CORPSE.add_reaction('LABILE_POM','SOIL_MICROBES',CUE=0.6,turnover_val=0.1,turnover_units='y',comment='Labile POM microbial decomposition',
                                monod_terms=[{'species':'SOIL_MICROBES','k':0.1e-3}])
decomp_network_CORPSE.add_reaction('RESISTANT_POM','SOIL_MICROBES',CUE=0.1,turnover_val=0.5,turnover_units='y',comment='Resistant POM microbial decomposition',
                                monod_terms=[{'species':'SOIL_MICROBES','k':0.1e-3}])
decomp_network_CORPSE.add_reaction('NECROMASS','SOIL_MICROBES',CUE=0.6,turnover_val=0.1,turnover_units='y',comment='Necromass microbial decomposition',
                                monod_terms=[{'species':'SOIL_MICROBES','k':0.1e-3}])
                                
# Microbial biomass turnover. Here, CUE would determine maintenance respiration loss as a fraction of turnover
decomp_network_CORPSE.add_reaction('SOIL_MICROBES','NECROMASS',CUE=0.6,turnover_val=0.1,turnover_units='y',
                                comment='Microbial biomass turnover. Inhibition  prevents biomass from declining to zero',
                                inhibition_terms=[{'species':'SOIL_MICROBES','k':0.1e-11,'type':'INVERSE_MONOD'}])
                                
# Physical/mineral protection of necromass. CUE=1.0 because it is an abiotic conservative reaction
decomp_network_CORPSE.add_reaction('NECROMASS','MAOM',CUE=1.0,turnover_val=0.01,turnover_units='y',comment='Physical/mineral protection of necromass',
                                inhibition_terms=[{'species':'MAOM','k':1e-7,'type':'MONOD'}])
# First order turnover of MAOM back to necromass
decomp_network_CORPSE.add_reaction('MAOM','NECROMASS',CUE=1.0,turnover_val=50.0,turnover_units='y',comment='First order turnover of MAOM back to necromass') 

# Dissolution into DOM
decomp_network_CORPSE.add_reaction('LABILE_POM','DOM1',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of labile POM into DOM',
                                inhibition_terms=[{'species':'DOM1','k':1e-12,'type':'MONOD'}])
                                
decomp_network_CORPSE.add_reaction('RESISTANT_POM','DOM2',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of resistant POM into DOM',
                                inhibition_terms=[{'species':'DOM1','k':1e-12,'type':'MONOD'}])
decomp_network_CORPSE.add_reaction('NECROMASS','DOM3',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of necromass into DOM',
                                inhibition_terms=[{'species':'DOM1','k':1e-12,'type':'MONOD'}])      


# CTC decomposition network
decomp_network_CTC=decomp_network.decomp_network()

decomp_network_CTC.add_pool('SOIL1',CN=  12.,initval=1e-10) 
decomp_network_CTC.add_pool('SOIL2',CN=  12.,initval=1e-10)
decomp_network_CTC.add_pool('SOIL3',CN=  10.,initval=1e-10)
decomp_network_CTC.add_pool('SOIL4',CN=  10.,initval=1e-10)
decomp_network_CTC.add_pool('LITR1',initval=1e3,initCN=20)
decomp_network_CTC.add_pool('LITR2',initval=1e-10,initCN=20)    
decomp_network_CTC.add_pool('LITR3',initval=1e-10,initCN=20)
decomp_network_CTC.add_pool('CWD',initval=1e-10,initCN=20)
decomp_network_CTC.add_pool('CO2',immobile=False)

# CWD decomposition to  litter
decomp_network_CTC.add_reaction('CWD','LITR2',CUE=1.0,turnover_val=0.00010*0.76,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION',
                            comment='Writing lossless CTC reactions  (fragmentation) as two reactions with CUE=1.0 and weighted decomposition rates')
decomp_network_CTC.add_reaction('CWD','LITR3',CUE=1.0,turnover_val=0.00010*0.24,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION',
                            comment='Writing lossless CTC reactions  (fragmentation) as two reactions with CUE=1.0 and weighted decomposition rates')

# Litter decomposition
decomp_network_CTC.add_reaction('LITR1','SOIL1',CUE=0.61,turnover_val=1.204,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION')
decomp_network_CTC.add_reaction('LITR2','SOIL2',CUE=0.45,turnover_val=0.0726,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION')
decomp_network_CTC.add_reaction('LITR3','SOIL3',CUE=0.71,turnover_val=0.0141,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION')

# SOM decomposition
decomp_network_CTC.add_reaction('SOIL1','SOIL2',CUE=0.72,turnover_val=0.0726,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION')
decomp_network_CTC.add_reaction('SOIL2','SOIL3',CUE=0.54,turnover_val=0.0141,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION')
decomp_network_CTC.add_reaction('SOIL3','SOIL4',CUE=0.45,turnover_val=0.00141,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION',RATE_AD_FACTOR=10.0)
decomp_network_CTC.add_reaction('SOIL4','CO2',CUE=0.0,turnover_val=0.0001,turnover_units='1/d',turnover_name='RATE_DECOMPOSITION',RATE_AD_FACTOR=100.0)



