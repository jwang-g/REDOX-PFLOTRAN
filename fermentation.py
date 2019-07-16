import decomp_network
import plot_pf_output

pflotran_exe='../pflotran-interface/src/pflotran/pflotran'

cellulose = decomp_network.decomp_pool(name='cellulose',CN=50,initval=1e3)
aqueous_waste = decomp_network.decomp_pool(name='DOM1',CN=50,initval=1e-15,immobile=False)
CO2=decomp_network.decomp_pool(name='CO2(aq)',immobile=False)
Hplus=decomp_network.decomp_pool(name='H+',immobile=False)
O2aq=decomp_network.decomp_pool(name='O2(aq)',immobile=False,initval=1e-4)
HCO3=decomp_network.decomp_pool(name='HCO3-',immobile=False)

ferm_hydrolysis = decomp_network.reaction(name='Fermentation hydrolysis',reactant_pools={'cellulose':1.0},product_pools={aqueous_waste['name']:0.67,'HCO3-':0.33},
                                        rate_constant=1e-1,rate_units='y', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                    inhibition_terms=[decomp_network.inhibition(species=aqueous_waste['name'],type='MONOD',k=1e-5)])

DOM_resp = decomp_network.reaction(name='DOM aerobic respiration',reactant_pools={aqueous_waste['name']:1.0,'O2(aq)':1.0},product_pools={'HCO3-':1.0,'H+':1.0},
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12),decomp_network.monod(species=aqueous_waste['name'],k=1e-8)],
                                    rate_constant=1.0e-5,reactiontype='MICROBIAL')

fermentation_network =  decomp_network.decomp_network(pools=[cellulose,aqueous_waste,Hplus,O2aq,HCO3],reactions=[ferm_hydrolysis,DOM_resp])

result,units=decomp_network.PF_network_writer(fermentation_network).run_simulation('SOMdecomp_template.txt','fermentation',pflotran_exe,print_output=False)

from pylab import *
result,units=plot_pf_output.convert_units(result,units,'M')
figure('Network diagram');clf()
pos=decomp_network.nx.drawing.nx_agraph.graphviz_layout(fermentation_network,prog='dot')
ax=subplot(212)
decomp_network.nx.draw_networkx(fermentation_network,pos=pos,with_labels=True,ax=ax,nodes=fermentation_network.nodes,node_color=decomp_network.make_nodecolors(fermentation_network.nodes),arrowsize=20)

subplot(211)
plot(result['cellulose'])
plot(result['Total DOM1'])
plot(result['Total HCO3-'])
legend()
show()