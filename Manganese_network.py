import decomp_network,plot_pf_output
from numpy import *
from matplotlib import pyplot

pools = [
decomp_network.decomp_pool(name='Cellulose',CN=50,constraints={'initial':1e2},kind='immobile'),
decomp_network.decomp_pool(name='Lignin',CN=50,constraints={'initial':1e2},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='DOM2',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'7.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':1e1}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
decomp_network.decomp_pool(name='Mn+++',kind='primary',constraints={'initial':'1.0e-30'}),
decomp_network.decomp_pool(name='Mn++',kind='primary',constraints={'initial':'1.0e-30'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
# decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='MnO4--',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
decomp_network.decomp_pool(name='MnIIIDOM2(aq)',kind='secondary'),


decomp_network.decomp_pool(name='Birnessite',rate='1.d-8 mol/m^2-sec',constraints={'initial':'1.0d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Mn(OH)3',rate='1.d-8 mol/m^2-sec',constraints={'initial':'1.0d-4  1.d2 m^2/m^3'},kind='mineral'),

decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=1.0e4,complexes=['>Carboxylic_acid']),
]

reactions = [
        decomp_network.reaction(name='Hydrolysis',reactant_pools={'Cellulose':1.0},product_pools={'DOM1':1.0},
                                        rate_constant=1e-1,rate_units='y', 
                                    inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=1e-5)]),
                                    
# Need a way to make a soluble reactive compound from lignin. Treating it like hydrolysis, but allowing only a small amount to be generated
        decomp_network.reaction(name='Lignin activation',reactant_pools={'Lignin':1.0},product_pools={'DOM2':1.0},
                                        rate_constant=1e-1,rate_units='y', 
                                    inhibition_terms=[decomp_network.inhibition(species='DOM2',type='MONOD',k=1e-10)]),
                                    
# Mn reduction turns DOM2 (lignin-based) into DOM1 (odixized, decomposable)
        decomp_network.reaction(name='Lignin Mn reduction',reactant_pools={'DOM2':1.0,'Mn+++':1.0},product_pools={'Mn++':1.0,'H+':1.0,'DOM1':1.0},
                                        monod_terms=[decomp_network.monod(species='DOM2',k=2e-11,threshold=1.1e-20),decomp_network.monod(species='Mn+++',k=1.3e-18,threshold=1.1e-25)],
                                        # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD')],
                                        rate_constant=2e-10,reactiontype='MICROBIAL'),
                                        
# Mn oxidation (Mn peroxidase?)
        decomp_network.reaction(name='Mn oxidation',reactant_pools={'O2(aq)':1.0,'H+':4.0,'Mn++':4.0},product_pools={'Mn+++':4.0},
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-5,threshold=1.1e-15),decomp_network.monod(species='Mn++',k=1.3e-12,threshold=1.1e-15)],
                                        # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD')],
                                        rate_constant=8e-10,reactiontype='MICROBIAL'),
                                    
# Calculating these as per unit carbon, so dividing by the 6 carbons in a glucose
# C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2
        decomp_network.reaction(name='fermentation',reactant_pools={'DOM1':6/6},product_pools={'Acetate-':2/6,'HCO3-':2/6,'H+':4/6+4*2/6,'Tracer':2/6}, # balancing pH of FeIII release requires an extra 5.5 H+ to be released here
                                        rate_constant=1e-10,reactiontype='MICROBIAL', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                    inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD'),decomp_network.inhibition(species='Acetate-',k=6.25e-5,type='MONOD')],
                                    monod_terms=[decomp_network.monod(species='DOM1',k=1e-5,threshold=1.1e-15)]),

# CH2O + H2O -> CO2 + 4H+ + 4 e-
# O2   + 4H+ + 4 e- -> 2H2O
        decomp_network.reaction(name='DOM aerobic respiration',reactant_pools={'DOM1':1.0,'O2(aq)':1.0},product_pools={'HCO3-':1.0,'H+':1.0,'Tracer':1.0},
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-5,threshold=1.1e-12),decomp_network.monod(species='DOM1',k=1e-8,threshold=1.1e-14)],
                                    rate_constant=1.0e-9,reactiontype='MICROBIAL'),

# C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
# 2 O2    + 8 H+ + 8 e- -> 4 H2O
        decomp_network.reaction(name='Acetate aerobic respiration',reactant_pools={'Acetate-':1.0,'O2(aq)':2.0},product_pools={'HCO3-':2.0,'H+':2.0,'Tracer':2.0},
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-5,threshold=1.1e-12),decomp_network.monod(species='Acetate-',k=1e-8,threshold=1.1e-14)],
                                    rate_constant=1.0e-9,reactiontype='MICROBIAL'),
                                    
] # End of reactions list

reaction_network =  decomp_network.decomp_network(pools=pools,reactions=reactions)

import copy
reaction_network_noMn = copy.deepcopy(reaction_network)
reaction_network_noMn.nodes['Mn(OH)3']['constraints']['initial']='0.0d-4  1.d2 m^2/m^3'
reaction_network_noMn.nodes['Birnessite']['constraints']['initial']='0.0d-4  1.d2 m^2/m^3'

networkfig=pyplot.figure('Reaction network');networkfig.clf()
drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','surf_complex','secondary','H+'],
        arrowstyle='-|>',font_size='medium',node_size=1000,font_color='k')


pflotran_exe='../pflotran-interface/src/pflotran/pflotran'
simlength=365

result,units=decomp_network.PF_network_writer(reaction_network).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
result_noMn,units_noMn=decomp_network.PF_network_writer(reaction_network_noMn).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)

def plot_result(result,SOM_ax=None,pH_ax=None,mineral_ax=None,gasflux_ax=None,porewater_ax=None,do_legend=False):

    if SOM_ax is not None:
        l=SOM_ax.plot(result['Cellulose']*1e-3,label='Cellulose')[0]
        SOM_ax.plot(result['Lignin']*1e-3,label='Lignin')[0]

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')
        SOM_ax.set_xlabel('Time (days)')
        if do_legend:
            SOM_ax.legend()

    if pH_ax is not None:
        pH_ax.plot(-log10(result['Free H+']))

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')
        
    if mineral_ax is not None:
        molar_volume=251.1700 # From database. cm3/mol
        molar_weight = 753.5724
        l=mineral_ax.plot(result['Birnessite VF']/molar_volume*1e6   ,label='Birnessite (Mn++)')[0]
        l=mineral_ax.plot(result['Mn(OH)3 VF']/molar_volume*1e6   ,label='Mn(OH)3 (Mn+++)')[0]
        
        # l=mineral_ax.plot(result['Total Mn+++']*result['Porosity']*1e3   ,label='Mn+++',ls='--')[0]
        
        # l=mineral_ax.plot(result['Total Mn++']*result['Porosity']*1e3 ,ls=':'  ,label='Mn++')[0]
        
        mineral_ax.set_title('Mn species')
        mineral_ax.set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
        mineral_ax.set_xlabel('Time (days)')
        if do_legend:
            mineral_ax.legend()
    
    if gasflux_ax is not None:
        gasflux_ax.set_yscale('log')
        
        # l=gasflux_ax.plot(result.index.values[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result.index.values)*1e3,label='CH4')[0]
        
        l=gasflux_ax.plot(result.index.values[:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result.index.values)*1e3,label='CO2',ls='-',c='C5')[0]

        gasflux_ax.set_title('Gas fluxes')
        gasflux_ax.set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')
        gasflux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     gasflux_ax.legend(fontsize='small')
        
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        porewater_ax.plot(result['Total DOM1'],label='DOM (oxidized)')
        porewater_ax.plot(result['Total DOM2'],label='DOM (lignin)')
        porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        porewater_ax.plot(result['Total O2(aq)'],'--',label='O2',c='C4')
        porewater_ax.plot(result['Total Mn+++'],'--',label='Mn+++',c='C1')
        porewater_ax.plot(result['Total Mn++'],':',label='Mn++',c='C2')
        
        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')
        
        if do_legend:
            porewater_ax.legend(fontsize='small',ncol=2)
    
# result,units=plot_pf_output.convert_units(result,units,'mol/m^3')

resultsfig,axes=pyplot.subplots(4,2,num='Results',clear=True,figsize=(4,10))
plot_result(result,SOM_ax=axes[0,0],pH_ax=axes[1,0],mineral_ax=axes[2,0],porewater_ax=axes[3,0],do_legend=True)
plot_result(result_noMn,SOM_ax=axes[0,1],pH_ax=axes[1,1],mineral_ax=axes[2,1],porewater_ax=axes[3,1],do_legend=False)

pyplot.show()
