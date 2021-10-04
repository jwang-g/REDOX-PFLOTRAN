import decomp_network
import copy


pools = [
decomp_network.decomp_pool(name='SOM',initCN=25,constraints={'initial':1e1},kind='immobile'),
# decomp_network.decomp_pool(name='Lignin',CN=50,constraints={'initial':1e2},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=25,constraints={'initial':1e-15},kind='primary'),
# decomp_network.decomp_pool(name='DOM2',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'6.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'0.2 G O2(g)'}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
# decomp_network.decomp_pool(name='Mn+++',kind='primary',constraints={'initial':'1.0e-30 M Manganite'}),
# decomp_network.decomp_pool(name='Mn++',kind='primary',constraints={'initial':'1.0e-30'}),
decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-10 M Fe(OH)3'}),
decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-15'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-5}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='NO3-',kind='primary',constraints={'initial':1e-5}), 
decomp_network.decomp_pool(name='SO4--',kind='primary',constraints={'initial':1e-5}), 
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='H2S(aq)',kind='secondary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='H2(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='N2(aq)',kind='primary',constraints={'initial':1e-15}),

decomp_network.decomp_pool(name='Na+',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Cl-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='HS-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='N2O(aq)',kind='primary',constraints={'initial':1e-15}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),
decomp_network.decomp_pool(name='N2(g)*',kind='gas'),
# decomp_network.decomp_pool(name='N2O(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
# decomp_network.decomp_pool(name='MnO4--',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='MnIIIDOM2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),

# decomp_network.decomp_pool(name='NH3(aq)',kind='secondary'),
decomp_network.decomp_pool(name='CO3--',kind='secondary'),
# decomp_network.decomp_pool(name='NH4SO4-',kind='secondary'),
# decomp_network.decomp_pool(name='Urea(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='HSO4-',kind='secondary'),
# decomp_network.decomp_pool(name='H2SO4(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='HNO3(aq)    ',kind='secondary'),
# decomp_network.decomp_pool(name='NaNO3(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='NaCl(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='NaSO4-',kind='secondary'),
# decomp_network.decomp_pool(name='NaCO3-',kind='secondary'),
# decomp_network.decomp_pool(name='NaHCO3(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='HCl(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CaCO3(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CaCl+',kind='secondary'),
# decomp_network.decomp_pool(name='CaCl2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CaHCO3+',kind='secondary'),
# decomp_network.decomp_pool(name='CaSO4(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CO(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='S--',kind='secondary'),

# Hui/Beth say pyrolusite probably not precipitating in soils. Look for delta-MnO2? Maybe try 'Mn(OH)2(am)'?
# See Roberts Earth Sci Rev article for some discussion of minerals? More Fe focused though
# decomp_network.decomp_pool(name='Birnessite',rate='0.d-16 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Mn(OH)2(am)',rate='0.d-16 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Manganite',rate='1.d-12 mol/m^2-sec',constraints={'initial':'0.0d-5  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)3',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.1d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Pyrite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Calcite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.875d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Pyrrhotite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),

# This is a workaround to allow pH buffering associated with solid substrate (Rock is a placeholder)
# decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=0.0e4,complexes=['>Carboxylic_acid']),

decomp_network.decomp_pool(name='H2O',kind='implicit'),
]



conc_scales={
    'DOM1':1e-1,
    'Acetate-':1e-1,
    'Fe+++':1e-10,
    'NH4+':1e-6,
    'HCO3-':1e-6,
    'NO3-':1e-6,
    'SO4--':1e-6,
    'CH4(aq)':1e-6,
    'H2(aq)':1e-6,
    'O2(aq)':1e-4,
}

rate_scale=10e-8
thresh=0.0
reactions = [
    # decomp_network.reaction(name='Aerobic decomposition',reactant_pools={'SOM':1.0},product_pools={'HCO3-':1.0},reactiontype='SOMDECOMP',
    #                                         rate_constant=1e-8,rate_units='1/sec',turnover_name='RATE_CONSTANT', 
    #                                         # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='THRESHOLD 1.0d20')]),
    #                                     monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)]),
                                        
                                        
    decomp_network.reaction(name='Hydrolysis',stoich='1.0 SOM -> 0.9 DOM1',reactiontype='SOMDECOMP',turnover_name='RATE_CONSTANT',
                                            rate_constant=rate_scale,rate_units='1/sec', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=conc_scales['DOM1']),
                                                          # decomp_network.inhibition(species='O2(aq)',type='MONOD',k=1e-11),
                                                          ]),
    
    # Calculating these as per unit carbon, so dividing by the 6 carbons in a glucose
    # C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2
    # Should it be inhibited by H2?
    decomp_network.reaction(name='fermentation',reactant_pools={'DOM1':6/6},product_pools={'Acetate-':2/6,'HCO3-':2/6,'H+':4/6,'H2(aq)':4*2/6,'Tracer':2/6}, # balancing pH of FeIII release requires an extra 5.5 H+ to be released here
                                            rate_constant=rate_scale,reactiontype='MICROBIAL', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='Acetate-',k=conc_scales['Acetate-'],type='MONOD')],
                                        monod_terms=[decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=thresh)]),

    # CH2O + H2O -> CO2 + 4H+ + 4 e-
    # O2   + 4H+ + 4 e- -> 2H2O
    decomp_network.reaction(name='DOM oxidation (O2)',stoich='1.0 DOM1 + 1.0 O2(aq) -> 1.0 HCO3- + 1.0 H+ + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),

    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 2 O2    + 8 H+ + 8 e- -> 4 H2O
    decomp_network.reaction(name='Acetate oxidation (O2)',stoich='1.0 Acetate-  + 2.0 O2(aq)  -> 2.0 HCO3-  + 2.0 H+  + 2.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),


    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 8 Fe+++ + 8 e- -> 8 Fe++ 
    decomp_network.reaction(name='Fe(III) reduction',stoich='1.0 Acetate- + 8.0 Fe+++ -> 2.0 HCO3- + 8.0 Fe++ + 9.0 H+ + 2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='Fe+++',k=conc_scales['Fe+++'],threshold=thresh)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD')],
                                            rate_constant=rate_scale*0.5,reactiontype='MICROBIAL'),
                                            
    # Oxidation of Fe++
    # Currently assuming backward reaction rate is zero
    # decomp_network.reaction(name='Fe(II) oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O',
    #                                         rate_constant=1e-2,backward_rate_constant=0.0,reactiontype='GENERAL'),

    # Acetoclastic methanogenesis
    # C2H3O2- + H+ -> CH4 + HCO3- + H+
    decomp_network.reaction(name='Acetoclastic methanogenesis',stoich='1.0 Acetate- -> 1.0 CH4(aq) + 1.0 HCO3- + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
                                            decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD'),
                                            decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                            decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')],
                                            rate_constant=rate_scale*0.1,reactiontype='MICROBIAL'),
                                            
    # Hydrogenotrophic methanogenesis
    decomp_network.reaction(name='Hydrogenotrophic methanogenesis',stoich='4.0 H2(aq) + 1.0 HCO3- + 1.0 H+ -> 1.0 CH4(aq) + 3.0 H2O',
                                            monod_terms=[decomp_network.monod(species='H2(aq)',k=conc_scales['H2(aq)'],threshold=thresh),decomp_network.monod(species='HCO3-',k=conc_scales['HCO3-'],threshold=thresh)],
                                            # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
                                            decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD'),
                                            decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                            decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')],
                                            rate_constant=rate_scale*0.1,reactiontype='MICROBIAL'),
    
    # H2 oxidation if oxygen available

    # Nitrification
    decomp_network.reaction(name='Nitrification',stoich='2.0 NH4+ + 4.0 O2(aq) -> 2.0 NO3- + 2.0 H2O + 4.0 H+',
                        monod_terms=[decomp_network.monod(species='NH4+',k=conc_scales['NH4+'],threshold=thresh),decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)],
                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
                        
    # Denitrification C2H3O2- + 2 NO3- + H+ -> 2 HCO3- + N2 + 2 H2O
    decomp_network.reaction(name='Denitrification',stoich='1.0 Acetate- + 2.0 NO3- + 1.0 H+ -> 2.0 HCO3- + 1.0 N2(aq) + 0.0 N2O(aq) + 2.0 H2O + 2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD')],
                                            rate_constant=rate_scale*0.8,reactiontype='MICROBIAL'),
                                            
    # Sulfate reduction C2H3O2- + SO4--  -> 2 HCO3- + HS-
    decomp_network.reaction(name='Sulfate reduction',stoich='1.0 Acetate- + 1.0 SO4--  -> 2.0 HCO3- + 2.0 HS- +  2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='SO4--',k=conc_scales['SO4--'],threshold=thresh)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
                                                                decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
                                                                decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                            rate_constant=rate_scale*0.3,reactiontype='MICROBIAL'),
                                            
                                            
    # Methane oxidation (O2)
    decomp_network.reaction(name='Methane oxidation (O2)',stoich='1.0 CH4(aq)  + 1.0 O2(aq)  -> 1.0 HCO3-  + 1.0 H+ + 1.0 H2O + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
    
    # Methane oxidation (NO3)
    decomp_network.reaction(name='Methane oxidation (NO3)',stoich='1.0 CH4(aq)  + 1.0 NO3-  -> 1.0 HCO3-  + 1.0 NH4+ + 1.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
    
    # Methane oxidation (SO4)
    decomp_network.reaction(name='Methane oxidation (SO4)',stoich='1.0 CH4(aq)  + 1.0 SO4-- -> 1.0 HCO3-  + 1.0 HS- + 1.0 H2O + 1.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='SO4--',k=conc_scales['SO4--'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
    
    # Methane oxidation (Fe)
    decomp_network.reaction(name='Methane oxidation (Fe)',stoich='1.0 CH4(aq)  + 8.0 Fe+++ + 3.0 H2O -> 1.0 HCO3-  + 8.0 Fe++ + 9.0 H+ + 1.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='Fe+++',k=conc_scales['Fe+++'],threshold=thresh),decomp_network.monod(species='CH4(aq)',k=conc_scales['CH4(aq)'],threshold=thresh)],
                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),

]

def yield_calculation(substrates,products):
    return 0.1

class microbe:
    def __init__(self,genes,size,CN,name):
        self.size=size
        self.genes=genes
        self.CN=CN
        self.name=name

    def make_reaction(self,react):
        reaction=copy.deepcopy(react)
        if 'HCO3-' in reaction['product_pools']:
            C_out='HCO3-'
        elif 'DOM1' in reaction['product_pools']:
            C_out='DOM1'
        else:
            raise ValueError('Products must include HCO3- or DOM1 to substract C from')
        # Not taking DOM C:N into account currently
        # This needs to check if the reaction is SOMDECOMP in which case representation of microbial products will be different
        microbe_yield = yield_calculation(reaction['reactant_pools'],reaction['product_pools'])
        if react['reactiontype'] == 'MICROBIAL':
            reaction['product_pools'][C_out]=reaction['product_pools'][C_out]-microbe_yield # Subtract microbe C content
            reaction['product_pools']['NH4+']=reaction['product_pools'].get('NH4+',0)+microbe_yield/self.CN
            # reaction['rate_constant']= 1.0 # Rate constant depending on microbial size, growth rate, gene copy, ...???
            reaction['biomass']=self.name
            reaction['biomass_yield']=microbe_yield
        elif react['reactiontype'] == 'SOMDECOMP':
            reaction['product_pools'][C_out]=reaction['product_pools'][C_out]-microbe_yield # Subtract microbe C content
            reaction['product_pools'][self.name]=microbe_yield
            # reaction['rate_constant']= 1.0 # Rate constant depending on microbial size, growth rate, gene copy, ...???
        else:
            raise ValueError('Microbe genes only defined for microbial reactions currently')
        reaction['name']=reaction['name']+f' ({self.name})'

        return reaction




# See Smeaton, Christina M., and Philippe Van Cappellen. 2018. “Gibbs Energy Dynamic Yield Method (GEDYM): Predicting Microbial Growth Yields under Energy-Limiting Conditions.” Geochimica et Cosmochimica Acta 241 (November): 1–16. https://doi.org/10.1016/j.gca.2018.08.023.
microbes=[
    microbe(genes=[reactions[0],reactions[1]],size=1.0,CN=10.0,name='microbe1'),
    microbe(genes=[reactions[1],reactions[2]],size=2.0,CN=12.0,name='microbe2'),
    microbe(genes=[reactions[5],reactions[4]],size=2.0,CN=12.0,name='microbe3')
]

initial_biomass=1e-3
microbe_reactions=[]
for m in microbes:
    for r in m.genes:
        microbe_reactions.append(m.make_reaction(r))
    # Microbial biomass pools need to be included in immobile pools
    pools.append(decomp_network.decomp_pool(name=m.name,CN=m.CN,constraints={'initial':initial_biomass},kind='immobile'))
    # Add microbial biomass decay
    microbe_reactions.append( decomp_network.reaction(name=f'{m.name} turnover',stoich=f'1.0 {m.name} -> 0.5 DOM1',reactiontype='SOMDECOMP',turnover_name='RATE_CONSTANT',
            rate_constant=1e-6,rate_units='1/sec',
            monod_terms=[decomp_network.monod(species=m.name,type='MONOD',k=1e-8)], # Monod dependence on biomass prevents it from exponentially decaying without limit
                                                          ))

network=decomp_network.decomp_network(pools,microbe_reactions)

from matplotlib import pyplot
networkfig=pyplot.figure('Reaction network',clear=True,figsize=(7.1,8.9))


drawn,pos=decomp_network.draw_network_with_reactions(network,
        omit=['gas','secondary','H+','>Carboxylate-','Carboxylic_acid','H2O','Na+','Cl-','Calcite','Ca++','Nitrification','NH4+','mineral','HCO3-'],
        font_size='xx-large',node_size=1800,arrowstyle='-|>',arrowsize=20.0,width=1.5,edge_color='k',node_alpha=0.9,markers={'Reaction':'h','mineral':'8'},
        namechanges={'SOM':'SOM','DOM1':'DOM','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
                     'H2(aq)':'H$_2$','DOM oxidation (O2)':'Aerobic\nresp','Acetate oxidation (O2)':'Aerobic\nresp','Aerobic decomposition':'Aerobic\nresp',
                     'NH4+':'NH$_4^+$','NO3-':'NO$_3^-$','N2O(aq)':'N$_2$O','N2(aq)':'N$_2$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$',
                     'SO4--':'SO$_4^{--}$','HS-':'H$_2$S','Methane oxidation (O2)':'Methane\noxidation','Methane oxidation (NO3)':'Methane\noxidation',
                     'Methane oxidation (Fe)':'Methane\noxidation','Methane oxidation (SO4)':'Methane\noxidation','Fe(III) reduction':'Fe(III)\nreduction','Sulfate reduction':'Sulfate\nreduction',
                     'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Acetoclastic methanogenesis':'Acetoclastic\nmethanogenesis',
                     'SOMdecomp Reaction':'SOM Reaction','General Reaction':'Abiotic Reaction','fermentation':'Fermentation','Primary aqueous':'Dissolved ion','Gas':'Dissolved gas'},connectionstyle='arc3, rad=0.2')

decomp_network.PF_network_writer(network).write_into_input_deck('SOMdecomp_template.txt','microbial_test_network.in',log_formulation=True,CO2name='HCO3-',truncate_concentration=1e-25,
                    database='/home/b0u/models/PFLOTRAN/REDOX-PFLOTRAN/hanford.dat',verbose=True,length_days=100)

import run_alquimia
result_alquimia,units_alq=run_alquimia.run_simulation('microbial_test_network.in',hands_off=True,simlength_days=10,dt=3600,initcond=pools,rateconstants={},truncate_concentration=0.0)


pyplot.show()