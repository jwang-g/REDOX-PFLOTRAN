## this code is to build reaction network for delta marsh site, and can be called by deltmarsh_profile code.

import decomp_network,plot_pf_output
from numpy import *
from matplotlib import pyplot
import matplotlib

pools = [
decomp_network.decomp_pool(name='SOM',CN=13,constraints={'initial':2.8},kind='immobile'),  #ref 60 mg/cm3 organic matter to C-mol/L at Wax lake McClellan et al., 2021
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),
#decomp_network.decomp_pool(name='Root_biomass',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=13,constraints={'initial':1e-15,'sed_air_interface':6.5e-4,'excretion':1e-20,'tide':1e-2,'sed_water_interface':1e-3},kind='primary'),
#decomp_network.decomp_pool(name='DOM2',CN=50,constraints={'initial':1e-30},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'6.0 P','sed_air_interface':'6.0 P','excretion':'6.0 P','tide':'7.5 P','sed_water_interface':'6.5 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'1.0 G O2(g)','sed_air_interface':1e-3,'excretion':1e-20,'tide':1e-2,'sed_water_interface':1e-5}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)','sed_air_interface':'200e-7 G CO2(g)','excretion':'100e-7 G CO2(g)','tide':'300e-7 G CO2(g)','sed_water_interface':'250e-7 G CO2(g)'}),
decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-10 M Fe(OH)3','sed_air_interface':1e-6,'excretion':1e-20,'tide':1e-20,'sed_water_interface':1e-15}),
decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-3','sed_air_interface':1e-8,'excretion':1e-20,'tide':1e-20,'sed_water_interface':1e-10}),

decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-10,'sed_air_interface':1e-10,'excretion':1e-20,'tide':1e-10,'sed_water_interface':1e-10}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='NO3-',kind='primary',constraints={'initial':3e-10,'sed_air_interface':1e-10,'excretion':1e-20,'tide':1e-10,'sed_water_interface':1e-10}),
decomp_network.decomp_pool(name='SO4--',kind='primary',constraints={'initial':1e-4,'sed_air_interface':1e-3,'excretion':1e-20,'tide':2e-3,'sed_water_interface':1.5e-3}),
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':7e-4,'sed_air_interface':6.5e-4,'excretion':1e-20,'tide':1e-2,'sed_water_interface':1e-3}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-4,'sed_air_interface':1e-4,'excretion':1e-20,'tide':1e-15,'sed_water_interface':1e-14}),
decomp_network.decomp_pool(name='H2S(aq)',kind='secondary',constraints={'initial':1e-35,'sed_air_interface':1e-36,'excretion':1e-20,'tide':1e-8,'sed_water_interface':1.5e-6}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15,'sed_air_interface':6.5e-4,'excretion':1e-20,'tide':1e-2,'sed_water_interface':1e-3}),
decomp_network.decomp_pool(name='H2(aq)',kind='primary',constraints={'initial':1e-6,'sed_air_interface':1e-4,'excretion':1e-20,'tide':1e-4,'sed_water_interface':1e-4}),
decomp_network.decomp_pool(name='N2(aq)',kind='primary',constraints={'initial':1e-6,'sed_air_interface':1e-4,'excretion':1e-20,'tide':1e-4,'sed_water_interface':1e-4}),
decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':0.5e-3,'sed_air_interface':0.5e-3,'excretion':1e-20,'tide':1e-1,'sed_water_interface':1e-2}),
decomp_network.decomp_pool(name='Na+',kind='primary',constraints={'initial':2e-3,'sed_air_interface':1.5e-3,'excretion':1e-20,'tide':2e-4,'sed_water_interface':1e-4}),
decomp_network.decomp_pool(name='Cl-',kind='primary',constraints={'initial':2e-3,'sed_air_interface':1.5e-3,'excretion':1e-20,'tide':2e-4,'sed_water_interface':1e-4}),
decomp_network.decomp_pool(name='HS-',kind='primary',constraints={'initial':1e-35,'sed_air_interface':1e-35,'excretion':1e-20,'tide':1e-15,'sed_water_interface':1e-15}),
decomp_network.decomp_pool(name='N2O(aq)',kind='primary',constraints={'initial':1e-15,'sed_air_interface':1e-15,'excretion':1e-20,'tide':1e-15,'sed_water_interface':1e-15}),


decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),
decomp_network.decomp_pool(name='N2(g)*',kind='gas'),
decomp_network.decomp_pool(name='H2O',kind='implicit'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='CO3--',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3+',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)4-',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeCH3COO+',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3(aq)',kind='secondary'),


decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-6 mol/m^2-sec',constraints={'initial':'0.1d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Pyrite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),
#decomp_network.decomp_pool(name='Calcite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.875d-3  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Pyrrhotite',rate='0.d-3 mol/m^2-sec',constraints={'initial':'0.0d-3  1.d2 m^2/m^3'},kind='mineral'),

]


# Fraction of Mn+++ that is not recycled in Mn-Peroxidase enzyme loop
##cjw make network reaction list below for marsh biogeochemistry. conc_scales are the half saturations for the primary species.
def make_network(change_constraints={},change_rate={}):
    ##cjw: below is the original setting for half saturation
    #conc_scales={
    #	'DOM1':1e-2,
    #	'Acetate-':1e-3,
    #	'Fe+++':1e-10,
    #    'Fe++':1e-1,
    #	'NH4+':1e-5,
    #	'HCO3-':1e-2,
    #	'NO3-':1e-6,
    #	'SO4--':1e-5,#1e-7,#(3166value),     #change here to test methane flux
    #    'HS-':1e-4,
    #	'CH4(aq)':1e-5,
    #	'H2(aq)':1e-5,
    #	'O2(aq)':1e-4,
    #}
    ##cjw new half saturation values for each reaction from Egger et al., 2016, Biogeosciences. Reed et al., gene-centric model, PNNL.
###half saturation unit is M.
    conc_scales={
    	'DOM1':1e-2,#1e-2,#7e-7,       #reed et al, PNAS
    	'Acetate-':1e-3,#1e-3,#7e-6,
    	'Fe+++':1e-10,    ##
        'Fe++':1e-1,
    	'NH4+':1e-5,      ##range is 57.1-107 uM
    	'HCO3-':1e-2,     #
    	'NO3-':2e-6,      # range is 0.002-0.04 mM
    	'SO4--':1e-3,#1e-5,#1e-7,#(3166value),     #range is 0.8 mM
        'HS-':1e-4,       #anaerobic oxidation is 121nM
    	'CH4(aq)':1e-5,
    	'H2(aq)':1e-5,
    	'O2(aq)':1e-4,    ## literature range is 0.0005-0.015 mM
    }
    rate_scale=2e-10
    truncate_conc=1e-30
    thresh=truncate_conc*1.01

    reactions = [
            decomp_network.reaction(name='Aerobic decomposition',reactant_pools={'SOM':1.0},product_pools={'HCO3-':1.0},reactiontype='SOMDECOMP',
	                                            rate_constant=1e-6,rate_units='1/s',turnover_name='RATE_CONSTANT',
                                            # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='THRESHOLD 1.0d20')]),
	                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)]),
            decomp_network.reaction(name='Hydrolysis',stoich='1.0 SOM -> 1.0 DOM1',reactiontype='SOMDECOMP',turnover_name='RATE_CONSTANT',
	                                            rate_constant=1e-7,rate_units='1/s', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
	                                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=conc_scales['DOM1']),
	                                                          decomp_network.inhibition(species='O2(aq)',k=6.25e-11,type='THRESHOLD -1.0d10')
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
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),
        # Oxidation of Fe++
        # Currently assuming backward reaction rate is zero
        decomp_network.reaction(name='Fe(II) abiotic oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O',
                                                rate_constant=1e-2,backward_rate_constant=0.0,reactiontype='GENERAL'),

        decomp_network.reaction(name='Fe(II) microbial oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ -> 1.0 Fe+++ + 0.5 H2O',
                                                monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh),decomp_network.monod(species='Fe++',k=conc_scales['Fe++'],threshold=thresh)],
                                                rate_constant=1e-8,reactiontype='MICROBIAL'),
	    # Oxidation of Fe++
	    # Currently assuming backward reaction rate is zero
#	    decomp_network.reaction(name='Fe(II) oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ -> 1.0 Fe+++ + 0.5 H2O',
#                                                monod_terms=[decomp_network.monod(species='Fe++',k=conc_scales['Fe++'],threshold=thresh),decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)],
#	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),

	    # Acetoclastic methanogenesis
	    # C2H3O2- + H+ -> CH4 + HCO3- + H+
	    decomp_network.reaction(name='Acetoclastic methanogenesis',stoich='1.0 Acetate- -> 1.0 CH4(aq) + 1.0 HCO3- + 1.0 Tracer',
	                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),],
	                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
	                                            decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD'),
	                                            decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
	                                            decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')],
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),

	    # Hydrogenotrophic methanogenesis
	    decomp_network.reaction(name='Hydrogenotrophic methanogenesis',stoich='4.0 H2(aq) + 1.0 HCO3- + 1.0 H+ -> 1.0 CH4(aq) + 3.0 H2O',
	                                            monod_terms=[decomp_network.monod(species='H2(aq)',k=conc_scales['H2(aq)']*10,threshold=thresh),decomp_network.monod(species='HCO3-',k=conc_scales['HCO3-']*10,threshold=thresh)],
	                                            # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
	                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
	                                            decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD'),
	                                            decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
	                                            decomp_network.inhibition(species='SO4--',k=conc_scales['SO4--'],type='MONOD')],
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),

	    # H2 oxidation if oxygen available
	    # Nitrification
	    decomp_network.reaction(name='Nitrification',stoich='2.0 NH4+ + 4.0 O2(aq) -> 2.0 NO3- + 2.0 H2O + 4.0 H+',
	                        monod_terms=[decomp_network.monod(species='NH4+',k=conc_scales['NH4+'],threshold=thresh),decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=thresh)],
	                        rate_constant=rate_scale,reactiontype='MICROBIAL'),

	    # Denitrification C2H3O2- + 2 NO3- + H+ -> 2 HCO3- + N2 + 2 H2O
	    decomp_network.reaction(name='Denitrification',stoich='1.0 Acetate- + 2.0 NO3- + 1.0 H+ -> 2.0 HCO3- + 1.0 N2(aq) + 0.0 N2O(aq) + 2.0 H2O + 2.0 Tracer',
	                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh)],
	                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD')],
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),
        ##added by Jiaze MMMM
        # Fermentative DNRA  CH2O + 0.5 NO3- +0.5 H2O= 0.5NH4+ + HCO3- or C2H3O2- +  NO3- + H2O + H+ -> 2 HCO3- + NH4+
        #decomp_network.reaction(name='Fermentative DNRA abiotic',stoich='1.0 DOM1 + 0.5 NO3- + 0.5 H2O <-> 0.5 NH4+ + 1.0 HCO3-',
        #                                        rate_constant=1e-2,backward_rate_constant=0.0,reactiontype='GENERAL'),

	    decomp_network.reaction(name='Fermentative DNRA',stoich='1.0 DOM1 + 0.5 NO3- + 0.5 H2O -> 0.5 NH4+ + 1.0 HCO3-',
                                                monod_terms=[decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh),decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=thresh)],
                                                inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD')],
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),
        #decomp_network.reaction(name='Fermentative DNRA',stoich='1.0 Acetate- + 1.0 NO3- + 1.0 H2O + 1.0 H+ -> 1.0 NH4+ + 2.0 HCO3- + 2.0 Tracer',
        #                                        monod_terms=[decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh),decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh)],
        #                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD')],
	    #                                        rate_constant=rate_scale,reactiontype='MICROBIAL'),
        # Autotrophic DNRA  HS- + NO3- + H+ + H2O = SO4-- + NH4+
        #decomp_network.reaction(name='Autotrophic DNRA abiotic',stoich='1.0 HS- + 1.0 NO3- + 1.0 H+ + 1.0 H2O <-> 1.0 SO4-- + 1.0 NH4+',
        #                                        rate_constant=1e-2,backward_rate_constant=0.0,reactiontype='GENERAL'),

	    decomp_network.reaction(name='Autotrophic DNRA',stoich='1.0 HS- + 1.0 NO3- + 1.0 H+ + 1.0 H2O -> 1.0 SO4-- + 1.0 NH4+',
                                                monod_terms=[decomp_network.monod(species='NO3-',k=conc_scales['NO3-'],threshold=thresh),decomp_network.monod(species='HS-',k=conc_scales['HS-'],threshold=thresh)],
                                                inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),decomp_network.inhibition(species='CH4(aq)',k=conc_scales['CH4(aq)'],type='MONOD')],
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),
        ##added by Jiaze wwww

	    # Sulfate reduction C2H3O2- + SO4--  -> 2 HCO3- + HS-
	    decomp_network.reaction(name='Sulfate reduction',stoich='1.0 Acetate- + 1.0 SO4--  -> 2.0 HCO3- + 2.0 HS- +  2.0 Tracer',
	                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=thresh),decomp_network.monod(species='SO4--',k=conc_scales['SO4--'],threshold=thresh)],
	                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=conc_scales['O2(aq)'],type='MONOD'),
	                                                                decomp_network.inhibition(species='NO3-',k=conc_scales['NO3-'],type='MONOD'),
	                                                                decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
	                                            rate_constant=rate_scale,reactiontype='MICROBIAL'),
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


    ] # End of reactions list

    import copy
    pools_copy=copy.deepcopy(pools)
    for n,p in enumerate(pools_copy):
        if p['name'] in change_rate:
            pools_copy[n]['rate']=change_rate[p['name']]
    return decomp_network.decomp_network(pools=decomp_network.change_constraints(pools_copy, change_constraints),reactions=reactions)


def plot_result(result,SOM_ax=None,pH_ax=None,mineral_ax=None,gasflux_ax=None,porewater_ax=None,do_legend=False):

    if SOM_ax is not None:
#        l=SOM_ax.plot(result['Cellulose']*1e-3,label='Cellulose')[0]
#        SOM_ax.plot(result['Lignin']*1e-3,label='Lignin')[0]
        l=SOM_ax.plot(result['SOM']*1e-3,label='SOM')[0]

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')
        SOM_ax.set_xlabel('Time (days)')
#        if do_legend:
#            SOM_ax.legend()

    if pH_ax is not None:
        pH_ax.plot(-log10(result['Free H+']))

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')

#    if mineral_ax is not None:
#        molar_volume_birnessite=251.1700 # From database. cm3/mol (first number)
#        molar_weight_birnessite = 753.5724 # (last number)
#        molar_volume_manganite = 24.45
#        molar_volume_pyrolusite = 500.0 # From database. Seems fishy
#        l=mineral_ax.plot(result['Birnessite2 VF']/molar_volume_birnessite*1e6   ,label='Birnessite')[0]
#        # l=mineral_ax.plot(result['Manganite VF']/molar_volume_manganite*1e6   ,label='Manganite (Mn+++)')[0]

#        # l=mineral_ax.plot(result['Total Mn+++']*result['Porosity']*1e3   ,label='Mn+++',ls='--')[0]

#        # l=mineral_ax.plot(result['Total Mn++']*result['Porosity']*1e3 ,ls=':'  ,label='Mn++')[0]

#        mineral_ax.set_title('Mn minerals')
#        mineral_ax.set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
#        mineral_ax.set_xlabel('Time (days)')
#        if do_legend:
#            mineral_ax.legend()
    if Fe_ax is not None:
        molar_volume=34.3600 # From database. cm3/mol
        molar_weight = 106.8690
        l=Fe_ax.plot(result['Fe(OH)3 VF']/molar_volume*1e6   ,label='Fe(OH)3')[0]

        # M/L to umol/cm3: 1e6/1e3=1e3
        l=Fe_ax.plot(result['Total Fe+++']*result['Porosity']*1e3   ,label='Fe+++',ls='--')[0]

        l=Fe_ax.plot(result['Total Fe++']*result['Porosity']*1e3 ,ls=':'  ,label='Fe++')[0]

        Fe_ax.set_title('Fe species')
        Fe_ax.set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
        Fe_ax.set_xlabel('Time (days)')
        if do_legend:
            Fe_ax.legend(fontsize='small')

    if gasflux_ax is not None:
        gasflux_ax.set_yscale('log')

        l=gasflux_ax.plot(result.index.values[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result.index.values)*1e3,label='CH4')[0]

        l=gasflux_ax.plot(result.index.values[:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result.index.values)*1e3,label='CO2',ls='-',c='C5')[0]

        gasflux_ax.set_title('Gas fluxes')
        gasflux_ax.set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')
        gasflux_ax.set_xlabel('Time (days)')
        if do_legend:
            gasflux_ax.legend(fontsize='small')

    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        # porewater_ax.plot(result['Total DOM1'],label='DOM (oxidized)')
        # porewater_ax.plot(result['Total DOM2'],label='DOM (lignin)')
        # porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        # porewater_ax.plot(result['Total O2(aq)'],'--',label='O2',c='C4')

        # porewater_ax.plot(result['Total DOM1'],label='DOM')
        # porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        porewater_ax.plot(result['Total O2(aq)'],'-',label='O2',c='C0')
        porewater_ax.plot(result['Total NO3-'],'-',c='C1',label='NO3-')
        porewater_ax.plot(result['Total N2(aq)'],'--',c='C1',label='N2')
        porewater_ax.plot(result['Total Fe+++'],'-',label='Fe+++',c='C2')
        # porewater_ax.plot(result['Free Fe+++'],':',label='Fe+++',c='C1')
        porewater_ax.plot(result['Total Fe++'],'--',label='Fe++',c='C2')

        porewater_ax.plot(result['Total SO4--'],'-',c='C3',label='SO4--')
        porewater_ax.plot(result['Total HS-'],'--',c='C3',label='H2S')
        porewater_ax.plot(result['Total CH4(aq)'],'--',c='C4',label='CH4')

        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')

        if do_legend:
            porewater_ax.legend(fontsize='small')

# All of this stuff is for the plot of the network diagram
mineral_col=265.0
TEA_col=175.0
reduction_col=150.0
oxidation_col=200.0
substrateox_col=reduction_col
substrate_col=130.0
gas_col=TEA_col

redox_seq={
    'O2':1250,
    'NO3-':900,
    # 'Mn4+':700,
    'HS-':700,
    'Fe3+':600,
    'SO4--':400,
    'CO2':200,
    'CH4':100,
    'SOM':1400,
    'DOM':1200,
    'Acetate':900,
    'H2':600
    }

pos={'Pyrrhotite': (mineral_col, (redox_seq['Fe3+']+redox_seq['SO4--'])/2),
 'Pyrite': (mineral_col, (redox_seq['Fe3+']+redox_seq['SO4--']+redox_seq['SO4--']-100)/3),
 'SO4--': (TEA_col, redox_seq['SO4--']),
 'HS-': (TEA_col, redox_seq['SO4--']-50),
 'H2S(aq)': (gas_col, redox_seq['SO4--']),
 'Sulfate reduction': (reduction_col, (redox_seq['SO4--'])-50),
 'Methane oxidation (SO4)': (oxidation_col, redox_seq['Fe3+']),

 'Fe(OH)3': (mineral_col, redox_seq['Fe3+']),
 'Fe+++': (TEA_col, redox_seq['Fe3+']),
 'Fe++': (TEA_col, redox_seq['Fe3+']-75),
 'Fe(II) abiotic oxidation': (oxidation_col, ((redox_seq['Fe3+']+redox_seq['O2'])+400)/2),
 'Fe(II) microbial oxidation':(oxidation_col, (redox_seq['Fe3+']+redox_seq['O2'])/2),
 'Fe(III) reduction': (reduction_col, redox_seq['Fe3+']-50),
 'Methane oxidation (Fe)': (oxidation_col, redox_seq['Fe3+']),

 'NH4+': (TEA_col, redox_seq['NO3-']+50),
 'NH4+': (TEA_col, redox_seq['NO3-']+50),
 'NO3-': (TEA_col, redox_seq['NO3-']),
 'Denitrification': (reduction_col, redox_seq['NO3-']-100),
 'Nitrification': (substrateox_col, redox_seq['NO3-']+75/2+55),
 'Fermentative DNRA': (substrateox_col, (redox_seq['NO3-']+redox_seq['DOM'])/2+190),
 'Autotrophic DNRA':(substrateox_col, (redox_seq['NO3-']+redox_seq['HS-'])/2+110),
 'N2(aq)': (gas_col, redox_seq['NO3-']-150),
 'N2O(aq)': (gas_col, redox_seq['NO3-']-75),
 'Methane oxidation (NO3)': (oxidation_col, (redox_seq['Fe3+'])),

 'SOM': (substrate_col, redox_seq['SOM']),
 'Hydrolysis': (reduction_col, (redox_seq['SOM'])+100),
 'DOM1': (substrate_col, redox_seq['DOM']),
 'fermentation': (reduction_col, (redox_seq['DOM']+redox_seq['Acetate'])/2+50),
 'Acetate-': (substrate_col, redox_seq['Acetate']),
 'H2(aq)': (substrate_col, redox_seq['H2']),

 'O2(aq)': (TEA_col, redox_seq['O2']),
 'HCO3-': (175, redox_seq['CO2']),

 'CH4(aq)': (TEA_col, redox_seq['CH4']),


 'DOM oxidation (O2)': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+300),
 'Aerobic decomposition': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+300),
 'Acetate oxidation (O2)': (substrateox_col, (redox_seq['O2']+redox_seq['Acetate'])/2+300),
 'Methane oxidation (O2)': (oxidation_col, (redox_seq['Fe3+'])),

 'Hydrogenotrophic methanogenesis': (reduction_col, (redox_seq['CO2'])-100),
 'Acetoclastic methanogenesis': (reduction_col, (redox_seq['CO2']+100-100))}

node_colors={
    'POM':'C0',
    'Microbe':'C1',
    'DOM':'C2',
    'MAOM':'C3',
    'Litter':'g',
    'CWD':'brown',
    'Mineral':'C4',
    'Primary aqueous':'C5',
    'Secondary aqueous':'C7',
    'Gas':'C9',
    'Metal ion':'orange',
    'Surface complex':'C8',
    'Unknown':'C5',
    'General Reaction':'#825708',
    'Microbial Reaction':'#e6990b',
    'SOMdecomp Reaction':'#c7a82c',
}

# result,units=plot_pf_output.convert_units(result,units,'mol/m^3')
if __name__ == '__main__':
    reaction_network=make_network()
    networkfig=pyplot.figure('Reaction network',clear=True,figsize=(11,13));networkfig.clf()
#    drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid'],
#            font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
#            namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'HCO$_3^-$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
#                         'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',})
    drawn,pos=decomp_network.draw_network_with_reactions(reaction_network,
            omit=['gas','secondary','H+','H2O','Na+','Cl-','Calcite','Ca++','HCO3-','mineral'],
            font_size='x-large',node_size=500,font_color='w',arrowstyle='-|>',arrowsize=20.0,width=1.5,edge_color='k',node_alpha=0.9,node_colors=node_colors,markers={'Reaction':'h','mineral':'8'},
            namechanges={'SOM':'SOM','DOM1':'DOM','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
                         'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
                         'H2(aq)':'H$_2$','DOM oxidation (O2)':'Aerobic\nresp','Acetate oxidation (O2)':'Aerobic\nresp','Aerobic decomposition':'Aerobic\nresp',
                         'NH4+':'NH$_4^+$','NO3-':'NO$_3^-$','N2O(aq)':'N$_2$O','N2(aq)':'N$_2$','Fermentative DNRA':'Fermentative\nDNRA','Autotrophic DNRA':'Autotrophic\nDNRA','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$','Fe(II) abiotic oxidation':'Fe(II)\nabiotic\noxidation',
                         'Fe(II) microbial oxidation':'Fe(II)\nmicrobial\noxidation','SO4--':'SO$_4^{--}$','HS-':'H$_2$S','Methane oxidation (O2)':'Methane\noxidation','Methane oxidation (NO3)':'Methane\noxidation',
                         'Methane oxidation (Fe)':'Methane\noxidation','Methane oxidation (SO4)':'Methane\noxidation','Fe(III) reduction':'Fe(III)\nreduction','Sulfate reduction':'Sulfate\nreduction',
                         'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Acetoclastic methanogenesis':'Acetoclastic\nmethanogenesis',
                         'SOMdecomp Reaction':'SOM Reaction','General Reaction':'Abiotic Reaction','fermentation':'Fermentation','Primary aqueous':'Dissolved ion','Gas':'Dissolved gas'},pos=pos,connectionstyle='arc3, rad=0.2')

    networkfig.axes[0].set_facecolor('#8a9ebf')

    for p in networkfig.axes[0].patches:
        ex=p.get_path().get_extents()
        if (ex.max[0]<155 or ex.max[0]>188) and ex.min[1]>300:
            p.set_connectionstyle('arc3,rad=-0.1')


#    pflotran_exe='../pflotran-elm-interface/src/pflotran/pflotran'
#    simlength=60
#    simlength=41*30 #Sun et al was 41 months

#    # Example results
#    results_all,units=decomp_network.PF_network_writer(make_network(leaf_Mn_mgkg=500.0,change_constraints={'Birnessite':'1.0d-5 1.d2 m^2/m^3','Manganite':'0.0973d-5 1.d2 m^2/m^3'})).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
#    resultsfig,axes=pyplot.subplots(4,1,num='Results',clear=True,figsize=(6.3,8),squeeze=False)
#    plot_result(results_all,SOM_ax=axes[0,0],pH_ax=axes[1,0],mineral_ax=axes[2,0],porewater_ax=axes[3,0],do_legend=True)


#    cm=pyplot.get_cmap('coolwarm')

#    import pandas

#    fig,axes=pyplot.subplots(3,1,num='Leaf Mn concentrations',clear=True)
#    axes[0].set_title('Lignin remaining')
#    axes[0].set_ylabel('Fraction of initial')
#    axes[1].set_title('Mn++ concentration')
#    axes[1].set_ylabel('Concentration (M)')
#    axes[2].set_title('Manganite concentration')
#    axes[2].set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
#    axes[2].set_xlabel('Time (days)')
#    # Range of leaf Mn concentrations
#    norm=matplotlib.colors.Normalize(1,5)
#    results_all_leafMn={}
#    for leafMn in arange(1.0,5.1) :
#        result,units=decomp_network.PF_network_writer(make_network(leaf_Mn_mgkg=10.0**leafMn)).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
#        axes[0].plot(result['Lignin']/result['Lignin'].iloc[0],label='Leaf Mn concentration = 10$^{%d}$ mg/kg'%int(leafMn),c=cm(norm(leafMn)))
#        axes[1].plot(result['Total Mn++'],c=cm(norm(leafMn)))
#        # axes[2].plot(result['Manganite VF']/24.45*1e6,c=cm(norm(leafMn)))
#        birnessite=axes[2].plot(result['Birnessite2 VF']/251.17*1e6,c=cm(norm(leafMn)),ls='--',label='Birnessite')[0]
#        results_all_leafMn[leafMn]=result

#    axes[0].legend()


#    nyears=40
#    results_long_leafMn=[]
#    nums=linspace(1,35,10)
#    nums=logspace(-1,1.5,10)

#    cellulose_eq=zeros(len(nums))
#    lignin_eq=zeros(len(nums))
#    for num in range(len(nums)):
#        result,units=decomp_network.PF_network_writer(make_network(leaf_Mn_mgkg=num,Mn2_scale=0.25e-2)).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=nyears*365)
#        results_long_leafMn.append(result)

#    norm=matplotlib.colors.LogNorm(nums.min(),nums.max())
#    fig,axes=pyplot.subplots(ncols=3,num='Multiple years',clear=True)
#    volume_factor=2/100*12e-3 # Converting from mol/m3 to kgC/m2 assuming 1 cm depth
#    for num in range(len(nums)):
#        cellulose=zeros(nyears*365)
#        lignin=zeros(nyears*365)
#        for yr in range(nyears):
#            cellulose[365*yr:]+=interp(arange(nyears*365)/365,results_long_leafMn[num].index/365,results_long_leafMn[num]['Cellulose'])[:365*(nyears-yr)]
#            lignin[365*yr:]+=interp(arange(nyears*365)/365,results_long_leafMn[num].index/365,results_long_leafMn[num]['Lignin'])[:365*(nyears-yr)]

#        axes[1].plot(arange(nyears*365)/365,lignin*volume_factor,ls='--',c=cm(norm(nums[num]))) # Assumes 1 cm thick organic horizon
#        cellulose_eq[num]=cellulose[-365:].mean()*volume_factor
#        lignin_eq[num]=lignin[-365:].mean()*volume_factor

#        axes[0].plot(results_long_leafMn[num].index/365,results_long_leafMn[num]['Lignin']*volume_factor,ls='--',c=cm(norm(nums[num])),label='Lignin: %1.1f mg/g Mn'%nums[num])
#        axes[2].plot(nums[num],cellulose_eq[num]+lignin_eq[num],'o',c=cm(norm(nums[num])),ms=8.0)

#    axes[2].plot(nums,cellulose_eq+lignin_eq,'k-')
#    axes[0].plot(results_long_leafMn[num].index/365,results_long_leafMn[-1]['Cellulose']*volume_factor,'k-',label='Cellulose')
#    axes[1].plot(arange(nyears*365)/365,cellulose*volume_factor,'k-')
#    axes[0].legend(loc='upper right')

#    axes[0].set(title='One leaf litter cohort',xlabel='Time (years)',ylabel='C stock (kg m$^{-2}$)')
#    axes[1].set(title='Cumulative litter layer',xlabel='Time (years)',ylabel='C stock (kg m$^{-2}$)')
#    axes[2].set(title='Litter layer vs Mn concentration',xlabel='Leaf Mn concentration (mg g$^{-1}$)',ylabel='C stock (kg m$^{-2}$)')



#    fig,axes=pyplot.subplots(3,1,num='Mn leakage amounts',clear=True)
#    axes[0].set_title('Lignin remaining')
#    axes[0].set_ylabel('Fraction of initial')
#    axes[1].set_title('Mn++ concentration')
#    axes[1].set_ylabel('Concentration (M)')
#    axes[2].set_title('Manganite concentration')
#    axes[2].set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
#    axes[2].set_xlabel('Time (days)')
#    # Range of leaf Mn concentrations
#    norm=matplotlib.colors.Normalize(-8,-3)
#    for Mnleakage in arange(-8.0,-3.0) :
#        result,units=decomp_network.PF_network_writer(make_network(Mn_peroxidase_Mn3_leakage=10.0**Mnleakage,leaf_Mn_mgkg=10.0)).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
#        axes[0].plot(result['Lignin']/result['Lignin'].iloc[0],label='Mn+++ leakage = 10$^{%d}$'%int(Mnleakage),c=cm(norm(Mnleakage)))
#        axes[1].plot(result['Total Mn++'],c=cm(norm(Mnleakage)))
#        axes[2].plot(result['Manganite VF']/24.45*1e6,c=cm(norm(Mnleakage)))
#        results_all=pandas.concat((results_all,result))

#    axes[0].legend()


#    fig,axes=pyplot.subplots(3,1,num='Birnessite rates',clear=True)
#    axes[0].set_title('Lignin remaining')
#    axes[0].set_ylabel('Fraction of initial')
#    axes[1].set_title('Mn++ concentration')
#    axes[1].set_ylabel('Concentration (M)')
#    axes[2].set_title('Manganite concentration')
#    axes[2].set_ylabel('Concentration\n($\mu$mol/cm$^{-3}$)')
#    axes[2].set_xlabel('Time (days)')
#    # Range of leaf Mn concentrations
#    norm=matplotlib.colors.Normalize(-16,-10)
#    for rate in arange(-16.0,-9+0.1,2) :
#        network=make_network(leaf_Mn_mgkg=1.0,Mn_peroxidase_Mn3_leakage=1e-3,change_rate={'Birnessite':'1.d%d mol/m^2-sec'%rate},change_constraints={'Birnessite':'1.0d-5 1.d2 m^2/m^3','Manganite':'0.0973d-5 1.d2 m^2/m^3'})
#        result,units=decomp_network.PF_network_writer(network).run_simulation('SOMdecomp_template.txt','manganese',pflotran_exe,print_output=False,length_days=simlength)
#        axes[0].plot(result['Lignin']/result['Lignin'].iloc[0],label='Dissolution rate = 10$^{%d}$ mol m$^{-2}$ s$^{-1}$'%int(rate),c=cm(norm(rate)))
#        axes[1].plot(result['Total Mn++'],c=cm(norm(rate)))
#        manganite=axes[2].plot(result['Manganite VF']/24.45*1e6,c=cm(norm(rate)),label='Manganite')[0]
#        # birnessite=axes[2].plot(result['Birnessite VF']/251.17*1e6,c=cm(norm(rate)),ls='--',label='Birnessite')[0]
#        results_all=pandas.concat((results_all,result))

#    axes[0].legend()
#    axes[2].legend(handles=[manganite])

#    # fig,axes=pyplot.subplots(1,1,num='Mn2 effect',clear=True)
#    # dlignin=results_all['Lignin'].diff()/(results_all.index[1]-results_all.index[0])
#    # dlignin=dlignin.mask(dlignin>0)
#    # axes.plot(results_all['Free Mn++'],-dlignin,'o')
#    # axes.set_title('Lignin degradation rate vs Mn++')
#    # axes.set_ylabel('Lignin degradataion rate (day$^{-1}$)')
#    # axes.set_xlabel('Mn++ concentration (M)')
#    # axes.set_xscale('linear')
    pyplot.savefig('reactionnetworkv.pdf')
    pyplot.show()
