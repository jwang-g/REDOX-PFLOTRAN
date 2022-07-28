import decomp_network

def make_network(
        tinyval = 1.0e-30,
        natomw  = 14.0067,
        catomw  = 12.0110,
        N_imm_lim=0.1, # Model crashes on too many time step cuts if this is too low
        N_uptake_lim=1e-4, # Plants don't take up enough N if this is too high
        thresh=0.0,
        adfactor_soil4=100.0,
        adfactor_soil3=10.0,
        oxygen_consuming=False,
                                ):
        # CTC decomposition network

        pools = [
                decomp_network.decomp_pool(name='SOIL1',CN=  12.0 ,constraints={'initial':tinyval/catomw},kind='immobile') ,
                decomp_network.decomp_pool(name='SOIL2',CN=  12.0 ,constraints={'initial':tinyval/catomw},kind='immobile'),
                decomp_network.decomp_pool(name='SOIL3',CN=  10.0 ,constraints={'initial':tinyval/catomw},kind='immobile'),
                decomp_network.decomp_pool(name='SOIL4',CN=  10.0 ,constraints={'initial':tinyval/catomw},kind='immobile'),
                decomp_network.decomp_pool(name='LITR1',constraints={'initial':1e3/catomw},initCN=20*natomw/catomw ,kind='immobile'),
                decomp_network.decomp_pool(name='LITR2',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile')  ,  
                decomp_network.decomp_pool(name='LITR3',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile'),
                decomp_network.decomp_pool(name='CWD',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile'),
                decomp_network.decomp_pool(name='CO2(aq)',kind='primary',constraints={'initial':'400e-6 G CO2(g)*'}),
                decomp_network.decomp_pool(name='NH4+',constraints={'initial':1e-10},kind='primary'),
                decomp_network.decomp_pool(name='NO3-',constraints={'initial':1e-2},kind='primary'),
                decomp_network.decomp_pool(name='HRimm',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Nimm',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Nimp',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Nmin',constraints={'initial':tinyval},kind='immobile'),

                decomp_network.decomp_pool(name='Plant_NH4_demand',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Plant_NO3_demand',constraints={'initial':tinyval},kind='immobile'),

                # These are for plant N uptake. Defining as aqueous pools so they can be represented as Microbial reaction
                decomp_network.decomp_pool(name='Tracer',constraints={'initial':tinyval},kind='primary'),
                decomp_network.decomp_pool(name='Tracer2',constraints={'initial':tinyval},kind='primary'),

                decomp_network.decomp_pool(name='CO2(g)*',kind='gas'),
        ]

        if oxygen_consuming:
                pools.append(decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'0.2 G O2(g)'}))
                pools.append(decomp_network.decomp_pool(name='O2(g)',kind='gas'))
        

        reactions = decomp_network.pools_list_to_dict([
                # CWD decomposition to  litter
                decomp_network.reaction(reactant_pools={'CWD':1.0},product_pools={'LITR2':0.76,'LITR3':0.24},
                                        rate_constant=0.00010,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',reactiontype='SOMDECOMP',
                                        name='CWD fragmentation'),

                # Litter decomposition
                # Monod dependence on NO3 and NH4 allows N limitation of immobilization to work without crashing ELM
                decomp_network.reaction(reactant_pools={'LITR1':1.0},product_pools={'SOIL1':0.61,'CO2(aq)':1-0.61},
                                rate_constant=1.204,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR1 decomp',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*1.204/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*1.204/24,threshold=thresh,pool_normalized=True) ]),
                decomp_network.reaction(reactant_pools={'LITR2':1.0},product_pools={'SOIL2':0.45,'CO2(aq)':1-0.45},
                                rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR2 decomp',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*0.0726/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*0.0726/24,threshold=thresh,pool_normalized=True) ]),
                decomp_network.reaction(reactant_pools={'LITR3':1.0},product_pools={'SOIL3':0.71,'CO2(aq)':1-0.71},
                                rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR3 decomp',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*0.0141/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*0.0141/24,threshold=thresh,pool_normalized=True) ]),

                # SOM decomposition
                decomp_network.reaction(reactant_pools={'SOIL1':1.0},product_pools={'SOIL2':0.72,'CO2(aq)':1-0.72},
                        rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL1 decomp',reactiontype='SOMDECOMP'),
                decomp_network.reaction(reactant_pools={'SOIL2':1.0},product_pools={'SOIL3':0.54,'CO2(aq)':1-0.54},
                        rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL2 decomp',reactiontype='SOMDECOMP'),
                decomp_network.reaction(reactant_pools={'SOIL3':1.0},product_pools={'SOIL4':0.45,'CO2(aq)':1-0.45},
                        rate_constant=0.00141*adfactor_soil3,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL3 decomp',reactiontype='SOMDECOMP'),
                decomp_network.reaction(reactant_pools={'SOIL4':1.0},product_pools={'CO2(aq)':1.0},
                        rate_constant=0.0001*adfactor_soil4,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL4 decomp',reactiontype='SOMDECOMP'),

                # Plant N uptake. Using monod reactions instead of sandbox for simplicity/flexibility
                # But this requires plant N uptake to be defined as an aqueous pool
                # Rate constant will be modified by plant N demand for each layer in ELM
                # Could add inhibition of NH4 uptake by NO3 if we want preferential uptake
                # Problem: If plant N demand sets rate constant, then we need to stop the sum of the two uptakes from exceeding plant N demand
                # Current approach is modifying rate constant in ELM/EMI for each reaction by relative amount of NO3 and NH4 so combined max rate is plant demand
                # Alternate way would be using inhibition here, or leaking excess N back out. Might allow more flexibility in how N uptake is defined
                decomp_network.reaction(stoich='1.0 NH4+ -> 1.0 Tracer2',reactiontype='MICROBIAL',
                        name='Plant NH4 uptake',monod_terms=[decomp_network.monod(species='NH4+',k=N_uptake_lim,threshold=thresh)],
                        rate_constant=1.0,biomass='Plant_NH4_demand',threshold=thresh),
                decomp_network.reaction(stoich='1.0 NO3- -> 1.0 Tracer',reactiontype='MICROBIAL',
                        name='Plant NO3 uptake',monod_terms=[decomp_network.monod(species='NO3-',k=N_uptake_lim,threshold=thresh)],
                        rate_constant=1.0,biomass='Plant_NO3_demand',threshold=thresh),

                # Nitrification. Simple version for CN reaction network without oxygen, H+, etc
                decomp_network.reaction(stoich='1.0 NH4+ -> 1.0 NO3-',reactiontype='MICROBIAL',
                        name='Nitrification',monod_terms=[decomp_network.monod(species='NH4+',k=1e-3)],threshold=thresh,rate_constant=1e-9),
        ])

        if oxygen_consuming:
                # Runs slow when oxygen limited with normal monod as litter and SOM pools get large. 
                # But pool-normalized monod slows down decomposition too much for large pools because oxygen transfers might have to happen faster than 1 hour
                O2lim=[decomp_network.monod(species='O2(aq)',k=1e-4)]
                xO2=0.025 # Runs through year 11 when xO2 is 0.1, but crashes at 0.05. And litter decomposition is very slow at 0.1
                reactions['LITR1 decomp']['monod_terms'].extend([decomp_network.monod(species='O2(aq)',k=xO2*1.204/24,threshold=thresh,pool_normalized=True)]*3)
                reactions['LITR2 decomp']['monod_terms'].extend([decomp_network.monod(species='O2(aq)',k=xO2*0.0726/24,threshold=thresh,pool_normalized=True)]*3)
                reactions['LITR3 decomp']['monod_terms'].extend([decomp_network.monod(species='O2(aq)',k=xO2*0.0141/24,threshold=thresh,pool_normalized=True)]*3)
                reactions['LITR1 decomp']['monod_terms'].extend(O2lim)
                reactions['LITR2 decomp']['monod_terms'].extend(O2lim)
                reactions['LITR3 decomp']['monod_terms'].extend(O2lim)
                # Monod terms were empty for soil pools before
                reactions['SOIL1 decomp']['monod_terms']=[decomp_network.monod(species='O2(aq)',k=xO2*0.0726/24,threshold=thresh,pool_normalized=True)]*3
                reactions['SOIL2 decomp']['monod_terms']=[decomp_network.monod(species='O2(aq)',k=xO2*0.0141/24,threshold=thresh,pool_normalized=True)]*3
                reactions['SOIL3 decomp']['monod_terms']=[decomp_network.monod(species='O2(aq)',k=xO2*0.00141*adfactor_soil3/24,threshold=thresh,pool_normalized=True)]*3
                reactions['SOIL4 decomp']['monod_terms']=[decomp_network.monod(species='O2(aq)',k=xO2*0.0001*adfactor_soil4/24,threshold=thresh,pool_normalized=True)]*3
                reactions['SOIL1 decomp']['monod_terms'].extend(O2lim)
                reactions['SOIL2 decomp']['monod_terms'].extend(O2lim)
                reactions['SOIL3 decomp']['monod_terms'].extend(O2lim)
                reactions['SOIL4 decomp']['monod_terms'].extend(O2lim)
                
        return decomp_network.decomp_network(pools,decomp_network.pools_dict_to_list(reactions))


def make_aqueous_network(
        tinyval = 1.0e-30,
        natomw  = 14.0067,
        catomw  = 12.0110,
        N_imm_lim=0.1, # Model crashes on too many time step cuts if this is too low
        N_uptake_lim=1e-4, # Plants don't take up enough N if this is too high
        thresh=0.0,
        adfactor_soil4=100.0,
        adfactor_soil3=10.0,
        DOM_scale=1e-1,O2_scale=1e-4,
        Fe=True,
        calcite=True,
                                ):
        # CTC decomposition network

        pools = [
                decomp_network.decomp_pool(name='SOIL1',CN=  12.0 ,constraints={'initial':tinyval/catomw},kind='immobile') ,
                decomp_network.decomp_pool(name='SOIL2',CN=  12.0 ,constraints={'initial':tinyval/catomw},kind='immobile'),
                decomp_network.decomp_pool(name='SOIL3',CN=  10.0 ,constraints={'initial':tinyval/catomw},kind='immobile'),
                decomp_network.decomp_pool(name='SOIL4',CN=  10.0 ,constraints={'initial':tinyval/catomw},kind='immobile'),
                decomp_network.decomp_pool(name='LITR1',constraints={'initial':1e3/catomw},initCN=20*natomw/catomw ,kind='immobile'),
                decomp_network.decomp_pool(name='LITR2',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile')  ,  
                decomp_network.decomp_pool(name='LITR3',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile'),
                decomp_network.decomp_pool(name='CWD',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile'),
                decomp_network.decomp_pool(name='CO2(aq)',kind='primary',constraints={'initial':'400e-6 G CO2(g)*'}),
                decomp_network.decomp_pool(name='NH4+',constraints={'initial':1e-10},kind='primary'),
                decomp_network.decomp_pool(name='NO3-',constraints={'initial':1e-10},kind='primary'),
                decomp_network.decomp_pool(name='HRimm',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Nimm',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Nimp',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Nmin',constraints={'initial':tinyval},kind='immobile'),

                decomp_network.decomp_pool(name='Plant_NH4_demand',constraints={'initial':tinyval},kind='immobile'),
                decomp_network.decomp_pool(name='Plant_NO3_demand',constraints={'initial':tinyval},kind='immobile'),

                # These are for plant N uptake. Defining as aqueous pools so they can be represented as Microbial reaction
                decomp_network.decomp_pool(name='Tracer',constraints={'initial':tinyval},kind='primary'),
                decomp_network.decomp_pool(name='Tracer2',constraints={'initial':tinyval},kind='primary'),

                decomp_network.decomp_pool(name='CO2(g)*',kind='gas'),
                decomp_network.decomp_pool(name='HCO3-',kind='secondary'),
                decomp_network.decomp_pool(name='CO3--',kind='secondary'),

                decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'0.2 G O2(g)'}),
                decomp_network.decomp_pool(name='O2(g)',kind='gas'),

                # DOM pools corresponding to the litter and soil pools. Litter DOM is flexible CN and there are two fixed CN ones for soil pools
                decomp_network.decomp_pool(name='DOM1',kind='primary',constraints={'initial':tinyval/catomw},CN=100.0),
                # decomp_network.decomp_pool(name='DOM2',kind='primary',constraints={'initial':tinyval/catomw},CN=10.0),
                # This accounts for stoichiometry of SOIL2 decomposing to SOIL3 which has a higher C:N
                # decomp_network.decomp_pool(name='DOM3',kind='primary',constraints={'initial':tinyval/catomw},CN=16.0),

                decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'6.0 P'}),

        ]
        # Pools for iron reduction
        if Fe:
                pools.extend([
                        decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-5 M Fe(OH)3'}),
                        decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':1e-6}),
                        decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
                        decomp_network.decomp_pool(name='FeIIDOM1(aq)',kind='secondary'),
                        decomp_network.decomp_pool(name='OH-',kind='secondary'),
                        # Specifying SSA on mass basis messes up repeated constraint equilibration for some reason
                        decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-10 mol/m^2-sec',constraints={'initial':'0.01  80.0 m^2/m^3'},kind='mineral'),
                        decomp_network.decomp_pool(name='H2O',kind='implicit'),        
                ])
        if calcite:
                pools.extend([
                        decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':'1.16e-3'}),
                        decomp_network.decomp_pool(name='CaCO3(aq)',kind='secondary'),
                        decomp_network.decomp_pool(name='CaHCO3+',kind='secondary'),
                        decomp_network.decomp_pool(name='CaOH+',kind='secondary'),
                        # I don't think PFLOTRAN will update SSA when volume changes when run via alquimia, might need to do that on driver side
                        decomp_network.decomp_pool(name='Calcite',rate='1.d-6 mol/m^2-sec',constraints={'initial':'0.d-6  80.0'},kind='mineral'),
                ])

        reactions = [
                # CWD decomposition to  litter
                decomp_network.reaction(reactant_pools={'CWD':1.0},product_pools={'LITR2':0.76,'LITR3':0.24},
                                        rate_constant=0.00010,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',reactiontype='SOMDECOMP',
                                        name='CWD fragmentation'),

                # Litter decomposition
                # Monod dependence on NO3 and NH4 allows N limitation of immobilization to work without crashing ELM
                # NOTE: Changing this to DOM with N content might mess up N immobilization
                decomp_network.reaction(reactant_pools={'LITR1':1.0},product_pools={'SOIL1':0.61,'DOM1':1-0.61},
                                rate_constant=1.204,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR1 decomp',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*1.204/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*1.204/24,threshold=thresh,pool_normalized=True) ],
                                inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),
                decomp_network.reaction(reactant_pools={'LITR2':1.0},product_pools={'SOIL2':0.45,'DOM1':1-0.45},
                                rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR2 decomp',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*0.0726/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*0.0726/24,threshold=thresh,pool_normalized=True) ],
                                inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),
                decomp_network.reaction(reactant_pools={'LITR3':1.0},product_pools={'SOIL3':0.71,'DOM1':1-0.71},
                                rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR3 decomp',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*0.0141/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*0.0141/24,threshold=thresh,pool_normalized=True) ],
                                inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),

                # SOM decomposition
                decomp_network.reaction(reactant_pools={'SOIL1':1.0},product_pools={'SOIL2':0.72,'DOM1':1-0.72},
                        rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL1 decomp',reactiontype='SOMDECOMP',
                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),
                decomp_network.reaction(reactant_pools={'SOIL2':1.0},product_pools={'SOIL3':0.54,'DOM1':1-0.54},
                        rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL2 decomp',reactiontype='SOMDECOMP',
                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),
                decomp_network.reaction(reactant_pools={'SOIL3':1.0},product_pools={'SOIL4':0.45,'DOM1':1-0.45},
                        rate_constant=0.00141*adfactor_soil3,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL3 decomp',reactiontype='SOMDECOMP',
                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),
                decomp_network.reaction(reactant_pools={'SOIL4':1.0},product_pools={'DOM1':1.0},
                        rate_constant=0.0001*adfactor_soil4,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL4 decomp',reactiontype='SOMDECOMP',
                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=DOM_scale)]),

                # decomp_network.reaction(name='DOM2 respiration',stoich = '10 DOM2 + 10 O2(aq) -> 10 CO2(aq) + NH4+',
                #         rate_constant=5e-7/10,reactiontype='MICROBIAL',
                #         monod_terms=[decomp_network.monod(species='DOM2',k=DOM_scale),decomp_network.monod(species='O2(aq)',k=O2_scale)]),
                # decomp_network.reaction(name='DOM3 respiration',stoich = '16.0 DOM2 + 16.0 O2(aq) -> 16.0 CO2(aq) + NH4+',
                #         rate_constant=5e-7/16.0,reactiontype='MICROBIAL',
                #         monod_terms=[decomp_network.monod(species='DOM2',k=DOM_scale),decomp_network.monod(species='O2(aq)',k=O2_scale)]),

                # decomp_network.reaction(name='DOM1 respiration',stoich = 'DOM1C + O2(aq) -> CO2(aq)',
                #         rate_constant=5e-7,reactiontype='MICROBIAL',
                #         monod_terms=[decomp_network.monod(species='DOM1C',k=DOM_scale),decomp_network.monod(species='O2(aq)',k=O2_scale)]),
                # decomp_network.reaction(name='DOM1 N mineralization',stoich = 'DOM1N -> NH4+',
                #         rate_constant=5e-7,reactiontype='MICROBIAL',
                #         monod_terms=[decomp_network.monod(species='DOM1N',k=DOM_scale),decomp_network.monod(species='O2(aq)',k=O2_scale)]),

                # Plant N uptake. Using monod reactions instead of sandbox for simplicity/flexibility
                # But this requires plant N uptake to be defined as an aqueous pool
                # Rate constant will be modified by plant N demand for each layer in ELM
                # Could add inhibition of NH4 uptake by NO3 if we want preferential uptake
                # Problem: If plant N demand sets rate constant, then we need to stop the sum of the two uptakes from exceeding plant N demand
                # Current approach is modifying rate constant in ELM/EMI for each reaction by relative amount of NO3 and NH4 so combined max rate is plant demand
                # Alternate way would be using inhibition here, or leaking excess N back out. Might allow more flexibility in how N uptake is defined
                decomp_network.reaction(stoich='1.0 NH4+ -> 1.0 Tracer2',reactiontype='MICROBIAL',
                        name='Plant NH4 uptake',monod_terms=[decomp_network.monod(species='NH4+',k=N_uptake_lim,threshold=thresh)],
                        rate_constant=1.0,biomass='Plant_NH4_demand',threshold=thresh),
                decomp_network.reaction(stoich='1.0 NO3- -> 1.0 Tracer',reactiontype='MICROBIAL',
                        name='Plant NO3 uptake',monod_terms=[decomp_network.monod(species='NO3-',k=N_uptake_lim,threshold=thresh)],
                        rate_constant=1.0,biomass='Plant_NO3_demand',threshold=thresh),

                # Nitrification. Simple version for CN reaction network without oxygen, H+, etc
                decomp_network.reaction(stoich='1.0 NH4+ -> 1.0 NO3-',reactiontype='MICROBIAL',
                        name='Nitrification',monod_terms=[decomp_network.monod(species='NH4+',k=1e-3)],threshold=thresh,rate_constant=1e-9),
        ]
        if Fe:
                reactions.extend([
                        # Iron reduction reactions
                        # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
                        # 8 Fe+++ + 8 e- -> 8 Fe++ 
                        # pH might not be balanced right with how we're handling DOM decomposition above, and buffering is probably messed up
                        decomp_network.reaction(name='Fe(III) reduction',reactant_pools={'DOM1':1.0,'Fe+++':4.0,'H2O':2.0},
                                                                        product_pools={'CO2(aq)':1.0,'Fe++':4.0,'H+':4.0,'NH4+':1.0/100*catomw/natomw},
                                                # stoich='2.0 DOM1 + 8.0 Fe+++ + 4.0 H2O -> 2.0 CO2(aq) + 8.0 Fe++ + 9.0 H + NH4+',
                                                monod_terms=[decomp_network.monod(species='DOM1',k=DOM_scale,threshold=1.1e-15),decomp_network.monod(species='Fe+++',k=1e-9)],
                                                inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=O2_scale,type='MONOD')],
                                                rate_constant=3e-8,reactiontype='MICROBIAL'),
                        decomp_network.reaction(name='Fe(II) microbial oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ -> 1.0 Fe+++ + 0.5 H2O',
                                                monod_terms=[decomp_network.monod(species='O2(aq)',k=O2_scale,threshold=0.0),decomp_network.monod(species='Fe++',k=0.1),
                                                                decomp_network.monod(species='H+',k=1e-5)],
                                                rate_constant=100e-8,reactiontype='MICROBIAL'),
                        # Assume as DOM is decomposed part per thousand Fe per unit C is released
                        decomp_network.reaction(name='DOM1 respiration',reactant_pools={'DOM1':1.0,'O2(aq)':1.0},product_pools={'CO2(aq)':1.0,'NH4+':1.0/100*catomw/natomw,'Fe+++':1e-3},
                                                rate_constant=5e-7,reactiontype='MICROBIAL',
                                                monod_terms=[decomp_network.monod(species='DOM1',k=0.1),decomp_network.monod(species='O2(aq)',k=O2_scale,threshold=0.0)])
                ])
        else:
                # DOM respiration. Assume stoichiometry of C6H12O6 + N0.5 + 6 O2 -> 6 CO2 + 6 H2O + 0.5 NH4+ (H and charge not balanced)
                # NOTE: Sandbox CN is mass units, so needs to be converted for molar units here (or use another approach)  
                reactions.extend([decomp_network.reaction(name='DOM1 respiration',reactant_pools={'DOM1':1.0,'O2(aq)':1.0},product_pools={'CO2(aq)':1.0,'NH4+':1.0/100*catomw/natomw},
                        rate_constant=5e-7,reactiontype='MICROBIAL',
                        monod_terms=[decomp_network.monod(species='DOM1',k=0.1),decomp_network.monod(species='O2(aq)',k=O2_scale,threshold=0.0)]),
                ])


                
        return decomp_network.decomp_network(pools,reactions)

def oxygen_demand(txt):
        out=0.0
        for line in txt.split('\n'):
                if line.strip().startswith('LITR1'):
                        out = out + float(line.split()[-1])*1.204/24*(1-.61)
                elif line.strip().startswith('LITR2'):
                        out = out + float(line.split()[-1])*0.0726/24*(1-.45)
                elif line.strip().startswith('LITR3'):
                        out = out + float(line.split()[-1])*0.0141/24*(1-.71)
                elif line.strip().startswith('SOIL1'):
                        out = out + float(line.split()[-1])*0.0726/24*(1-.72)
                elif line.strip().startswith('SOIL2'):
                        out = out + float(line.split()[-1])*0.0141/24*(1-.54)
                elif line.strip().startswith('SOIL3'):
                        out = out + float(line.split()[-1])*0.00141/24*(1-.45)*10
                elif line.strip().startswith('SOIL4'):
                        out = out + float(line.split()[-1])*0.0001/24*100
        return out
                

decomp_network_ad=make_network()
decomp_network_notad=make_network(adfactor_soil3=1.0,adfactor_soil4=1.0)
decomp_network_O2_ad=make_aqueous_network(calcite=False,Fe=True)
decomp_network_O2_notad=make_aqueous_network(adfactor_soil3=1.0,adfactor_soil4=1.0,calcite=False,Fe=True)


def load_state_from_logfile(filename):
        lines=[]
        with open(filename) as f:
                for line in f:
                        # if line.strip() == 'Alquimia primary species (mol/m3 bulk):':
                        if line.strip() == 'Alquimia aux doubles:':
                                lines=f.readlines()
                                break
                else:
                        raise ValueError('No alquimia state printout found')
        return lines


def str_to_alquimia(inputfile,state_text,temperature=20.0):
        import run_alquimia
        import numpy as np
        ffi=run_alquimia.ffi

        (chem,data,sizes,status)=run_alquimia.init_alquimia(inputfile,hands_off=True)
        pools=run_alquimia.get_alquimiavector(data.meta_data.primary_names)
        minerals=run_alquimia.get_alquimiavector(data.meta_data.mineral_names)

        data.state.temperature=temperature
        data.properties.volume=1.0
        data.state.water_density=1000.0
        data.state.aqueous_pressure=101325.0
        data.state.porosity=0.5


        init_cond=run_alquimia.convert_condition_to_alquimia(None,'initial')
        chem.ProcessCondition(ffi.new('void **',data.engine_state),init_cond,ffi.addressof(data.properties),ffi.addressof(data.state),ffi.addressof(data.aux_data),status)
        run_alquimia.check_status(status,False)

        data.state.temperature=temperature
  
        if isinstance(state_text,str):
                state_text=state_text.split('\n')
        for line in state_text:
                if line.startswith('Porosity'):
                        data.state.porosity=float(line.split('=')[-1])
                elif line.startswith('Saturation'):
                        data.properties.saturation=float(line.split('=')[-1])
                elif line.startswith('Temperature'):
                        data.state.temperature=float(line.split('=')[-1])
                elif line.startswith('Pressure'):
                        data.state.aqueous_pressure=float(line.split('=')[-1])
                elif line.startswith('Water density'):
                        data.state.water_density=float(line.split('=')[-1])
                elif line.startswith('Volume'):
                        data.properties.volume=float(line.split('=')[-1])
                elif line.startswith('Aux double'):
                        l=line.split()
                        if len(l)>2:
                                data.aux_data.aux_doubles.data[int(l[2])-1]=float(l[3])
                else:
                        l=line.split()
                        if len(l)>0 and l[1] in pools:
                                poolnum=pools.index(l[1])
                                data.state.total_mobile.data[poolnum]=float(l[3])
                                data.state.total_immobile.data[poolnum]=float(l[2])
                        elif len(l)>0 and l[1] in minerals:
                                poolnum=minerals.index(l[1])
                                data.state.mineral_volume_fraction.data[poolnum]=float(l[2])
                                data.state.mineral_specific_surface_area.data[poolnum]=float(l[3])

        return chem,data,status


def test_solve(chem,data,status,dt=3600):
        import run_alquimia
        import numpy as np
        
        # run_alquimia.print_alquimia_object(data.state)
        mobile_before=run_alquimia.get_alquimiavector(data.state.total_mobile)
        immobile_before=run_alquimia.get_alquimiavector(data.state.total_immobile)
        aux_before=run_alquimia.get_alquimiavector(data.aux_data.aux_doubles)
        mineral_before=run_alquimia.get_alquimiavector(data.state.mineral_volume_fraction)
        mineral_SSA_before=run_alquimia.get_alquimiavector(data.state.mineral_specific_surface_area)
        free_before=np.array(run_alquimia.get_alquimiavector(data.aux_output.primary_free_ion_concentration))*data.state.water_density/1000
        pools=run_alquimia.get_alquimiavector(data.meta_data.primary_names)
        mineralnames=run_alquimia.get_alquimiavector(data.meta_data.mineral_names)
        try:
                max_cuts=run_alquimia.run_onestep(chem,data,dt,status,hands_off=True)
                print(f'Converged at {max_cuts:d} cuts, {status.num_newton_iterations} Newton iterations')
        except RuntimeError as err:
                print(f'****** Solve failed: {status.num_newton_iterations} Newton iterations ******')
                print(err)
                # run_alquimia.print_alquimia_object(status)
        mobile_after=run_alquimia.get_alquimiavector(data.state.total_mobile)
        immobile_after=run_alquimia.get_alquimiavector(data.state.total_immobile)
        aux_after=run_alquimia.get_alquimiavector(data.aux_data.aux_doubles)
        mineral_after=run_alquimia.get_alquimiavector(data.state.mineral_volume_fraction)
        mineral_SSA_after=run_alquimia.get_alquimiavector(data.state.mineral_specific_surface_area)
        free_after=np.array(run_alquimia.get_alquimiavector(data.aux_output.primary_free_ion_concentration))*data.state.water_density/1000
        import pandas
        df=pandas.DataFrame({'Mobile_before':mobile_before,'Mobile_after':mobile_after,
                           'Immobile_before':immobile_before,'Immobile_after':immobile_after,
                           'Free after':free_after},index=pools)
        df['Mobile_change']=df['Mobile_after']-df['Mobile_before']
        df['Immobile_change']=df['Immobile_after']-df['Immobile_before']
        df2=pandas.DataFrame({'aux_before':aux_before,'aux_after':aux_after})
        df2['aux_change']=df2['aux_after']-df2['aux_before']
        minerals=pandas.DataFrame({'Mineral VF before':mineral_before,'Mineral VF after':mineral_after,
                                   'Mineral SSA before':mineral_SSA_before,'Mineral SSA after':mineral_SSA_after,
                                },index=mineralnames)
        minerals['Mineral VF change']=minerals['Mineral VF after']-minerals['Mineral VF before']
        minerals['Mineral SSA change']=minerals['Mineral SSA after']-minerals['Mineral SSA before']
        return df,df2,minerals


if __name__=='__main__':
        # Write out input deck
        decomp_network.PF_network_writer(decomp_network_ad,precision=8).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM_adspinup.in',
                CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
                chem_args={'MAX_RESIDUAL_TOLERANCE':.5e-13,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-10}
                )

        decomp_network.PF_network_writer(decomp_network_notad,precision=8).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM.in',
                CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
                chem_args={'MAX_RESIDUAL_TOLERANCE':.5e-13,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-10}
                )


        # Write out input deck
        decomp_network.PF_network_writer(decomp_network_O2_ad,precision=8).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM_O2consuming_adspinup.in',
                CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
                chem_args={'MAX_RESIDUAL_TOLERANCE':.5e-13,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-10}
                )

        decomp_network.PF_network_writer(decomp_network_O2_notad,precision=8).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM_O2consuming.in',
                CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
                chem_args={'MAX_RESIDUAL_TOLERANCE':.5e-13,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-10}
                )

        import sys
        if 'check_cons' in sys.argv:
                from plot_pf_output import read_tecfile
                import matplotlib.pyplot as plt
                result,units=decomp_network.PF_network_writer(decomp_network_ad,precision=4).run_simulation('SOMdecomp_template.txt','CTC_alquimia_forELM_adspinup','../pflotran-interface/src/pflotran/pflotran',
                        print_output=False,length_days=20,log_formulation=True,CO2name='CO2(aq)',SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
                        chem_args={'MAX_RESIDUAL_TOLERANCE':1e-15,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-15})
                d,u=read_tecfile('CTC_alquimia_forELM_adspinup_generated-mas.dat')

                result2,units=decomp_network.PF_network_writer(decomp_network_O2_ad,precision=4).run_simulation('SOMdecomp_template.txt','CTC_alquimia_forELM_O2consuming_adspinup','../pflotran-interface/src/pflotran/pflotran',
                        print_output=False,length_days=20,log_formulation=True,CO2name='CO2(aq)',SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
                        chem_args={'MAX_RESIDUAL_TOLERANCE':1e-15,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-15})
                d2,u2=read_tecfile('CTC_alquimia_forELM_O2consuming_adspinup_generated-mas.dat')

                totalC1=(d['Global LITR1C']+d['Global LITR2C']+d['Global LITR3C']+d['Global SOIL1']+d['Global SOIL2']+d['Global SOIL3']+d['Global SOIL4']+d['Global CWDC']+d['Global CO2(aq)'])
                totalC1_2=(result['LITR1C']+result['LITR2C']+result['LITR3C']+result['SOIL1']+result['SOIL2']+result['SOIL3']+result['SOIL4']+result['CWDC']+result['HRimm'])
                totalC2=(d2['Global LITR1C']+d2['Global LITR2C']+d2['Global LITR3C']+d2['Global SOIL1']+d2['Global SOIL2']+d2['Global SOIL3']+d2['Global SOIL4']+d2['Global CWDC']+d2['Global CO2(aq)']+d2['Global DOM1'])
                totalC2_2=(result2['LITR1C']+result2['LITR2C']+result2['LITR3C']+result2['SOIL1']+result2['SOIL2']+result2['SOIL3']+result2['SOIL4']+result2['CWDC']+result2['Total CO2(aq)']*.25*1000+result2['Total DOM1']*0.25*1000)

                totalN1=(d['Global LITR1N']+d['Global LITR2N']+d['Global LITR3N']+d['Global CWDN']+(d['Global SOIL1']/12+d['Global SOIL2']/12+d['Global SOIL3']/10+d['Global SOIL4']/10)*12.0110/14.0067+d['Global NO3-']+d['Global NH4+']+d['Global Tracer']+d['Global Tracer2'])
                totalN2=(d2['Global LITR1N']+d2['Global LITR2N']+d2['Global LITR3N']+d2['Global CWDN']+(d2['Global SOIL1']/12+d2['Global SOIL2']/12+d2['Global SOIL3']/10+d2['Global SOIL4']/10)*12.0110/14.0067+d2['Global NO3-']+d2['Global NH4+']+d2['Global Tracer']+d2['Global Tracer2']+
                                (d2['Global DOM1']/20)*12.0110/14.0067)

                totalN1_2=(result['LITR1N']+result['LITR2N']+result['LITR3N']+result['CWDN']+(result['SOIL1']/12+result['SOIL2']/12+result['SOIL3']/10+result['SOIL4']/10)*12.0110/14.0067+(result['Total NO3-']+result['Total NH4+']+result['Total Tracer']+result['Total Tracer2'])*.25*1000)
                totalN2_2=(result2['LITR1N']+result2['LITR2N']+result2['LITR3N']+result2['CWDN']+(result2['SOIL1']/12+result2['SOIL2']/12+result2['SOIL3']/10+result2['SOIL4']/10)*12.0110/14.0067+(result2['Total NO3-']+result2['Total NH4+']+result2['Total Tracer']+result2['Total Tracer2']+
                                (result2['Total DOM1']/20)*12.0110/14.0067)*.25*1000)

                f,a=plt.subplots(num=1,clear=True,nrows=2)
                (totalC1-totalC1.iloc[0]).plot(ax=a[0],c='C0')
                (totalC1_2-totalC1_2.iloc[0]).plot(ax=a[0],c='C0',ls='--')
                (totalC2-totalC2.iloc[0]).plot(ax=a[0],c='C1')
                (totalC2_2-totalC2_2.iloc[0]).plot(ax=a[0],c='C1',ls='--')

                (totalN1-totalN1.iloc[0]).plot(ax=a[1],c='C0')
                (totalN1_2-totalN1_2.iloc[0]).plot(ax=a[1],c='C0',ls='--')
                (totalN2-totalN2.iloc[0]).plot(ax=a[1],c='C1')
                (totalN2_2-totalN2_2.iloc[0]).plot(ax=a[1],c='C1',ls='--')

        if 'show_networks' in sys.argv:
                import matplotlib.pyplot as plt
                f,a=plt.subplots(ncols=3,clear=True,num='Reaction networks',figsize=(15,6))

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
                opts=dict(
                omit=['gas','secondary','H2O','Na+','Cl-','Calcite','Ca++','Tracer','Tracer2',
                        'HRimm','Nimm','Nimp','Nmin','Plant_NH4_demand','Plant_NO3_demand'],

                font_size='medium',node_size=800,font_color='w',arrowstyle='-|>',arrowsize=10.0,width=1.0,edge_color='k',node_alpha=0.9,node_colors=node_colors,markers={'Reaction':'h','mineral':'8'},
                namechanges={'SOM':'SOM','DOM1':'DOM','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
                     'H2(aq)':'H$_2$','DOM oxidation (O2)':'Aerobic\nresp','Acetate oxidation (O2)':'Aerobic\nresp','Aerobic decomposition':'Aerobic\nresp',
                     'NH4+':'NH$_4^+$','NO3-':'NO$_3^-$','N2O(aq)':'N$_2$O','N2(aq)':'N$_2$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$',
                     'SO4--':'SO$_4^{--}$','HS-':'H$_2$S','Methane oxidation (O2)':'Methane\noxidation','Methane oxidation (NO3)':'Methane\noxidation',
                     'Methane oxidation (Fe)':'Methane\noxidation','Methane oxidation (SO4)':'Methane\noxidation','Fe(III) reduction':'Fe(III)\nreduction','Sulfate reduction':'Sulfate\nreduction',
                     'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Acetoclastic methanogenesis':'Acetoclastic\nmethanogenesis',
                     'SOMdecomp Reaction':'SOM Reaction','General Reaction':'Abiotic Reaction','fermentation':'Fermentation','Primary aqueous':'Dissolved ion','Gas':'Dissolved gas',
                     'CO2(aq)':'CO$_2$','CWD fragmentation':'CWD\nfragmentation','Plant NH4 uptake':'Plant NH$_4$\nuptake','Plant NO3 uptake':'Plant NO$_3$ uptake',
                     'DOM1 respiration':'DOM\nrespiration','Fe(II) microbial oxidation':'Fe(II)\noxidation'},connectionstyle='arc3, rad=0.2')
                
                decomp_network.draw_network_with_reactions(decomp_network_notad,ax=a[0],do_legend=False,**opts)
                decomp_network.draw_network_with_reactions(make_aqueous_network(Fe=False),do_legend=False,ax=a[1],**opts)
                decomp_network.draw_network_with_reactions(decomp_network_O2_notad,ax=a[2],do_legend=False,**opts)

                a[0].set(facecolor='#8a9ebf',title='ELM default')
                a[1].set(facecolor='#8a9ebf',title='ELM pools with oxygen')
                a[2].set(facecolor='#8a9ebf',title='ELM pools with oxygen and Fe reduction')

                plt.show()

        if 'test_solve' in sys.argv:
                if len(sys.argv)>sys.argv.index('test_solve')+1:
                        state=load_state_from_logfile(sys.argv[sys.argv.index('test_solve')+1])
                else:
                        state=''
                mobile,aux,miner=test_solve(*str_to_alquimia('CTC_alquimia_forELM_O2consuming_adspinup.in',state),dt=3600*0.5**12)
                print(aux)
                print(mobile)
                print(miner)