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
                                rate_constant=1.204,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR1 decomposition',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*1.204/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*1.204/24,threshold=thresh,pool_normalized=True) ]),
                decomp_network.reaction(reactant_pools={'LITR2':1.0},product_pools={'SOIL2':0.45,'CO2(aq)':1-0.45},
                                rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR2 decomposition',reactiontype='SOMDECOMP',
                                monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim*0.0726/24,threshold=thresh,pool_normalized=True) ,
                                        decomp_network.monod(species='NO3-',k=N_imm_lim*0.0726/24,threshold=thresh,pool_normalized=True) ]),
                decomp_network.reaction(reactant_pools={'LITR3':1.0},product_pools={'SOIL3':0.71,'CO2(aq)':1-0.71},
                                rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR3 decomposition',reactiontype='SOMDECOMP',
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
                reactions['LITR1 decomposition']['monod_terms'].extend([decomp_network.monod(species='O2(aq)',k=xO2*1.204/24,threshold=thresh,pool_normalized=True)]*3)
                reactions['LITR2 decomposition']['monod_terms'].extend([decomp_network.monod(species='O2(aq)',k=xO2*0.0726/24,threshold=thresh,pool_normalized=True)]*3)
                reactions['LITR3 decomposition']['monod_terms'].extend([decomp_network.monod(species='O2(aq)',k=xO2*0.0141/24,threshold=thresh,pool_normalized=True)]*3)
                reactions['LITR1 decomposition']['monod_terms'].extend(O2lim)
                reactions['LITR2 decomposition']['monod_terms'].extend(O2lim)
                reactions['LITR3 decomposition']['monod_terms'].extend(O2lim)
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
# Write out input deck
decomp_network.PF_network_writer(decomp_network_ad,precision=4).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM_adspinup.in',
        CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
        chem_args={'MAX_RESIDUAL_TOLERANCE':1e-15,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-15}
        )

decomp_network.PF_network_writer(decomp_network_notad,precision=4).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM.in',
        CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
        chem_args={'MAX_RESIDUAL_TOLERANCE':1e-15,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-15}
        )

decomp_network_O2_ad=make_network(oxygen_consuming=True)
decomp_network_O2_notad=make_network(oxygen_consuming=True,adfactor_soil3=1.0,adfactor_soil4=1.0)
# Write out input deck
decomp_network.PF_network_writer(decomp_network_O2_ad,precision=4).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM_O2consuming_adspinup.in',
        CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
        chem_args={'MAX_RESIDUAL_TOLERANCE':1e-15,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-15}
        )

decomp_network.PF_network_writer(decomp_network_O2_notad,precision=4).write_into_input_deck('SOMdecomp_template.txt','CTC_alquimia_forELM_O2consuming.in',
        CO2name='CO2(aq)',log_formulation=True,SOMdecomp_Q10=1.5,moisturefunc='LOGTHETA',
        chem_args={'MAX_RESIDUAL_TOLERANCE':1e-15,'MAX_RELATIVE_CHANGE_TOLERANCE':1e-15}
        )