import decomp_network
import plot_pf_output

def make_network(
        tinyval = 1.0e-30,
        natomw  = 14.0067,
        catomw  = 12.0110,
        N_imm_lim=1e-2,
        thresh=0.0,
        adfactor_soil4=100.0,
        adfactor_soil3=10.0,
                                ):
        # CTC decomposition network
        decomp_network_CTC=decomp_network.decomp_network()

        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL1',CN=  12.0 ,constraints={'initial':tinyval/catomw},kind='immobile') )
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL2',CN=  12.0 ,constraints={'initial':tinyval/catomw},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL3',CN=  10.0 ,constraints={'initial':tinyval/catomw},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL4',CN=  10.0 ,constraints={'initial':tinyval/catomw},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR1',constraints={'initial':1e3/catomw},initCN=20*natomw/catomw ,kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR2',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile')    )
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR3',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='CWD',constraints={'initial':tinyval/catomw},initCN=20*natomw/catomw ,kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='CO2(aq)',kind='primary'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='NH4+',constraints={'initial':1e-10},kind='primary'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='NO3-',constraints={'initial':1e-10},kind='primary'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='HRimm',constraints={'initial':tinyval},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Nimm',constraints={'initial':tinyval},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Nimp',constraints={'initial':tinyval},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Nmin',constraints={'initial':tinyval},kind='immobile'))

        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Plant_NH4_demand',constraints={'initial':tinyval},kind='immobile'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Plant_NO3_demand',constraints={'initial':tinyval},kind='immobile'))

        # These are for plant N uptake. Defining as aqueous pools so they can be represented as Microbial reaction
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Tracer',constraints={'initial':tinyval},kind='primary'))
        decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='Tracer2',constraints={'initial':tinyval},kind='primary'))

        # CWD decomposition to  litter
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'CWD':1.0},product_pools={'LITR2':0.76,'LITR3':0.24},
                                rate_constant=0.00010,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',reactiontype='SOMDECOMP',
                                name='CWD fragmentation'))

        # Litter decomposition
        # Monod dependence on NO3 and NH4 allows N limitation of immobilization to work without crashing ELM
        # Although this would limit decomposition if NH4 is high but NO3 is low?

        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR1':1.0},product_pools={'SOIL1':0.61,'CO2(aq)':1-0.61},
                        rate_constant=1.204,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR1 decomposition',reactiontype='SOMDECOMP',
                        monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim,threshold=thresh),decomp_network.monod(species='NO3-',k=N_imm_lim,threshold=thresh)]))
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR2':1.0},product_pools={'SOIL2':0.45,'CO2(aq)':1-0.45},
                        rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR2 decomposition',reactiontype='SOMDECOMP',
                        monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim,threshold=thresh),decomp_network.monod(species='NO3-',k=N_imm_lim,threshold=thresh)]))
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR3':1.0},product_pools={'SOIL3':0.71,'CO2(aq)':1-0.71},
                        rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR3 decomposition',reactiontype='SOMDECOMP',
                        monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim,threshold=thresh),decomp_network.monod(species='NO3-',k=N_imm_lim,threshold=thresh)]))

        # SOM decomposition
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL1':1.0},product_pools={'SOIL2':0.72,'CO2(aq)':1-0.72},
                rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL1 decomp',reactiontype='SOMDECOMP'))
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL2':1.0},product_pools={'SOIL3':0.54,'CO2(aq)':1-0.54},
                rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL2 decomp',reactiontype='SOMDECOMP'))
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL3':1.0},product_pools={'SOIL4':0.45,'CO2(aq)':1-0.45},
                rate_constant=0.00141*adfactor_soil3,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL3 decomp',reactiontype='SOMDECOMP'))
        decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL4':1.0},product_pools={'CO2(aq)':1.0},
                rate_constant=0.0001*adfactor_soil4,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL4 decomp',reactiontype='SOMDECOMP'))

        # Plant N uptake. Using monod reactions instead of sandbox for simplicity/flexibility
        # But this requires plant N uptake to be defined as an aqueous pool
        # Rate constant will be modified by plant N demand for each layer in ELM
        # Could add inhibition of NH4 uptake by NO3 if we want preferential uptake
        # Problem: If plant N demand sets rate constant, then we need to stop the sum of the two uptakes from exceeding plant N demand
        # Current approach is modifying rate constant in ELM/EMI for each reaction by relative amount of NO3 and NH4 so combined max rate is plant demand
        # Alternate way would be using inhibition here, or leaking excess N back out. Might allow more flexibility in how N uptake is defined
        decomp_network_CTC.add_reaction(decomp_network.reaction(stoich='1.0 NH4+ -> 1.0 Tracer2',reactiontype='MICROBIAL',
                name='Plant NH4 uptake',monod_terms=[decomp_network.monod(species='NH4+',k=N_imm_lim,threshold=thresh)],
                rate_constant=1.0,biomass='Plant_NH4_demand',threshold=thresh))
        decomp_network_CTC.add_reaction(decomp_network.reaction(stoich='1.0 NO3- -> 1.0 Tracer',reactiontype='MICROBIAL',
                name='Plant NO3 uptake',monod_terms=[decomp_network.monod(species='NO3-',k=N_imm_lim,threshold=thresh)],
                rate_constant=1.0,biomass='Plant_NO3_demand',threshold=thresh))

        # Nitrification. Simple version for CN reaction network without oxygen, H+, etc
        decomp_network_CTC.add_reaction(decomp_network.reaction(stoich='1.0 NH4+ -> 1.0 NO3-',reactiontype='MICROBIAL',
                name='Nitrification',monod_terms=[decomp_network.monod(species='NH4+',k=1e-3)],threshold=thresh,rate_constant=1e-9))

        return decomp_network_CTC

# Not yet implemented: Add temperature response function to sandbox as option in network writer and make sure it works with PFLOTRAN code
'''
      Abiotic_Factors
        TEMPERATURE_RESPONSE_FUNCTION # Might need to implement changes in PFLOTRAN code?
          Q10 1.5d0
        /
      /
'''

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