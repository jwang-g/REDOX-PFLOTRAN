import run_alquimia
import decomp_network
from matplotlib import pyplot
import numpy
import copy

pools = [
decomp_network.decomp_pool(name='cellulose',CN=50,constraints={'initial':8e3},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':0.5e-1},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'5.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':'1.0 G O2(g)'}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-10 M Fe(OH)3'}),
decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-3'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='H2(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Mg++',kind='primary',constraints={'initial':0.5e-3}),
decomp_network.decomp_pool(name='Ca++',kind='primary',constraints={'initial':0.5e-3}),
decomp_network.decomp_pool(name='Na+',kind='primary',constraints={'initial':2e-3}),
decomp_network.decomp_pool(name='K+',kind='primary',constraints={'initial':2e-5}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3+',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)4-',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeCH3COO+',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='Fe(OH)2+',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3(aq)',kind='secondary'),

# See Roberts Earth Sci Rev article for discussion of iron minerals, including some time scales
decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-6 mol/m^2-sec',constraints={'initial':'0.875d-3  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Goethite',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d1 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Fe',rate='1.d-7 mol/m^2-sec',constraints={'initial':'1.0e-6  1. m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

# Maybe this should be a cation exchange reaction instead?
# Probably makes sense to include both
# There should be some Fe(II) sorption or complexation with OM
decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=2000.0*5.0,complexes=['>Carboxylic_acid']), # site density in mol/m3

decomp_network.decomp_pool(name='H2O',kind='implicit'),
]
             
# Herndon et al 2015, BGC: Fe(III) average concentration 0.37 mmol/L (same as mM). SO4- was 0.07, NO3- was 0.03.
# Fe(III) was 60% of dissolved Fe
# ** What form was it in? Does it just dissolve or was it released/complexed by microbial activity?
# DOC 5-15 mmol/L
# Fe probably present as iron oxide precipitates or Fe-OM complexes, not primary minerals
#
# Herndon et al 2015, JGR-B: Organic acids 12% of WEOC in organic layer. formate 48.6 umol/ g SOC; acetate 66, propionate 0.76
# Each mol of acetate produces 1 mol of methane or 8 mol of Fe(II)
# Ferrihydrite specific surface area ~250 m2/g (Kaiser and Guggenberger 2008 https://doi.org/10.1046/j.1365-2389.2003.00544.x)
# 

# Beth: There should be a lot of solid OM functional groups, probably carboxylate/carboxylic acid that exchange protons and buffer pH in these soils
# Fe(II) is probably sorbing onto OM, or staying dissolved
# Fe(III) should precipitate in this system, so if it's dissolved it would be complexed with something

conc_scales={
    'DOM1':5e-1,
    'Acetate-':1e-1,
    'Fe+++':1e-10,
    'O2(aq)':1e-4,
    'Fe++':1e-1,
}
anox_inhib_conc = 1e-5

truncate_conc=1e-15

reactions = [
    decomp_network.reaction(name='Aerobic decomposition',reactant_pools={'cellulose':1.0},product_pools={'HCO3-':1.0},reactiontype='SOMDECOMP',
                                            rate_constant=1e-1,rate_units='y', 
                                            # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='THRESHOLD 1.0d20')]),
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=1.1e-9)]),
    
    decomp_network.reaction(name='Hydrolysis',stoich='1.0 cellulose -> 1.0 DOM1',reactiontype='SOMDECOMP',
                                            rate_constant=1e-1,rate_units='y', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=conc_scales['DOM1']),
                                                          # decomp_network.inhibition(species='O2(aq)',k=6.25e-11,type='THRESHOLD -1.0d10')
                                                          decomp_network.inhibition(species='O2(aq)',type='MONOD',k=1e-11),
                                                          ]),
    
    # Calculating these as per unit carbon, so dividing by the 6 carbons in a glucose
    # C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2
    # Should it be inhibited by H2?
    decomp_network.reaction(name='fermentation',reactant_pools={'DOM1':6/6},product_pools={'Acetate-':2/6,'HCO3-':2/6,'H+':4/6,'H2(aq)':4*2/6,'Tracer':2/6}, # balancing pH of FeIII release requires an extra 5.5 H+ to be released here
                                            rate_constant=1e-10,reactiontype='MICROBIAL', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD'),decomp_network.inhibition(species='Acetate-',k=conc_scales['Acetate-'],type='MONOD')],
                                        monod_terms=[decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=1.1e-15)]),

    # CH2O + H2O -> CO2 + 4H+ + 4 e-
    # O2   + 4H+ + 4 e- -> 2H2O
    decomp_network.reaction(name='DOM aerobic respiration',stoich='1.0 DOM1 + 1.0 O2(aq) -> 1.0 HCO3- + 1.0 H+ + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=1.1e-10),decomp_network.monod(species='DOM1',k=conc_scales['DOM1'],threshold=1.1e-16)],
                                        rate_constant=1.0e-9,reactiontype='MICROBIAL'),

    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 2 O2    + 8 H+ + 8 e- -> 4 H2O
    decomp_network.reaction(name='Acetate aerobic respiration',stoich='1.0 Acetate-  + 2.0 O2(aq)  -> 2.0 HCO3-  + 2.0 H+  + 2.0 Tracer ',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=1.1e-10),decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-16)],
                                        rate_constant=1.0e-9,reactiontype='MICROBIAL'),


    # C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
    # 8 Fe+++ + 8 e- -> 8 Fe++ 
    decomp_network.reaction(name='Fe(III) reduction',stoich='1.0 Acetate- + 8.0 Fe+++ -> 2.0 HCO3- + 8.0 Fe++ + 9.0 H+ + 2.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-15),decomp_network.monod(species='Fe+++',k=conc_scales['Fe+++'],threshold=1.1e-15)],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD')],
                                            rate_constant=2e-10,reactiontype='MICROBIAL'),
                                            
    # Oxidation of Fe++
    # Currently assuming backward reaction rate is zero
    decomp_network.reaction(name='Fe(II) abiotic oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O',
                                            rate_constant=1e-2,backward_rate_constant=0.0,reactiontype='GENERAL'),
                                            
    decomp_network.reaction(name='Fe(II) microbial oxidation',stoich='1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ -> 1.0 Fe+++ + 0.5 H2O',
                                            monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=1.1e-10),decomp_network.monod(species='Fe++',k=conc_scales['Fe++'],threshold=1.1e-15)],
                                            rate_constant=1e-8,reactiontype='MICROBIAL'),

    # Acetoclastic methanogenesis
    # C2H3O2- + H+ -> CH4 + HCO3- + H+
    decomp_network.reaction(name='Methanogenesis',stoich='1.0 Acetate- -> 1.0 CH4(aq) + 1.0 HCO3- + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-15),
                                                         ],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD'),decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                            rate_constant=1e-10,reactiontype='MICROBIAL'),

    # Hydrogenotrophic methanogenesis
    
    # H2 oxidation if oxygen available
    
    # Cation exchange
    decomp_network.ion_exchange(name='Cation exchange',cations={'Fe++':0.1,'Fe+++':0.3,'Mg++':1.1,'Ca++':4.1,'Na+':1.0,'K+':0.9,'H+':1.1},CEC=2000.0,mineral=None)

]

rate_scale=1e-8
rateconstants={
       '1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  <-> 1.00e+00 Fe+++  + 5.00e-01 H2O'                      : 1.0e0*1.0e1*0,
       '1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  -> 1.00e+00 Fe+++  + 5.00e-01 H2O'                      : rate_scale*100,
       '1.00e+00 DOM1  -> 3.33e-01 Acetate-  + 3.33e-01 HCO3-  + 6.67e-01 H+  + 1.33e+00 H2(aq)  + 3.33e-01 Tracer'  : rate_scale*120,
       "1.00e+00 DOM1  + 1.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 Tracer"                    : rate_scale*3,
       "1.00e+00 Acetate-  + 2.00e+00 O2(aq)  -> 2.00e+00 HCO3-  + 2.00e+00 H+  + 2.00e+00 Tracer"                : rate_scale*3,
       "1.00e+00 Acetate-  + 8.00e+00 Fe+++  -> 2.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 2.00e+00 Tracer" : rate_scale*2,
       "1.00e+00 Acetate-  -> 1.00e+00 CH4(aq)  + 1.00e+00 HCO3-  + 1.00e+00 Tracer"                             : rate_scale*2,
       "cellulose decay to CO2 (SOMDEC sandbox)"                                                              : rate_scale*15.85, #1.0/(365*24*3600)*1.0
       "cellulose decay to DOM1 (SOMDEC sandbox)"                                                              : rate_scale*15.85 #1.0/(365*24*3600)*1.0
}

reaction_network =  decomp_network.decomp_network(pools=pools,reactions=reactions)


# Generate PFLOTRAN input file with correct reactions
decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt','Arctic_redox.in',log_formulation=True,CO2name='Tracer',truncate_concentration=1e-25)
            
dq=0.01 # Diffusion coefficient when aerobic
O2_const=numpy.zeros(365*24)+dq

simlength=150
initfrac=0.0
timestep=3600


O2_initial=numpy.zeros(simlength*24)
O2_initial[:int(simlength*24*initfrac)]=dq

# Incubation data
import pandas,glob
datadir='/Users/b0u/Documents/NGEE-Arctic/Redox_sims/Data'
Barrow_synthesis_data=pandas.concat([pandas.read_csv(fname,na_values=-9999,header=8,skiprows=[9]) for fname in glob.glob(datadir+'/Barrow_soil_geochem_synthesis/*.csv')])  
units=pandas.read_csv(glob.glob(datadir+'/Barrow_soil_geochem_synthesis/*.csv')[0],na_values=-9999,header=8).iloc[0]
SOC_layermean=Barrow_synthesis_data.groupby(Barrow_synthesis_data['Soil_layer'].str.capitalize())['SOC'].mean()
SOC_layer_microtopo_mean=Barrow_synthesis_data.groupby([Barrow_synthesis_data['Soil_layer'].str.capitalize(),Barrow_synthesis_data['Microtopography'].str.capitalize()=='Trough'])['SOC'].mean()


# SOC units (%) need to be converted to mol C/m3
# Bulk density reported for some cores. Best bet maybe to find relationship between SOC% and BD and assume it's constant accross these incubations
# Although, might be worth checking if they have some bulk density that I haven't found
corephysdata=pandas.read_csv(datadir+'/Barrow_soilcores2/core_physical_chem_data_20180321.csv',header=5,na_values=-9999,skiprows=[6])
corephysunits=pandas.read_csv(datadir+'/Barrow_soilcores2/core_physical_chem_data_20180321.csv',header=5,na_values=-9999).iloc[0]
from scipy.stats import linregress
BD_SOM_fit=linregress(corephysdata[['Organic_Matter_Content','Dry_Bulk_Density']].mask(corephysdata['Organic_Matter_Content']<6.0).dropna().to_numpy()) 
BD_layerest=BD_SOM_fit.intercept + BD_SOM_fit.slope*SOC_layermean/0.6 # Guessing at conversion from SOC to SOM weight
# Bockheim et al 2003, looks like a broader dataset:
#  They use a fit C[%] = -9.7872 ln(BD[g/cm3])+8.2432
# BD_layerest2=numpy.exp((8.2432-SOC_layermean)/9.7872)
BD_layerest2_trough=numpy.exp((8.2432-SOC_layer_microtopo_mean)/9.7872)

waterchemistry=pandas.read_csv(datadir+'/Barrow_porewater_chem/BGC_BarrowWaterChemistry_2013_2014_v1.csv',header=8,skiprows=[9,10],na_values=-9999)
waterchemistry_units=pandas.read_csv(datadir+'/Barrow_porewater_chem/BGC_BarrowWaterChemistry_2013_2014_v1.csv',header=8,na_values=-9999).iloc[0]

cellulosefrac=0.05
molar_volume_FeOH3=34.3600 #cm3/mol
porosity=0.5

def get_layer(Core_ID,layer,minT=4.0):
    if isinstance(Core_ID,int):
        Core_ID='NGADG%04d'%Core_ID
    data_layer=Barrow_synthesis_data[(Barrow_synthesis_data['Soil_layer'].str.capitalize()==layer.capitalize())&(Barrow_synthesis_data['Core_ID']==Core_ID)&(Barrow_synthesis_data['Incubation_Temperature']>minT)]
    if len(data_layer)==0:
        raise ValueError('No data found for %s layer in core %s'%(layer,Core_ID))
    return data_layer


# CEC numbers estimated from http://www.soilquality.org.au/factsheets/cation-exchange-capacity
def make_initcond(Core_ID,layer,cellulosefrac=0.05,porosity=None,minT=4.0,FeII_FeVF_factor=1.5,CEC_OM=200,CEC_mineral=25,oxic=False,Fe_SSA=1e2,carboxylate_conc=0.1,BD_factor=1.0,otherconstraints={}):
    if isinstance(Core_ID,int):
        Core_ID='NGADG%04d'%Core_ID
    data_layer=Barrow_synthesis_data[(Barrow_synthesis_data['Soil_layer'].str.capitalize()==layer.capitalize())&(Barrow_synthesis_data['Core_ID']==Core_ID)&(Barrow_synthesis_data['Incubation_Temperature']>minT)]
    if len(data_layer)==0:
        raise ValueError('No data found for %s layer in core %s'%(layer,Core_ID))
    
    SOC_layer=data_layer['SOC'].mean()/100 # In fraction, not percent
    
    # Based on Bockheim et al 2003 fit. Units are g/cm3
    BD_layer=numpy.exp((8.2432-SOC_layer*100)/9.7872)*BD_factor
    
    if porosity is None:
        if oxic:
            raise ValueError('Cannot calculate porosity for oxic incubation')
        porosity=1.0-data_layer['Moisture'].max()*BD_layer
    
    # Set initial iron concentration based on maximum dissolved Fe(II), multiplied by a correction factor to account for Fe(II) sorption and/or incomplete Fe(III) reduction by end of incubation
    Fe_VF=(data_layer['Fe_II'].max()-data_layer['Fe_II'].min())*FeII_FeVF_factor*1e-6*BD_layer*molar_volume_FeOH3
    FeII_t0=data_layer['Fe_II'].dropna().iloc[0]*1e-6*BD_layer/porosity*1e3
    
    # CEC affects sorption of cations as well as buffering capacity. Calculate as weighted average of CEC for organic and mineral components
    # Units of input CEC are meq/100g which is standard in soil science measurements. Convert to eq/m3 for PFLOTRAN
    CEC=(CEC_OM*SOC_layer + CEC_mineral*(1-SOC_layer))*1e-2*(BD_layer*1e3)
    
    pH_t0=data_layer['pH'].dropna().iloc[0]
    
    newconstraints=otherconstraints.copy()
    newconstraints['cellulose']=SOC_layer/12*BD_layer*100**3*cellulosefrac
    newconstraints['Fe(OH)3']='%1.4e  %1.2e m^2/m^3'%(Fe_VF,Fe_SSA)
    newconstraints['Fe++']=FeII_t0
    newconstraints['H+']='%1.2f P'%pH_t0
    
    # Cations. Using mean for now, could change to be depth-aware
    newconstraints['Mg++']=waterchemistry['Magnesium'].mean()*1e-3
    newconstraints['Ca++']=waterchemistry['Calcium'].mean()*1e-3
    newconstraints['Na+']=waterchemistry['Sodium'].mean()*1e-3
    newconstraints['K+']=waterchemistry['Potassium'].mean()*1e-3
                      
    if oxic:
        newconstraints['O2(aq)']='1.0 G O2(g)'
    
    new_initcond=decomp_network.change_constraints(pools,newconstraints)
    # if oxic:
    #     new_initcond=decomp_network.change_site_density(new_initcond, '>Carboxylate-', SOC_layer*100*0.5e1*BD_layer)
    # else:
    new_initcond=decomp_network.change_site_density(new_initcond, '>Carboxylate-', SOC_layer*carboxylate_conc*BD_layer*100**3)
    
    return new_initcond,BD_layer,SOC_layer,CEC,porosity
    

cores_trough=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')]['Core_ID'].unique()
cores_rim=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='rim')]['Core_ID'].unique()
cores_center=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='center')]['Core_ID'].unique()
cores_oxic=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)]['Core_ID'].unique()

good_cores_trough_anoxic=9
good_cores_rim_anoxic=5
good_cores_center_anoxic=[17,3]
good_cores_oxic=3 # Note, Core 3 is oxic in organic horizon and anoxic in mineral. No crossover.



pools_atmoO2_organic,BD_atmoO2_organic,SOC_atmoO2_organic,CEC_atmoO2_organic,porosity_atmoO2_organic=make_initcond(3,'Organic',cellulosefrac=.05,oxic=True,porosity=porosity)
pools_organic_trough,BD_organic_trough,SOC_organic_trough,CEC_organic_trough,porosity_organic_trough=make_initcond(9,'Organic')   
pools_organic_nottrough,BD_organic_nottrough,SOC_organic_nottrough,CEC_organic_nottrough,porosity_organic_nottrough=make_initcond(5,'Organic',BD_factor=2.0)   
# result_lowFe_organic,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_lowFe_organic,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc)

# No oxic mineral horizon
# pools_atmoO2_mineral,BD_atmoO2_mineral,SOC_atmoO2_mineral,CEC_atmoO2_mineral=make_initcond(3,'Mineral',cellulosefrac=.05*.04,oxic=True)
# result_highO2_mineral,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_atmoO2_mineral,bc=pools_atmoO2_mineral,diffquo={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_atmoO2_mineral)
pools_mineral_trough,BD_mineral_trough,SOC_mineral_trough,CEC_mineral_trough,porosity_mineral_trough=make_initcond(9,'Mineral')
pools_mineral_nottrough,BD_mineral_nottrough,SOC_mineral_nottrough,CEC_mineral_nottrough,porosity_mineral_nottrough=make_initcond(5,'Mineral',BD_factor=2.0)   

result_highO2_organic,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_atmoO2_organic,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_atmoO2_organic,porosity=porosity_atmoO2_organic)
result_organic_trough,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_organic_trough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_organic_trough,porosity=porosity_organic_trough) 
result_organic_nottrough,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_organic_nottrough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_organic_nottrough,porosity=porosity_organic_nottrough) 
result_mineral_trough,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_mineral_trough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_mineral_trough,porosity=porosity_mineral_trough)    
result_mineral_nottrough,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_mineral_nottrough,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc,CEC=CEC_mineral_nottrough,porosity=porosity_mineral_nottrough) 

oxicfrac=0.1
nperiodic=3

O2_periodic=numpy.zeros(simlength*24)
O2_periodic[:int(simlength*24*oxicfrac)]=dq
result_periodicO2,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength*nperiodic,timestep,initcond=pools_organic_trough,bc=pools_atmoO2_organic,
    diffquo={'O2(aq)':O2_periodic},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_organic_trough)


nperiodic=1
pH_obs=float(decomp_network.pools_list_to_dict(pools_organic_trough)['H+']['constraints']['initial'].split()[0])  
periodic_out={}
for ndryperiods in range(1,5):
    periodic_out[ndryperiods]={}
    for pH in pH_obs+numpy.arange(-.5,1.1,0.5):

        O2_periodic=numpy.zeros(int(simlength/ndryperiods)*24)
        O2_periodic[:int(simlength/ndryperiods*24*oxicfrac)]=dq
        result_periodicO2,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength*nperiodic,timestep,initcond=decomp_network.change_constraint(pools_organic_trough,'H+','%1.2f P'%pH),bc=pools_atmoO2_organic,
            diffquo={'O2(aq)':O2_periodic},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc,CEC=CEC_organic_trough)
        periodic_out[ndryperiods][pH]=result_periodicO2




import plot_pf_output
from pylab import *

networkfig2=pyplot.figure('Reaction network (with reactions)',clear=True)
drawn=decomp_network.draw_network_with_reactions(reaction_network,node_size=700,arrowstyle='->',arrowsize=8.5,#edge_color='gray',
            omit=['NH4+','Rock(s)','gas','surf_complex','secondary','H+','implicit','H2(aq)'],
            namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'CO2(aq)',
                         'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate'})


def plot_result(result,SOM_ax=None,pH_ax=None,Fe_ax=None,CO2flux_ax=None,CH4flux_ax=None,porewater_ax=None,cation_ax=None,
                SOM_args={},pH_args={},Fe_args={},CO2flux_args={},CH4flux_args={},porewater_args={},cation_args={},
                do_legend=False,gdrywt=False,BD=None,SOC_pct=None,cellulose_SOC_frac=1.0,obs=None,obs_marker=None):
    if gdrywt:
        if SOC_pct is None and BD is None:
            raise TypeError('SOC_pct or BD must be a number if gdrywt is True')
        SOC_mol_m3=result['Total Sorbed cellulose'].iloc[0]/cellulose_SOC_frac # mol SOC/m3
        SOC_gC_cm3=SOC_mol_m3*12/100**3
        SOC_gC_gdwt=SOC_pct/100
        cm3_to_dwt=SOC_gC_gdwt/SOC_gC_cm3 # Conversion from /cm3 to /gdwt. This is 1/ bulk density in g/cm3
        if BD is not None:
            cm3_to_dwt=1.0/BD
        print('cm3_to_dwt = %1.3f'%cm3_to_dwt)
    if SOM_ax is not None:
        if gdrywt:
            # Plot in %, i.e. gC/gdrywt *100.
            l=SOM_ax.plot(result['Total Sorbed cellulose']/SOC_mol_m3*SOC_pct+SOC_pct*(1-cellulose_SOC_frac),label='SOM',**SOM_args)[0]
            SOM_ax.set_ylabel('Concentration\n(SOC %)')
        else:
            l=SOM_ax.plot(result['Total Sorbed cellulose']*1e-3,label='SOM',**SOM_args)[0]
            SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_xlabel('Time (days)')

    if pH_ax is not None:
        pH_ax.plot(-numpy.log10(result['Free H+']))

        pH_ax.set_title('pH')
        pH_ax.set_ylabel('pH')
        pH_ax.set_xlabel('Time (days)')
        
        if obs is not None:
            if obs_marker is not None:
                markers=['o','s','^','+','<','x']
                obs_cats=unique(obs[obs_marker])
                print(obs_cats)
                for x in range(len(obs_cats)):
                    xx=obs[obs_marker]==obs_cats[x]
                    pH_ax.plot(obs.Incubation_Time[xx],obs['pH'][xx],marker=markers[x],linestyle='None',label='Measured',color='C0')
            else:
                pH_ax.plot(obs.Incubation_Time,obs['pH'],marker='o',linestyle='None',label='Measured',color='C0',**pH_args)
        
    if Fe_ax is not None:
        molar_volume=34.3600 # From database. cm3/mol
        molar_weight = 106.8690
        # Add sorbed Fe to this
        if gdrywt:
            # Assume we can use SOC % to convert from volume to dry weight
            l=Fe_ax.plot(result['Fe(OH)3 VF']/molar_volume*1e6*cm3_to_dwt   ,label='Iron oxide',**Fe_args)[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol g dwt$^{-1}$)')
            l=Fe_ax.plot(result['Total Fe+++']*result['Porosity']*1e3*cm3_to_dwt   ,label='Fe$^{3\!\!+}$',ls='--',**Fe_args)[0]
            
            l=Fe_ax.plot(result['Total Fe++']*result['Porosity']*1e3*cm3_to_dwt ,ls=':'  ,label='Fe$^{2\!\!+}$',**Fe_args)[0]
            Fe_ax.plot((result['Total Sorbed Fe++']+result['Total Sorbed Fe+++'])*1e6/100**3*cm3_to_dwt,**Fe_args,label='Sorbed Fe')
            if obs is not None:
                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        Fe_ax.plot(obs.Incubation_Time[xx],obs['Fe_II'][xx],marker=markers[x],linestyle='None',label='Measured Fe$^{2\!\!+}$',color='C2')
                else:
                    Fe_ax.plot(obs.Incubation_Time,obs['Fe_II'],marker='o',linestyle='None',label='Measured Fe$^{2\!\!+}$',color='C2')
        else:
            l=Fe_ax.plot(result['Fe(OH)3 VF']/molar_volume*1e6   ,label='Iron oxide (solid)')[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol cm$^{-3}$)')
        
            # M/L to umol/cm3: 1e6/1e3=1e3
            l=Fe_ax.plot(result['Total Fe+++']*result['Porosity']*1e3   ,label='Fe+++',ls='--')[0]
            
            l=Fe_ax.plot(result['Total Fe++']*result['Porosity']*1e3 ,ls=':'  ,label='Fe++')[0]
        
        Fe_ax.set_title('Fe species')
        
        Fe_ax.set_xlabel('Time (days)')
        if do_legend:
            Fe_ax.legend(fontsize='small')
    
    if CO2flux_ax is not None:
        # gasflux_ax.set_yscale('log')
        if gdrywt:
            
            l=CO2flux_ax.plot(result.index.values[:-1],numpy.diff(result['Total Tracer']*result['Porosity'])/numpy.diff(result.index.values)*1e3*cm3_to_dwt,label='CO2',**CO2flux_args)[0]
            CO2flux_ax.set_ylabel('CO$_2$ flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
            
            if obs is not None:
                CO2mean=obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff()
                CO2std=obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /obs['Incubation_Time'].diff()
                

                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        CO2flux_ax.errorbar(obs.Incubation_Time[xx],CO2mean[xx],yerr=CO2std[xx],marker=markers[x],linestyle='None',label='Measured CO$_2$',color='C0')
                else:
                    CO2flux_ax.errorbar(obs.Incubation_Time,CO2mean,yerr=CO2std,marker='o',linestyle='None',label='Measured CO$_2$',color='C0')
                
        else:
            
            l=CO2flux_ax.plot(result.index.values[:-1],numpy.diff(result['Total Tracer']*result['Porosity'])/numpy.diff(result.index.values)*1e3,label='CO2',**CO2flux_args)[0]
            CO2flux_ax.set_ylabel('CO$_2$ flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        CO2flux_ax.set_title('CO$_2$ flux')
        
        CO2flux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     CO2flux_ax.legend(fontsize='small')
    if CH4flux_ax is not None:
        if gdrywt:
            l=CH4flux_ax.plot(result.index.values[:-1],numpy.diff(result['Total CH4(aq)']*result['Porosity'])/numpy.diff(result.index.values)*1e3*cm3_to_dwt,label='CH4',**CH4flux_args)[0]
            CH4flux_ax.set_ylabel('CH$_4$ Flux rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
            
            if obs is not None:
                CH4mean=obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff() 
                CH4std=obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().std(axis=1) /obs['Incubation_Time'].diff() 
                
                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        CH4flux_ax.errorbar(obs.Incubation_Time[xx],CH4mean[xx],yerr=CH4std[xx],marker=markers[x],linestyle='None',label='Measured CH$_4$',color='C0')
                else:
                    CH4flux_ax.errorbar(obs.Incubation_Time,CH4mean,yerr=CH4std,marker='o',linestyle='None',label='Measured CH$_4$',color='C0')
                
                
        else:
            l=CH4flux_ax.plot(result.index.values[:-1],numpy.diff(result['Total CH4(aq)']*result['Porosity'])/numpy.diff(result.index.values)*1e3,label='CH4',**CH4flux_args)[0]
            
            CH4flux_ax.set_ylabel('CH$_4$ flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        CH4flux_ax.set_title('CH$_4$ flux')
        
        CH4flux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     CH4flux_ax.legend(fontsize='small')
        
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        porewater_ax.plot(result['Total DOM1'],label='DOM')
        porewater_ax.plot(result['Total Acetate-'],label='Acetate',c='C3')
        porewater_ax.plot(result['Total O2(aq)'],'--',label='O2',c='C4')
        porewater_ax.plot(result['Total Fe+++'],'--',label='Fe+++',c='C1')
        porewater_ax.plot(result['Free Fe+++'],':',label='Fe+++',c='C1')
        porewater_ax.plot(result['Total Fe++'],':',label='Fe++',c='C2')
        
        porewater_ax.set_title('Porewater concentrations')
        porewater_ax.set_ylabel('Concentration (M)')
        porewater_ax.set_xlabel('Time (days)')
        
        
        if obs is not None:
            if obs_marker is not None:
                markers=['o','s','^','+','<','x']
                obs_cats=unique(obs[obs_marker])
                print(obs_cats)
                for x in range(len(obs_cats)):
                    xx=obs[obs_marker]==obs_cats[x]
                    porewater_ax.plot(obs.Incubation_Time[xx],obs['TOAC'][xx]*1e-6/BD*1000/porosity,marker=markers[x],linestyle='None',label='Measured TOAC',color='C3',**porewater_args)
                    porewater_ax.plot(obs.Incubation_Time[xx],obs['WEOC'][xx]*1e-6/BD*1000/porosity,marker=markers[x],linestyle='None',label='Measured WEOC',color='C0',**porewater_args)
            else:
                porewater_ax.plot(obs.Incubation_Time,obs['TOAC']*1e-6/BD*1000/porosity,marker='o',linestyle='None',label='Measured TOAC',color='C3')
                porewater_ax.plot(obs.Incubation_Time,obs['WEOC']*1e-6/BD*1000/porosity,marker='o',linestyle='None',label='Measured WEOC',color='C0')


        
        if do_legend:
            porewater_ax.legend(fontsize='small')
            
    if cation_ax is not None:
        for n,cation in enumerate(decomp_network.pools_list_to_dict(reactions)['Cation exchange']['cations'].keys() ):
            if cation == 'H+':
                continue
            cation_ax.plot(result['Total Sorbed '+cation]/(BD*1e-3*100**3)*1000,c='C'+str(n),label=cation,**cation_args)
        cation_ax.plot(result['CEC H+']/(BD*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),**cation_args)
        
        cation_ax.set_ylabel('Exch conc\n(mmol/kg)')
        
        if do_legend:
            cation_ax.legend(fontsize='small')
        
    


colors={'Anaerobic':'C1','Periodic':'C2','Low Fe':'C3','Aerobic':'C0'}


# Figures for talk
# fig,axes=subplots(3,1,num='Just oxygen',figsize=(6,7.4),clear=True)
# plot_result(result_lowFe,SOM_ax=axes[0])
# axes[0].set_ylim(bottom=-0.01)
# 
# axes[1].plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total Tracer']*result_lowFe['Porosity'])/diff(result_lowFe.index.values)*1e3,label='CO2',ls='--')
# # axes[2].plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total CH4(aq)']*result_lowFe['Porosity'])/diff(result_lowFe.index.values)*1e6,label='CH4')
# methane_simple=zeros(len(result_lowFe.index.values[:-1]))
# methane_simple[24*100:]=0.02
# axes[2].plot(result_lowFe.index.values[:-1],methane_simple)
# 
# axes[1].set_title('CO$_2$ flux')
# axes[2].set_title('CH$_4$ flux')
# axes[1].set_ylabel('Flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')
# axes[2].set_ylabel('Flux rate\n(nmol cm$^{-3}$ day$^{-1}$)')
# axes[2].set_xlabel('Time (days)')
# axes[2].set_ylim(top=0.031)
# 
# xmax=axes[0].get_xlim()[1]
# for ax in axes:
#     ax.axvspan(100,xmax,color='b',alpha=0.1)
#     ax.set_xlim(right=xmax)
#     ax.text(90,ax.get_ylim()[1]*0.9,'Aerobic',ha='right')
#     ax.text(110,ax.get_ylim()[1]*0.9,'Inundated',ha='left')


# With Fe reduction
# fig,axes=subplots(3,1,num='Organic horizon Anoxic',figsize=(6,8.4),clear=True)
f,cation_axes=subplots(ncols=5,num='Cations',clear=True,figsize=(13,4))
fig,axs=subplots(nrows=4,ncols=5,num='Time series plots',clear=True,figsize=(16,8),sharex=False)
from string import ascii_lowercase
for x in range(5):
    for y in range(4):
        axs[y,x].text(0.05,1.05,'('+ascii_lowercase[x*4+y]+')',transform=axs[y,x].transAxes)

axes=axs[:,1]
obs=get_layer(9,'Organic')
plot_result(result_organic_trough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],cation_ax=cation_axes[1],do_legend=False,BD=BD_organic_trough,
            gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,#BD=BD_layerest2_trough['Organic',True],
            # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')
                    # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')])
                    obs=obs)
# axes[0].set_ylim(bottom=-0.01)

# axes[0].set_ylim(bottom=0.9e-11,top=15)
axes[0].set_yscale('linear')
axes[0].text(0.5,1.2,'%s (%s)\nOrganic horizon (Anoxic)'%(obs['Core_ID'].iloc[0],obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')


# xmax=axes[0].get_xlim()[1]
xmax=max(simlength,90)
for ax in axes:
    ax.axvspan(initfrac*simlength,xmax,color='b',alpha=0.1)
    ax.set_xlim(left=-5,right=xmax)
    
axes=axs[:,2]
obs=get_layer(5,'Organic')
plot_result(result_organic_nottrough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],cation_ax=cation_axes[2],do_legend=False,
            gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,BD=BD_organic_nottrough,
            # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')
                    # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()!='trough')])
                    obs=obs)
# axes[0].set_ylim(bottom=-0.01)

# axes[0].set_ylim(bottom=0.9e-11,top=15)
axes[0].set_yscale('linear')
axes[0].text(0.5,1.2,'%s (%s)\nOrganic horizon (Anoxic)'%(obs['Core_ID'].iloc[0],obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')


# xmax=axes[0].get_xlim()[1]
# xmax=90
for ax in axes:
    ax.axvspan(initfrac*simlength,xmax,color='b',alpha=0.1)
    ax.set_xlim(left=-5,right=xmax)


# fig,axes=subplots(3,1,num='Mineral horizon Anoxic',figsize=(6,8.4),clear=True)
axes=axs[:,3]
obs=get_layer(9,'Mineral')
plot_result(result_mineral_trough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],cation_ax=cation_axes[3],do_legend=False,
    gdrywt=True,SOC_pct=SOC_layermean['Mineral'],cellulose_SOC_frac=cellulosefrac,BD=BD_mineral_trough,
    # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')
        # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')])
        obs=obs)
# axes[0].set_ylim(bottom=-0.01)

# axes[0].set_ylim(bottom=0.9e-11,top=3.0)
axes[0].set_yscale('linear')
axes[0].text(0.5,1.2,'%s (%s)\nMineral horizon (Anoxic)'%(obs['Core_ID'].iloc[0],obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

# xmax=axes[0].get_xlim()[1]
# xmax=90
for ax in axes:
    ax.axvspan(initfrac*simlength,xmax,color='b',alpha=0.1)
    ax.set_xlim(left=-5,right=xmax)

axes=axs[:,4]
obs=get_layer(5,'Mineral')
plot_result(result_mineral_nottrough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],cation_ax=cation_axes[4],do_legend=False,
    gdrywt=True,SOC_pct=SOC_layermean['Mineral'],cellulose_SOC_frac=cellulosefrac,BD=BD_mineral_nottrough,
    # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')&
        # (Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()!='trough')])
        obs=obs)
# axes[0].set_ylim(bottom=-0.01)

# axes[0].set_ylim(bottom=0.9e-11,top=3.0)
axes[0].set_yscale('linear')
axes[0].text(0.5,1.2,'%s (%s)\nMineral horizon (Anoxic)'%(obs['Core_ID'].iloc[0],obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

# xmax=axes[0].get_xlim()[1]
# xmax=90
for ax in axes:
    ax.axvspan(initfrac*simlength,xmax,color='b',alpha=0.1)
    ax.set_xlim(left=-5,right=xmax)


# Oxic
# fig,axes=subplots(3,1,num='Organic horizon Oxic',figsize=(6,8.4),clear=True)
axes=axs[:,0]
obs=get_layer(3,'Organic')
plot_result(result_highO2_organic,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],cation_ax=cation_axes[0],do_legend=False,
        gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,BD=BD_atmoO2_organic,
        # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')&(Barrow_synthesis_data['Incubation_Temperature']>4)])
        obs=obs)
# axes[0].set_ylim(bottom=-0.01,top=30)
axes[1].set_ylim(bottom=-0.01,top=0.2)

axes[0].set_yscale('linear')
# axes[0].set_ylim(-5e-2,155e-2) # This leaves out Core NGADG0003 which has 10x higher flux for some reason
# axes[1].set_ylim(bottom=0.9e-11)
axes[0].text(0.5,1.2,'%s (%s)\nOrganic horizon (Oxic)'%(obs['Core_ID'].iloc[0],obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

# axes[0].set_xlim(right=71)
for ax in fig.axes:
    ax.set_xlim(left=-5,right=xmax)
for ax in axs[3,:]:
    ax.set_ylim(3.8,6.1)

leg=axs[2,4].legend(loc=(0.03,0.45),ncol=2,edgecolor='k')
cation_axes[0].legend()

# Measurements did not include oxic mineral horizon incubation
# 
# # Aerobic
# fig,axes=subplots(3,1,num='Mineral horizon Oxic',figsize=(6,8.4),clear=True)
# plot_result(result_highO2_mineral,Fe_ax=axes[1],CO2flux_ax=axes[0],pH_ax=axes[2],do_legend=True,
#             gdrywt=True,SOC_pct=SOC_layermean['Mineral'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2_trough['Mineral',False],
#             obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')&(Barrow_synthesis_data['Incubation_Temperature']>4)])
# # axes[0].set_ylim(bottom=-0.01,top=30)
# 
# axes[0].set_yscale('linear')
# axes[0].set_ylim(-5e-2,155e-2) # This leaves out Core NGADG0003 which has 10x higher flux for some reason
# # axes[1].set_ylim(bottom=0.9e-11)
# 
# for ax in axes:
#     axes[0].set_xlim(left=-5,right=71)
#     axes[1].set_xlim(left=-5,right=71)
#     axes[2].set_xlim(left=-5,right=71)



# Periodic inundation
# rateconstants['1.0e+00 Fe++  + 2.5e-01 O2(aq)  + 1.0e+00 H+  <-> 1.0e+00 Fe+++  + 5.0e-01 H2O'                 ] =  20.0e0
# rateconstants["cellulose decay to CO2 (SOMDEC sandbox)" ]=                                                      3.17e-8*5.0
# rateconstants["cellulose decay to DOM1 (SOMDEC sandbox)" ]=                                                      3.17e-8*1.0
# rateconstants["1.0e+00 Acetate-  + 8.0e+00 Fe+++  -> 2.0e+00 HCO3-  + 8.0e+00 Fe++  + 9.0e+00 H+  + 2.0e+00 Tracer"] = 1.0e-8*5.0
# pools_atmoO2_organic2=decomp_network.change_constraints(pools_atmoO2_organic,{'O2(aq)':'1.0 G O2(g)'})


fig,axes=subplots(4,1,num='Periodic inundation',figsize=(6,8.4),clear=True)
# CH4ax=axes[0].twinx()
plot_result(result_periodicO2,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],#porewater_ax=axes[3],
        gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2_trough['Organic',False])
# axes[0].set_ylim(bottom=-0.01)
axes[2].legend(labels=['Fe(OH)$_3$','Fe$^{3\!\!+}$','Fe$^{2\!\!+}$','Sorbed Fe'],loc=(0.83,0.2),fontsize='small')

axes[0].set_ylim(bottom=0.9e-11)
# axes[0].legend()
# CH4ax.legend()
# axes[3].legend(loc='lower left')
t=fig.axes[0].text(0.01,1.05,'(a)',transform=fig.axes[0].transAxes)  
t=fig.axes[1].text(0.01,1.05,'(b)',transform=fig.axes[1].transAxes)  
t=fig.axes[2].text(0.01,1.05,'(c)',transform=fig.axes[2].transAxes)  
t=fig.axes[3].text(0.01,1.05,'(d)',transform=fig.axes[3].transAxes)  


xmax=axes[0].get_xlim()[1]
O2_periodic=numpy.zeros(simlength*24)
O2_periodic[:int(simlength*24*oxicfrac)]=dq
for ax in axes:
    for num in range(3):
        ax.axvspan(num*simlength+nonzero(diff(O2_periodic))[0]/24,(num+1)*simlength,color='b',alpha=0.1)
    ax.set_xlim(right=xmax,left=0)







cm=pyplot.get_cmap('plasma')
norm=matplotlib.colors.Normalize(min(periodic_out[1].keys()),max(periodic_out[1].keys()))
fig,axs=subplots(ncols=2,nrows=1,num='Periodic by dry periods',clear=True,squeeze=False,figsize=(7,4))
for pH in periodic_out[1]:
    axs[0,0].plot(array([1,2,3,4]),array([periodic_out[n][pH]['Total Tracer'].iloc[simlength*24-1] for n in range(1,5)])*result_periodicO2['Porosity'].iloc[1]/BD_layerest2_trough['Organic',False],'o',label=pH,c=cm(norm(pH)))
    axs[0,1].plot(array([1,2,3,4]),array([periodic_out[n][pH]['Total CH4(aq)'].iloc[simlength*24-1] for n in range(1,5)])*result_periodicO2['Porosity'].iloc[1]/BD_layerest2_trough['Organic',False],'o',label=pH,c=cm(norm(pH)))

axs[0,0].set_xlabel('Number of oxic-anoxic cycles')
axs[0,1].set_xlabel('Number of oxic-anoxic cycles')
axs[0,0].set_ylabel('Cumulative CO$_2$ flux\n(mmol g dwt$^{-1}$)')
axs[0,1].set_ylabel('Cumulative CH$_4$ flux\n(mmol g dwt$^{-1}$)')
axs[0,0].set_title('Cumulative CO$_2$ flux')
axs[0,1].set_title('Cumulative CH$_4$ flux')
axs[0,0].set_ylim(bottom=0)
axs[0,1].set_ylim(bottom=0)
axs[0,1].legend(title='Initial pH')

axs[0,0].text(0.01,1.02,'(a)',transform=axs[0,0].transAxes)
axs[0,1].text(0.01,1.02,'(b)',transform=axs[0,1].transAxes)


ndry_plotted=[3]
fig,axs=subplots(ncols=len(ndry_plotted),nrows=5,num='Periodic comparison',clear=True,figsize=(10,8.8),squeeze=False)


for nn,ndry in enumerate(ndry_plotted):
    c='k'
    plot_result(periodic_out[ndry][pH_obs-1],CH4flux_ax=axs[1,nn],gdrywt=True,BD=BD_layerest2_trough['Organic',False],SOC_pct=SOC_layermean['Organic'],CH4flux_args={'color':c,'ls':'--'})
    plot_result(periodic_out[ndry][pH_obs],CH4flux_ax=axs[1,nn],gdrywt=True,BD=BD_layerest2_trough['Organic',False],SOC_pct=SOC_layermean['Organic'],CH4flux_args={'color':c})
    plot_result(periodic_out[ndry][pH_obs+1],CH4flux_ax=axs[1,nn],gdrywt=True,BD=BD_layerest2_trough['Organic',False],SOC_pct=SOC_layermean['Organic'],CH4flux_args={'color':c,'ls':':'})
    
    plot_result(periodic_out[ndry][pH_obs-1],CO2flux_ax=axs[0,nn],gdrywt=True,BD=BD_layerest2_trough['Organic',False],SOC_pct=SOC_layermean['Organic'],CO2flux_args={'color':c,'ls':'--'})
    plot_result(periodic_out[ndry][pH_obs],CO2flux_ax=axs[0,nn],gdrywt=True,BD=BD_layerest2_trough['Organic',False],SOC_pct=SOC_layermean['Organic'],CO2flux_args={'color':c})
    plot_result(periodic_out[ndry][pH_obs+1],CO2flux_ax=axs[0,nn],gdrywt=True,BD=BD_layerest2_trough['Organic',False],SOC_pct=SOC_layermean['Organic'],CO2flux_args={'color':c,'ls':':'})

    axs[0,0].legend(labels=[pH_obs-1,pH_obs,pH_obs+1],title='Initial pH')
    axs[0,nn].set_xlim(0,simlength)
    axs[1,nn].set_xlim(0,simlength)
    axs[0,nn].set_title('CO$_2$ flux rate')
    axs[1,nn].set_title('CH$_4$ flux rate')

    axs[2,nn].set_ylim(bottom=0,top=50)
    axs[2,nn].set_xlabel('Time (days)')
    axs[2,nn].set_ylabel('Fe(II) production rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
    axs[2,nn].set_title('Fe(II) production rate')
    
    
    axs[2,nn].plot((periodic_out[ndry][pH_obs]['Total Fe++'].diff()*24*1e3*result_periodicO2['Porosity'].iloc[1]+periodic_out[ndry][pH_obs]['Total Sorbed Fe++'].diff()*24*1e6/100**3)/BD_layerest2_trough['Organic',False],'-',c=c,label=str(ndry))
    axs[2,nn].plot((periodic_out[ndry][pH_obs-1]['Total Fe++'].diff()*24*1e3*result_periodicO2['Porosity'].iloc[1]++periodic_out[ndry][pH_obs-1]['Total Sorbed Fe++'].diff()*24*1e6/100**3)/BD_layerest2_trough['Organic',False],'--',c=c)
    axs[2,nn].plot((periodic_out[ndry][pH_obs+1]['Total Fe++'].diff()*24*1e3*result_periodicO2['Porosity'].iloc[1]++periodic_out[ndry][pH_obs+1]['Total Sorbed Fe++'].diff()*24*1e6/100**3)/BD_layerest2_trough['Organic',False],':',c=c)
        
    axs[3,nn].plot(periodic_out[ndry][pH_obs]['Fe(OH)3 VF']/molar_volume_FeOH3*1e3/BD_layerest2_trough['Organic',False],'-',c=c)
    axs[3,nn].plot(periodic_out[ndry][pH_obs-1]['Fe(OH)3 VF']/molar_volume_FeOH3*1e3/BD_layerest2_trough['Organic',False],'--',c=c)
    axs[3,nn].plot(periodic_out[ndry][pH_obs+1]['Fe(OH)3 VF']/molar_volume_FeOH3*1e3/BD_layerest2_trough['Organic',False],':',c=c)
    axs[3,nn].set_xlabel('Time (days)')
    axs[3,nn].set_ylabel('Fe oxides \n(mmol g dwt$^{-1}$)')
    axs[3,nn].set_title('Fe oxide minerals ')
    # axs[3,nn].set_ylim(-0.1,3.4)

    axs[4,nn].plot(-log10(periodic_out[ndry][pH_obs]['Free H+']),'-',c=c)
    axs[4,nn].plot(-log10(periodic_out[ndry][pH_obs-1]['Free H+']),'--',c=c)
    axs[4,nn].plot(-log10(periodic_out[ndry][pH_obs+1]['Free H+']),':',c=c)


    # axs[2,1].set_ylim(bottom=0,top=50)
    axs[4,nn].set_xlabel('Time (days)')
    axs[4,nn].set_ylabel('pH')
    axs[4,nn].set_title('pH')
    
    for x in range(5):
        for num in range(ndry):
            axs[x,nn].axvspan(simlength*oxicfrac/ndry+num*simlength/ndry,(num+1)*simlength/ndry,color='b',alpha=0.1,zorder=-1)
        axs[x,nn].set_xlim(0,simlength)
        axs[x,nn].text(0.01,1.05,'('+ascii_lowercase[x*len(ndry_plotted)+nn]+')',transform=axs[x,nn].transAxes)

plot_obs_data=False
if plot_obs_data:

    inc_time=Barrow_synthesis_data.Incubation_Time
    inc_T=Barrow_synthesis_data.Incubation_Temperature
    inc_anoxic=Barrow_synthesis_data.Headspace=='Anoxic'
    
    

    CO2mean=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/inc_time.diff()/(Barrow_synthesis_data['SOC'].astype(float)/100/12)   
    CO2std=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /inc_time.diff()/(Barrow_synthesis_data['SOC'].astype(float)/100/12)   
    CH4mean=Barrow_synthesis_data[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/inc_time.diff()/(Barrow_synthesis_data['SOC'].astype(float)/100/12)   
    CH4std=Barrow_synthesis_data[['CH4_1','CH4_2','CH4_3']].astype(float).diff().std(axis=1) /inc_time.diff()/(Barrow_synthesis_data['SOC'].astype(float)/100/12)   

    f,a=pyplot.subplots(nrows=6,ncols=10,squeeze=True,clear=True,num='Incubation data',figsize=(15,8))
    cm=pyplot.get_cmap('coolwarm')
    norm=matplotlib.colors.Normalize(inc_T.min(),inc_T.max())
    markers={'Organic':'o','Mineral':'s','Transition':'^','Permafrost':'x'}
    for num in range(10):
        column=Barrow_synthesis_data['Core_ID'].unique()[num]
        for temp in inc_T.unique():
            c=cm(norm(temp))
            for layer in ['Organic','Mineral','Transition','Permafrost']:
                marker=markers[layer]
                xx=(inc_T==temp)&(inc_anoxic)&(Barrow_synthesis_data['Core_ID']==column)&(Barrow_synthesis_data['Soil_layer'].str.capitalize()==layer)
                a[0,num].errorbar(inc_time[xx],CO2mean[xx],yerr=CO2std[xx],label='T = %d C, %s'%(temp,layer),marker=marker,linestyle='-',c=c)
                a[0,num].set_title(column)
                a[1,num].errorbar(inc_time[xx],CH4mean[xx],yerr=CH4std[xx],label='T = %d C, %s'%(temp,layer),marker=marker,linestyle='-',c=c)
                a[2,num].errorbar(inc_time[xx],Barrow_synthesis_data['Fe_II'].astype(float)[xx],marker=marker,linestyle='-',c=c)
                a[3,num].errorbar(inc_time[xx],Barrow_synthesis_data['TOAC'].astype(float)[xx],marker=marker,linestyle='-',c=c)
                a[3,num].errorbar(inc_time[xx],Barrow_synthesis_data['WEOC'].astype(float)[xx],marker='o',linestyle='--',c=c)
                a[4,num].errorbar(inc_time[xx],Barrow_synthesis_data['pH'].astype(float)[xx],marker=marker,linestyle='-',c=c)
                
                xx=(inc_T==temp)&(~inc_anoxic)&(Barrow_synthesis_data['Core_ID']==column)&(Barrow_synthesis_data['Soil_layer'].str.capitalize()==layer)
                a[5,num].errorbar(inc_time[xx],CO2mean[xx],yerr=CO2std[xx],label='T = %d C, %s'%(temp,layer),marker=marker,linestyle='-',c=c)


        
    a[0,0].set_ylabel('CO$_2$ flux (anoxic)\n($\mu$mol g$^{-1}$ dwt day$^{-1}$)')
    a[5,0].set_ylabel('CO$_2$ flux (oxic)\n($\mu$mol g$^{-1}$ dwt day$^{-1}$)')
    a[1,0].set_ylabel('CH$_4$ flux\n($\mu$mol g$^{-1}$ dwt day$^{-1}$)')
    a[2,0].set_ylabel('Fe(II)\n($\mu$mol g$^{-1}$ dwt)')
    a[3,0].set_ylabel('Org acids\n($\mu$mol g$^{-1}$ dwt)')
    a[4,0].set_ylabel('pH')

    # a[1,0].legend()

    for ax in a.ravel():
        ax.set_xlabel('Time (days)')

pyplot.show()
