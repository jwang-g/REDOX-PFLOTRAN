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

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3+',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)4-',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeCH3COO+',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIDOM1(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeIIIAcetate(aq)',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)2(aq)',kind='secondary'),
# decomp_network.decomp_pool(name='Fe(OH)2+',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3(aq)',kind='secondary'),

# See Roberts Earth Sci Rev article for discussion of iron minerals, including some time scales
decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-3 mol/m^2-sec',constraints={'initial':'0.875d-3  1.d2 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Goethite',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d1 m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Fe',rate='1.d-7 mol/m^2-sec',constraints={'initial':'1.0e-6  1. m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.d2 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

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
}
anox_inhib_conc = 1e-8

truncate_conc=1e-15

reactions = [
    decomp_network.reaction(name='Aerobic decomposition',reactant_pools={'cellulose':1.0},product_pools={'HCO3-':1.0},reactiontype='SOMDECOMP',
                                            rate_constant=1e-1,rate_units='y', 
                                            # inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='THRESHOLD 1.0d20')]),
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=conc_scales['O2(aq)'],threshold=1.1e-9)]),
    
    decomp_network.reaction(name='Hydrolysis',stoich='1.0 cellulose -> 1.0 DOM1',reactiontype='SOMDECOMP',
                                            rate_constant=1e-1,rate_units='y', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                        inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=conc_scales['DOM1']),
                                                          decomp_network.inhibition(species='O2(aq)',k=6.25e-11,type='THRESHOLD -1.0d10')
                                                          # decomp_network.inhibition(species='O2(aq)',type='MONOD',k=1e-11),
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

    # Acetoclastic methanogenesis
    # C2H3O2- + H+ -> CH4 + HCO3- + H+
    decomp_network.reaction(name='Methanogenesis',stoich='1.0 Acetate- -> 1.0 CH4(aq) + 1.0 HCO3- + 1.0 Tracer',
                                            monod_terms=[decomp_network.monod(species='Acetate-',k=conc_scales['Acetate-'],threshold=1.1e-15),
                                                         ],
                                            inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=anox_inhib_conc,type='MONOD'),decomp_network.inhibition(species='Fe+++',k=conc_scales['Fe+++'],type='MONOD')],
                                            rate_constant=1e-10,reactiontype='MICROBIAL'),

    # Hydrogenotrophic methanogenesis
    
    # H2 oxidation if oxygen available

]

rateconstants={
       '1.0e+00 Fe++  + 2.5e-01 O2(aq)  + 1.0e+00 H+  <-> 1.0e+00 Fe+++  + 5.0e-01 H2O'                      : 1.0e0*1.0e-1,
       '1.0e+00 DOM1  -> 3.3e-01 Acetate-  + 3.3e-01 HCO3-  + 6.7e-01 H+  + 1.3e+00 H2(aq)  + 3.3e-01 Tracer'  : 6.0e-8*20.0,
       "1.0e+00 DOM1  + 1.0e+00 O2(aq)  -> 1.0e+00 HCO3-  + 1.0e+00 H+  + 1.0e+00 Tracer"                    : 3.0e-8*1.0,
       "1.0e+00 Acetate-  + 2.0e+00 O2(aq)  -> 2.0e+00 HCO3-  + 2.0e+00 H+  + 2.0e+00 Tracer"                : 3.0e-8*1.0,
       "1.0e+00 Acetate-  + 8.0e+00 Fe+++  -> 2.0e+00 HCO3-  + 8.0e+00 Fe++  + 9.0e+00 H+  + 2.0e+00 Tracer" : 1.0e-8*5.0,
       "1.0e+00 Acetate-  -> 1.0e+00 CH4(aq)  + 1.0e+00 HCO3-  + 1.0e+00 Tracer"                             : 1.0e-9*10.0,
       "cellulose decay to CO2 (SOMDEC sandbox)"                                                              : 3.17e-8*5.0, #1.0/(365*24*3600)*1.0
       "cellulose decay to DOM1 (SOMDEC sandbox)"                                                              : 3.17e-8*1.0 #1.0/(365*24*3600)*1.0
}

reaction_network =  decomp_network.decomp_network(pools=pools,reactions=reactions)

networkfig2=pyplot.figure('Reaction network (with reactions)',clear=True)
drawn=decomp_network.draw_network_with_reactions(reaction_network,node_size=700,arrowstyle='->',arrowsize=8.5,edge_color='gray',
            omit=['NH4+','Rock(s)','gas','surf_complex','secondary','H+','implicit','H2(aq)'],
            namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'CO2(aq)',
                         'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate'})

# Generate PFLOTRAN input file with correct reactions
decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt','Arctic_redox.in',log_formulation=False,CO2name='Tracer',truncate_concentration=1e-25)
            
dq=0.01 # Diffusion coefficient when aerobic
O2_const=numpy.zeros(365*24)+dq

simlength=70
initfrac=0.0
timestep=3600
nperiodic=3

O2_initial=numpy.zeros(simlength*24)
O2_initial[:int(simlength*24*initfrac)]=dq

# Incubation data
import pandas,glob
datadir='/Users/b0u/Documents/NGEE-Arctic/Redox_sims/Data'
Barrow_synthesis_data=pandas.concat([pandas.read_csv(fname,na_values=-9999,header=8,skiprows=[9]) for fname in glob.glob(datadir+'/Barrow_soil_geochem_synthesis/*.csv')])  
units=pandas.read_csv(glob.glob(datadir+'/Barrow_soil_geochem_synthesis/*.csv')[0],na_values=-9999,header=8).iloc[0]
SOC_layermean=Barrow_synthesis_data.groupby(Barrow_synthesis_data['Soil_layer'].str.capitalize())['SOC'].mean()

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
BD_layerest2=numpy.exp((8.2432-SOC_layermean)/9.7872)

cellulosefrac=0.05
molar_volume_FeOH3=34.3600 #cm3/mol
porosity=0.25
# Assuming all iron is converted to FeII by end of incubation (?)
Fe_VF_organic=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Soil_layer']=='Organic')]['Fe_II'].max()*1e-6*BD_layerest2['Organic']*molar_volume_FeOH3
pools_organic=decomp_network.change_constraints(pools,{'cellulose':SOC_layermean['Organic']/100/12*BD_layerest2['Organic']*100**3*cellulosefrac,
                                                      'Fe(OH)3':'%1.4f  1.e2 m^2/m^3'%Fe_VF_organic}) #
pools_organic=decomp_network.change_site_density(pools_organic, '>Carboxylate-', SOC_layermean['Organic']*2e2)
pools_atmoO2_organic=decomp_network.change_constraints(pools_organic,{'O2(aq)':'1.0 G O2(g)','H+':'4.5 P','DOM1':1e-20,
                    'cellulose':SOC_layermean['Organic']/100/12*BD_layerest2['Organic']*100**3*cellulosefrac*0.04, # Accounts for soils being pretty dry
                    'Fe++':Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data.Incubation_Time==0)&(Barrow_synthesis_data['Incubation_Temperature']>4)]\
                    .drop_duplicates(['Core_ID','Moisture','Fe_II'])['Fe_II'].mean()*1e-6*BD_layerest2['Organic']/porosity*1e3},inplace=False)
pools_atmoO2_organic=decomp_network.change_site_density(pools_atmoO2_organic, '>Carboxylate-', SOC_layermean['Organic']*0.5e1)
pools_lowFe_organic=decomp_network.change_constraints(pools_organic,{'Fe(OH)3':'0.0  1.e2 m^2/m^3','Fe+++':1.0e-10,'Fe++':1.0e-10},inplace=False)


Fe_VF_mineral=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')]['Fe_II'].max()*1e-6*BD_layerest2['Mineral']*molar_volume_FeOH3
Fe_VF_mineral=100*1e-6*BD_layerest2['Mineral']*molar_volume_FeOH3

pools_mineral=decomp_network.change_constraints(pools,{'cellulose':SOC_layermean['Mineral']/100/12*BD_layerest2['Mineral']*100**3*cellulosefrac,
                                                      'Fe(OH)3':'%1.4f  1.e2 m^2/m^3'%Fe_VF_mineral}) #
pools_mineral=decomp_network.change_site_density(pools_mineral, '>Carboxylate-', SOC_layermean['Mineral']*2e2)

pools_atmoO2_mineral=decomp_network.change_constraints(pools_mineral,{'O2(aq)':'1.0 G O2(g)','H+':'4.5 P','DOM':1e-20,
                    'cellulose':SOC_layermean['Mineral']/100/12*BD_layerest2['Mineral']*100**3*cellulosefrac*0.04, # Accounts for soils being pretty dry
                    'Fe++':Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data.Incubation_Time==0)&(Barrow_synthesis_data['Incubation_Temperature']>4)]\
                    .drop_duplicates(['Core_ID','Moisture','Fe_II'])['Fe_II'].mean()*1e-6*BD_layerest2['Mineral']/porosity*1e3},inplace=False)
pools_atmoO2_mineral=decomp_network.change_site_density(pools_atmoO2_mineral, '>Carboxylate-', SOC_layermean['Mineral']*0.5e1)

pools_lowFe_mineral=decomp_network.change_constraints(pools_mineral,{'Fe(OH)3':'0.0  1.e2 m^2/m^3','Fe+++':1.0e-10,'Fe++':1.0e-10},inplace=False)

pools_permafrost=decomp_network.change_constraint(pools,'cellulose',SOC_layermean['Permafrost']/100/12*BD_layerest2['Permafrost']*100**3)



result_highO2_organic,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_atmoO2_organic,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc)
result_organic,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_organic,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc)    
# result_lowFe_organic,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_lowFe_organic,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc)

result_highO2_mineral,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_atmoO2_mineral,bc=pools_atmoO2_mineral,diffquo={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc)
result_mineral,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength,timestep,initcond=pools_mineral,hands_off=False,rateconstants=rateconstants,bc=pools_atmoO2_mineral,diffquo={'O2(aq)':O2_initial},truncate_concentration=truncate_conc)    


import plot_pf_output
from pylab import *
    
from run_alquimia import plot_result


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
fig,axes=subplots(4,1,num='Organic horizon Anoxic',figsize=(6,8.4),clear=True)
plot_result(result_organic,Fe_ax=axes[1],gasflux_ax=axes[0],pH_ax=axes[2],porewater_ax=axes[3],gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2['Organic'])
# axes[0].set_ylim(bottom=-0.01)
axes[1].legend(loc='right',labels=['Iron oxide (solid)','Fe$^{3\!\!+}$ (dissolved)','Fe$^{2\!\!+}$ (dissolved)'])

axes[0].set_ylim(bottom=0.9e-11,top=15)
axes[0].set_yscale('linear')
# axes[0].text(90,0.0,'Aerobic',ha='right')
# axes[0].text(110,0.0,'Inundated',ha='left')
axes[0].legend(labels=['CH$_4$','CO$_2$'])
axes[3].legend(loc='lower left')

xx=(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')&(Barrow_synthesis_data['Incubation_Temperature']>4)
CO2mean=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/Barrow_synthesis_data['Incubation_Time'].diff()
CO2std=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /Barrow_synthesis_data['Incubation_Time'].diff()
CH4mean=Barrow_synthesis_data[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/Barrow_synthesis_data['Incubation_Time'].diff() 
CH4std=Barrow_synthesis_data[['CH4_1','CH4_2','CH4_3']].astype(float).diff().std(axis=1) /Barrow_synthesis_data['Incubation_Time'].diff() 


axes[0].errorbar(Barrow_synthesis_data.Incubation_Time[xx],CO2mean[xx],yerr=CO2std[xx],marker='o',linestyle='None',label='Measured',color='C5')
axes[0].errorbar(Barrow_synthesis_data.Incubation_Time[xx],CH4mean[xx],yerr=CH4std[xx],marker='o',linestyle='None',label='Measured',color='C0')
axes[1].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['Fe_II'][xx],marker='o',linestyle='None',label='Measured',color='C2')

axes[2].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['pH'][xx],marker='o',linestyle='None',label='Measured',color='C0')
axes[3].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['TOAC'][xx]*1e-6/BD_layerest2['Organic']*1000/porosity,marker='o',linestyle='None',label='Measured TOAC',color='C3')
axes[3].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['WEOC'][xx]*1e-6/BD_layerest2['Organic']*1000/porosity,marker='o',linestyle='None',label='Measured WEOC',color='C0')


# xmax=axes[0].get_xlim()[1]
xmax=90
for ax in axes:
    ax.axvspan(initfrac*simlength,xmax,color='b',alpha=0.1)
    ax.set_xlim(right=xmax)


fig,axes=subplots(4,1,num='Mineral horizon Anoxic',figsize=(6,8.4),clear=True)
plot_result(result_mineral,Fe_ax=axes[1],gasflux_ax=axes[0],pH_ax=axes[2],porewater_ax=axes[3],gdrywt=True,SOC_pct=SOC_layermean['Mineral'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2['Mineral'])
# axes[0].set_ylim(bottom=-0.01)
axes[1].legend(loc='right',labels=['Iron oxide (solid)','Fe$^{3\!\!+}$ (dissolved)','Fe$^{2\!\!+}$ (dissolved)'])

axes[0].set_ylim(bottom=0.9e-11,top=3.0)
axes[0].set_yscale('linear')
# axes[0].text(90,0.0,'Aerobic',ha='right')
# axes[0].text(110,0.0,'Inundated',ha='left')
axes[0].legend(labels=['CH$_4$','CO$_2$'])
axes[3].legend(loc='lower left')

xx=(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')&(Barrow_synthesis_data['Incubation_Temperature']>4)
CO2mean=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/Barrow_synthesis_data['Incubation_Time'].diff()
CO2std=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /Barrow_synthesis_data['Incubation_Time'].diff()
CH4mean=Barrow_synthesis_data[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/Barrow_synthesis_data['Incubation_Time'].diff() 
CH4std=Barrow_synthesis_data[['CH4_1','CH4_2','CH4_3']].astype(float).diff().std(axis=1) /Barrow_synthesis_data['Incubation_Time'].diff() 


axes[0].errorbar(Barrow_synthesis_data.Incubation_Time[xx],CO2mean[xx],yerr=CO2std[xx],marker='o',linestyle='None',label='Measured',color='C5')
axes[0].errorbar(Barrow_synthesis_data.Incubation_Time[xx],CH4mean[xx],yerr=CH4std[xx],marker='o',linestyle='None',label='Measured',color='C0')
axes[1].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['Fe_II'][xx],marker='o',linestyle='None',label='Measured',color='C2')

axes[2].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['pH'][xx],marker='o',linestyle='None',label='Measured',color='C0')
axes[3].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['TOAC'][xx]*1e-6/BD_layerest2['Mineral']*1000/porosity,marker='o',linestyle='None',label='Measured TOAC',color='C3')
axes[3].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['WEOC'][xx]*1e-6/BD_layerest2['Mineral']*1000/porosity,marker='o',linestyle='None',label='Measured WEOC',color='C0')


# xmax=axes[0].get_xlim()[1]
xmax=90
for ax in axes:
    ax.axvspan(initfrac*simlength,xmax,color='b',alpha=0.1)
    ax.set_xlim(right=xmax)


# Aerobic
fig,axes=subplots(3,1,num='Organic horizon Oxic',figsize=(6,8.4),clear=True)
plot_result(result_highO2_organic,Fe_ax=axes[1],gasflux_ax=axes[0],pH_ax=axes[2],gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2['Organic'])
# axes[0].set_ylim(bottom=-0.01,top=30)
axes[1].legend(loc='right',labels=['Iron oxide (solid)','Fe$^{3\!\!+}$ (dissolved)','Fe$^{2\!\!+}$ (dissolved)'])

xx=(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')&(Barrow_synthesis_data['Incubation_Temperature']>4)
CO2mean=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/Barrow_synthesis_data['Incubation_Time'].diff() 
CO2std=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /Barrow_synthesis_data['Incubation_Time'].diff()

axes[0].errorbar(Barrow_synthesis_data.Incubation_Time[xx],CO2mean[xx],yerr=CO2std[xx],marker='o',linestyle='None',label='Measured',color='C5')
axes[0].set_yscale('linear')
axes[0].set_ylim(-5e-2,155e-2) # This leaves out Core NGADG0003 which has 10x higher flux for some reason
# axes[1].set_ylim(bottom=0.9e-11)
axes[0].legend(labels=['CH$_4$','CO$_2$'])

axes[1].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['Fe_II'][xx],marker='o',linestyle='None',label='Measured',color='C2')
axes[2].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['pH'][xx],marker='o',linestyle='None',label='Measured',color='C0')

# axes[0].set_xlim(right=71)
axes[0].set_xlim(right=71)
axes[1].set_xlim(right=71)
axes[2].set_xlim(right=71)


# Aerobic
fig,axes=subplots(3,1,num='Mineral horizon Oxic',figsize=(6,8.4),clear=True)
plot_result(result_highO2_mineral,Fe_ax=axes[1],gasflux_ax=axes[0],pH_ax=axes[2],gdrywt=True,SOC_pct=SOC_layermean['Mineral'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2['Mineral'])
# axes[0].set_ylim(bottom=-0.01,top=30)
axes[1].legend(loc='right',labels=['Iron oxide (solid)','Fe$^{3\!\!+}$ (dissolved)','Fe$^{2\!\!+}$ (dissolved)'])

xx=(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')&(Barrow_synthesis_data['Incubation_Temperature']>4)
CO2mean=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/Barrow_synthesis_data['Incubation_Time'].diff() 
CO2std=Barrow_synthesis_data[['CO2_1','CO2_2','CO2_3']].astype(float).diff().std(axis=1) /Barrow_synthesis_data['Incubation_Time'].diff()

axes[0].errorbar(Barrow_synthesis_data.Incubation_Time[xx],CO2mean[xx],yerr=CO2std[xx],marker='o',linestyle='None',label='Measured',color='C5')
axes[0].set_yscale('linear')
axes[0].set_ylim(-5e-2,155e-2) # This leaves out Core NGADG0003 which has 10x higher flux for some reason
# axes[1].set_ylim(bottom=0.9e-11)
axes[0].legend(labels=['CH$_4$','CO$_2$'])

axes[1].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['Fe_II'][xx],marker='o',linestyle='None',label='Measured',color='C2')
axes[2].plot(Barrow_synthesis_data.Incubation_Time[xx],Barrow_synthesis_data['pH'][xx],marker='o',linestyle='None',label='Measured',color='C0')

# axes[0].set_xlim(right=71)
axes[0].set_xlim(right=71)
axes[1].set_xlim(right=71)
axes[2].set_xlim(right=71)



# Periodic inundation
oxicfrac=0.1
O2_periodic=numpy.zeros(simlength*24)
O2_periodic[:int(simlength*24*oxicfrac)]=dq
result_periodicO2,output_units=run_alquimia.run_simulation('Arctic_redox.in',simlength*nperiodic,timestep,initcond=pools_atmoO2_organic,bc=pools_atmoO2_organic,diffquo={'O2(aq)':O2_periodic},hands_off=False,rateconstants=rateconstants,truncate_concentration=truncate_conc)


fig,axes=subplots(4,1,num='Periodic inundation',figsize=(6,8.4),clear=True)
plot_result(result_periodicO2,Fe_ax=axes[1],gasflux_ax=axes[0],pH_ax=axes[2],porewater_ax=axes[3],gdrywt=True,SOC_pct=SOC_layermean['Organic'],cellulose_SOC_frac=cellulosefrac,BD=BD_layerest2['Organic'])
# axes[0].set_ylim(bottom=-0.01)
axes[1].legend(labels=['Fe(OH)$_3$','Fe$^{3\!\!+}$','Fe$^{2\!\!+}$'],loc=(0.83,0.2),fontsize='small')

axes[0].set_ylim(bottom=0.9e-11)
axes[0].legend(labels=['CH$_4$','CO$_2$'])
axes[3].legend(loc='lower left')


xmax=axes[0].get_xlim()[1]
for ax in axes:
    for num in range(3):
        ax.axvspan(num*simlength+nonzero(diff(O2_periodic))[0]/24,(num+1)*simlength,color='b',alpha=0.1)
    ax.set_xlim(right=xmax)

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