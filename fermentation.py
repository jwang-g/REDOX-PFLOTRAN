import decomp_network
import plot_pf_output

pflotran_exe='../pflotran-interface/src/pflotran/pflotran'

simlength=365

pools = [
decomp_network.decomp_pool(name='cellulose',CN=50,constraints={'initial':1e2},kind='immobile'),
decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':1e-3},kind='primary'),
decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'4.0 P'}),
decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':1e-12}),
decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-3 M Fe(OH)3'}),
decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-3'}),
decomp_network.decomp_pool(name='NH4+',kind='primary',constraints={'initial':1e-15}), # SOMDecomp sandbox requires this
decomp_network.decomp_pool(name='Tracer',kind='primary',constraints={'initial':1e-15}), # Just to accumulate CO2 loss
decomp_network.decomp_pool(name='CH4(aq)',kind='primary',constraints={'initial':1e-15}),
decomp_network.decomp_pool(name='Acetate-',kind='primary',constraints={'initial':1e-15}),

decomp_network.decomp_pool(name='CO2(g)',kind='gas'),
decomp_network.decomp_pool(name='O2(g)',kind='gas'),

decomp_network.decomp_pool(name='CO2(aq)',kind='secondary'),
decomp_network.decomp_pool(name='OH-',kind='secondary'),
decomp_network.decomp_pool(name='FeCO3+',kind='secondary'),
decomp_network.decomp_pool(name='Fe(OH)4-',kind='secondary'),
decomp_network.decomp_pool(name='Acetic_acid(aq)',kind='secondary'),
decomp_network.decomp_pool(name='FeCH3COO+',kind='secondary'),
# decomp_network.decomp_pool(name='FeCO3(aq)',kind='secondary'),

decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-1  1.d2 m^2/m^3','bc':'1.75d-2  1. m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Goethite',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d1 m^2/m^3','bc':'1.75d-2  1. m^2/m^3'},kind='mineral'),
# decomp_network.decomp_pool(name='Fe',rate='1.d-7 mol/m^2-sec',constraints={'initial':'1.0e-6  1. m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.e-10 m^2/m^3'},kind='mineral'),
decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=2000.0,complexes=['>Carboxylic_acid']),
]

pools_lowFe=pools.copy()
for n,p in enumerate(pools_lowFe):
    if p['name']=='Fe(OH)3':
        pools_lowFe[n]=pools_lowFe[n].copy()
        pools_lowFe[n].update(constraints={'initial':'0.0  1.e2 m^2/m^3'})
        
pools_highO2=pools.copy()
for n,p in enumerate(pools_highO2):
    if p['name']=='O2(aq)':
        pools_highO2[n]=pools_highO2[n].copy()
        pools_highO2[n].update(constraints={'initial':'1.0e1'})
        
    # if p['name']=='Fe+++':
    #     pools_lowFe[n]=pools_lowFe[n].copy()
    #     pools_lowFe[n].update(constraints={'initial':'1.d-20'})        
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
# Major issue so far is balancing pH: Fe(III) release from mineral sucks up a lot of H+, and I'm not sure whether it's being replenished from
# some reaction involving Fe++

# Beth: There should be a lot of solid OM functional groups, probably carboxylate/carboxylic acid that exchange protons and buffer pH in these soils
# Fe(II) is probably sorbing onto OM, or staying dissolved
# Fe(III) should precipitate in this system, so if it's dissolved it would be complexed with something
# Per Persson has papers on Fe partitioning in soils as a function of pH
# Derek Lovely has done a lot of Fe reduction lab work to look at reduction rates
# Kinetics of water-rock interaction (textbook) has lots of rate constants. Look for chapter by Eric Rodin
# Aaron Thompson does work on Fe oxides and redox fluctuations

hydrolysis = decomp_network.reaction(name='Hydrolysis',reactant_pools={'cellulose':1.0},product_pools={'DOM1':1.0},
                                        rate_constant=1e-1,rate_units='y', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                    inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=1e-5)])
                                    
# Calculating these as per unit carbon, so dividing by the 6 carbons in a glucose
# C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2
fermentation = decomp_network.reaction(name='fermentation',reactant_pools={'DOM1':6/6},product_pools={'Acetate-':2/6,'HCO3-':2/6,'H+':4/6+4*2/6,'Tracer':2/6}, # balancing pH of FeIII release requires an extra 5.5 H+ to be released here
                                        rate_constant=1e-6,reactiontype='MICROBIAL', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                    inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD'),decomp_network.inhibition(species='Acetate-',k=6.25e-8,type='MONOD')],
                                    monod_terms=[decomp_network.monod(species='DOM1',k=1e-5,threshold=1e-15)])

# CH2O + H2O -> CO2 + 4H+ + 4 e-
# O2   + 4H+ + 4 e- -> 2H2O
DOM_resp = decomp_network.reaction(name='DOM aerobic respiration',reactant_pools={'DOM1':1.0,'O2(aq)':1.0},product_pools={'HCO3-':1.0,'H+':1.0,'Tracer':1.0},
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12),decomp_network.monod(species='DOM1',k=1e-8,threshold=1e-14)],
                                    rate_constant=1.0e-5,reactiontype='MICROBIAL')

# C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
# 2 O2    + 8 H+ + 8 e- -> 4 H2O
acetate_resp = decomp_network.reaction(name='Acetate aerobic respiration',reactant_pools={'Acetate-':1.0,'O2(aq)':2.0},product_pools={'HCO3-':2.0,'H+':2.0,'Tracer':2.0},
                                        monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12),decomp_network.monod(species='Acetate-',k=1e-8,threshold=1e-14)],
                                    rate_constant=1.0e-5,reactiontype='MICROBIAL')


# C2H3O2- + 2 H2O -> 2 CO2 + 7 H+ + 8 e-
# 8 Fe+++ + 8 e- -> 8 Fe++ 
Fe_reduction = decomp_network.reaction(name='Fe reduction',reactant_pools={'Acetate-':1.0,'Fe+++':8.0},product_pools={'HCO3-':2.0,'Fe++':8.0,'H+':9.0,'Tracer':2.0},
                                        monod_terms=[decomp_network.monod(species='Acetate-',k=2e-3,threshold=1e-15),decomp_network.monod(species='Fe+++',k=1.3e-12,threshold=1e-15)],
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD')],
                                        rate_constant=2e-7,reactiontype='MICROBIAL')
                                        
# C2H3O2- + H+ -> CH4 + HCO3- + H+
Methanogenesis = decomp_network.reaction(name='Methanogenesis',reactant_pools={'Acetate-':1.0},product_pools={'CH4(aq)':1.0,'HCO3-':1.0, 'Tracer':1.0},
                                        monod_terms=[decomp_network.monod(species='Acetate-',k=2e-3,threshold=1e-15),
                                                     ],
                                        inhibition_terms=[decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD'),decomp_network.inhibition(species='Fe+++',k=6.25e-9,type='MONOD')],
                                        rate_constant=1e-8,reactiontype='MICROBIAL')

fermentation_network =  decomp_network.decomp_network(pools=pools_lowFe,reactions=[hydrolysis,fermentation,DOM_resp,acetate_resp,Fe_reduction])
fermentation_network_Fe =  decomp_network.decomp_network(pools=pools,reactions=[hydrolysis,fermentation,DOM_resp,acetate_resp,Fe_reduction])
fermentation_network_Fe_CH4 =  decomp_network.decomp_network(pools=pools_lowFe,reactions=[hydrolysis,fermentation,DOM_resp,acetate_resp,Fe_reduction,Methanogenesis])
fermentation_network_highO2 =  decomp_network.decomp_network(pools=pools_highO2,reactions=[hydrolysis,fermentation,DOM_resp,acetate_resp,Fe_reduction,Methanogenesis])

result,units=decomp_network.PF_network_writer(fermentation_network).run_simulation('SOMdecomp_template.txt','fermentation',pflotran_exe,print_output=False,length_days=simlength)
result_Fe,units_Fe=decomp_network.PF_network_writer(fermentation_network_Fe).run_simulation('SOMdecomp_template.txt','fermentation',pflotran_exe,print_output=False,length_days=simlength)
result_Fe_CH4,units_Fe_CH4=decomp_network.PF_network_writer(fermentation_network_Fe_CH4).run_simulation('SOMdecomp_template.txt','fermentation',pflotran_exe,print_output=False,length_days=simlength)
result_highO2,units_highO2=decomp_network.PF_network_writer(fermentation_network_highO2).run_simulation('SOMdecomp_template.txt','fermentation',pflotran_exe,print_output=False,length_days=simlength)

from pylab import *
result,units=plot_pf_output.convert_units(result,units,'M')
result_Fe,units_Fe=plot_pf_output.convert_units(result_Fe,units_Fe,'M')
result_Fe_CH4,units_Fe_CH4=plot_pf_output.convert_units(result_Fe_CH4,units_Fe_CH4,'M')
result_highO2,units_highO2=plot_pf_output.convert_units(result_highO2,units_highO2,'M')
figure('Network diagram',figsize=(11.8,4.8));clf()
ax=subplot(223)
to_draw=fermentation_network_Fe_CH4.copy()
for p in pools:
    if p['kind'] in ['secondary','surf_complex'] or p['name'] in ['HRimm','Tracer']:
        to_draw.remove_node(p['name'])
    elif p['kind'] in ['mineral','gas']:
        to_draw=decomp_network.nx.compose(decomp_network.get_reaction_from_database('hanford.dat',p),to_draw)
pos=decomp_network.nx.drawing.nx_agraph.graphviz_layout(to_draw,prog='dot')
decomp_network.nx.draw_networkx(to_draw,pos=pos,with_labels=True,ax=ax,nodes=to_draw.nodes,node_color=decomp_network.make_nodecolors(to_draw.nodes),arrowsize=15,font_size='small',arrowstyle='->')
title('Decomposition network diagram (without complexes)')

ax=subplot(224)
to_draw=fermentation_network_Fe_CH4.copy()
for p in pools:
    if p['name'] in ['HRimm','Tracer']:
        to_draw.remove_node(p['name'])
    if p['kind'] in ['mineral','gas','secondary','surf_complex']:
        to_draw=decomp_network.nx.compose(decomp_network.get_reaction_from_database('hanford.dat',p),to_draw)
pos=decomp_network.nx.drawing.nx_agraph.graphviz_layout(to_draw,prog='dot')
decomp_network.nx.draw_networkx(to_draw,pos=pos,with_labels=True,ax=ax,nodes=to_draw.nodes,node_color=decomp_network.make_nodecolors(to_draw.nodes),arrowsize=15,font_size='small',arrowstyle='->')
title('Decomposition network diagram (with complexes)')

ax=subplot(221)
# ax.set_yscale('log')
for pool in ['cellulose','Free DOM1','Free Acetate-','Total CH4(aq)','Total Fe+++','Total Fe++']:
    l=plot(result[pool],label=pool)[0]
    plot(result_Fe[pool],ls='--',c=l.get_color())
    plot(result_Fe_CH4[pool],ls=':',c=l.get_color())
    plot(result_highO2[pool],ls='-.',c=l.get_color())
l=plot(result['Total Tracer']+result['HRimm'],label='CO2 produced')[0]
plot(result_Fe['Total Tracer']+result_Fe['HRimm'],c=l.get_color(),ls='--')
plot(result_Fe_CH4['Total Tracer']+result_Fe_CH4['HRimm'],c=l.get_color(),ls=':')
plot(result_highO2['Total Tracer']+result_highO2['HRimm'],c=l.get_color(),ls='-.')
l=plot(-log10(result['Free H+']),label='pH')[0]
plot(-log10(result_Fe['Free H+']),c=l.get_color(),ls='--')
plot(-log10(result_Fe_CH4['Free H+']),c=l.get_color(),ls=':')
plot(-log10(result_highO2['Free H+']),c=l.get_color(),ls='-.')
legend(handles=[Line2D([0],[0],ls='-',color='k',label='No Fe(III)'),Line2D([0],[0],color='k',ls='--',label='Fe(III)'),Line2D([0],[0],color='k',ls=':',label='Fe(III) and CH4')],
        fontsize='small',loc='upper right')

title('Concentrations')
ylabel('Concentration (M)')
xlabel('Time (days)')

ax=subplot(222)
ax.set_yscale('log')
for pool in ['cellulose','Free DOM1','Free Acetate-','Total CH4(aq)','Total Fe+++','Total Fe++']:
    l=plot(result[pool],label=pool)[0]
    plot(result_Fe[pool],ls='--',c=l.get_color())
    plot(result_Fe_CH4[pool],ls=':',c=l.get_color())
    plot(result_highO2[pool],ls='-.',c=l.get_color())
l=plot(result['Total Tracer']+result['HRimm'],label='CO2 produced')[0]
plot(result_Fe['Total Tracer']+result_Fe['HRimm'],c=l.get_color(),ls='--')
plot(result_Fe_CH4['Total Tracer']+result_Fe_CH4['HRimm'],c=l.get_color(),ls=':')
plot(result_highO2['Total Tracer']+result_highO2['HRimm'],c=l.get_color(),ls='-.')
l=plot(result['Free H+'],label='H+')[0]
plot(result_Fe['Free H+'],c=l.get_color(),ls='--')
plot(result_Fe_CH4['Free H+'],c=l.get_color(),ls=':')
plot(result_highO2['Free H+'],c=l.get_color(),ls='-.')

title('Concentrations (log scale)')
ylabel('Concentration (M)')
xlabel('Time (days)')
legend(fontsize='small',ncol=2,loc='lower right')

tight_layout()


# CTC network with fermentation
simlength=365*10

# CTC decomposition network
decomp_network_CTC=decomp_network.decomp_network()

decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL1',CN=  12.,constraints={'initial':1e-10},kind='immobile') )
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL2',CN=  12.,constraints={'initial':1e-10},kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL3',CN=  10.,constraints={'initial':1e-10},kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='SOIL4',CN=  10.,constraints={'initial':1e-10},kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR1',constraints={'initial':1e2},initCN=20,kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR2',constraints={'initial':1e-10},initCN=20,kind='immobile')    )
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='LITR3',constraints={'initial':1e-10},initCN=20,kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='CWD',constraints={'initial':1e-10},initCN=20,kind='immobile'))
decomp_network_CTC.add_pool(decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-10},kind='immobile'))

for pool in pools[1:]:
    decomp_network_CTC.add_pool(pool)

# CWD decomposition to  litter
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'CWD':1.0},product_pools={'LITR2':0.76,'LITR3':0.24},
                            rate_constant=0.00010,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',
                            monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)],
                            name='CWD fragmentation'))

# Litter decomposition
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR1':1.0},product_pools={'SOIL1':0.61,'HCO3-':1-0.61},
                rate_constant=1.204,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR1 decomposition',
                monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR2':1.0},product_pools={'SOIL2':0.45,'HCO3-':1-0.45},
                rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR2 decomposition',
                monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'LITR3':1.0},product_pools={'SOIL3':0.71,'HCO3-':1-0.71},
                rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='LITR3 decomposition',
                monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))

# SOM decomposition
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL1':1.0},product_pools={'SOIL2':0.72,'HCO3-':1-0.72},
            rate_constant=0.0726,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL1 decomp',
            monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL2':1.0},product_pools={'SOIL3':0.54,'HCO3-':1-0.54},
            rate_constant=0.0141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL2 decomp',
            monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL3':1.0},product_pools={'SOIL4':0.45,'HCO3-':1-0.45},
            rate_constant=0.00141,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL3 decomp',
            monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))
decomp_network_CTC.add_reaction(decomp_network.reaction(reactant_pools={'SOIL4':1.0},product_pools={'HCO3-':1.0},
            rate_constant=0.0001,rate_units='1/d',turnover_name='RATE_DECOMPOSITION',name='SOIL4 decomp',
            monod_terms=[decomp_network.monod(species='O2(aq)',k=1e-8,threshold=1e-12)]))

decomp_network_CTC.add_reaction(DOM_resp)
decomp_network_CTC.add_reaction(Fe_reduction)
decomp_network_CTC.add_reaction(acetate_resp)
decomp_network_CTC.add_reaction(fermentation)
decomp_network_CTC.add_reaction( decomp_network.reaction(name='hydrolysis',reactant_pools={'LITR1':1.0},product_pools={'DOM1':1.0},
                                        rate_constant=1e-1,rate_units='y', #  Jianqiu Zheng et al., 2019: One third of fermented C is converted to CO2
                                    inhibition_terms=[decomp_network.inhibition(species='DOM1',type='MONOD',k=1e-5),decomp_network.inhibition(species='O2(aq)',k=6.25e-8,type='MONOD')]))


CTC_result,CTC_units=decomp_network.PF_network_writer(decomp_network_CTC).run_simulation('SOMdecomp_template.txt','CTC',pflotran_exe,length_days=simlength)

# Run with low Fe mineral concentration to cut off Fe reduction pathway
decomp_network_CTC_lowFe=decomp_network_CTC.copy()
decomp_network_CTC_lowFe.nodes['Fe(OH)3']['constraints']={'initial':'0.0d-5  1. m^2/m^3','bc':'1.0d-5  1. m^2/m^3'}
CTC_result_lowFe,CTC_units_lowFe=decomp_network.PF_network_writer(decomp_network_CTC_lowFe).run_simulation('SOMdecomp_template.txt','CTC',pflotran_exe,length_days=simlength)

# Run with abundant oxygen so aerobic CTC reactions will proceed
decomp_network_CTC_highO2=decomp_network_CTC.copy()
decomp_network_CTC_highO2.nodes['O2(aq)']['constraints']={'initial':1.0e1,'bc':'0.2 G O2(g)'}
CTC_result_highO2,CTC_units_highO2=decomp_network.PF_network_writer(decomp_network_CTC_highO2).run_simulation('SOMdecomp_template.txt','CTC',pflotran_exe,length_days=simlength)

figure('CTC network');clf()
subplot(122)
to_draw=decomp_network_CTC.copy()
for p in pools:
    if p['kind'] in ['secondary'] or p['name'] in ['Tracer','HRimm']:
        to_draw.remove_node(p['name'])
    elif p['kind'] in ['mineral','gas','surf_complex']:
        to_draw=decomp_network.nx.compose(decomp_network.get_reaction_from_database('hanford.dat',p),to_draw)
pos=decomp_network.nx.drawing.nx_agraph.graphviz_layout(to_draw,prog='dot')
decomp_network.nx.draw_networkx(to_draw,pos=pos,with_labels=True,nodes=to_draw.nodes,node_color=decomp_network.make_nodecolors(to_draw.nodes),arrowsize=15,font_size='small',arrowstyle='->')

CTC_result,CTC_units=plot_pf_output.convert_units(CTC_result,CTC_units,'M')
CTC_result_highO2,CTC_units_highO2=plot_pf_output.convert_units(CTC_result_highO2,CTC_units_highO2,'M')
CTC_result_lowFe,CTC_units_lowFe=plot_pf_output.convert_units(CTC_result_lowFe,CTC_units_lowFe,'M')
subplot(221)
handles=[]
for pool in ['LITR1C','SOIL1','SOIL2','SOIL3','SOIL4']:
    l=plot(CTC_result_highO2[pool],label=pool)[0]
    plot(CTC_result[pool],ls='--',c=l.get_color())
    plot(CTC_result_lowFe[pool],ls=':',c=l.get_color())
    handles.append(l)
l=plot(CTC_result_highO2['HRimm']+CTC_result_highO2['Total Tracer'],label='CO2 prod')[0]
plot(CTC_result['HRimm']+CTC_result['Total Tracer'],c=l.get_color(),ls='--')
plot(CTC_result_lowFe['HRimm']+CTC_result_lowFe['Total Tracer'],c=l.get_color(),ls=':')
handles.append(l)
handles.append(Line2D([0],[0],c='k',ls='-',label='High O2'))
handles.append(Line2D([0],[0],c='k',ls='--',label='Low O2, high Fe'))
handles.append(Line2D([0],[0],c='k',ls=':',label='Low O2, low Fe'))
legend(handles=handles,fontsize='small',ncol=2,loc='center right')
ylabel('Concentration (M)')
xlabel('Time (days)')

ax=subplot(223)
ax.set_yscale('log')
handles=[]
for pool in ['Free DOM1','Free Acetate-','Free Fe+++','Free H+','Free O2(aq)']:
    l=plot(CTC_result_highO2[pool],label=pool)[0]
    plot(CTC_result[pool],ls='--',c=l.get_color())
    plot(CTC_result_lowFe[pool],ls=':',c=l.get_color())
    handles.append(l)
# handles.append(Line2D([0],[0],c='k',ls='-',label='High O2'))
# handles.append(Line2D([0],[0],c='k',ls='--',label='Low O2, high Fe'))
# handles.append(Line2D([0],[0],c='k',ls=':',label='Low O2, low Fe'))
legend(handles=handles,fontsize='small',loc='lower right',ncol=1)
ylabel('Concentration (M)')
xlabel('Time (days)')


tight_layout()

show()
