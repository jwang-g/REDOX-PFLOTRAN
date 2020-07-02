from run_alquimia import get_alquimiavector,ffi,lib,check_status,init_alquimia,convert_condition_to_alquimia,print_alquimia_object,convert_output,run_onestep
import decomp_network
import Manganese_network as Mn
from matplotlib import pyplot
import numpy


class layer:
    def __init__(self,volume,saturation=1.0,temperature=20.0,water_density=1000.0,porosity=0.25,pressure=101325.0,rateconstants={},diffquo={}):
        self.volume=volume # m3
        self.saturation=saturation
        self.temperature=temperature
        self.water_density=water_density
        self.porosity=porosity
        self.pressure=pressure
        self.mineral_specific_surface_area={}
        self.mineral_volume_fraction={}
        self.surface_site_density={}
        self.total_immobile={}
        self.total_mobile={}
        self.free_mobile={}
        self.rateconstants=rateconstants
        self.mineral_rate_cnst={}
        self.aux_ints=[]
        self.aux_doubles=[]
        self.mineral_reaction_rate={}
        self.secondary_free_ion_concentration={}
        self.diffquo=diffquo
        
        # Some things that I'm not using are skipped
        
    def copy_from_alquimia(self,data):
        self.volume=data.properties.volume
        self.saturation=data.properties.saturation
        self.temperature=data.state.temperature
        self.water_density=data.state.water_density
        self.porosity=data.state.porosity
        self.pressure=data.state.aqueous_pressure
        
        self.mineral_names=get_alquimiavector(data.meta_data.mineral_names)
        for num in range(data.meta_data.mineral_names.size):
            self.mineral_rate_cnst[self.mineral_names[num]]=data.properties.mineral_rate_cnst.data[num]
            self.mineral_specific_surface_area[self.mineral_names[num]]=data.state.mineral_specific_surface_area.data[num]
            self.mineral_volume_fraction[self.mineral_names[num]]=data.state.mineral_volume_fraction.data[num]
            self.mineral_reaction_rate[self.mineral_names[num]]=data.aux_output.mineral_reaction_rate.data[num]
            
        self.primary_names=get_alquimiavector(data.meta_data.primary_names)
        for num in range(len(self.primary_names)):
            self.total_immobile[self.primary_names[num]]=data.state.total_immobile.data[num]
            self.total_mobile[self.primary_names[num]]=data.state.total_mobile.data[num]
            self.free_mobile[self.primary_names[num]]=data.aux_output.primary_free_ion_concentration.data[num]
            
        for num in range(len(self.secondary_names)):
            self.secondary_free_ion_concentration[self.secondary_names[num]]=data.aux_output.secondary_free_ion_concentration.data[num]
            
        self.surface_site_names=get_alquimiavector(data.meta_data.surface_site_names)
        for num in range(len(self.surface_site_names)):
            self.surface_site_density[self.surface_site_names[num]]=data.state.surface_site_density.data[num]
            
        self.aux_ints=get_alquimiavector(data.aux_data.aux_ints)
        self.aux_doubles=get_alquimiavector(data.aux_data.aux_doubles)
        
        
    def copy_to_alquimia(self,data):
        data.properties.volume=self.volume
        data.properties.saturation=self.saturation
        data.state.temperature=self.temperature
        data.state.water_density=self.water_density
        data.state.porosity=self.porosity
        data.state.aqueous_pressure=self.pressure
        
        
        for num in range(len(self.mineral_names)):
            data.properties.mineral_rate_cnst.data[num]=self.mineral_rate_cnst[self.mineral_names[num]]
            data.state.mineral_specific_surface_area.data[num]=self.mineral_specific_surface_area[self.mineral_names[num]]
            data.state.mineral_volume_fraction.data[num]=self.mineral_volume_fraction[self.mineral_names[num]]
            
        for num in range(len(self.surface_site_names)):
            data.state.surface_site_density.data[num]=self.surface_site_density[self.surface_site_names[num]]
            
        for num in range(len(self.primary_names)):
            data.state.total_immobile.data[num]=self.total_immobile[self.primary_names[num]]
            data.state.total_mobile.data[num]=self.total_mobile[self.primary_names[num]]
            
        for num in range(len(self.aux_ints)):
            data.aux_data.aux_ints.data[num]=self.aux_ints[num]
        for num in range(len(self.aux_doubles)):
            data.aux_data.aux_doubles.data[num]=self.aux_doubles[num]
        

    def setup_output(self,nsteps,dt):
        self.output={
            'total_mobile':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'free':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'immobile':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'aq_complex':numpy.ma.masked_all((nsteps,len(self.secondary_names)),dtype=float),
            'mineral_VF':numpy.ma.masked_all((nsteps,len(self.mineral_rate_cnst)),dtype=float),
            'mineral_rate':numpy.ma.masked_all((nsteps,len(self.mineral_rate_cnst)),dtype=float),
            'time':numpy.arange(nsteps,dtype=float)*dt,
            'actual_dt':numpy.ma.masked_all(nsteps,dtype=float),'ncuts':numpy.ma.masked_all(nsteps,dtype=int),
            'porosity':numpy.ma.masked_all(nsteps,dtype=float)
        }
        
    def write_output(self,data,step,dt,num_cuts=0):
        self.output['total_mobile'][step,:]=numpy.array([self.total_mobile[name] for name in self.primary_names])
        self.output['free'][step,:]=numpy.array([self.free_mobile[name] for name in self.primary_names])
        self.output['aq_complex'][step,:]=numpy.array([self.secondary_free_ion_concentration[name] for name in self.secondary_names])
        self.output['immobile'][step,:]=numpy.array([self.total_immobile[name] for name in self.primary_names])
        self.output['mineral_VF'][step,:]=numpy.array([self.mineral_volume_fraction[name] for name in self.mineral_names])
        self.output['time'][step]=step*dt
        self.output['mineral_rate'][step,:]=numpy.array([self.mineral_reaction_rate[name] for name in self.mineral_names])
        self.output['porosity'][step]=self.porosity
        self.output['actual_dt'][step]=dt/2**num_cuts
        self.output['ncuts'][step]=num_cuts
        
    def convert_output(self):
        import pandas
        output_DF=pandas.DataFrame(index=self.output['time'])
        output_DF['Porosity']=self.output['porosity']
        output_DF['ncuts']=self.output['ncuts']
        output_DF['actual_dt']=self.output['actual_dt']
        output_units={'time':'s','Porosity':'NA','ncuts':'NA','actual_dt':'s'}
        total=pandas.DataFrame(self.output['total_mobile'],columns=self.primary_names,index=self.output['time']).add_prefix('Total ')
        for col in total.columns:
            output_DF[col]=total[col]
        output_units.update([('Total '+s,'M') for s in self.primary_names])
        free=pandas.DataFrame(self.output['free'],columns=self.primary_names,index=self.output['time']).add_prefix('Free ')
        for col in free.columns:
            output_DF[col]=free[col]
        output_units.update([('Free '+s,'M') for s in self.primary_names])
        sorbed=pandas.DataFrame(self.output['immobile'],columns=self.primary_names,index=self.output['time']).add_prefix('Total Sorbed ')
        for col in sorbed.columns:
            output_DF[col]=sorbed[col]
        output_units.update([('Total Sorbed '+s,'mol/m^3') for s in self.primary_names])
        mineral_VF=pandas.DataFrame(self.output['mineral_VF'],columns=self.mineral_names,index=self.output['time']).add_suffix(' VF')
        for col in mineral_VF.columns:
            output_DF[col]=mineral_VF[col]
        output_units.update([(s+' VF','m^3 mnrl/m^3 bulk') for s in self.mineral_names])
        mineral_rate=pandas.DataFrame(self.output['mineral_rate'],columns=self.mineral_names,index=self.output['time']).add_suffix(' Rate')
        for col in mineral_rate.columns:
            output_DF[col]=mineral_rate[col]
        output_units.update([(s+' Rate','mol/m^3/sec') for s in self.mineral_names])
        secondary=pandas.DataFrame(self.output['aq_complex'],columns=secondary_names,index=self.output['time'])
        for col in secondary.columns:
            output_DF[col]=secondary[col]
        output_units.update([(s,'M') for s in self.secondary_names])
        
        self.output_DF=output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600))
        self.output_units=output_units
        

reaction_network=Mn.make_network(leaf_Mn_mgkg=0.0,Mn2_scale=0.25e-2,Mn_peroxidase_Mn3_leakage=2e-4,Mn3_scale=1e-12) # We will add the Mn along with leaf litter manually instead of generating it through decomposition

rateconstants={
'1.0e+00 DOM1  + 1.0e+00 O2(aq)  -> 1.0e+00 HCO3-  + 1.0e+00 H+  + 1.0e+00 Tracer  + 0.0e+00 Mn++':1e-6,
 '1.0e+00 DOM2  + 2.0e-04 Mn++  + 2.0e-04 H+  -> 1.0e+00 DOM1  + 2.0e-04 Mn+++':0.5e-6,
 'Cellulose decay to DOM1 (SOMDEC sandbox)':2.0/(365*24*3600),
 'Lignin decay to DOM2 (SOMDEC sandbox)':1.0/(365*24*3600),
 '1.0e+00 Mn++  -> 1.0e+00 Tracer2':1.0e-8/100**3, # This is root uptake
 '1.0e+00 DOM1  + 4.0e+00 Mn+++  -> 1.0e+00 HCO3-  + 4.0e+00 Mn++  + 5.0e+00 H+':1e-7 # Manganese reduction
}

input_file='manganese.in'

decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=False,CO2name='Tracer',truncate_concentration=1e-25)


chem,data,sizes,status=init_alquimia(input_file,hands_off=False)



# Read secondary complex names from input file since Alquimia does not provide them
with open(input_file,'r') as infile:
    secondary_names=[]
    for line in infile:
        if 'SECONDARY_SPECIES' in line.split('#')[0]:
            break
    for line in infile:
        l=line.strip().split('#')[0]
        if l.startswith('END') or l.startswith('/'):
            break
        if len(l)>0:
            secondary_names.append(l)

# Set up layers
# layers=[layer(0.1,rateconstants=rateconstants),layer(0.1,rateconstants=rateconstants),layer(0.3,rateconstants=rateconstants)]
layers=[layer(0.1,rateconstants=rateconstants) for num in range(5)]

for l in layers:
    l.secondary_names=secondary_names

# Herndon et al 2014 (BGC): Deep bulk soil Mn concentration ~1500 ug/g = 54 umol/cm3 assuming bulk density=2 g/cm3
Mn_VF=1500e-6/54.94*2 *22.36
initcond=decomp_network.change_constraints(Mn.pools,{'Cellulose':1.0e2,'Lignin':1.0e-8,'H+':'6.0 P',
                                                    # 'Manganite':'1.0d-7  1.d2 m^2/m^3','Mn(OH)2(am)':'%1.2g 1.d2 m^2/m^3'%Mn_VF,
                                                    'Birnessite2':'%1.2g 1.d2 m^2/m^3'%Mn_VF,
                                                    'O2(aq)':'1.0 G O2(g)'})
initcond=decomp_network.change_site_density(initcond, '>DOM1', 1e4)
bc=initcond
dt=3600*12
nsteps=365*24//(dt//3600)*10

# Set up initial condition
for l in layers:
    # Initialize state data
    data.properties.volume=l.volume
    data.properties.saturation=l.saturation

    data.state.temperature=l.temperature
    data.state.water_density=l.water_density
    data.state.porosity=l.porosity
    data.state.aqueous_pressure=l.pressure


    # Set properties: surface site density and mineral rate constants
    # This is necessary when running in hands-on mode
    if initcond is not None:
        for constraint in initcond:
            if constraint['kind']=='surf_complex':
                sitename=constraint['name']
                sitenum=get_alquimiavector(data.meta_data.surface_site_names).index(sitename)
                print('Applying site density: {name:s} (position={pos:d}): {dens:1.1f}'.format(name=sitename,pos=sitenum,dens=constraint['site_density']))
                data.state.surface_site_density.data[sitenum]=constraint['site_density']
            elif constraint['kind']=='mineral':
                name=constraint['name']
                num=get_alquimiavector(data.meta_data.mineral_names).index(name)
                data.properties.mineral_rate_cnst.data[num]=float(constraint['rate'].split()[0].replace('d','e'))
                
    # Set up boundary condition if applicable
    if bc is not None:
        bc_cond=convert_condition_to_alquimia(bc,'initial')
        bc_state=ffi.new('AlquimiaState *')
        bc_state.temperature=layers[0].temperature
        bc_state.water_density=layers[0].water_density
        bc_state.porosity=layers[0].porosity
        bc_state.aqueous_pressure=layers[0].pressure
        bc_auxdata=ffi.new('AlquimiaAuxiliaryData *')
        lib.AllocateAlquimiaState(sizes,bc_state)
        for num in range(data.state.surface_site_density.size):
            bc_state.surface_site_density.data[num]=data.state.surface_site_density.data[num]
        lib.AllocateAlquimiaAuxiliaryData(sizes,bc_auxdata)
        chem.ProcessCondition(ffi.new('void **',data.engine_state),bc_cond,ffi.addressof(data.properties),bc_state,bc_auxdata,status)
        check_status(status,False)
    else:
        bc_state=None
                
    # Aqueous kinetic rate constants also need to be specified in hands-off mode
    # Alquimia interface always sets backward rate to zero
    for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
        data.properties.aqueous_kinetic_rate_cnst.data[num]=l.rateconstants[reactname]

    init_cond=convert_condition_to_alquimia(initcond,'initial')
    chem.ProcessCondition(ffi.new('void **',data.engine_state),init_cond,ffi.addressof(data.properties),ffi.addressof(data.state),ffi.addressof(data.aux_data),status)
    check_status(status,False)

    # Pflotran sets porosity based on minerals or something? Needs to be reset
    data.state.porosity=l.porosity
    
    l.copy_from_alquimia(data)
    l.setup_output(nsteps+1,dt)

    ##### At this point, the model should be initialized ##########
    print('''

    *****************************************************
    Successfully initialized alquimia geochemical engine
    *****************************************************

    ''')


flow_rate=numpy.zeros(len(layers)-1)+1e-7 # cm/s = 10 L/m2/s, positive is downward
flow_rate=numpy.linspace(1e-7,0,len(layers)) # Rate declinine linearly with depth, assumes removal or accumulation in lower layers
min_dt=0.1
truncate_concentration=1e-20


litter_mass=1.0 # kg C/m2

Mn_molarmass=54.94        #g/mol
C_molarmass=12.01         #g/mol
litter_ligninfrac=0.5
litter_Cfrac_mass=0.4       #g/g
litter_Mn_mg_g_initial=2.0  #mg/g


root_efolding=0.2 # e-folding depth of root biomass in m
root_biomass_top=0.3 # gC/cm3
root_Cfrac=0.4

z=numpy.array([0]+[l.volume for l in layers]).cumsum()
z_mid=(z[:-1]+z[1:])/2 
for l in range(len(layers)):
    layers[l].total_immobile['Root_biomass']=root_biomass_top*numpy.exp(-z_mid[l]/root_efolding)/(root_Cfrac*12)*100**3
    layers[l].mineral_rate_cnst['Mn(OH)2(am)']=2e-11
    layers[l].mineral_rate_cnst['Birnessite2']=2e-10
    layers[l].surface_site_density['>DOM1']=1e3
    layers[l].diffquo={'O2(aq)':0.01**(l+1)}

# Treat top layer as O horizon with less root biomass and mineral Mn
layers[0].total_immobile['Root_biomass']=root_biomass_top/(root_Cfrac*12)*100**3*1e-3
# layers[0].mineral_volume_fraction['Mn(OH)2(am)']=1e-7
layers[0].mineral_volume_fraction['Birnessite2']=1e-7
layers[0].surface_site_density['>DOM1']=1e2

for l in layers:
    l.write_output(data, 0,dt)

flow_in=numpy.zeros(len(layers),dtype=float)
flow_out=numpy.zeros(len(layers),dtype=float)
immobile_specs=['DOM2','Tracer2']

import time
t0=time.time()
tprev=t0

success=True
for step in range(nsteps):
    # First deposit leaf litter, including last year worth of root Mn uptake
    if (step*dt)%(365*24*3600) == 0: # End of year. Needs dt to be even divisor of one day to work
        layers[0].total_immobile['Cellulose']=layers[0].total_immobile['Cellulose'] + litter_mass*(1-litter_ligninfrac)/layers[0].volume/C_molarmass*1000
        layers[0].total_immobile['Lignin']=layers[0].total_immobile['Lignin'] + litter_mass*litter_ligninfrac/layers[0].volume/C_molarmass*1000
        
        # Add up all Mn uptake from previous year to add to leaf litter, and reset those pools
        Mn_uptake_total=0.0 # mol Mn
        for l in layers:
            Mn_uptake_total=Mn_uptake_total+l.total_mobile['Tracer2']*l.volume*1000*l.porosity*l.saturation
            l.total_mobile['Tracer2']=0.0

        layers[0].total_mobile['Mn++']=layers[0].total_mobile['Mn++'] + Mn_uptake_total/(layers[0].porosity*layers[0].saturation*layers[0].volume*1000)
        print('Leaf litter Mn concentration = %1.1f mmol/kg'%(Mn_uptake_total*1e3/litter_mass))
    
    # Next advect mobile species. This assumes we can separate advective transport and reactions at this time step length
    for spec in layers[0].primary_names:
        if spec in immobile_specs:
            continue # Skip this one because it's not really mobile, more a workaround for Lignin decomp
        flow_in[:]=0.0
        flow_out[:]=0.0
        for num in range(len(layers)-1):
            # First calculate net flow for each layer
            flow_out[num]=flow_rate[num]*layers[num].total_mobile[spec]*10*dt # mol/L*cm/s*(10L/m2)/cm -> mol/s
            flow_in[num+1]=flow_rate[num]*layers[num].total_mobile[spec]*10*dt
        for num in range(len(layers)):
            # Then calculate change in concentration using amount of water in the layer
            # Units of mobile species are M (mol/L water)
            spec_mol_layer=layers[num].total_mobile[spec]*layers[num].volume*1000*layers[num].porosity*layers[num].saturation # mol of the species in the layer initially
            spec_mol_layer=spec_mol_layer+flow_in[num]-flow_out[num]
            layers[num].total_mobile[spec]=spec_mol_layer/(layers[num].volume*1000*layers[num].porosity*layers[num].saturation)
            
            
    # Next do chemistry.
    try:
        for n,l in enumerate(layers):
            l.copy_to_alquimia(data)
            dq={}
            if bc is not None:
                primarynames=get_alquimiavector(data.meta_data.primary_names)
                for spec in primarynames:
                    if spec in l.diffquo.keys():
                        if numpy.iterable(l.diffquo[spec]):
                            dq[spec]=l.diffquo[spec][step%len(l.diffquo[spec])]
                        else:
                            dq[spec]=l.diffquo[spec]
                    else:
                        dq[spec]=0.0
                        # print('O2 before: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
            num_cuts=run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=rateconstants)
            l.copy_from_alquimia(data)
            # Write output
            l.write_output(data,step+1,dt,num_cuts)
            # print('O2 after: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
    except RuntimeError as err:
        print('ERROR on timestep %d, layer %d: %s'%(step,n,err))
        print('Returning output so far')
        l.write_output(data,step,dt,num_cuts)
        success=False
        break
    except KeyboardInterrupt as err:
        print('INTERRUPTED on timestep %d'%(step))
        print('Returning output so far')
        l.write_output(data,step,dt,num_cuts)
        success=False
        break
    

            
    
    if step%100==0:
        t1=time.time()
        mean_cuts=numpy.mean([l.output['ncuts'][step-100:step].mean() for l in layers])
        mean_dt=numpy.mean([l.output['actual_dt'][step-100:step].mean() for l in layers])
        print('*** Step {step:d} of {nsteps:d}. Time elapsed: {t:d} s ({tperstep:1.1f} s per {steplength:1.1f} hour timestep). Mean cuts: {meancuts:1.1f} Mean dt: {meandt:1.1f} s ***'.format(
                step=step,nsteps=nsteps,t=int(t1-t0),tperstep=(t1-tprev)/25,meancuts=mean_cuts,meandt=mean_dt,steplength=dt/3600))
        tprev=t1


networkfig=pyplot.figure('Reaction network',clear=True)
drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid'],
        font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
        namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'HCO$_3^-$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',})


for l in layers:
    l.convert_output()

molar_volume_manganite = 24.45 # cm3/mol
molar_volume_MnOH2am = 22.3600
molar_volume_birnessite = 251.1700
BD=2.0 # Approximate soil bulk density

f,axs=pyplot.subplots(ncols=3,nrows=len(layers),sharex=True,clear=True,num='Simulation results',figsize=(7,6.5))
t=layers[0].output_DF.index/365
for num in range(len(layers)):
    axs[num,0].plot(t,layers[num].output_DF['Total Sorbed Cellulose']*12/100**3,label='Cellulose')
    axs[num,0].plot(t,layers[num].output_DF['Total Sorbed Lignin']*12/100**3,label='Lignin')
    axs[num,0].plot(t,layers[num].output_DF['Total DOM1']*12/1000*layers[num].porosity*layers[num].saturation,label='DOM')
    axs[num,0].plot(t,layers[num].output_DF['Total Sorbed DOM1']*12/100**3,label='Sorbed DOM')
    
    axs[num,0].set_ylabel('C density\n(g C cm$^{-3}$)')
    
    axs[num,1].plot(t,layers[num].output_DF['Total Mn++']*1e6/1000*layers[num].porosity*layers[num].saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+}$')
    axs[num,1].plot(t,layers[num].output_DF['Total Mn+++']*1e6/1000*layers[num].porosity*layers[num].saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+\!\!+}$')
    # axs[num,1].plot(t,layers[num].output_DF['Manganite VF']/molar_volume_manganite*1e6*Mn_molarmass/BD,label='Manganite')
    axs[num,2].plot(t,layers[num].output_DF['Birnessite2 VF']*7/molar_volume_birnessite*1e6*Mn_molarmass/BD,label='Birnessite')
    
    axs[num,1].plot(t,layers[num].output_DF['Total Tracer2']*1e6/1000*layers[num].porosity*layers[num].saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+}$ root uptake')
    
    # ax=axs[num,1].twinx()
    # axs[num,2].plot(t,layers[num].output_DF['Mn(OH)2(am) VF']/molar_volume_MnOH2am*1e6*Mn_molarmass/BD,label='Mn(OH)$_2$(am)',ls='-')
    axs[num,2].plot(t,(layers[num].output_DF['Total Mn++']+layers[num].output_DF['Total Mn+++'])*1e6/1000*layers[num].porosity*layers[num].saturation*Mn_molarmass/BD + 
                        (layers[num].output_DF['Birnessite2 VF']*7/molar_volume_birnessite)*1e6*Mn_molarmass/BD,'k-',label='Total Mn')
    
    axs[num,1].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
    axs[num,2].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
    # axs[num,1].set_ylim(*axs[0,1].get_ylim())

    
axs[0,0].legend()
axs[-1,0].set_xlabel('Time (years)')
axs[0,1].legend()
axs[0,2].legend()

axs[-1,1].set_xlabel('Time (years)')
axs[-1,2].set_xlabel('Time (years)')
axs[0,0].set_title('Organic Carbon')
axs[0,1].set_title('Manganese')
