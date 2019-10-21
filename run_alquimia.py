from _alquimia import ffi,lib

import sys

def check_status(status,printmessage=True):
    if status.error != 0:
        message=ffi.string(status.message)
        lib.FreeAlquimiaEngineStatus(status)
        raise RuntimeError(message)
    else:
        if printmessage:
            print(str(ffi.string(status.message),encoding='utf8'))

def print_metadata(data):
    for datatype in ['aqueous_kinetic_names','ion_exchange_names','isotherm_species_names','mineral_names','primary_names','surface_site_names']:
        field=getattr(data.meta_data,datatype)
        print('%d %s:'%(field.size,datatype))
        for num in range(field.size):
            print('   '+str(ffi.string(field.data[num]),encoding='utf8'))

def get_alquimiavector(vec):
    s=vec.size
    if s==0:
        return []
    elif isinstance(vec.data[0],(float,int)):
        return [vec.data[n] for n in range(s)]
    else: # Assume it's a string
        return [str(ffi.string(vec.data[n]),encoding='utf8') for n in range(s)]

def print_alquimia_object(obj,indent=0):
    if hasattr(obj,'size'):
        # print(' '*indent+str(get_alquimiavector(obj)))
        for element in range(obj.size):
            print(' '*indent+'---')
            print_alquimia_object(obj.data[element],indent+2)
            print(' '*indent+'---')
    elif isinstance(obj,(int,float)):
        print(' '*indent+str(obj))
    elif ffi.typeof(obj) == ffi.typeof(ffi.new('char *')):
        print(' '*indent+str(ffi.string(obj),encoding='utf8'))
    elif isinstance(obj,ffi.CData):
        for field in dir(obj):
            print(' '*indent+field)
            print_alquimia_object(getattr(obj,field),indent+2)
    else:
        print(' '*indent+obj)


def convert_output(output,meta_data):
    import pandas
    output_DF=pandas.DataFrame(index=output['time'])
    output_DF['Porosity']=output['porosity']
    output_DF['ncuts']=output['ncuts']
    output_DF['actual_dt']=output['actual_dt']
    output_units={'time':'s','Porosity':'NA','ncuts':'NA','actual_dt':'s'}
    total=pandas.DataFrame(output['total_mobile'],columns=get_alquimiavector(meta_data.primary_names),index=output['time']).add_prefix('Total ')
    for col in total.columns:
        output_DF[col]=total[col]
    output_units.update([('Total '+s,'M') for s in get_alquimiavector(meta_data.primary_names)])
    free=pandas.DataFrame(output['free'],columns=get_alquimiavector(meta_data.primary_names),index=output['time']).add_prefix('Free ')
    for col in free.columns:
        output_DF[col]=free[col]
    output_units.update([('Free '+s,'M') for s in get_alquimiavector(meta_data.primary_names)])
    sorbed=pandas.DataFrame(output['immobile'],columns=get_alquimiavector(meta_data.primary_names),index=output['time']).add_prefix('Total Sorbed ')
    for col in sorbed.columns:
        output_DF[col]=sorbed[col]
    output_units.update([('Total Sorbed '+s,'mol/m^3') for s in get_alquimiavector(meta_data.primary_names)])
    mineral_VF=pandas.DataFrame(output['mineral_VF'],columns=get_alquimiavector(meta_data.mineral_names),index=output['time']).add_suffix(' VF')
    for col in mineral_VF.columns:
        output_DF[col]=mineral_VF[col]
    output_units.update([(s+' VF','m^3 mnrl/m^3 bulk') for s in get_alquimiavector(meta_data.mineral_names)])
    mineral_rate=pandas.DataFrame(output['mineral_rate'],columns=get_alquimiavector(meta_data.mineral_names),index=output['time']).add_suffix(' Rate')
    for col in mineral_rate.columns:
        output_DF[col]=mineral_rate[col]
    output_units.update([(s+' Rate','mol/m^3/sec') for s in get_alquimiavector(meta_data.mineral_names)])
    
    return output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600)),output_units

def run_onestep(chem,data,dt,status,min_dt=0.1,num_cuts=0,influx={}):
    converged=False
    porosity=data.state.porosity
    
    max_cuts=num_cuts
    actual_dt=dt/2**num_cuts
    if actual_dt<min_dt:
        raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
        
    for key in influx.keys():
        pos=get_alquimiavector(data.meta_data.primary_names).index(key)
        data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + influx[key]/actual_dt
        
    chem.ReactionStepOperatorSplit(ffi.new('void **',data.engine_state),actual_dt,ffi.addressof(data.properties),ffi.addressof(data.state),
                                    ffi.addressof(data.aux_data),status)
    data.state.porosity=porosity
    converged=status.converged
    if converged:
        check_status(status,False)
        
        chem.GetAuxiliaryOutput(ffi.new('void **',data.engine_state),ffi.addressof(data.properties),ffi.addressof(data.state),
                                        ffi.addressof(data.aux_data),ffi.addressof(data.aux_output),status)
        check_status(status,False)
        return max_cuts
    else:
        # Undo influx because we will do it again in the shorter run
        for key in influx.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(key)
            data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] - influx[key]/actual_dt
            
        # Run it twice, because we cut the time step in half
        ncuts=run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,influx=influx)
        if ncuts>max_cuts:
            max_cuts=ncuts
            
        # This starts from ncuts so it doesn't have to try all the ones that failed again
        for n in range(2**(ncuts-(num_cuts+1))):
            ncuts2=run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,influx=influx)
            if ncuts2>max_cuts:
                max_cuts=ncuts2

        return max_cuts
    

def convert_condition_to_alquimia(cond,name):
    primary=[]
    immobile=[]
    mineral=[]
    surf_complex=[]
    
    for constraint in cond:
        if constraint['kind']=='primary':
            primary.append(constraint)
        if constraint['kind']=='immobile':
            immobile.append(constraint)
        if constraint['kind']=='mineral':
            mineral.append(constraint)
        if constraint['kind']=='surf_complex':
            surf_complex.append(constraint)
    
    condition=ffi.new('AlquimiaGeochemicalCondition *')
    lib.AllocateAlquimiaGeochemicalCondition(len(name),len(primary)+len(immobile),len(mineral),condition)
    condition.name[0:len(name)]=name.encode()
    
    for num in range(len(primary)):
        constraint=primary[num]['constraints'][name]
        lib.AllocateAlquimiaAqueousConstraint(ffi.addressof(condition.aqueous_constraints.data[num]))
        condition.aqueous_constraints.data[num].primary_species_name[0:len(primary[num]['name'])]=primary[num]['name'].encode()
        if isinstance(constraint,float):
            condition.aqueous_constraints.data[num].constraint_type[0:len('total_aqueous')]='total_aqueous'.encode()
            condition.aqueous_constraints.data[num].value=constraint
        elif len(constraint.split())==1:
            condition.aqueous_constraints.data[num].constraint_type[0:len('total_aqueous')]='total_aqueous'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.replace('d','e'))
        elif constraint.split()[1]=='P':
            condition.aqueous_constraints.data[num].constraint_type[0:2]='pH'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
        elif constraint.split()[1]=='G':
            condition.aqueous_constraints.data[num].constraint_type[0:3]='gas'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
            condition.aqueous_constraints.data[num].associated_species[0:len(constraint.split()[2])]=constraint.split()[2].encode()
        elif constraint.split()[1]=='M':
            condition.aqueous_constraints.data[num].constraint_type[0:len('mineral')]='mineral'.encode()
            condition.aqueous_constraints.data[num].value=float(constraint.split()[0].replace('d','e'))
            condition.aqueous_constraints.data[num].associated_species[0:len(constraint.split()[2])]=constraint.split()[2].encode()
        else:
            raise ValueError('Unrecognized constraint type: %s'%constraint)

    for num in range(len(immobile)):
        constraint=immobile[num]['constraints'][name]
        lib.AllocateAlquimiaAqueousConstraint(ffi.addressof(condition.aqueous_constraints.data[num+len(primary)]))
        condition.aqueous_constraints.data[num+len(primary)].primary_species_name[0:len(immobile[num]['name'])]=immobile[num]['name'].encode()
        # Need to check within alquimia to make this an immobile constraint, I guess
        condition.aqueous_constraints.data[num+len(primary)].constraint_type[0:len('total_sorbed')]='total_sorbed'.encode()
        condition.aqueous_constraints.data[num+len(primary)].value=constraint
        
    for num in range(len(mineral)):
        constraint=mineral[num]['constraints'][name]
        lib.AllocateAlquimiaMineralConstraint(ffi.addressof(condition.mineral_constraints.data[num]))
        condition.mineral_constraints.data[num].mineral_name[0:len(mineral[num]['name'])]=mineral[num]['name'].encode()
        condition.mineral_constraints.data[num].volume_fraction=float(constraint.split()[0].replace('d','e'))
        condition.mineral_constraints.data[num].specific_surface_area=float(constraint.split()[1].replace('d','e'))
        
        
    return condition

def run_simulation(input_file,simlength_days,dt=3600*12,min_dt=0.1,volume=1.0,saturation=1.0,hands_off=True,
                    water_density=1000.0,temperature=25.0,porosity=0.25,pressure=101325.0,initcond=None,influx={},rateconstants={}):
    
    import time
    t0=time.time()
    nsteps=int(simlength_days/(dt/(3600*24)))

    engine_name='pflotran'
    status=ffi.new('AlquimiaEngineStatus *') 
    lib.AllocateAlquimiaEngineStatus(status)

    chem=ffi.new('AlquimiaInterface *')

    lib.CreateAlquimiaInterface(engine_name.encode(),chem,status)
    check_status(status)


    data=ffi.new('AlquimiaData *')
    sizes=ffi.addressof(data.sizes)
    functionality=ffi.addressof(data.functionality)
    engine_state=ffi.new('void **',data.engine_state)

    # input_file='../../alquimia-dev/benchmarks/batch_chem/ABCD_microbial.in'
    chem.Setup(input_file.encode(),ffi.cast('_Bool',hands_off),engine_state,sizes,functionality,status)
    check_status(status)
    data.engine_state=engine_state[0]

    # Allocates memory for all components of data structure
    lib.AllocateAlquimiaData(data)

    # Get metadata
    chem.GetProblemMetaData(engine_state,ffi.addressof(data.meta_data),status)
    check_status(status,False)

    print_metadata(data)
            

    # Set up initial condition
    init_cond_name='initial'
    if initcond is None:
        chem_cond=ffi.new("AlquimiaGeochemicalCondition *")
        lib.AllocateAlquimiaGeochemicalCondition(len(init_cond_name), 0, 0, chem_cond);
        chem_cond.name[0:len(init_cond_name)]=init_cond_name.encode()
    else:
        chem_cond=convert_condition_to_alquimia(initcond,init_cond_name)
        # print_alquimia_object(chem_cond)

    # Initialize state data
    data.properties.volume=volume
    data.properties.saturation=saturation

    data.state.temperature=temperature
    data.state.water_density=water_density
    data.state.porosity=porosity
    data.state.aqueous_pressure=pressure
    
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
                
    # Aqueous kinetic rate constants also need to be specified in hands-off mode
    # but it's wonky because alquimia doesn't store names for the reactions. So order might be undefined?
    # Also, Alquimia interface always sets backward rate to zero
    # Will be a problem with more than one reaction
    # Also only allows these to be specified for GENERAL_REACTION, not microbial, sandbox, etc
    for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
        data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]

    chem.ProcessCondition(ffi.new('void **',data.engine_state),chem_cond,ffi.addressof(data.properties),ffi.addressof(data.state),ffi.addressof(data.aux_data),status)
    check_status(status,False)
    
    # print_alquimia_object(data.state)
    print_alquimia_object(data.properties)

    # Pflotran sets porosity based on minerals or something? Needs to be reset
    data.state.porosity=porosity


    ##### At this point, the model should be initialized ##########
    print('''

    *****************************************************
    Successfully initialized alquimia geochemical engine
    *****************************************************

    ''')


    import numpy
    output={
        'total_mobile':numpy.ma.masked_all((nsteps,data.state.total_mobile.size),dtype=float),
        'free':numpy.ma.masked_all((nsteps,data.state.total_mobile.size),dtype=float),
        'immobile':numpy.ma.masked_all((nsteps,data.state.total_mobile.size),dtype=float),
        'aq_complex':numpy.ma.masked_all((nsteps,data.aux_output.secondary_free_ion_concentration.size),dtype=float),
        'mineral_VF':numpy.ma.masked_all((nsteps,data.state.mineral_volume_fraction.size),dtype=float),
        'mineral_rate':numpy.ma.masked_all((nsteps,data.state.mineral_volume_fraction.size),dtype=float),
        'time':numpy.arange(nsteps,dtype=float)*dt,
        'actual_dt':numpy.ma.masked_all(nsteps,dtype=float),'ncuts':numpy.ma.masked_all(nsteps,dtype=int),
        'porosity':numpy.ma.masked_all(nsteps,dtype=float)
        # There is some other stuff to add. Secondary species, mineral reaction rates, etc
    }

    output['total_mobile'][0,:]=get_alquimiavector(data.state.total_mobile)
    output['free'][0,:]=get_alquimiavector(data.aux_output.primary_free_ion_concentration)
    output['aq_complex'][0,:]=get_alquimiavector(data.aux_output.secondary_free_ion_concentration)
    output['immobile'][0,:]=get_alquimiavector(data.state.total_immobile)
    output['mineral_VF'][0,:]=get_alquimiavector(data.state.mineral_volume_fraction)
    output['time'][0]=0.0

    for step in range(1,nsteps):
        try:
            influx_step=dict([(key,influx[key][step%len(influx[key])]) for key in influx.keys()])
            # print('O2 before: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
            num_cuts=run_onestep(chem,data,dt,status,min_dt=min_dt,influx=influx_step)
            # print('O2 after: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
        except RuntimeError as err:
            print('ERROR on timestep %d: %s'%(step,err))
            print('Returning output so far')
            break
                

        if step%25==0:
            t1=time.time()
            print('*** Step {step:d} of {nsteps:d}. Time elapsed: {t:d} s ({tperstep:1.1g} s per timestep). Mean cuts: {meancuts:1.1f} Mean dt: {meandt:1.1f} s ***'.format(
                    step=step,nsteps=nsteps,t=int(t1-t0),tperstep=(t1-t0)/step,meancuts=output['ncuts'][step-10:step].mean(),meandt=output['actual_dt'][step-10:step].mean()))
        output['total_mobile'][step,:]=get_alquimiavector(data.state.total_mobile)
        output['immobile'][step,:]=get_alquimiavector(data.state.total_immobile)
        output['mineral_VF'][step,:]=get_alquimiavector(data.state.mineral_volume_fraction)
        output['time'][step]=step*dt
        output['actual_dt'][step]=dt/2**num_cuts
        output['ncuts'][step]=num_cuts
        output['mineral_rate'][step,:]=get_alquimiavector(data.aux_output.mineral_reaction_rate)
        output['free'][step,:]=get_alquimiavector(data.aux_output.primary_free_ion_concentration)
        output['aq_complex'][step,:]=get_alquimiavector(data.aux_output.secondary_free_ion_concentration)
        output['porosity'][step]=data.state.porosity
        
    
        
    import pandas
    # Put into a dataframe in the same format as reading it out of a tecfile
    output_DF,output_units=convert_output(output,data.meta_data)
    
    # Need to free all the Alquimia arrays?
    lib.FreeAlquimiaData(data)
    return output_DF,output_units

if __name__ == '__main__':
    
    
    
    import decomp_network
    pools = [
    decomp_network.decomp_pool(name='cellulose',CN=50,constraints={'initial':1e2},kind='immobile'),
    decomp_network.decomp_pool(name='HRimm',constraints={'initial':1e-20},kind='immobile'),

    decomp_network.decomp_pool(name='DOM1',CN=50,constraints={'initial':1e-30},kind='primary'),
    decomp_network.decomp_pool(name='H+',kind='primary',constraints={'initial':'4.0 P'}),
    decomp_network.decomp_pool(name='O2(aq)',kind='primary',constraints={'initial':1e-12}),
    decomp_network.decomp_pool(name='HCO3-',kind='primary',constraints={'initial':'400e-6 G CO2(g)'}),
    decomp_network.decomp_pool(name='Fe+++',kind='primary',constraints={'initial':'.37e-3 M Fe(OH)3'}),
    decomp_network.decomp_pool(name='Fe++',kind='primary',constraints={'initial':'0.37e-30'}),
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

    decomp_network.decomp_pool(name='Fe(OH)3',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d2 m^2/m^3'},kind='mineral'),
    # decomp_network.decomp_pool(name='Goethite',rate='1.d-5 mol/m^2-sec',constraints={'initial':'1.75d-2  1.d1 m^2/m^3'},kind='mineral'),
    # decomp_network.decomp_pool(name='Fe',rate='1.d-7 mol/m^2-sec',constraints={'initial':'1.0e-6  1. m^2/m^3'},kind='mineral'),
    decomp_network.decomp_pool(name='Fe(OH)2',rate='1.d-7 mol/m^2-sec',constraints={'initial':'0.0e-20  1.e-10 m^2/m^3'},kind='mineral'),
    decomp_network.decomp_pool(name='Rock(s)',rate='0.0 mol/m^2-sec',constraints={'initial':'0.5  5.0e3 m^2/m^3'},kind='mineral'),

    decomp_network.decomp_pool(name='>Carboxylate-',kind='surf_complex',mineral='Rock(s)',site_density=10.0,complexes=['>Carboxylic_acid']),
    ]
    
    # Initialize Petsc first
    initialized = ffi.new('PetscBool *')
    lib.PetscInitialized(initialized)   
    if not initialized[0]:
        err = lib.PetscInitializeNoArguments()
        assert err==0, 'Error in Petsc initialization'
            
            
    pools_atmoO2=pools.copy()
    for n,p in enumerate(pools_atmoO2):
        if p['name']=='O2(aq)':
            pools_atmoO2[n]=pools_atmoO2[n].copy()
            pools_atmoO2[n].update(constraints={'initial':'0.2 G O2(g)'})
            
    pools_highO2=pools.copy()
    for n,p in enumerate(pools_highO2):
        if p['name']=='O2(aq)':
            pools_highO2[n]=pools_highO2[n].copy()
            pools_highO2[n].update(constraints={'initial':'0.8e2 G O2(g)'})
    
    pools_lowFe=pools_atmoO2.copy()
    for n,p in enumerate(pools_lowFe):
        if p['name']=='Fe(OH)3':
            pools_lowFe[n]=pools_lowFe[n].copy()
            pools_lowFe[n].update(constraints={'initial':'0.0  1.e2 m^2/m^3'})
    
    # Using names assigned in alquimia. Not the best solution...
    rateconstants={
           '1.0 Fe++ + 0.25 O2(aq) + 1.0 H+ <-> 1.0 Fe+++ + 0.5 H2O'                                             : 1.0e2*1.0,
           "1.0e+00 DOM1  -> 3.3e-01 Acetate-  + 3.3e-01 HCO3-  + 2.0e+00 H+  + 3.3e-01 Tracer"                  : 1.0e-9*1.0,
           "1.0e+00 DOM1  + 1.0e+00 O2(aq)  -> 1.0e+00 HCO3-  + 1.0e+00 H+  + 1.0e+00 Tracer"                    : 1.0e-8*1.0,
           "1.0e+00 Acetate-  + 2.0e+00 O2(aq)  -> 2.0e+00 HCO3-  + 2.0e+00 H+  + 2.0e+00 Tracer"                : 1.0e-8*1.0,
           "1.0e+00 Acetate-  + 8.0e+00 Fe+++  -> 2.0e+00 HCO3-  + 8.0e+00 Fe++  + 9.0e+00 H+  + 2.0e+00 Tracer" : 5.0e-10*1.0,
           "1.0e+00 Acetate-  -> 1.0e+00 CH4(aq)  + 1.0e+00 HCO3-  + 1.0e+00 Tracer"                             : 1.0e-11*1.0,
           "cellulose decay (SOMDEC sandbox)"                                                                    : 10.0/(365*24*3600)
    }
    
    from numpy import zeros
    O2_const=zeros(365*24)+0.01/3600
    O2_periodic=zeros(365*24)
    O2_periodic[0::60*24]=1.0e2/3600
    output_highO2,output_units=run_simulation('fermentation_generated.in',365,3600,initcond=pools_atmoO2,influx={'O2(aq)':O2_const},hands_off=False,rateconstants=rateconstants,saturation=1.0)
    output,output_units=run_simulation('fermentation_generated.in',365,3600,initcond=pools_atmoO2,hands_off=False,rateconstants=rateconstants)    
    output_lowFe,output_units=run_simulation('fermentation_generated.in',365,3600,initcond=pools_lowFe,hands_off=False,rateconstants=rateconstants)
    output_periodicO2,output_units=run_simulation('fermentation_generated.in',365,3600,initcond=pools_atmoO2,influx={'O2(aq)':O2_periodic},hands_off=False,rateconstants=rateconstants)
    
    # err = lib.PetscFinalize()
    
    import decomp_network
    import plot_pf_output
    from pylab import *
    result,units=plot_pf_output.convert_units(output,output_units,'M')
    result_lowFe,units=plot_pf_output.convert_units(output_lowFe,output_units,'M')
    result_highO2,units=plot_pf_output.convert_units(output_highO2,output_units,'M')
    result_periodicO2,units=plot_pf_output.convert_units(output_periodicO2,output_units,'M')
    # figure('Network diagram',figsize=(11.8,4.8));clf()
    # ax=subplot(121)
    # decomp_network.draw_network(fermentation_network_Fe_CH4,omit=['secondary','surf_complex','NH4+','Rock(s)'],arrowstyle='-|>')
    # title('Decomposition network diagram (without complexes)')
    # 
    # ax=subplot(122)
    # decomp_network.draw_network(fermentation_network_Fe_CH4,omit=['NH4+','Rock(s)'],arrowstyle='-|>')
    # title('Decomposition network diagram (with aqueous complexes)')

    figure('Cellulose sim (alquimia)');clf()

    subplot(221)
    # ax.set_yscale('log')
    handles=[]
    for pool in ['Total Sorbed cellulose','Total CH4(aq)']:
        l=plot(result[pool],label=pool)[0]
        plot(result_lowFe[pool],ls='--',c=l.get_color())
        plot(result_highO2[pool],ls='-.',c=l.get_color())
        plot(result_periodicO2[pool],ls=':',c=l.get_color())
        handles.append(l)
    l=plot(result['Total Tracer']+result['Total Sorbed HRimm'],label='CO2 produced')[0]
    handles.append(l)
    plot(result_lowFe['Total Tracer']+result_lowFe['Total Sorbed HRimm'],c=l.get_color(),ls='--')
    plot(result_highO2['Total Tracer']+result_highO2['Total Sorbed HRimm'],c=l.get_color(),ls='-.')
    plot(result_periodicO2['Total Tracer']+result_periodicO2['Total Sorbed HRimm'],c=l.get_color(),ls=':')
    # handles.append(Line2D([0,0],[0,0],color='k',ls='-',label='Aerobic'))
    # handles.append(Line2D([0,0],[0,0],color='k',ls='--',label='With Fe(III)'))
    # handles.append(Line2D([0,0],[0,0],color='k',ls=':',label='With methanogenesis'))
    # handles.append(Line2D([0,0],[0,0],color='k',ls='-.',label='Anaerobic'))
    legend(handles=handles,fontsize='small',ncol=2)
    title('Concentrations')
    ylabel('Concentration (M)')
    xlabel('Time (days)')

    # figure('Cellulose sim pH and log',figsize=(6,8));clf()

    subplot(223)
    l=plot(-log10(result['Free H+']),label='Anaerobic')[0]
    plot(-log10(result_highO2['Free H+']),c=l.get_color(),ls='-.',label='Aerobic')
    plot(-log10(result_periodicO2['Free H+']),c=l.get_color(),ls=':',label='Periodically aerobic')
    plot(-log10(result_lowFe['Free H+']),c=l.get_color(),ls='--',label='Low Fe(III)')
    legend(fontsize='small')

    title('pH')
    ylabel('pH')
    xlabel('Time (days)')

    ax=subplot(222)
    ax.set_yscale('log')
    # for pool in ['Free DOM1','Free Acetate-']:
    #     l=plot(result[pool],label=pool)[0]
    #     # plot(result_Fe[pool],ls='--',c=l.get_color())
    #     # plot(result_Fe_CH4[pool],ls=':',c=l.get_color())
    #     plot(result_highO2[pool],ls='-.',c=l.get_color())
    # 
    # title('Concentrations (log scale)')
    # ylabel('Concentration (M)')
    # xlabel('Time (days)')
    # legend(fontsize='small',ncol=1)
    # ylim(1e-15,1e-1)

    l=plot(result.index.values[:-1],diff(result['Total CH4(aq)'])/diff(result.index.values),label='CH4')[0]
    plot(result_highO2.index.values[:-1],diff(result_highO2['Total CH4(aq)'])/diff(result_highO2.index.values),ls='-.',c=l.get_color())
    plot(result_periodicO2.index.values[:-1],diff(result_periodicO2['Total CH4(aq)'])/diff(result_periodicO2.index.values),ls=':',c=l.get_color())
    plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total CH4(aq)'])/diff(result_lowFe.index.values),ls='--',c=l.get_color())

    # l=plot(result.index.values[:-1],diff(result['Total Fe++'])/diff(result.index.values),label='Fe(II)')[0]
    l=plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total Fe++'])/diff(result_lowFe.index.values),ls='--',label='Fe++')[0]
    plot(result_highO2.index.values[:-1],diff(result_highO2['Total Fe++'])/diff(result_highO2.index.values),ls='-.',c=l.get_color())
    plot(result_periodicO2.index.values[:-1],diff(result_periodicO2['Total Fe++'])/diff(result_periodicO2.index.values),ls=':',c=l.get_color())
    # plot(result_Fe_CH4.index.values[:-1],diff(result_Fe_CH4['Total Fe++'])/diff(result_Fe_CH4.index.values),ls=':',c=l.get_color())

    l=plot(result.index.values[:-1],diff(result['Total Tracer']+result['Total Sorbed HRimm'])/diff(result.index.values),label='CO2')[0]
    # plot(result_highO2.index.values[:-1],diff(result_highO2['Total Tracer']+result_highO2['HRimm'])/diff(result_highO2.index.values),ls='-.',c=l.get_color())
    plot(result_lowFe.index.values[:-1],diff(result_lowFe['Total Tracer']+result_lowFe['Total Sorbed HRimm'])/diff(result_lowFe.index.values),ls='--',c=l.get_color())
    

    title('Production rates')
    ylabel('Rate (M/day)')
    xlabel('Time (days)')
    legend(fontsize='small')



    ax=subplot(224)
    ax.set_yscale('log')
    for pool in ['Free DOM1','Free Acetate-','Free Fe+++','Free Fe++','Free O2(aq)']:
        l=plot(result[pool],label=pool)[0]
        plot(result_lowFe[pool],ls='--',c=l.get_color())
        plot(result_periodicO2[pool],ls=':',c=l.get_color())
        plot(result_highO2[pool],ls='-.',c=l.get_color())
    
    title('Concentrations (log scale)')
    ylabel('Concentration (M)')
    xlabel('Time (days)')
    legend(fontsize='small',ncol=1)
    ylim(1e-15,1e-1)



    title('Log concentrations')
    ylabel('Concentration (M)')
    xlabel('Time (days)')
    legend(fontsize='small')


    tight_layout()
    show()

