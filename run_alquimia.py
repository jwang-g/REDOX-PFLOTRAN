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
        print(' '*indent+str(get_alquimiavector(obj)))
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
    output_units={'time':'s','Porosity':'NA'}
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
    
    return output_DF,output_units

def run_onestep(chem,data,dt,status,min_dt=0.1):
    converged=False
    porosity=data.state.porosity
    # Need to set this up so it keeps going up to the actual time step
    num_cuts=0
    while True:
        actual_dt=dt/2**num_cuts
        chem.ReactionStepOperatorSplit(ffi.new('void **',data.engine_state),actual_dt,ffi.addressof(data.properties),ffi.addressof(data.state),
                                        ffi.addressof(data.aux_data),status)
        data.state.porosity=porosity
        converged=status.converged
        if converged:
            break
        elif actual_dt>min_dt:
            num_cuts += 1
        else:
            raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
        
    # Run the remaining steps until the end of dt, assuming we can continue the same time step length
    for n in range(2**num_cuts-1):
        chem.ReactionStepOperatorSplit(ffi.new('void **',data.engine_state),actual_dt,ffi.addressof(data.properties),ffi.addressof(data.state),
                                        ffi.addressof(data.aux_data),status)
        data.state.porosity=porosity
    check_status(status,False)
    
    chem.GetAuxiliaryOutput(ffi.new('void **',data.engine_state),ffi.addressof(data.properties),ffi.addressof(data.state),
                                    ffi.addressof(data.aux_data),ffi.addressof(data.aux_output),status)
    check_status(status,False)
    return num_cuts
    

def run_simulation(input_file,simlength_days,dt=3600*12,min_dt=0.1,volume=1.0,saturation=1.0,
                    water_density=1000.0,temperature=25.0,porosity=0.25,pressure=101325.0):
    
    
    nsteps=int(365/(dt/(3600*24)))

    # Initialize Petsc first
    err = lib.PetscInitializeNoArguments()
    assert err==0, 'Error in Petsc initialization'

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
    chem.Setup(input_file.encode(),ffi.cast('_Bool',True),engine_state,sizes,functionality,status)
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
    chem_cond=ffi.new("AlquimiaGeochemicalCondition *")
    lib.AllocateAlquimiaGeochemicalCondition(len(init_cond_name), 0, 0, chem_cond);
    chem_cond.name[0:len(init_cond_name)]=init_cond_name.encode()

    # In alquimia batchChemDriver, surface site densities and isotherms are copied in from input file at this point

    # Initialize state data
    data.properties.volume=volume
    data.properties.saturation=saturation

    data.state.temperature=temperature
    data.state.water_density=water_density
    data.state.porosity=porosity
    data.state.aqueous_pressure=pressure

    chem.ProcessCondition(ffi.new('void **',data.engine_state),chem_cond,ffi.addressof(data.properties),ffi.addressof(data.state),ffi.addressof(data.aux_data),status)
    check_status(status,False)

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
        'actual_dt':numpy.ma.masked_all(nsteps,dtype=float),'ncuts':numpy.ma.masked_all(nsteps,dtype=float),
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
            num_cuts=run_onestep(chem,data,dt,status,min_dt=min_dt)
        except RuntimeError as err:
            print('ERROR on timestep %d: %s'%(step,err))
            print('Returning output so far')
            break
                
        if num_cuts>0:
            print('Timestep %d: Converged after %d cuts to dt = %1.1f s'%(step,num_cuts,dt/2**num_cuts))
        elif step%100==0:
            print('*** Step %d of %d ***'%(step,nsteps))
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
    return output_DF,output_units
