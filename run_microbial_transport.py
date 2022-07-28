#Created by Ben Sulman for the Manganese soil profile, and revised by Jiaze Wang to simulate delta marsh soil profile,02/18/2021

from run_alquimia import get_alquimiavector,ffi,lib,check_status,init_alquimia,convert_condition_to_alquimia,print_metadata,convert_rateconstants
import decomp_network
#import delmar_network_tidev2 as Mar
import microbes as Mar
from matplotlib import pyplot
import numpy
import xarray
import pdb


class layer:
    def __init__(self,volume,saturation=1.0,temperature=20.0,water_density=1000.0,porosity=0.25,pressure=101325.0,BD=1.5,CEC=50.0,rateconstants={},diffquo={}):
        self.volume=volume # m3
        self.saturation=saturation   # out of 1.0
        self.temperature=temperature # degrees C
        self.water_density=water_density # kg/m3
        self.porosity=porosity # out of 1.0
        self.pressure=pressure # Pa
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
        self.flow_in={}
        self.flow_out={}
        self.diffquo=diffquo
        self.BD=BD   # g/cm3
        self.CEC=CEC # meq/kg
        
        # Some things that I'm not using are skipped

    def __str__(self):
        out='Layer with volume = %1.2e m, porosity = %1.2f, CEC = %1.2f'%(self.volume,self.porosity,self.CEC)
        out=out + '\n                    Mobile    Immobile\n'
        for spec in self.total_mobile:
            out=out + spec.ljust(20)+'%1.2e   %1.2e\n'%(self.total_mobile[spec],self.total_immobile[spec])
        out=out+'Rateconstants:\n'
        for key in sorted(self.rateconstants):
            out=out+key.ljust(30)+': %1.2e\n'%self.rateconstants[key]
        
        return out
        
    def copy_from_alquimia(self,data):
        self.volume=data.properties.volume
        self.saturation=data.properties.saturation
        self.temperature=data.state.temperature
        self.water_density=data.state.water_density
        self.porosity=data.state.porosity
        self.pressure=data.state.aqueous_pressure
        
        self.mineral_names=get_alquimiavector(data.meta_data.mineral_names)
        for num in range(data.meta_data.mineral_names.size):
            self.mineral_rate_cnst[self.mineral_names[num]]=data.properties.mineral_rate_cnst.data[num] # Alquimia doesn't actually update this
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
            'porosity':numpy.ma.masked_all(nsteps,dtype=float),
            #'CEC H+':numpy.ma.masked_all(nsteps,dtype=float),
            'flow_in':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            'flow_out':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
        }
        
    def write_output(self,step,dt,num_cuts=0):
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
        # Keep track of how much H+ is in the carboxylate buffer rather than the CEC site, so we can plot CEC-exchangeable H+ later. 
        # Assumes there is one ion exchange site, and that the buffering sorption site is on Rock(s)
        # Aux doubles in this spot stores free surface site density (probably best not to rely on aux_doubles in general though)
        #self.output['CEC H+'][step]=self.total_immobile['H+']-(self.surface_site_density['>Carboxylate-']*self.mineral_volume_fraction['Rock(s)']-self.aux_doubles[len(self.total_mobile)*2+len(self.secondary_free_ion_concentration)+1+self.surface_site_names.index('>Carboxylate-')])
        self.output['flow_in'][step,:]=numpy.array([self.flow_in.get(name,0.0) for name in self.primary_names])
        self.output['flow_out'][step,:]=numpy.array([self.flow_out.get(name,0.0) for name in self.primary_names])
        
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
        secondary=pandas.DataFrame(self.output['aq_complex'],columns=self.secondary_names,index=self.output['time'])
        for col in secondary.columns:
            output_DF[col]=secondary[col]
        output_units.update([(s,'M') for s in self.secondary_names])
        
        flow_in=pandas.DataFrame(self.output['flow_in'],columns=self.primary_names,index=self.output['time']).add_suffix(' inflow')
        for col in flow_in.columns:
            output_DF[col]=flow_in[col]
        output_units.update([(s+' inflow','mol/m2/sec') for s in self.primary_names])
        flow_out=pandas.DataFrame(self.output['flow_out'],columns=self.primary_names,index=self.output['time']).add_suffix(' outflow')
        for col in flow_out.columns:
            output_DF[col]=flow_out[col]
        output_units.update([(s+' outflow','mol/m2/sec') for s in self.primary_names])
        
        #output_DF['CEC H+']=pandas.DataFrame(self.output['CEC H+'],index=self.output['time'])
        #output_units['CEC H+']='mol/m^3'
        
        self.output_DF=output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600))
        self.output_units=output_units
        
        
    def run_onestep(self,chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=10,truncate_concentration=0.0,rateconstants={},min_cuts=0,TRc={},Ebc={}):
        self.copy_to_alquimia(data)
        bc_tmp=bc
        converged=False
        porosity=data.state.porosity
        saturation=data.properties.saturation
        max_cuts=num_cuts
        actual_dt=dt/2**num_cuts
        Eb={}
        imbspec=['SOMC','SOMN']
        #print(data.properties.saturation)
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]

        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            if data.state.total_mobile.data[pos] < truncate_concentration and data.state.total_mobile.data[pos] != 0.0:
                for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                    if '->' in reactname and spec in reactname.split('->')[0].split():  # reactants
                        data.properties.aqueous_kinetic_rate_cnst.data[num]=0.0
        cut_for_flux=False

        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)

            bc_conc=bc.total_mobile.data[pos]
            local_conc=data.state.total_mobile.data[pos]

            if local_conc<0:
                raise RuntimeError('Initial concentration of %s < 0'%spec)

        if cut_for_flux or num_cuts<min_cuts:
            converged=False
        else:
            chem.ReactionStepOperatorSplit(ffi.new('void **',data.engine_state),actual_dt,ffi.addressof(data.properties),ffi.addressof(data.state),
                                            ffi.addressof(data.aux_data),status)
            data.state.porosity=porosity
            converged=status.converged
            # Check for negative concentrations
            # if (numpy.array(get_alquimiavector(data.state.total_mobile))<0).any() or (numpy.array(get_alquimiavector(data.state.total_immobile))<0).any():
            #     converged=False

        if converged:
            check_status(status,False)

            chem.GetAuxiliaryOutput(ffi.new('void **',data.engine_state),ffi.addressof(data.properties),ffi.addressof(data.state),
                                            ffi.addressof(data.aux_data),ffi.addressof(data.aux_output),status)
            check_status(status,False)

        if converged:
            for spec in diffquo.keys():
                if spec in TRc.keys():
                    pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                    #if spec != 'SOMC':
                    if spec not in imbspec:
                        data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + TRc[spec]*actual_dt
                        #print('before ebull',data.state.total_mobile.data[pos])
                        #ctmp=data.state.total_mobile.data[pos] + TRc[spec]*actual_dt + Ebc[spec]*actual_dt
                        if Ebc[spec] < 0.0:
                            ctmp=data.state.total_mobile.data[pos] + Ebc[spec]*actual_dt
                            #if spec=='CH4(aq)':
                                #Ebc[spec]=0.0
                                #print('check ebull or not ctmp,data,Ebc*dt',ctmp,data.state.total_mobile.data[pos],Ebc[spec]*actual_dt)
                            if ctmp > 0:
                                data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + Ebc[spec]*actual_dt
                                #print('ebull',data.state.total_mobile.data[pos])
                                Eb[spec]=Ebc[spec]
                                if spec == 'HCO3-':
                                    Hpos=get_alquimiavector(data.meta_data.primary_names).index('H+')
                                    data.state.total_immobile.data[Hpos]=data.state.total_immobile.data[Hpos]+Ebc[spec]*1000*porosity*saturation*actual_dt
                                elif spec == 'HS-':
                                    Hpos=get_alquimiavector(data.meta_data.primary_names).index('H+')
                                    data.state.total_immobile.data[Hpos]=data.state.total_immobile.data[Hpos]+Ebc[spec]*1000*porosity*saturation*actual_dt
                            else:
                                Eb[spec]=Ebc[spec]*0
                                #print('non ebull',data.state.total_mobile.data[pos])
                        else:
                            Eb[spec]=Ebc[spec]
                    else:
                        data.state.total_immobile.data[pos] = data.state.total_immobile.data[pos] + TRc[spec]*actual_dt
                        Eb[spec]=Ebc[spec]
                    #if spec=='CH4(aq)':
                    #    print('transport induced change %3.20f,%2.30f,%d'%(TRc[spec]*actual_dt,data.state.total_mobile.data[pos],actual_dt))
                else:
                    pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                    data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos]
                    data.state.total_immobile.data[pos] = data.state.total_immobile.data[pos]
            self.copy_from_alquimia(data)
            return max_cuts,Eb
        else:
            if actual_dt/2<min_dt:
                if cut_for_flux:
                    raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
                else:
                    raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
            #flux_tmp=numpy.zeros(data.meta_data.primary_names.size)
            # data will be reset to layer contents at beginning of next call
            # Run it twice, because we cut the time step in half
            ncuts,Ebgas=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc_tmp,truncate_concentration=truncate_concentration,rateconstants=rateconstants,TRc=TRc,Ebc=Ebc)
            # If this completes, it means that run_onestep was successful
            self.copy_from_alquimia(data)
            bc_tmp=bc

            if ncuts>max_cuts:
                max_cuts=ncuts
            Eb=Ebgas

            # This starts from ncuts so it doesn't have to try all the ones that failed again
            for n in range(2**(ncuts-(num_cuts+1))):
                ncuts2,Ebgas=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc_tmp,truncate_concentration=truncate_concentration,rateconstants=rateconstants,TRc=TRc,Ebc=Ebc)

                self.copy_from_alquimia(data)
                bc_tmp=bc
                if ncuts2>max_cuts:
                    max_cuts=ncuts2
                Eb=Ebgas
            return max_cuts,Eb

def convert_to_xarray(layers,t0=0.0,drop_nas=True,convert_output=True):
    for l in layers:
        if convert_output or not hasattr(l,'output_DF'):
            l.convert_output()
    data_array = xarray.concat([xarray.Dataset.from_dataframe(layer.output_DF) for layer in layers],dim='layer').rename({'index':'time','layer':'depth'})
    data_array['dz']=xarray.DataArray([layer.volume for layer in layers],dims='depth',attrs={'units':'cm'})*100     #cjw why multiply by 100
    data_array['z_bottom']=data_array['dz'].cumsum()
    data_array['z_top']=data_array['z_bottom']-data_array['dz']
    data_array['z_middle']=data_array['z_bottom']-data_array['dz']/2
    data_array['z_bottom'].attrs['units']='cm'
    data_array['z_middle'].attrs['units']='cm'
    data_array['z_top'].attrs['units']='cm'

    data_array['depth']=data_array['z_middle']

    for var in layers[0].output_units:
        data_array[var].attrs['units']=layers[0].output_units[var]
# cjw: revise needed to construct layers for marsh site.
    # To do: Add other layer properties/attributes and also leaf Mn concentration
    data_array['saturation']=xarray.DataArray([layer.saturation for layer in layers],dims='depth',attrs={'units':'fraction'})
    data_array['BD']=xarray.DataArray([layer.BD for layer in layers],dims='depth',attrs={'units':'g cm-3'})
#    data_array['CEC']=xarray.DataArray([layer.CEC for layer in layers],dims='depth',attrs={'units':'meq kg-1'})

    data_array['time']=data_array['time']+t0
    data_array['time'].attrs['units']='days'

#cjw    if leaf_Mn is not None:
#cjw        data_array['litter_Mn']=xarray.DataArray(leaf_Mn,dims='litter_year',attrs={'units':'mmol/kg dry mass'})

    if drop_nas:
        data_array=data_array.dropna(dim='time')

    return data_array

def copy_to_layers(data_xarray,layers):
    for var in data_xarray.variables:
        for depth in range(len(data_xarray.depth)):
            if var.startswith('Total Sorbed'):
                layers[depth].total_immobile[var[len('Total Sorbed '):]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()
            elif var.startswith('Total '):
                layers[depth].total_mobile[var[len('Total '):]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()
            elif var.startswith('Free '):
                layers[depth].free_mobile[var[len('Free '):]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()
            elif var.endswith('VF'):
                layers[depth].mineral_volume_fraction[var[:-3]]=data_xarray[var].dropna(dim='time').isel(time=-1,depth=depth).item()

def Ebullition(Cw,volume,porosity,saturation,Press,Temp,tstep,Flag):
    k=1/dt        #compute Ebflux per second
    R=8.3144621  ## gas constant in unit of m3Pa/K/mol
    Tstd=273.15
    T0=Tstd+25    ## kalvin temp K  # standard T
    T=Temp+Tstd   ## convert from degree C to K
    Psum=0.0        #sum of partial pressure for all gas species in Pa
    Pgas={}         #partial pressure for each gas
    por=porosity
    sat=saturation
    vol_watr=volume*por*sat      #unit is m3
    alpha={}        ## Bunsen coefficient, dimensioness
    Eclim={}        #Ebulision concentration limit
    Ebflx={}
    #Hcp/Hgas={}     #Henry's law constant unit is mol/(m3*Pa) Hcc=Hcp * RT, Hcc is dimensionless which is Conc_aq/Concgas
    Hcp={'O2(aq)':1.2*1e-5*np.exp(-1700*(1/T-1/T0)),
    'CH4(aq)':1.4*1e-5*np.exp(-1900*(1/T-1/T0)),# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'HCO3-':3.3*1e-4*np.exp(-2400*(1/T-1/T0)),
    'HS-':1.0*1e-3*np.exp(-2100*(1/T-1/T0)),
    'H2(aq)':7.7*1e-6*np.exp(-530*(1/T-1/T0)),#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'N2(aq)':6.4*1e-6*np.exp(-1600*(1/T-1/T0)),#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'N2O(aq)':2.4*1e-4*np.exp(-2700*(1/T-1/T0)),
    }    ##unit in mol/(m3*Pa) # Ref Sander R., Compilation of Henry's law constants(version 4.0) for water as solvent, 2015

    ## Method 0 using concentration threshold for ebullision
    if Flag==0:
        gas_frac={'O2(aq)':0,#
        'CH4(aq)':0.5,# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'HCO3-':0.1,
        'HS-':0.05,
        'H2(aq)':0.05,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2(aq)':0.27,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2O(aq)':0.03,}    ## bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        for gidx in Hcp.keys():
            alpha[gidx]=Hcp[gidx]*Tstd*R #Bunsen coefficient or solubility of gas Ref. Tang 2010, Biogeosciences; Sander R., 2015
            Eclim[gidx]=Press*alpha[gidx]/(R*T*1000)   ##mol/L
        for Gspec in Cw.keys():
            if Gspec in Hcp.keys():
                if Cw[Gspec] > Eclim[Gspec]:
                    Ebflx[Gspec]=(Eclim[Gspec]-Cw[Gspec])*gas_frac[Gspec]*k
                else:
                    Ebflx[Gspec]=0.0
            else:
                Ebflx[Gspec]=0.0
    ## Method 1 using pressure threshold for ebullision, which only happens when pores are filled with water instead of unsatureated environment
    elif Flag==1:
        #trunc_limit=1e-30    #mol/L
        for Gspec in Hcp.keys():
            Pgas[Gspec]=Cw[Gspec]*1000/Hcp[Gspec]   #Cw is in unit of mol?l need to convert to mol/m3
            Psum=Psum+Pgas[Gspec]
        if Psum > Press:          #ebullison occurs when Psum exceed the hydrostatic pressure + Patm
            Ebfrac=(Press-Psum)/Psum
        else:
            Ebfrac=0.0
        for Gspec in Cw.keys():
            if Gspec in Hcp.keys():
                Ebflx[Gspec]=k*Ebfrac*Pgas[Gspec]/(R*T*1000)     # Ebflux is positive when ebullision occurs Eb mol/M3 to mol/L
            else:
                Ebflx[Gspec]=0.0
    return Ebflx

rate_scale=1e-6
truncate_conc=1e-30
thresh=truncate_conc*1.01
#reaction_network=Mar.make_network()
#network=decomp_network.decomp_network(pools,microbe_reactions)
# See Smeaton, Christina M., and Philippe Van Cappellen. 2018. “Gibbs Energy Dynamic Yield Method (GEDYM): Predicting Microbial Growth Yields under Energy-Limiting Conditions.” Geochimica et Cosmochimica Acta 241 (November): 1–16. https://doi.org/10.1016/j.gca.2018.08.023.

##cjw update half saturation term for each substrates
for idx in range(len(Mar.microbe_reactions)):             # a list with different reaction dict
    #rxnnm=Mar.microbe_reactions[idx]['name'].split(" ")                                  # alist of words in names of reaction seperated by space
    rxnnm=Mar.microbe_reactions[idx]['name']        #cjw try a new way to include cost-benefit for Ks
    for m in Mar.microbes:
        mxgrow,gnk=Mar.growth_rate(m)
        for igene in range(len(m.genes)):
            micnm=m.genes[igene]['name']+ f' ({m.name})'
            cmplx=Mar.gnpnlty[m.genes[igene]['name']]  #cost-benefit for each reaction or feature gene function
            if micnm == rxnnm:
                if 'inhibition_terms' in Mar.microbe_reactions[idx].keys():
                    inhterm=Mar.microbe_reactions[idx]['inhibition_terms']      # a list of dictionary of inhibition informariton
                    for inh in range(len(inhterm)):
                        inhspec=inhterm[inh]['species']
                        if inhspec in gnk.keys():
                            inhterm[inh]['k']=gnk[inhspec]*cmplx
                if 'monod_terms' in Mar.microbe_reactions[idx].keys():
                    modterm=Mar.microbe_reactions[idx]['monod_terms']
                    for mod in range(len(modterm)):
                        modspec=modterm[mod]['species']
                        if modspec in gnk.keys():
                            modterm[mod]['k']=gnk[modspec]*cmplx
#        micnm=f'({m.name})'
#        for micnm in rxnnm:
#            if 'inhibition_terms' in Mar.microbe_reactions[idx].keys():
#                inhterm=Mar.microbe_reactions[idx]['inhibition_terms']      # a list of dictionary of inhibition informariton
#                for inh in range(len(inhterm)):
#                    inhspec=inhterm[inh]['species']
#                    if inhspec in gnk.keys():
#                        inhterm[inh]['k']=gnk[inhspec]
#            if 'monod_terms' in Mar.microbe_reactions[idx].keys():
#                modterm=Mar.microbe_reactions[idx]['monod_terms'] 
#                for mod in range(len(modterm)):
#                    modspec=modterm[mod]['species']
#                    if modspec in gnk.keys():
#                        modterm[mod]['k']=gnk[modspec]
#print(Mar.microbe_reactions)
reaction_network=decomp_network.decomp_network(Mar.pools,Mar.microbe_reactions)

#print('reactions',Mar.microbe_reactions)
rateconstants={}
for r in range(len(Mar.microbe_reactions)):
    rateconstants[Mar.microbe_reactions[r]['name']]=rate_scale
precision=2


##cjw diffusion coefficient
##diffucoef={'O2(aq)':0.001**((l+1)*0.75),'CH4(aq)':0.001**((l+1)*0.75)}
#cjw
input_file='microbial_test_network.in'
decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=True,CO2name='HCO3-',truncate_concentration=1e-25,database='/home/46w/wetland/REDOX-PFLOTRAN/hanford.dat',verbose=True,length_days=100)

#decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=False,truncate_concentration=truncate_conc)
#cjwdecomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt','deltamarsh.in',length_days=30,log_formulation=False)

# Read secondary complex names from input file since Alquimia does not provide them

##cjw run_alquimia.run_simulation already include this

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

import time
starting_time=time.time()

rateconstants_warmed=rateconstants.copy()
#cjw pH for fresh, brackish and salt is from DeLaune, 1983 (average data from 0-50cm).
#bavg,poravg,satavg=READSOIL_CRMS.soils()
#for station in stnm:
nm=['CRMS2825']#,'CRMS2825','CRMS3166']

#for pH in numpy.arange(6.3,7.3):
for stnm in nm:
    chem,data,sizes,status=init_alquimia(input_file,hands_off=True)
        #cjw below is the codes for layer setting up.
        # Set up layers
        # Top (organic) layer should be thinner and have lower bulk density though
        # Low bulk density causes simulation to slow or crash though. Actually CEC being too low (<100 combined with BD<1) is the problem
#cjw        layers=[layer(0.05,rateconstants=rateconstants_warmed,BD=0.05,porosity=0.5)]+[layer(0.1,BD=0.25,rateconstants=rateconstants_warmed) for num in range(3)]
    #layers=[layer(0.1,rateconstants=rateconstants_warmed,BD=0.08,porosity=0.95,saturation=1)]
#cjw assume it is always saturated even when water table is low (precipitation could bring water or river diversion)
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=bavg[stnm][0],porosity=poravg[stnm][0],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][1],porosity=poravg[stnm][1],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][2],porosity=poravg[stnm][2],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][3],porosity=poravg[stnm][3],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][4],porosity=poravg[stnm][4],saturation=1)]+[layer(0.1,rateconstants=rateconstants_warmed,BD=bavg[stnm][5],porosity=poravg[stnm][5],saturation=1)]
#station 2825
    layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.15,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.15,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.85)]
    #layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.15,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.15,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.85)]+[layer(0.05,BD=0.13,rateconstants=rateconstants_warmed,porosity=0.95,saturation=0.88)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.89)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.96,saturation=0.91)]+[layer(0.1,BD=0.09,rateconstants=rateconstants_warmed,porosity=0.97,saturation=0.92)]    
    
    for l in layers:
        l.secondary_names=secondary_names
    #for l in layers:
#           l.initcond=Mar.pools.copy()
        l.initcond=decomp_network.change_constraints(Mar.pools,{'Ca++':'1e-30',
                                                        'Na+':'1e-30','Cl-':'1e-30','CH4(aq)':'2.37998572e-09',
                                                        'SO4--':'1e-30',#})     #'%1.8f'%((0.75-suldpth[sulidx])*1e-15/30)})
                                                        'O2(aq)':'0.2 G O2(g)'})               #cjw: set the sulfate profile so4 decrease with depth (0.75-suldpth[sulidx])/30)


    #for l in layers:
    #    l.initcond=decomp_network.change_constraints(Mar.pools,{'O2(aq)':0.002})
    bc=layers[0].initcond
    #print(bc,layers[0].initcond)

#cjw    bc=None
    #bio_dt=1
    #nyears=1
    #nsteps=365*24//(dt//3600)*nyears
    #nsteps=36000

    #cjw transport time step
    bio_dt=1
    dt=3600
    repyr=4            ##spin-up at least 6 cycle
    nyears=3*repyr      ## frequency of rerun the entire observation dataset
    nsteps=(365+366+365)*24//(dt//3600)*repyr
    skip_tm=dt/bio_dt    ## biology skip time regarding physical time step
        # Set up initial condition
    for l in layers:
            # Initialize state data
        data.properties.volume=l.volume
        data.properties.saturation=l.saturation

        data.state.temperature=l.temperature
        data.state.water_density=l.water_density
        data.state.porosity=l.porosity
        data.state.aqueous_pressure=l.pressure

        #print('layer porosity,BD: %1.4f,%1.4f'%(l.porosity,l.BD)) 

            # Set properties: surface site density and mineral rate constants
            # This is necessary when running in hands-on mode
        if l.initcond is not None:
            for constraint in l.initcond:
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
            bc_state.temperature=l.temperature #layers[0].temperature           #replace layers[0]
            bc_state.water_density=l.water_density #layers[0].water_density
            bc_state.porosity=l.porosity #layers[0].porosity
            bc_state.aqueous_pressure=l.pressure #layers[0].pressure
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
#cjw        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
#cjw            data.properties.aqueous_kinetic_rate_cnst.data[num]=l.rateconstants[reactname]
        # CEC also needs to be specified in hands-off mode
        # CEC in PFLOTRAN is in eq/m3, so it must be converted from normal units of meq/kg
#cjw        CEC_pf=l.CEC*1e-3*(l.BD*1e-3*100**3)
#cjw        print('Applying CEC: %1.2g'%l.CEC)
#cjw        data.state.cation_exchange_capacity.data[0]=CEC_pf        
        
        init_cond=convert_condition_to_alquimia(l.initcond,'initial')
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

    print('''


            *************************************
            * Starting simulation with station = %s, Ndep = %03d, Warming = %d *
            *************************************


        '''%(stnm,0,20))        #cjw remove warming



    initial_HCO3 = l.total_mobile['HCO3-']
    initial_O2 = l.total_mobile['O2(aq)']
    initial_H = l.total_mobile['H+']
        # Flow rate cm/s = 10 L/m2/s, positive is downward
    #flow_rate=numpy.linspace(1e-6,1e-7,len(layers)) # Rate declines linearly with depth, assumes removal or accumulation in lower layers
    #flow_rate=numpy.zeros(len(layers))*0
#    flow_rate=numpy.linspace(1e-7,1e-8,len(layers))
    min_dt=0.1
    truncate_concentration=1e-20
    
    ###cjw compute tide water level below to
    import scipy as sp
    import mpmath as mp
    import numpy as np
    import math
    import pandas as pd
    import operator
    import seaborn as sns

    ##cjw parameter needed for transport
    tide_t=range(0,nsteps)
    Nlyr=range(0,len(layers))
    tdt=3600
    #hypath=r"/Users/46w/Documents/estuary wetland/CRMS_data/10WaterLevel_3166.csv"
    hyfile1='10WaterLevel_'+stnm[4:8]+'.csv'
    hyfile2='10Salinity_'+stnm[4:8]+'.csv'
    hyfile3='10Temp_'+stnm[4:8]+'.csv'
    #wtl=pd.read_csv('./10WaterLevel_2825.csv')         ##
    #salinity=pd.read_csv('./10Salinity_2825.csv') 
    wtl=pd.read_csv('./'+hyfile1)         ##
    salinity=pd.read_csv('./'+hyfile2)
    wtemp=pd.read_csv('./'+hyfile3)
    wtemp=wtemp.interpolate()

    #wtl=pd.read_csv('./10WaterLevel_2825.csv')         ##
    #wtl=pd.read_csv('./10WaterLevel_3166.csv')         ##
    #salinity=pd.read_csv('./10Salinity_3166.csv')      ## when salinity is high, Fe+++ is tend to be negative or model is not converge
    #wtl=pd.read_csv('./10WaterLevel_2825.csv')         ##
    #salinity=pd.read_csv('./10Salinity_2825.csv')      ## when salinity is high, Fe+++ is tend to be negative or model is not converge

#c    abmsl=2.59 #units in meter
#c    alpha=[163,154.6,176.1,153.8,37.4,30.8,37.2,19.2]    #phase angle in degrees M2,S2,N2,K2,K1,O1,P1,Q1
#c    aM2=[0.013,0.007,0.005,0.002,0.114,0.114,0.036,0.025]       #amplitude in meters M2,S2,N2,K2,K1,O1,P1,Q1
#c    w=[28.984104,30.0,28.43973,30.082138,15.041069,13.943035,14.958931,13.398661] #angual speed degree/hour
#c    tide_t=range(0,nsteps)        #t in hours
#c    Zs=2.59+0.1   #marsh elevation above mean sea level in meters
    Zs=0.0
    Zt=[]
    sf=[]       #salinity in tide
    Tw=[]
#c    for i in tide_t:
#c        sz=abmsl
#c        for j in range(0,len(w)):
#c            z=aM2[j]*math.cos(w[j]*math.pi/(180*dt)*i*dt+alpha[j])
#c            sz=sz+z
#c        Zt+=[sz,]
##create salinity time series in tide, this is just an arbitrary pattern for testing
#c    for i in tide_t:
#c        sal = 36*math.cos(w[1]*i+alpha[1])
#c        sf+=[abs(sal),]
    Zt=wtl['WaterLevel'].values.tolist()
    sf=salinity['Salinity'].values.tolist()
    Tw=wtemp['Temperature'].values.tolist()
##replace NAN value with mean value
    rep_wtl=np.nanmean(Zt)
    rep_sal=np.nanmean(sf)
    rep_temp=np.nanmean(Tw)

    Zt=[rep_wtl if math.isnan(x) else x for x in Zt]
    sf=[rep_sal if math.isnan(x) else x for x in sf]
##subset data to match rainfall data, here subset data to year 2009. subset 2011-2013.
    Zt=[x for x in Zt[24*365*2:24*(365*4+366)-1]]
    sf=[x for x in sf[24*365*2:24*(365*4+366)-1]]
    Tw=[x for x in Tw[24*365*2:24*(365*4+366)-1]]

    if len(Zt) < nsteps//repyr:               #run twice of the observation data.
        Zt=Zt+[rep_wtl]*(nsteps//repyr-len(Zt))
    if len(sf) < nsteps//repyr:               #run twice of the observation data.
        sf=sf+[rep_sal]*(nsteps//repyr-len(sf))
    if len(Tw) < nsteps//repyr:               #run twice of the observation data.
        Tw=Tw+[rep_sal]*(nsteps//repyr-len(Tw))

    Zt=Zt*repyr
    sf=sf*repyr
    Tw=Tw*repyr
    #Zt=[i for i in Zt]
    ##uniform salinity
    #sf=[5+x*0 for x in sf]
#cjw create a tidal forcing concentration of compounds in tide for water-sediment interface boundary
##tconc_frac is the fraction of each primary species relative to salinity in water.
#    rstcl=0.14   #unit is mg/ml:mg/ml
    salfc=0.00180665
#ppt to sulfate concetration
#    sulfc=salfc/0.14 ## sulfate concetration is mg/L, and 1e-3*mg/L/molar_mass=M
    tconc_fc={
    #'DOM1':0.01,
    #'Acetate-':0.01,
    'HCO3-':0,#0.01,
    'O2(aq)':0,
    'NO3-':0,#1.29046e-10,
    'NH4+':0,#1.29046e-10,
    'SO4--':0.14*1e-3/salfc, #1.29046e-2,
    'Cl-':1e-3/salfc,
    'Na+':0.5769*1e-3/salfc,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':0.024*1e-3/salfc,
    'HS-':1e-9*1e-3/salfc,
    'Fe+++':0,#1e-15,
    'Fe++':0,
    'CH4(aq)':0,#0.1,
    'H2(aq)':0,#0.01,
    'N2(aq)':0,#0.01,
    'N2O(aq)':0,#0.001,
    'H+':0,#0.01,         #this could be related to the alklinity in porewater and tide
    #'Tracer':0,#0.01,
    }

    molar_mass={#'SOM':180.156,
    #'HRimm':59.052,
    #'DOM1':180.156,
    #'Acetate-':59.052,
    'HCO3-':61.0168,
    'O2(aq)':31.998,
    'NO3-':62.0049,
    'NH4+':18.039,
    'SO4--':96.06,
    'Cl-':35.454,
    'Na+':22.98977,
    'Ca++':40.0780,
    'HS-':33.1,
    'Fe+++':55.845,
    'Fe++':55.845,
    'CH4(aq)':16.04,
    'H2(aq)':2,
    'N2(aq)':28.0134,
    'N2O(aq)':44.013,
    'H+':1,         #this could be related to the alklinity in porewater and tide
    #'Tracer':1,
    }   ##molar mass is g/mol
    tide_conc={'salt':sf}
    for spec in tconc_fc.keys():
        if spec=='O2(aq)':
            tmp=[(i*tconc_fc[spec]+9*1e-3)/molar_mass[spec] for i in sf]
        elif spec=='NO3-':
            tmp=[i*0.0+45*1e-6 for i in sf]          #Upreti et al., 2020 et al., 45 umol/L
        elif spec=='Fe+++':
            tmp=[i*0.0+0.1*1e-6 for i in sf]         #Telfeyan et al., 2017, surface water at Myrtle Grove 0.1-0.7 umol/Kg -> 0.1-0.7 umol/L taking water density as 1024 kg/m3
        #elif spec=='Fe++':
        #    tmp=[i*0.0+0*1e-15 for i in sf]
        elif spec=='H+':
            tmp=[i*0.0+1*1e-7 for i in sf]
        elif spec=='HCO3-':
            tmp=[i*0.0+initial_HCO3 for i in sf]
        #elif spec=='N2(aq)':
        #    tmp=[i*0.0+0.00068 for i in sf]
        #elif spec=='CH4(aq)':
        #    tmp=[i*0.0+1e-10 for i in sf]
        #elif spec=='H2(aq)':
        #    tmp=[i*0.0+4.325e-10 for i in sf]
        else:
            tmp=[i*tconc_fc[spec]/molar_mass[spec] for i in sf]
        tide_conc[spec]=tmp

    #    soillyr=[l.volume for l in layers]
    zbio=0.1                #maximum bioturbation depth in meter below which bioturbation decrease exponentially
    thick=[l.volume for l in layers]
    mid_thick=[i/2 for i in thick]
    zdpth=np.array([l.volume for l in layers]).cumsum()
    mid_dpth=list(map(operator.sub,zdpth,mid_thick))

    #wInf=(3.3*10**(-0.87478367-0.00043512*zdpth[len(layers)-1]))*0.01/(3600*24*365)   ## advection rate or sedimentation rate in cm/yr -> m/s and zdpth in meter is the depth of the bottom layer
    porInf=layers[len(layers)-1].porosity-0.05#0.9#poravg[stnm][5]-0.05#0.9               ## porosity at the bottom layer or deepest depth of the column
    #print(porInf)
    por0=layers[0].porosity #0.95#poravg[stnm][0]#0.95

    coeffp=4*0.01                 ## unit in cm, coeffecient for exponential porosity change -> change unit to m
    #por=porInf+(por0-porInf)*np.exp([i*(-1/coeffp) for i in mid_dpth])  #porosity within layers
    por=[l.porosity for l in layers]
    #por=poravg[stnm]
    #for i in range(len(layers)-1):
    #    protmp=layers[i].porosity
    #    print(protmp,por[i])

    pori=porInf+(por0-porInf)*np.exp([i*(-1/coeffp) for i in zdpth])  #porosity at interfaces
    #print(pori)
    alpha=np.array([0]+[(thick[i]+thick[i+1])/2 for i in range(len(thick)-1)])   ##[0] create a 0 value for alpha to have an extra column or element ie without [0], alpha will have 3 element, with it will have 4 element

    Temp=0                   ## temperature in degree.
    coeffDb=1*0.01                # coefficient for exponential bioturbation decrease unit in cm ->m
    sig=0.3*1e-3             ## diffusion boundary layer thickness in meter (which means it is 0.3 mm)
    #Db0=15*(wInf*(365*24*3600)/0.01)**0.6*(1e-4/(24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/yr -> m^2/s) above zbio

    #ad=0.0336                     ## ion-specific coefficient with unit of cm2/d/degree C for diffusion coeff, ad is different among substrates -> change unit to m2/s/degree C
    adFc=0.0001/(24*3600)
    ad={'DOM1':0.0336,
    'Acetate-':0.0336,
    'HCO3-':0.0336,
    'O2(aq)':0.0386,
    'NO3-':0.0336,
    'NH4+':0.0336,
    'SO4--':0.0336, #1.29046e-2,
    'Cl-':0.0336,
    'Na+':0.0336,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':0.0336,
    'HS-':0.0242,
    'Fe+++':0.0336,
    'Fe++':0.0336,
    'CH4(aq)':0.0336,
    'H2(aq)':0.0336,
    'N2(aq)':0.0336,
    'N2O(aq)':0.0336,
    'H+':0.0336,         #this could be related to the alklinity in porewater and tide
    'Tracer':0.0336,
    }
    for nm in ad.keys():
        ad[nm]=ad[nm]*adFc
    #D0=1e-9                  ## zero-degree coefficient m2/s, change among substrates
    DFc=1.1574e-9             ##diffusion convert factor from cm2/day to m2/s; reference from Greenblatt and Rodney J. Sobey, 2001, Journal of Coastal Research
    D0={'DOM1':0.2,
    'Acetate-':0.2,
    'HCO3-':0.415589,
    'O2(aq)':0.955,
    'NO3-':0.845,
    'NH4+':0.847,
    'SO4--':0.432, #1.29046e-2,
    'Cl-':0.432,
    'Na+':0.432,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':0.432,
    'HS-':0.842,
    'Fe+++':0.432,
    'Fe++':0.432,
    'CH4(aq)':0.955,
    'H2(aq)':0.955,
    'N2(aq)':0.955,
    'N2O(aq)':0.955,
    'H+':0.415589,         #this could be related to the alklinity in porewater and tide
    'Tracer':0.415589,
    }
    for nm in ad.keys():
        D0[nm]=D0[nm]*DFc


    z=numpy.array([0]+[l.volume for l in layers]).cumsum()
    z_mid=(z[:-1]+z[1:])/2
#    for l in range(len(layers)):
#cjw            layers[l].total_immobile['Root_biomass']=root_biomass_top*numpy.exp(-z_mid[l]/root_efolding)/(root_Cfrac*12)*100**3
#cjw update following diffusion coefficient
#        layers[l].diffquo={'O2(aq)':0.01**((l+1)*0.75)}
        # Treat top layer as O horizon with less root biomass and mineral Mn
#cjw        layers[0].total_immobile['Root_biomass']=root_biomass_top/(root_Cfrac*12)*100**3*1e-3
        # layers[0].mineral_volume_fraction['Mn(OH)2(am)']=1e-7
#cjw        layers[0].mineral_volume_fraction['Birnessite2']=1e-7
#    layers[0].surface_site_density['>DOM1']=1e2

    for l in layers:
        l.write_output(0,dt)

    flow_in=numpy.zeros(len(layers),dtype=float)
    flow_out=numpy.zeros(len(layers),dtype=float)

    t0=time.time()
    tprev=t0

    tstart=0

        # Restart from existing state?
    restart_state=None
    if isinstance(restart_state,str):
        restart_state=xarray.open_dataset(restart_state)
    if isinstance(restart_state,xarray.Dataset):
        copy_to_layers(restart_state, layers)
        tstart=restart_state.dropna(dim='time')['time'].isel(time=-1).item()
    elif isinstance(restart_state,list) and isinstance(restart_state[0],layer):
        layers=restart_state

#cjw tmp for bc state to be reset as the way it is at dry condition; set dumpy layers for flux diffusion
    #bc_tmp=bc_state
    satur_initial=[l.saturation for l in layers]   #layer thickness
    gas_species={'O2(aq)':initial_O2,#0.00027235027,#0.000282, O2 partial pressure/kH=0.20950/769.23
    'CH4(aq)':2.37998572e-09,#0.00175,#1e-9, partial pressure/ henry constant=1.7*1e-6/714.29 under standard condition T=25, Wania R. et al., 2010
    'HCO3-':initial_HCO3,#1.32e-5,
    'HS-':0.0,
    'H2(aq)':4.325e-10,#4.325e-10,   partial pressure/ henry constant=0.78/1282.05
    'N2(aq)':0.00068,#0.000714, partial pressure/ henry constant=0.78/1282.05
    'N2O(aq)':0.0,}#16.4*1e-9,}    ## air equilibrium concentration CO2 is 0.00039/29.41=1.32607956e-05

    gas_frac={'O2(aq)':0.209,#
    'CH4(aq)':1.7*1e-6,# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'HCO3-':0.0004,
    'HS-':0.0,
    'H2(aq)':5e-7,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'N2(aq)':0.786,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'N2O(aq)':0.0,}    ## bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    #methane embulltion
    R=8.3144621  ## gas constant in unit of m3Pa/K/mol
    T0=273.15+25    ## kalvin temp K
    Tstd=273.15
    Patm=101325     ##unit is Pa
    #met_frac=0.4    #
    #ebS=0.05708-0.001545*T+0.00002069*T*T
    #ebV=ebS*layers[0].volume*layers[0].porosity*layers[0].saturation
    #press=Patm+ryo*g*H     #g in 9.81 m/s2; h is 1 m, water density in 1000 kg/m3, pressure in 9810 Pa
    #n_mol=Patm*ebV/(R*T)
    #Ebmet=n_mol*0.4/(layers[0].volume*layers[0].porosity*layers[0].saturation*1000)   ##mol/L
    #print(Ebmet)
    immobile_species=['SOMC','SOMN']#['SOM']
    Fimb={'SOMC':4e-9,
    'SOMN':4e-9*0.04}    #C/N ratio is 25;4e-9 #immobile species flux at sediment water interface unit is molC/m3. SOM to SOC is SOM*0.56;
    ### fresh 376/12=9.9e-7 mol-C/m2/s or kgram unit 0.376/0.56=0.6714 kg/m2/yr (2e-8 kg/m2/s); brackish 9.116e-7 mol-C/m2/s,0.345/0.56=0.616 kg m2/yr (2e-8 kg/m2/s); saline 1.149e-6 mol-c/m2/s, 435/0.56=0.776 kg m2/yr (2e-8); reference from Suir et al., 2019, Delaung et al., 2012. #mol/m3-bluk
#cjw define diffusion dict
    #Dt={}
    #Dsed={}
    #Dsedi={}
    ##cjw define dict for flux
    LtranR=numpy.zeros((len(D0),len(layers),nsteps))
    SOMtranR=numpy.zeros((len(layers),nsteps))
    Ebull=numpy.zeros((len(D0),len(layers),nsteps))
    gas_conc={}
##flux at sediment water surface
    Jswi=numpy.zeros((len(D0),nsteps))
    Ebflag=1   #choose ebullition method, flag=0 means concentration threshold, flag=1 means pressure threshold
#cjw above it is transport    

#cjw biology starts below
    micro_rate=numpy.zeros((len(rateconstants),nsteps))
    Bgrow=numpy.zeros((len(Mar.microbes),nsteps))
    grw_tmp=numpy.zeros((len(Mar.microbes),len(Mar.microbes[0].genes))) 
    C_subs=['SOM','DOM1','Acetate-','HCO3-','CH4(aq)']   ## carbon sources and also as energy source/electron donor for microbes to grow
    #mdeath=0.16/86400     #unit is per second, death rate
    #mcdeath=0.64/86400
    #rateconstants_tmp={}   #store temporary rateconstants for each reaction with microbes names
    rateconstants_bio={}

    success=True
    for step in range(nsteps):  
        #for rectnm in rateconstants.keys():
        #    nmtmp=rectnm.split(" ")
        #    newnm=' '.join(nmtmp[:-1])
        #    rateconstants_tmp[newnm]=rateconstants[rectnm]*0
        mic=True    
#cjw compute biomass growth
        if mic:
        #rateconst=rateconstants
            for num in range(len(layers)):
                micdx=-1
                for m in Mar.microbes:
                    micdx=micdx+1
                    #grw_tmp=numpy.zeros(len(Mar.microbes),len(m.genes))   #growth contribution from each gene effort
                    growth,gnks=Mar.growth_rate(m)   
                    B_old=layers[num].total_immobile[m.name]
                    #nelim=1.0   #energy and nutrient limitation term
                    #ihlim=1.0
                    subs_nm=[]
##compute maximum carbon limitation term and pass it to energy only limitation term
                    for idx in range(len(m.genes)):
                        substrates=m.genes[idx]['reactant_pools']      ##dict of substrates with stoic microbes work on/uptake
                        cmplx=Mar.gnpnlty[m.genes[idx]['name']]
                        subs_nm=subs_nm+list(substrates.keys())
                        C_tmp=list(set(subs_nm).intersection(C_subs))
                        celim=0.0                                      #carbon + energy limitation term for non-carbon reactions
                        for spec in C_tmp:
                            if spec == 'SOM':
                                soilc='SOMC'
                                #soiln='SOMN'  #cjw soilN change ???
                                c_tmp=layers[num].total_immobile[soilc]/(layers[num].total_immobile[soilc]+1)   #nutrient/energy limitation for growth; organic carbon species are both nutrient and engergy sources sources.
                            else:
                                c_tmp=layers[num].total_mobile[spec]/(layers[num].total_mobile[spec]+gnks[spec])
                            celim=max(c_tmp,celim)
                    for idx in range(len(m.genes)):
                        nelim=1.0   #energy and nutrient limitation term
                        ihlim=1.0
                        substrates=m.genes[idx]['reactant_pools']      ##dict of substrates with stoic microbes work on/uptake
                        cmplx=Mar.gnpnlty[m.genes[idx]['name']]
                        C_tmp=list(set(list(substrates.keys())).intersection(C_subs))
                        for spec in substrates.keys():
                            if spec == 'SOM':
                                soilc='SOMC'
                                #soiln='SOMN'  #cjw soilN change ???
                                tmp_lim=layers[num].total_immobile[soilc]/(layers[num].total_immobile[soilc]+1*cmplx)   #nutrient/energy limitation for growth; organic carbon species are both nutrient and engergy sources sources.
                            elif spec == 'H2O':
                                tmp_lim=1
                            else:
                                tmp_lim=layers[num].total_mobile[spec]/(layers[num].total_mobile[spec]+gnks[spec]*cmplx)
                                #celim=tmp_lim
                            nelim=nelim*tmp_lim
                        if 'inhibition_terms' in m.genes[idx].keys():
                            inhterm=m.genes[idx]['inhibition_terms']      # a list of dictionary of inhibition informariton
                            for inh in range(len(inhterm)):
                                inhspec=inhterm[inh]['species']
                                tmp_inh=gnks[inhspec]/(layers[num].total_mobile[inhspec]+gnks[inhspec]*cmplx)
                                ihlim=ihlim*tmp_inh
## carbon limitation for energy only reaction impact on growth
                        if len(C_tmp) == 0:                            ## no carbon source in the gene reaction, the reaction provides only energy
                            nelim=celim*nelim                                  #must use carbon gene function to grow
                        else:
                            nelim=nelim                                #can grow without other functional genes                                                
                        grw_tmp[micdx,idx]=growth*nelim*ihlim
                        #print(growth,nelim,ihlim,m.name,idx)
                    grwconst=numpy.amax(grw_tmp[micdx,:])
                    #print(numpy.amax(grw_tmp[2,:]))
                    Bgrow[micdx,step]=grwconst*B_old
                    B_new=Bgrow[micdx,step]*bio_dt*skip_tm+B_old
                    #print(B_new,grwconst,m.name)
            ##cjw biology growing: compute the growth of biomass for each microbes
                    layers[num].total_immobile[m.name]=B_new
                    #if m.name=='microbe3':
                    #    print(B_new,layers[num].total_immobile[m.name],m.name,grwconst,B_old)
            ##compute rateconstants for chemical reaction
                    for idx in range(len(m.genes)):  
                        microbe_yield = Mar.yield_calculation(m.genes[idx]['name'],m.genes[idx]['reactant_pools'],m.genes[idx]['product_pools'],Mar.Gib_std)                           
                        #rate_tmp=growth*B_old/microbe_yield                ##final rate constant is the sum of all microbes growth on the reaction   
                        gfrac=grw_tmp[micdx,idx]/numpy.sum(grw_tmp[micdx,:])
                        #rate_tmp=grwconst*B_old/microbe_yield
                        rate_tmp=gfrac*grwconst*B_old/microbe_yield
                        rectnm=m.genes[idx]['name']+f' ({m.name})'
                        if rectnm in rateconstants.keys():
                            rateconstants[rectnm]=rate_tmp
                #print('start run with microbes',rate_tmp,rectnm)
                #for rectnm in rateconstants.keys():
                #    if rateconstants[rectnm] == 0:
                #        rateconstants[rectnm] = 1e-7
                #print(rateconstants)
                #rateconstants_bio=rateconstants.copy()   
                #rateconstants_microbe=convert_rateconstants(rateconstants_bio,reaction_network,precision=precision)
                #print(rateconstants.keys())
                for name in rateconstants:
                    rxnidx=decomp_network.get_stoich_from_name(name,Mar.microbe_reactions,precision)
                    rateconstants_bio[rxnidx]=rateconstants[name]
                for nmidx in rateconstants.keys():
                    ridx=list(rateconstants.keys()).index(nmidx)
                    micro_rate[ridx,step]=rateconstants[nmidx]
                #print('rate_network_bio',rateconstants_microbe)
#cjw define diffusion dict
        Dt={}
        Dsed={}
        Dsedi={}
        ##cjw define dict for flux
        Rc={}                     #flux between layers
        Ebgas={}
        dzw=Zt[step]-Zs
        #dzw=abs(dzw)
        #Temp=Tw[step]
        T=Tw[step]+273.15

        Hcp={'O2(aq)':1.2*1e-5*np.exp(-1700*(1/T-1/T0)),
        'CH4(aq)':1.4*1e-5*np.exp(-1900*(1/T-1/T0)),# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        #'HCO3-':3.3*1e-4*np.exp(-2400*(1/T-1/T0)),
        'HS-':1.0*1e-3*np.exp(-2100*(1/T-1/T0)),
        'H2(aq)':7.7*1e-6*np.exp(-530*(1/T-1/T0)),#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2(aq)':6.4*1e-6*np.exp(-1600*(1/T-1/T0)),#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2O(aq)':2.4*1e-4*np.exp(-2700*(1/T-1/T0)),
        }
#cjw diffusion changes with water depth
        if dzw >0:
            Fir=max(1,15.9*dzw**(-0.1))                     ## change -0.43 to -0.1 to reduce Fir. depth is in meter, diffusion enhancement factor to represent bio irrigation impact within layers
            wInf=(3.3*10**(-0.87478367-0.00043512*dzw))*0.01/(3600*24*365)   ## advection rate might be too high with 1.4e-10, which is 4mm per year #advection rate or sedimentation rate in cm/yr -> m/s and zdpth in meter is the depth of the bottom layer
            Db0=5.2*10**(0.76241122-0.00039724*dzw)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s
            #wInf=wInf*0.5    ##reduce advection by half
            #for num in range(len(layers)):
            #    layers[num].saturation=1.0                  #saturated with water when it is flooded
        else:
            Fir=1                                            ## diffusion enhancement factor caused by bio irrigation at layer interface
            wInf=3.9e-11
            Db0=5.2*10**(0.76241122-0.00039724*0.0)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s, when dzw <=0 or when wetland is under dry condition
            #wInf=wInf*0.5    ##reduce advection by half
            #for num in range(len(layers)):
            #    if mid_dpth[num]+dzw >0:
            #        layers[num].saturation=1.0
            #    else:
            #        layers[num].saturation=satur_initial[num]   # use the same inital saturation level and unsaturated with water when it is not flooded
        #Db0=15*(wInf*(365*24*3600)/0.01)**0.6/365*(1e-4/(24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/d -> m^2/s) above zbio
        #Db0=5.2*10**(0.76241122-0.00039724*dzw)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s
        Dbi=[Db0 if zdpth[i] <= zbio else Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(zdpth))]
        for spec in D0.keys():
            Dt[spec] = D0[spec]*0.01 + ad[spec]*Temp       ## free solution diffusion coefficient at ambient temperature T
            tmp=[Dt[spec]*por[i]**2*Fir+Db0 if mid_dpth[i] <= zbio else Dt[spec]*por[i]**2*Fir+Db0*np.exp(-1*(mid_dpth[i]-zbio)/coeffDb) for i in range(len(por))]
            tmpi=[Dt[spec]*pori[i]**2*Fir+Db0 if zdpth[i] <= zbio else Dt[spec]*pori[i]**2*Fir+Db0*np.exp(-1*(zdpth[i]-zbio)/coeffDb) for i in range(len(pori))]
            Dsed[spec]=tmp                            ## enhanced difussion coeff within layer with/without bioturbation
            Dsedi[spec]=tmpi                          ## enhanced diffusion coeff at layer interface below bioturbation depth

#cjw compute diffusion for each timpstep for each layer
        for l in range(len(layers)):
            for spec in layers[l].primary_names:
                if spec in Dsed.keys():
                    layers[l].diffquo[spec]=Dsed[spec][l]
            #print('Dsedi %2.20f at layer l %2d'%(Dsedi['SO4--'][l],l))
        for num in range(len(layers)):
            if num < len(layers)-1:
                dz=(thick[num+1]+thick[num])/2
            if num > 0:
                dz2=(thick[num]+thick[num-1])/2
            #layers[num].porosity=por[num]
            if num == 0:
                if dzw >0 :
                    if bc is not None:
                        primarynames=layers[num].primary_names
                        for spec in primarynames:
                            if spec in D0.keys():
                                #Ebflux=0.0
                                if spec in tconc_fc.keys():            ##the upper boundary of sediment column or SWI is imposed to be the tide_conc.
                                    Cw_tmp=tide_conc[spec][step]              ##concentration in overlaying water supposed is uniform within water column, otherwise this value should be bottomwater concentration
                                    #Cswi=(layers[num].total_mobile[spec]-Cw_tmp)*sig/(sig+thick[0]/2)+Cw_tmp    ##concentration at sediment-water interface
                                    Cswi=Cw_tmp
                                    #Jflx=layers[num].porosity**2*Dsedi[spec][num]*(Cw_tmp-layers[num].total_mobile[spec])*1e3/(layers[0].volume)   #(lyrdpth[0]+lyrdpth[1]) mol/m2/s ; Jflx positive means downward flux or flux into soil and negative means upward flux outside of soil
                                else:
                                    #Jflx=0
                                    Cswi=layers[num].total_mobile[spec]
##cjw ebullition for gas
                                #cje compute flux from underlie layer to top layer
                                rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                                rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                                rtmp3=por[num]*Dsed[spec][num]*(layers[num].total_mobile[spec]-Cswi)/(layers[num].porosity*(thick[num]/2+sig)*thick[num])
                                rtmp4=wInf*porInf*(0.5*Cswi+(1-0.5)*layers[num].total_mobile[spec])/(layers[num].porosity*thick[num])
                                Rtmp=rtmp1-rtmp2-rtmp3+rtmp4#-Ebflux
                                #print(rtmp1,rtmp3,por[num],Dsed[spec][num],wInf)
                                #if spec =='NO3-':
                                #    print('Rtmp of NO3 %2.30f'%Rtmp)
                                ##cjw SWI flux -- first order fick's law
                                Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num]#+Ebflux*1000/thick[num]                         ##unit is mol/m2/s
                                Jidx=list(D0.keys()).index(spec)
                                Jswi[Jidx,step]=Jflx
                                #Ebull[num,step]=Ebflux*1000/thick[num]
                                #print('before Ebflux %2.9f,%2.9f at surface wet',Jswi[Jidx,step],Ebflux*1000/thick[num] )
                                #print('Jflx %2.30f'%Rtmp)
                            elif spec in immobile_species:
                                rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                                rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                                rtmp3=Fimb[spec]*1/((1-layers[num].porosity)*thick[num])    #when flooded, immobile flux reduce to 0.1*Fimb
                                rtmp4=0.0
                                Rtmp=rtmp3-rtmp4+rtmp1-rtmp2
                            else:
                                Rtmp=0.0
                            if spec not in Rc.keys():
                                Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                            else:
                                Rc[spec].append(Rtmp)
                                #print('Jflx %2.30f'%Rtmp)
                elif dzw <=0 :
                    #bc_state=bc_tmp
                    if bc is not None:
                        primarynames=layers[num].primary_names
                        for spec in primarynames:
                            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)  ##this should work because data is copied to alquimia#is empty now
                            if spec in D0.keys():
                                #Ebflux=0.0
                                if spec in gas_species.keys():
                                    if spec in Hcp.keys():
                                        Cw_tmp=Hcp[spec]*Patm*gas_frac[spec]*1e-3#*Tstd/(T*1000)   ##mol/L#gas_species[spec]
                                        #print(Cw_tmp,spec,initial_HCO3)
                                    else:
                                        Cw_tmp=gas_species[spec]
                                        #print(Hcp[spec],Cw_tmp)
                                    Cswi=Cw_tmp
                                else:
                                    Cswi=layers[num].total_mobile[spec]
                                    #Ebflux=0.0
                                rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                                rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                                rtmp3=por[num]*Dsed[spec][num]*100*(layers[num].total_mobile[spec]-Cswi)/(layers[num].porosity*thick[num]*(sig+thick[num]/2))  ##surface o2 diffusion increase 2 order (100)
                                rtmp4=wInf*porInf*(0.5*Cswi+(1-0.5)*layers[num].total_mobile[spec])/(layers[num].porosity*thick[num])
                                if spec == 'HCO3-':
                                    dCO2=rtmp3-rtmp4
                                    layers[num].total_immobile['H+']=layers[num].total_immobile['H+']-dCO2*1000*layers[num].porosity*layers[num].saturation*dt
                                if spec == 'HS-':
                                    dH2S=rtmp3-rtmp4
                                    layers[num].total_immobile['H+']=layers[num].total_immobile['H+']-dH2S*1000*layers[num].porosity*layers[num].saturation*dt
                                Rtmp=rtmp1-rtmp2-rtmp3+rtmp4#-Ebflux

                                Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num] #+ Ebflux*1000/thick[num]     ##positive means flux is upward, negative means flux is in soil

                                Jidx=list(D0.keys()).index(spec)

                                Jswi[Jidx,step]=Jflx

                            elif spec in immobile_species:
                                rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                                rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                                rtmp3=Fimb[spec]/((1-layers[num].porosity)*thick[num])
                                rtmp4=0.0
                                Rtmp=rtmp3+rtmp4+rtmp1-rtmp2
                            else:
                                Rtmp=0.0
                            if spec not in Rc.keys():
                                Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                            else:
                                Rc[spec].append(Rtmp)
            elif num < len(layers)-1 and num >0:
                #layers[num].copy_to_alquimia(data)
                #layers[num+1].porosity=por[num+1]
                primarynames=layers[num].primary_names
                for spec in primarynames:
                    if spec in D0.keys():
                        rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                        rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        rtmp3=pori[num-1]*Dsedi[spec][num-1]*(layers[num].total_mobile[spec]-layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num]*dz2)
                        rtmp4=wInf*porInf*(alpha[num]*layers[num-1].total_mobile[spec]+(1-alpha[num])*layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                        Jidx=list(D0.keys()).index(spec)
                        Jswi[Jidx,step]=Jswi[Jidx,step]#+Ebflux*1000/thick[num]
                    elif spec in immobile_species:
                        rtmp1=(1-pori[num])*Dbi[num]*(layers[num+1].total_immobile[spec]-layers[num].total_immobile[spec])/((1-layers[num].porosity)*dz*thick[num])
                        rtmp2=wInf*(1-porInf)*(alpha[num+1]*layers[num].total_immobile[spec]+(1-alpha[num+1])*layers[num+1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                        rtmp3=(1-pori[num-1])*Dbi[num-1]*(layers[num].total_immobile[spec]-layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num]*dz2)
                        rtmp4=wInf*(1-porInf)*(alpha[num]*layers[num-1].total_immobile[spec]+(1-alpha[num])*layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                    else:
                        Rtmp=0.0
                    if spec not in Rc.keys():
                        Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                    else:
                        Rc[spec].append(Rtmp)
            elif num==len(layers)-1:
                #layers[num].copy_to_alquimia(data)
                primarynames=layers[num].primary_names
                for spec in primarynames:
                    Cbtm=layers[num].total_mobile[spec]      ##bottom boundary layer for mobile concentration
                    Sbtm=layers[num].total_immobile[spec]    ## bottom boundary for immobile concentration
                    if spec in D0.keys():
                        #Rtmp=0.0
                        rtmp1=0.0
                        rtmp2=wInf*porInf*(alpha[num]*layers[num].total_mobile[spec]+(1-alpha[num])*Cbtm)/(layers[num].porosity*thick[num])
                        rtmp3=pori[num-1]*Dsedi[spec][num-1]*(layers[num].total_mobile[spec]-layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num]*dz2)
                        rtmp4=wInf*porInf*(alpha[num]*layers[num-1].total_mobile[spec]+(1-alpha[num])*layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                        Jidx=list(D0.keys()).index(spec)
                        Jswi[Jidx,step]=Jswi[Jidx,step]#+Ebflux*0*1000/thick[num]            ##unit is mol/m2/s

                    elif spec in immobile_species:
                        rtmp1=0.0
                        rtmp2=wInf*(1-porInf)*(alpha[num]*layers[num].total_immobile[spec]+(1-alpha[num])*Cbtm)/((1-layers[num].porosity)*thick[num])
                        rtmp3=(1-pori[num-1])*Dbi[num-1]*(layers[num].total_immobile[spec]-layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num]*dz2)
                        rtmp4=wInf*(1-porInf)*(alpha[num]*layers[num-1].total_immobile[spec]+(1-alpha[num])*layers[num-1].total_immobile[spec])/((1-layers[num].porosity)*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                    else:
                        Rtmp=0.0
                    if spec not in Rc.keys():
                        Rc[spec]=[Rtmp]     ##suppose the flux happens at one unit m2 surface
                    else:
                        Rc[spec].append(Rtmp)
# Next advect mobile species. This assumes we can separate advective transport and reactions at this time step length
            # Alternate strategy would be to pass these rates to run_onestep and spread them over the variable timesteps        
        #dCO2=layers[0].total_mobile['HCO3-']-initial_HCO3
        #layers[0].total_mobile['HCO3-']=layers[0].total_mobile['HCO3-']-dCO2
        #layers[0].total_immobile['H+']=layers[0].total_immobile['H+']-dCO2*1000*layers[0].porosity*layers[0].saturation
        try:
#cjw                for n,l in enumerate(layers+[incubation_layer]):
            for n,l in enumerate(layers):
                l.copy_to_alquimia(data)
                ##cjw: update rateconstant each timestep to represent microbial activities.
                #rateconstants_microbe=rateconst.copy()
                rateconstants_microbe=rateconstants_bio.copy()   
                #print('rateconst before alquimia',rateconstants_bio)
                #rateconstants_microbe=convert_rateconstants(rateconstants_bio,reaction_network,precision=precision)
                #print('rateconst in alquimia',rateconstants_microbe)
                l.rateconstants=rateconstants_microbe
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
            #cjw get Transport Rate or flux for run_onestep
                Rc_lyr={}
                for spec in layers[0].primary_names:
                    if spec in Rc.keys():
                        Rc_lyr[spec]=Rc[spec][n]
                        if spec in D0.keys():
                            Ridx=list(D0.keys()).index(spec)
                            LtranR[Ridx,n,step]=Rc[spec][n]
                            #print('Rc_lyr %2.30f for %s at layer %d'%(LtranR[Ridx,n,step],spec,n))
                        SOMtranR[n,step]=Rc['SOMC'][n]
                    else:
                        Rc_lyr[spec]=0.0

            ##cjw ebullition limit for each layer
                if mid_dpth[n]+dzw >0:
                    press=Patm+9.81*1000*(mid_dpth[n]+dzw)     #g in 9.81 m/s2; h is 1 m, water density in 1000 kg/m3, pressure in 9810 Pa
                else:
                    press=Patm
                #for gidx in gas_species.keys():
                #    if gas_frac[gidx]>0.0:
                #        gas_conc[gidx]=layers[num].total_mobile[gidx]
                Eb_tmp=Ebullition(l.total_mobile,l.volume,l.porosity,l.saturation,press,Tw[step],dt,Flag=Ebflag)

                num_cuts,Ebrate=l.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=l.rateconstants,TRc=Rc_lyr,Ebc=Eb_tmp)

                for spec in D0.keys():
                    Ridx=list(D0.keys()).index(spec)
                    Ebull[Ridx,n,step]=Ebrate[spec]*l.volume*l.porosity*l.saturation*1000      ##mol/L to mol/m3 to mol

                l.copy_from_alquimia(data)
                # Write output
                l.write_output(step+1,dt,num_cuts)
                    #layers[lyrcnt].write_output(step+1,dt,num_cuts)
                    # print('O2 after: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
                    #print('flux exchange %1.9f'%(excflx[8]))
        except RuntimeError as err:
            print('ERROR on timestep %d, layer %d: %s'%(step,n,err))
            print('Returning output so far')
            l.write_output(step,dt,num_cuts)
            success=False
            break
        except KeyboardInterrupt as err:
            print('INTERRUPTED on timestep %d'%(step))
            print('Returning output so far')
            l.write_output(step,dt,num_cuts)
            success=False
            break


        if step%100==0 and step>0:
            t1=time.time()
            cuts=[l.output['ncuts'][step-100:step].mean() for l in layers]
            mean_dt=numpy.mean([l.output['actual_dt'][step-100:step].mean() for l in layers])
            print('*** Step {step:d} of {nsteps:d} ({nyears:1.1f} of {totalyears:d} years). Time elapsed: {t:1.1f} min ({tperstep:1.1f} s per {steplength:1.1f} hour timestep). Mean cuts: {meancuts:s} Mean dt: {meandt:1.1f} s ***'.format(
                    step=step,nsteps=nsteps,t=(t1-t0)/60,tperstep=(t1-tprev)/25,meancuts=str(cuts),meandt=mean_dt,steplength=dt/3600,nyears=step*dt/(3600*24*365),totalyears=int(nsteps*dt/(3600*24*365))))
            tprev=t1


#cjw not sure what is in output
    output=convert_to_xarray(layers,t0=tstart)
    if isinstance(restart_state,xarray.Dataset):
        output=xarray.concat([restart_state,output.isel(time=slice(1,None))],dim='time')


    import datetime
    today=datetime.datetime.today()
#cjw        output.to_netcdf('Mn_output/Mn_pH{ph:1.1f}_Ndep{Ndep:03d}_warming{warming:d}_{year:04d}-{month:02d}-{day:02d}.nc'.format(ph=pH,Ndep=int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600) )),year=today.year,month=today.month,day=today.day,warming=warming))
#cjw        convert_to_xarray([incubation_layer],leaf_Mn=leaf_Mn_concs).to_netcdf('Mn_output/Mn_incubations_pH{ph:1.1f}_Ndep{Ndep:03d}_warming{warming:d}_{year:04d}-{month:02d}-{day:02d}.nc'.format(ph=pH,Ndep=int(Ndep/(1000/molar_mass['N']/100**2/(365*24*3600) )),year=today.year,month=today.month,day=today.day,warming=warming))


print('\n\n\n Simulation finished. Total time: %1.1f minutes\n'%((time.time()-starting_time)/60))

#cjw molar_volume_manganite = 24.45 # cm3/mol
#cjw molar_volume_MnOH2am = 22.3600
#cjw molar_volume_birnessite = 251.1700


def plot_output(output,axs,subsample=1,do_legend=True,**kwargs):
    for num in range(len(output.depth)):
        out=output.isel(depth=num,time=slice(None,None,subsample)).dropna(dim='time')
        t=out.time/(1)
        porosity=out['Porosity']
        saturation=out['saturation']
        BD=out['BD']
        #axs[num,0].plot(subt,out['Total DOM1'][std:etd]*1e6,c='C2',label='DOM',**kwargs)
        axs[num,0].plot(out['Total Sorbed microbe1'],label='microbe$_1$',c='C0',**kwargs)
        axs[num,1].plot(out['Total Sorbed microbe2'],label='microbe$_2$',c='C5',**kwargs)
        axs[num,2].plot(out['Total Sorbed microbe3'],label='microbe$_3$',c='C1',**kwargs)
        axs[num,0].set_ylabel('\u03BCM')


        ###cjw methane saturation conc in water is 1746 umol/L
        
    axs[-1,0].set_xlabel('Time (second)')
    if do_legend:
        axs[0,0].legend()
        axs[0,1].legend()
        axs[0,2].legend()
out=output.isel(depth=0,time=slice(None,None,1)).dropna(dim='time')
#pyplot.plot(out['Total Sorbed microbe1'])
microbe,(ax1,ax2,ax3,ax4,ax5)=pyplot.subplots(5, 1,figsize=(8,6))
ax1.plot(out['Total Sorbed microbe6'],label='microbe$_1$',linewidth=2)
ax2.plot(out['Total Sorbed microbe7'],label='microbe$_2$',linewidth=2)
ax3.plot(out['Total Sorbed microbe8'],label='microbe$_3$',linewidth=2)
ax4.plot(out['Total Sorbed microbe9'],label='SOMC',linewidth=2)
ax5.plot(out['Total Sorbed microbe10'],label='SOMC',linewidth=2)
pyplot.savefig('microbes.pdf')

geochem,(ax1,ax2,ax3,ax4)=pyplot.subplots(4, 1,figsize=(8,6))
ax1.plot(out['Total CH4(aq)'],label='Fe+++',linewidth=2)
ax2.plot(out['Total Acetate-'],label='Acetate-',linewidth=2)
ax3.plot(out['Total O2(aq)'],label='DOM1',linewidth=2)
ax4.plot(out['Total SO4--'],label='SO4--',linewidth=2)
pyplot.savefig('geochem.pdf')

#ratechem,rx=pyplot.subplots(1, 1,figsize=(8,6))
#rx.plot(micro_rate[:,:],linewidth=2)
#pyplot.savefig('ratechem.pdf')

pyplot.tight_layout()
pyplot.show()
#cjw networkfig=pyplot.figure('Reaction network',clear=True)
#cjw drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid'],
#c        font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
#c        namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'HCO$_3^-$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
#c                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',})
