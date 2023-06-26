#Created by Ben Sulman for the Manganese soil profile, and revised by Jiaze Wang to simulate delta marsh soil profile,02/18/2021

from run_alquimia import get_alquimiavector,ffi,lib,check_status,init_alquimia,convert_condition_to_alquimia,print_metadata
import decomp_network
import MRDwet_RedoxNet as Mar       #delmar_network_tidev3 is methanogenesis without inhibition
from matplotlib import pyplot
import numpy
import xarray
import pdb
#cjw
import READSOIL_CRMS

class layer:
    def __init__(self,volume,saturation=0.9,temperature=20.0,water_density=1000.0,porosity=0.25,pressure=101325.0,BD=1.5,CEC=None,rateconstants={},diffquo={}):
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
        self.CEC=CEC # meq/kg       #cjw-- remove CEC
#cjw above function defines some essential attributes for each instance
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
            'porosity':numpy.ma.masked_all(nsteps,dtype=float),
#            'CEC H+':numpy.ma.masked_all(nsteps,dtype=float),             #cjw-- remove CEC
            #'flow_in':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
            #'flow_out':numpy.ma.masked_all((nsteps,len(self.total_mobile)),dtype=float),
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
#        self.output['CEC H+'][step]=self.total_immobile['H+']-(self.surface_site_density['>Carboxylate-']*self.mineral_volume_fraction['Rock(s)']-self.aux_doubles[len(self.total_mobile)*2+len(self.secondary_free_ion_concentration)+1+self.surface_site_names.index('>Carboxylate-')])  ##cjw--remove CEC
        #self.output['flow_in'][step,:]=numpy.array([self.flow_in.get(name,0.0) for name in self.primary_names])
        #self.output['flow_out'][step,:]=numpy.array([self.flow_out.get(name,0.0) for name in self.primary_names])

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

        #flow_in=pandas.DataFrame(self.output['flow_in'],columns=self.primary_names,index=self.output['time']).add_suffix(' inflow')
        #for col in flow_in.columns:
        #    output_DF[col]=flow_in[col]
        #output_units.update([(s+' inflow','mol/m2/sec') for s in self.primary_names])
        #flow_out=pandas.DataFrame(self.output['flow_out'],columns=self.primary_names,index=self.output['time']).add_suffix(' outflow')
        #for col in flow_out.columns:
        #    output_DF[col]=flow_out[col]
        #output_units.update([(s+' outflow','mol/m2/sec') for s in self.primary_names])

#        output_DF['CEC H+']=pandas.DataFrame(self.output['CEC H+'],index=self.output['time'])
#        output_units['CEC H+']='mol/m^3'

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
        #print(data.properties.saturation)
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]

        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            if data.state.total_mobile.data[pos] < truncate_concentration and data.state.total_mobile.data[pos] != 0.0:
                for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                    if '->' in reactname and spec in reactname.split('->')[0].split():  # reactants
                        data.properties.aqueous_kinetic_rate_cnst.data[num]=0.0

#cjw update truncate_conc for mineral precipitation and dissolution
        #mineral_nm=get_alquimiavector(data.meta_data.mineral_names)
        #mineral_sp=['Pyrite']
        #pos1=get_alquimiavector(data.meta_data.primary_names).index('SO4--')
        #pos2=get_alquimiavector(data.meta_data.primary_names).index('HS-')
        #pos3=get_alquimiavector(data.meta_data.primary_names).index('Fe++')
        #thresh_sulf=1e-6      #turn off mineral precipitation when sulfate is smaller than 1e-5 M.
        #for num in range(data.meta_data.mineral_names.size):
            #if mineral_nm[num] in mineral_sp:
                #data.properties.mineral_rate_cnst.data[num]=rateconstants[mineral_nm[num]]
        #if data.state.total_mobile.data[pos1] < thresh_sulf or data.state.total_mobile.data[pos2] < thresh_sulf or data.state.total_mobile.data[pos3] < 1e-10:
            #for num in range(data.meta_data.mineral_names.size):
                #if mineral_nm[num] in mineral_sp:
                    #data.properties.mineral_rate_cnst.data[num]=0.0

        #mineral_nm=get_alquimiavector(data.meta_data.mineral_names)
        #mineral_sp=['Pyrrhotite']
        #pos=get_alquimiavector(data.meta_data.primary_names).index('HS-')
        #thresh=1e-6
        #for num in range(data.meta_data.mineral_names.size):
        #    if mineral_nm[num] in mineral_sp:
        #        data.properties.mineral_rate_cnst.data[num]=rateconstants[mineral_nm[num]]
        #if data.state.total_mobile.data[pos] <= thresh:
        #    for num in range(data.meta_data.mineral_names.size):
        #        if mineral_nm[num] in mineral_sp:
        #            data.properties.mineral_rate_cnst.data[num]=0.0

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
                    if spec != 'SOM':
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


##cjw remove leaf_Mn
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

rate_scale=1e-9
truncate_conc=1e-30
thresh=truncate_conc*1.01

reaction_network=Mar.make_network() # We will add the Mn along with leaf litter manually instead of generating it through decomposition

rateconstants={
'1.00e+00 DOM1  -> 3.33e-01 Acetate-  + 3.33e-01 HCO3-  + 6.67e-01 H+  + 1.33e+00 H2(aq)  + 3.33e-01 Tracer':rate_scale*1e1,
'1.00e+00 DOM1  + 1.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 Tracer':rate_scale*1e1,
'4.00e+00 H2(aq)  + 1.00e+00 HCO3-  + 1.00e+00 H+  -> 1.00e+00 CH4(aq)  + 3.00e+00 H2O':rate_scale,#*0.31,#1.157,#0.31,#5.26,#,#0.00295,#,       #Delaune et al., 1983[fresh marsh 440 mgC/m2/day, rate constant 1.157e-10; brackish marsh 200mgC/m2/d,rate constant 5.26e-10,salt marsh,11.8 mgC/m2/day, 3.10e-11 brackish open water 13 mgC/m2/d, rate constant 2.95e-13, salt water 3.6 mgC/m2/d rate cons 8.18e-14, fresh water 37 mg C/m2/day 8.17e-11]
#'1.25e+00 Acetate-  + 2.00e+00 NO3-  + 7.50e-01 H+  -> 2.50e+00 HCO3-  + 1.00e+00 N2(aq)  + 1.00e+00 H2O  + 2.50e+00 Tracer':rate_scale*0.0036,#0.00407,#0.00652,#   ##Kanchan 2020 paper, barataria bay lake direct rate at 20 degree convert from flux to rate constant unit.[6.52e-13 marsh, 3.6e-13 lake3166, 4.07e-13channel3169]
#'2.25e+00 Acetate-  + 4.00e+00 NO3-  + 1.75e+00 H+  -> 4.50e+00 HCO3-  + 1.00e+00 N2(aq)  + 1.00e+00 N2O(aq)  + 2.00e+00 H2O  + 4.50e+00 Tracer':rate_scale*0.0036,
'1.00e+00 Acetate-  + 2.00e+00 NO3-  + 1.00e+00 H+  -> 2.00e+00 HCO3-  + 1.00e+00 N2O(aq)  + 1.00e+00 H2O  + 2.00e+00 Tracer':rate_scale*0.0036,
'2.50e-01 Acetate-  + 1.00e+00 N2O(aq)  + 1.00e+00 H+  -> 5.00e-01 HCO3-  + 1.00e+00 N2(aq)  + 2.50e-01 H+  + 5.00e-01 Tracer':rate_scale*0.0036,
'1.00e+00 Acetate-  + 2.00e+00 O2(aq)  -> 2.00e+00 HCO3-  + 2.00e+00 H+  + 2.00e+00 Tracer':rate_scale*1e2,
'1.00e+00 CH4(aq)  + 2.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 H2O  + 1.00e+00 Tracer':rate_scale*2.31,    ## Roslev and King, 1996
'2.00e+00 NH4+  + 4.00e+00 O2(aq)  -> 2.00e+00 NO3-  + 2.00e+00 H2O  + 4.00e+00 H+':rate_scale*0.23, #0.358,          ## oxy uptake rate from Kanchan 2020 [3.58e-11 marsh, 2.3e-11 lake, 2.3e-11channel]
#'1.00e+00 Acetate- + 1.00e+00 NO3- + 1.00e+00 H2O + 1.00e+00 H+ -> 1.00e+00 NH4+ + 2.00e+00 HCO3- + 2.00e+00 Tracer':rate_scale,
#'1.00e+00 DOM1 + 5.00e-01 NO3- + 5.00e-01 H2O <-> 5.00e-01 NH4+ + 1.00e+00 HCO3-':rate_scale,
'1.00e+00 DOM1  + 5.00e-01 NO3-  + 5.00e-01 H2O  -> 5.00e-01 NH4+  + 1.00e+00 HCO3-':rate_scale,
#'1.00e+00 HS- + 1.00e+00 NO3- + 1.00e+00 H+ + 1.00e+00 H2O <-> 1.00e+00 SO4-- + 1.00e+00 NH4+':rate_scale*0.1,
#'1.00e+00 HS- + 1.00e+00 NO3- + 1.00e+00 H+ + 1.00e+00 H2O -> 1.00e+00 SO4-- + 1.00e+00 NH4+':rate_scale*0.1,
'1.00e+00 HS-  + 1.00e+00 NO3-  + 1.00e+00 H+  + 1.00e+00 H2O  -> 1.00e+00 SO4--  + 1.00e+00 NH4+':rate_scale*0.0+1e-11,
#'5.00e-01 HS-  + 1.00e+00 NO3-  + 5.50e+00 H+  -> 5.00e-01 SO4--  + 1.00e+00 NH4+  + 1.00e+00 H2O':rate_scale*0.0+1e-11,
#'5.00e-01 HS-  + 1.00e+00 NO3-  + 5.50e+00 H+  -> 5.00e-01 SO4--  + 1.00e+00 NH4+  + 1.00e+00 H2O':rate_scale*0.0+1e-11,
#'1.00e+00 HS-  + 2.00e+00 O2(aq)  -> 1.00e+00 SO4--  + 1.00e+00 H+':rate_scale,
'1.00e+00 Acetate-  + 8.00e+00 Fe+++  -> 2.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 2.00e+00 Tracer':rate_scale,   #from Schoepfer et al., 2014, JGR-B in scale of 1e-8 mol/L/s
'1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  <-> 1.00e+00 Fe+++  + 5.00e-01 H2O':rate_scale*0,
'1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  -> 1.00e+00 Fe+++  + 5.00e-01 H2O':1e-9,#rate_scale*100,#*10000,
'1.00e+00 CH4(aq)  + 8.00e+00 Fe+++  + 3.00e+00 H2O  -> 1.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 1.00e+00 Tracer':rate_scale*0.1,    ##Roslev and King, 1996
'5.00e+00 CH4(aq)  + 8.00e+00 NO3-  + 3.00e+00 H+  -> 5.00e+00 HCO3-  + 4.00e+00 N2(aq)  + 9.00e+00 H2O  + 1.00e+00 Tracer':rate_scale*0.1,
#'1.00e+00 CH4(aq)  + 1.00e+00 NO3-  + 1.00e+00 H+  -> 1.00e+00 HCO3-  + 1.00e+00 NH4+  + 1.00e+00 Tracer':rate_scale*0.1,
'1.00e+00 Acetate-  + 1.00e+00 SO4--  -> 2.00e+00 HCO3-  + 1.00e+00 HS-  + 2.00e+00 Tracer':rate_scale,   #Abdul M. Al-Raei et al., 2009, Ocean Dynamics 1e-8, rate_scale*1e2; sulfate reduction rate range from 1.157e-08 to 1.157e-12 mol/L/s tabel 1 and figure 1 from Pester et al., 2012, frontiers in microbiology; Kanchan suggests 6.9e-10 - 3.3e-9
'1.00e+00 CH4(aq)  + 1.00e+00 SO4--  -> 1.00e+00 HCO3-  + 1.00e+00 HS-  + 1.00e+00 H2O  + 1.00e+00 Tracer':rate_scale*0.1, #Abdul M. Al-Raei et al., 2009, Ocean Dynamics
'1.00e+00 Acetate-  + 1.00e+00 H2O  -> 1.00e+00 CH4(aq)  + 1.00e+00 HCO3-  + 1.00e+00 Tracer':rate_scale,#0.00295,#100,#,
'SOM decay to CO2 (SOMDEC sandbox)': 1e-6,#1e-6,
'SOM decay to DOM1 (SOMDEC sandbox)':1e-7,#1e-7,
#'Pyrite':1e-8,
#'Pyrrhotite':1e-11,
}

#cjw
input_file='deltamarsh.in'

decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=False,truncate_concentration=truncate_conc)

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
#cjw getting soil properties for all stations
bavg,poravg,satavg=READSOIL_CRMS.soils()
#wind station name
wndid='Gisl'        #wind data from grand isle
#for station in stnm:
nm=['CRMS4245','CRMS3166','CRMS2825','CRMS0220','CRMS0224']
for stnm in nm:#bavg.keys():
    chem,data,sizes,status=init_alquimia(input_file,hands_off=False)
    
#cjw assume it is always saturated even when water table is low (precipitation could bring water or river diversion)
    layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=bavg[stnm][0],porosity=poravg[stnm][0],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][1],porosity=poravg[stnm][1],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][2],porosity=poravg[stnm][2],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][3],porosity=poravg[stnm][3],saturation=1)]+[layer(0.05,rateconstants=rateconstants_warmed,BD=bavg[stnm][4],porosity=poravg[stnm][4],saturation=1)]+[layer(0.1,rateconstants=rateconstants_warmed,BD=bavg[stnm][5],porosity=poravg[stnm][5],saturation=1)]
    for l in layers:
        l.secondary_names=secondary_names
    for l in layers:
        l.initcond=decomp_network.change_constraints(Mar.pools,{'Ca++':'1e-30',
                                                        'Na+':'1e-30','Cl-':'1e-30','CH4(aq)':'2.37998572e-09',
                                                        'SO4--':'1e-30',#})     #'%1.8f'%((0.75-suldpth[sulidx])*1e-15/30)})
                                                        'O2(aq)':'0.2 G O2(g)'})               #cjw: set the sulfate profile so4 decrease with depth (0.75-suldpth[sulidx])/30)
       
    pools_tide=Mar.pools.copy()
    for n,p in enumerate(pools_tide):
        if p['name']=='O2(aq)':
            pools_tide[n]=pools_tide[n].copy()
            pools_tide[n].update(constraints={'initial':'0.2 G O2(g)'})
    bc=layers[0].initcond

    dt=3600
#cjw 10-yr run set up, change to the following setting for 10-year run including leap year and non-leap year
    repyr=3
    yrchain=[2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
    nyears=len(yrchain)*repyr      ## frequency of rerun the entire observation dataset 2009,2010,2011,2012(),2013,2014,2015,2016(),2017,2018,2019
    ##nsteps=365*24//(dt//3600)*nyears
    nsteps=0
    for i in yrchain:
        if i%4==0:
            yday=366
        else:
            yday=365
        nsteps=nsteps+yday*24//(dt//3600)
    nsteps=nsteps*repyr
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
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=l.rateconstants[reactname]

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

    min_dt=0.1
    truncate_concentration=1e-20

#cjw update wet-dry condition for dynamic boundary conditions
##if water level > marsh elevation, then bc=sed_water_interface, else bc=sed_air_interface
###cjw compute tide water level below to
    import scipy as sp
    import mpmath as mp
    import numpy as np
    import math
    import pandas as pd
    import operator
    import seaborn as sns

## cjw set flag for ebullition scheme flag=0 means concentration threshold, flag=1 means pressure threshold
    Ebul_flag=1   #choose ebullition method, flag=0 means concentration threshold, flag=1 means pressure threshold
##cjw parameter needed for transport
    tide_t=range(0,nsteps)
    Nlyr=range(0,len(layers))
    tdt=3600
    ##cjw path for hydrological and meterological hourly data.
    hypath='./forcing/'
    hyfile1='10WaterLevel_'+stnm[4:8]+'.csv'
    hyfile2='10Salinity_'+stnm[4:8]+'.csv'
    hyfile3='10temp_'+stnm[4:8]+'.csv'
    hyfile4='10wind_'+wndid+'.csv'
    wtl=pd.read_csv(hypath+hyfile1)         ##
    salinity=pd.read_csv(hypath+hyfile2)
    wtemp=pd.read_csv(hypath+hyfile3)
    wtemp=wtemp.interpolate()
    wnd=pd.read_csv(hypath+hyfile4)

    Zs=0.0
    Zt=[]
    sf=[]       #salinity in tide
    Tw=[]
    wdspd=[]
    Zt=wtl['WaterLevel'].values.tolist()
    sf=salinity['Salinity'].values.tolist()
    Tw=wtemp['Temperature'].values.tolist()
    wdspd=wnd['wind_spd'].values.tolist()
##replace NAN value with mean value
    rep_wtl=np.nanmean(Zt)
    rep_sal=np.nanmean(sf)
    rep_temp=np.nanmean(Tw)
    rep_wnd=np.nanmean(wdspd)   # already clean the data so no need to redo it again

    Zt=[rep_wtl if math.isnan(x) else x for x in Zt]
    sf=[rep_sal if math.isnan(x) else x for x in sf]
    Tw=[rep_temp if math.isnan(x) else x for x in Tw]
    #wdspd=[rep_wnd if math.isnan(x) else x for x in wdspd]

    if len(Zt) < nsteps//repyr:               #run twice of the observation data.
        Zt=Zt+[rep_wtl]*(nsteps//repyr-len(Zt))
    if len(sf) < nsteps//repyr:               #run twice of the observation data.
        sf=sf+[rep_sal]*(nsteps//repyr-len(sf))
    if len(Tw) < nsteps//repyr:               #run twice of the observation data.
        Tw=Tw+[rep_sal]*(nsteps//repyr-len(Tw))
    if len(wdspd) < nsteps//repyr:               #run twice of the observation data.
        wdspd=wdspd+[rep_wnd]*(nsteps//repyr-len(wdspd))

    Zt_tmp=Zt
    sf_tmp=sf
    Tw_tmp=Tw
    wdspd_tmp=wdspd

    Zt=Zt*(repyr-1)
    sf=sf*(repyr-1)
    Tw=Tw*(repyr-1)
    wdspd=wdspd*(repyr-1)

    ## set senario runs for flooding, the relative sea level rising rate is 10mm/yr in MR (Jankowski et al., 2017); assume elevated mean sea level does not impact hydrodynamic is the same but
    ## elevated water level after 20 yrs will be 20*10mm=200mm=20cm=0.2m; after 10 yrs will be 0.1m, after 40 yrs will be 0.4m
    MSL=0.0      ## elevated mean sea level in meters after 20 yrs with 10mm/yr rising rate; 10 yrs with 10mm/yr
    dsal=0.0       ## elevated salinity increase or elevated salt water intrusion
    Zt_tmp=[x+MSL for x in Zt_tmp]
    sf_tmp=[x+dsal for x in sf_tmp]
    Zt=Zt+Zt_tmp
    sf=sf+sf_tmp
    Tw=Tw+Tw_tmp
    wdspd=wdspd+wdspd_tmp

#cjw create a tidal forcing concentration of compounds in tide for water-sediment interface boundary
##tconc_frac is the fraction of each primary species relative to salinity in water.
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
            #tmp=[i*0.0+167.8188*1e-6 for i in sf]
            tmp=[i*0.0+initial_HCO3 for i in sf]
        #elif spec=='N2(aq)':
        #    tmp=[i*0.0+0.00068 for i in sf]
        elif spec=='CH4(aq)':
            tmp=[i*0.0+0.022*1e-6 for i in sf]          #0.022 uM is from Brooks, 1975; 7.48*1e-6M is from Songjie He et al., water measurement in Barataria Bay.
        #elif spec=='H2(aq)':
        #    tmp=[i*0.0+4.325e-10 for i in sf]
        else:
            tmp=[i*tconc_fc[spec]/molar_mass[spec] for i in sf]
        tide_conc[spec]=tmp                            ##all gas species are set to be equilibrium with atmosphere, so the above setting is changed in the below code

    zbio=0.1                #maximum bioturbation depth in meter below which bioturbation decrease exponentially
    thick=[l.volume for l in layers]
    mid_thick=[i/2 for i in thick]
    zdpth=np.array([l.volume for l in layers]).cumsum()
    mid_dpth=list(map(operator.sub,zdpth,mid_thick))

    #wInf=(3.3*10**(-0.87478367-0.00043512*zdpth[len(layers)-1]))*0.01/(3600*24*365)   ## advection rate or sedimentation rate in cm/yr -> m/s and zdpth in meter is the depth of the bottom layer
    porInf=layers[len(layers)-1].porosity-0.05#0.9#poravg[stnm][5]-0.05#0.9               ## porosity at the bottom layer or deepest depth of the column

    por0=layers[0].porosity #0.95#poravg[stnm][0]#0.95

    coeffp=4*0.01                 ## unit in cm, coeffecient for exponential porosity change -> change unit to m
    #por=porInf+(por0-porInf)*np.exp([i*(-1/coeffp) for i in mid_dpth])  #porosity within layers
    por=[l.porosity for l in layers]
    #por=poravg[stnm]
    #for i in range(len(layers)-1):
    #    protmp=layers[i].porosity
    #    print(protmp,por[i])

    pori=porInf+(por0-porInf)*np.exp([i*(-1/coeffp) for i in zdpth])  #porosity at interfaces

    alpha=np.array([0]+[(thick[i]+thick[i+1])/2 for i in range(len(thick)-1)])   ##[0] create a 0 value for alpha to have an extra column or element ie without [0], alpha will have 3 element, with it will have 4 element

    Temp=0                   ## temperature in degree.
    coeffDb=1*0.01                # coefficient for exponential bioturbation decrease unit in cm ->m
    sig=0.3*1e-3             ## diffusion boundary layer thickness in meter (which means it is 0.3 mm)
    #Db0=15*(wInf*(365*24*3600)/0.01)**0.6*(1e-4/(24*3600))     ## bioturbation diffusion coeff in unit m^2/s (cm^2/yr -> m^2/s) above zbio

    #ad=0.0336                     ## ion-specific coefficient with unit of cm2/d/degree C for diffusion coeff, ad is different among substrates -> change unit to m2/s/degree C
    adFc=0.00*0.0001/(24*3600)     ##0.01 is arbitrary chosen to make model works
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
### cjw Zt is tide water level
##update bc state by different layer property
    z=numpy.array([0]+[l.volume for l in layers]).cumsum()
    z_mid=(z[:-1]+z[1:])/2

    for l in layers:
        l.write_output(0,dt)

    t0=time.time()
    tprev=t0

    tstart=0

    restart_state=None
    if isinstance(restart_state,str):
        restart_state=xarray.open_dataset(restart_state)
    if isinstance(restart_state,xarray.Dataset):
        copy_to_layers(restart_state, layers)
        tstart=restart_state.dropna(dim='time')['time'].isel(time=-1).item()
    elif isinstance(restart_state,list) and isinstance(restart_state[0],layer):
        layers=restart_state

#cjw tmp for bc state to be reset as the way it is at dry condition; set dummy layers for flux diffusion
    satur_initial=[l.saturation for l in layers]   #layer thickness
    gas_species={'O2(aq)':initial_O2,#0.00027235027,#0.000282, O2 partial pressure/kH=0.20950/769.23
    'CH4(aq)':2.37998572e-09,#0.00175,#1e-9, partial pressure/ henry constant=1.7*1e-6/714.29 under standard condition T=25, Wania R. et al., 2010
    'HCO3-':initial_HCO3,#1.32e-5,
    'HS-':0.0,
    'H2(aq)':4.325e-10,#4.325e-10,   partial pressure/ henry constant=0.78/1282.05
    'N2(aq)':0.00068,#0.000714, partial pressure/ henry constant=0.78/1282.05
    'N2O(aq)':0.0,}#4.325e-10; 16.4*1e-9,}    ## air equilibrium concentration CO2 is 0.00039/29.41=1.32607956e-05

    gas_frac={'O2(aq)':0.209,#
    'CH4(aq)':1.7*1e-6,# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'HCO3-':0.0004,
    'HS-':0.0,
    'H2(aq)':5e-7,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'N2(aq)':0.786,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    'N2O(aq)':0*3e-7,}    ## bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
    #methane embulltion
    R=8.3144621  ## gas constant in unit of m3Pa/K/mol
    T0=273.15+25    ## kalvin temp K
    Tstd=273.15
    Patm=101325     ##unit is Pa

    immobile_species=['SOM']#['SOM']
    Fimb={'SOM':4e-9}    #4e-9 #immobile species flux at sediment water interface unit is molC/m3. SOM to SOC is SOM*0.56;
    ### fresh 376/12=9.9e-7 mol-C/m2/s or kgram unit 0.376/0.56=0.6714 kg/m2/yr (2e-8 kg/m2/s); brackish 9.116e-7 mol-C/m2/s,0.345/0.56=0.616 kg m2/yr (2e-8 kg/m2/s); saline 1.149e-6 mol-c/m2/s, 435/0.56=0.776 kg m2/yr (2e-8); reference from Suir et al., 2019, Delaung et al., 2012. #mol/m3-bluk

    ##cjw define dict for flux
    LtranR=numpy.zeros((len(D0),len(layers),nsteps))
    SOMtranR=numpy.zeros((len(layers),nsteps))
    Ebull=numpy.zeros((len(D0),len(layers),nsteps))
    gas_conc={}
##flux at sediment water surface
    Jswi=numpy.zeros((len(D0),nsteps))
    Ebflag=Ebul_flag   #choose ebullition method, flag=0 means concentration threshold, flag=1 means pressure threshold
#cjw question about the unit of total_mobile (mol/L?) and immobile species (immobile seems to have mol/m3-bulk)
    success=True
    for step in range(nsteps):
        #cjw define diffusion dict
        Dt={}
        Dsed={}
        Dsedi={}
        ##cjw define dict for flux
        Rc={}                     #flux between layers
        Ebgas={}
        dzw=Zt[step]-Zs
        T=Tw[step]+273.15

        Hcp={'O2(aq)':1.2*1e-5*np.exp(-1700*(1/T-1/T0)),
        'CH4(aq)':1.4*1e-5*np.exp(-1900*(1/T-1/T0)),# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        #'HCO3-':3.3*1e-4*np.exp(-2400*(1/T-1/T0)),
        'HS-':1.0*1e-3*np.exp(-2100*(1/T-1/T0)),
        'H2(aq)':7.7*1e-6*np.exp(-530*(1/T-1/T0)),#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2(aq)':6.4*1e-6*np.exp(-1600*(1/T-1/T0)),#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2O(aq)':2.4*1e-4*np.exp(-2700*(1/T-1/T0)),
        }
        ## Schmidt number for each gas in freshwater
        Sc_f={'O2(aq)':1800.6-120.1*Tw[step]+3.7818*Tw[step]**2-0.047608*Tw[step]**3,
        'CH4(aq)':1897.8-114.28*Tw[step]+3.2902*Tw[step]**2-0.039061*Tw[step]**3,# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'HCO3-':1911.1-118.11*Tw[step]+3.4527*Tw[step]**2-0.041320*Tw[step]**3,
        'HS-':2055.6-137.11*Tw[step]+4.3173*Tw[step]**2-0.054350*Tw[step]**3, ##same as N2O
        'H2(aq)':377.09-19.154*Tw[step]+0.50137*Tw[step]**2-0.005669*Tw[step]**3,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2(aq)':1970.7-131.45*Tw[step]+4.1390*Tw[step]**2-0.052106*Tw[step]**3,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2O(aq)':2055.6-137.11*Tw[step]+4.3173*Tw[step]**2-0.054350*Tw[step]**3,
        }
        ##Schmidt number for each gas in sea water
        Sc_s={'O2(aq)':1953.4-128.0*Tw[step]+3.9918*Tw[step]**2-0.050091*Tw[step]**3,
        'CH4(aq)':2039.2-120.31*Tw[step]+3.4209*Tw[step]**2-0.040437*Tw[step]**3,# bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'HCO3-':2073.1-125.62*Tw[step]+3.6276*Tw[step]**2-0.043219*Tw[step]**3,
        'HS-':2301.1-151.1*Tw[step]+4.7364*Tw[step]**2-0.059431*Tw[step]**3,  ##same as N2O
        'H2(aq)':410.14-20.503*Tw[step]+0.53175*Tw[step]**2-0.0060111*Tw[step]**3,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2(aq)':2206.1-144.86*Tw[step]+4.5413*Tw[step]**2-0.056988*Tw[step]**3,#bublle gas fraction, Tokida et al., 2005; Kellner et al., 2006
        'N2O(aq)':2301.1-151.1*Tw[step]+4.7364*Tw[step]**2-0.059431*Tw[step]**3,
        }
#cjw diffusion changes with water depth
        if dzw >0:
            Fir=max(1,15.9*dzw**(-0.1))                     ## change -0.43 to -0.1 to reduce Fir. depth is in meter, diffusion enhancement factor to represent bio irrigation impact within layers
            wInf=(3.3*10**(-0.87478367-0.00043512*dzw))*0.01/(3600*24*365)   ## advection rate might be too high with 1.4e-10, which is 4mm per year #advection rate or sedimentation rate in cm/yr -> m/s and zdpth in meter is the depth of the bottom layer
            drycoef=0.01                ##arbitory chosen to make model without negative
            Db0=5.2*10**(0.76241122-0.00039724*dzw)*(0.0001/(365*24*3600))             ##unit is cm2/yr -> m2/s
            #wInf=wInf*0.5    ##reduce advection by half
            #for num in range(len(layers)):
            #    layers[num].saturation=1.0                  #saturated with water when it is flooded
        else:
            Fir=1                                            ## diffusion enhancement factor caused by bio irrigation at layer interface
            wInf=3.9e-11
            drycoef=0.01
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
            Dt[spec] = D0[spec]*drycoef + ad[spec]*Tw[step]       ## free solution diffusion coefficient at ambient temperature T
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
                                    if spec in Hcp.keys() and spec != 'CH4(aq)':
                                        Cw_tmp=Hcp[spec]*(Patm+9.81*1000*(dzw/2))*gas_frac[spec]*1e-3#*Tstd/(T*1000)   ##mol/L#gas_species[spec]
                                        ##print(Cw_tmp,spec,initial_HCO3)
                                    else:
                                        Cw_tmp=tide_conc[spec][step]              ##concentration in overlaying water supposed is uniform within water column, otherwise this value should be bottomwater concentration
                                    Cswi=Cw_tmp
                                else:
                                    Cswi=layers[num].total_mobile[spec]
                                #cje compute flux from underlie layer to top layer
                                rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                                rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                                rtmp3=por[num]*Dsed[spec][num]*(layers[num].total_mobile[spec]-Cswi)/(layers[num].porosity*(thick[num]/2+sig)*thick[num])
                                rtmp4=wInf*porInf*(0.5*Cswi+(1-0.5)*layers[num].total_mobile[spec])/(layers[num].porosity*thick[num])
                                Rtmp=rtmp1-rtmp2-rtmp3+rtmp4#-Ebflux
                                ##cjw SWI flux -- first order fick's law
                                Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num]#+Ebflux*1000/thick[num]                         ##unit is mol/m2/s
                                gaswt=Cswi+Jflx/(dzw*1000)
                                if spec in Sc_s.keys():
                                    if spec in Hcp.keys():
                                        gaseq=Hcp[spec]*(Patm+9.81*1000*(dzw/2))*gas_frac[spec]*1e-3#*Tstd/(T*1000)   ##mol/L#gas_species[spec]
                                    else:
                                        gaseq=tide_conc[spec][step]
                                    if sf[step] > 1:
                                        Sc=Sc_s[spec]
                                    else:
                                        Sc=Sc_f[spec]
                                    Jawi=(gaswt-gaseq)*(1000*1e-2/3600)*1.91*numpy.exp(0.35*wdspd[step])*(Sc/600)**(-0.5)
                                    Jflx=Jawi
                                else:
                                    Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num]#+Ebflux*1000/thick[num]                         ##unit is mol/m2/s
                                Jidx=list(D0.keys()).index(spec)
                                Jswi[Jidx,step]=Jflx
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
                elif dzw <=0 :
                    #bc_state=bc_tmp
                    if bc is not None:
                        primarynames=layers[num].primary_names
                        for spec in primarynames:
                            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)  ##this should work because data is copied to alquimia#is empty now
                            if spec in D0.keys():
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
##comment out exposure wind effect, because it does not quite make sense
#                                if spec in Sc_s.keys():
#                                    gassed=Cswi+Jflx*1/(sig*1*1000)
#                                    if spec in Hcp.keys():
#                                        gaseq=Hcp[spec]*Patm*gas_frac[spec]*1e-3#*Tstd/(T*1000)   ##mol/L#gas_species[spec]
#                                    else:
#                                        gaseq=gas_species[spec]
#                                    if sf[step] > 1:
#                                        Sc=Sc_s[spec]
#                                    else:
#                                        Sc=Sc_f[spec]
#                                    #Jasi=(gassed-gaseq)*(1000*1e-2/3600)*(2.07+0.215*wdspd[step]**1.7)*(Sc/600)**(-0.5)
#                                    Jasi=(gassed-gaseq)*(1000*1e-2/3600)*1.91*numpy.exp(0.35*wdspd[step])*(Sc/600)**(-0.5)
#                                    Jflx=Jasi
#                                    Rtmp=rtmp1-rtmp2-Jasi/(layers[num].porosity*thick[num]*1000)+rtmp4
#                                else:
#                                    Jflx=(rtmp3*1000-rtmp4*1000)*layers[num].porosity*thick[num] #+ Ebflux*1000/thick[num]     ##positive means flux is upward, negative means flux is in soil
#                                Jidx=list(D0.keys()).index(spec)
#                                Jswi[Jidx,step]=Jflx
                                #if spec in Sc_s.keys():
                                    #if spec in Hcp.keys():
                                        #gaseq=Hcp[spec]*(Patm+9.81*1000*(dzw/2))*gas_frac[spec]*1e-3#*Tstd/(T*1000)   ##mol/L#gas_species[spec]
                                    #Jawi=(gaswt-gaseq)*(1e-2/3600)*1.91*exp(0.35*wdspd[step])*(Sc_s[spec]/600)**(-0.5)

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
                primarynames=layers[num].primary_names
                for spec in primarynames:
                    if spec in D0.keys():
                        rtmp1=pori[num]*Dsedi[spec][num]*(layers[num+1].total_mobile[spec]-layers[num].total_mobile[spec])/(layers[num].porosity*dz*thick[num])
                        rtmp2=wInf*porInf*(alpha[num+1]*layers[num].total_mobile[spec]+(1-alpha[num+1])*layers[num+1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        rtmp3=pori[num-1]*Dsedi[spec][num-1]*(layers[num].total_mobile[spec]-layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num]*dz2)
                        rtmp4=wInf*porInf*(alpha[num]*layers[num-1].total_mobile[spec]+(1-alpha[num])*layers[num-1].total_mobile[spec])/(layers[num].porosity*thick[num])
                        Rtmp=rtmp1-rtmp2-rtmp3+rtmp4
                        #Jidx=list(D0.keys()).index(spec)
                        #Jswi[Jidx,step]=Jswi[Jidx,step]#+Ebflux*1000/thick[num]
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
                        #Jidx=list(D0.keys()).index(spec)
                        #Jswi[Jidx,step]=Jswi[Jidx,step]#+Ebflux*0*1000/thick[num]            ##unit is mol/m2/s

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
##cjw above is the first try on migrating boundary condition setup into time step loop.
        try:
            for n,l in enumerate(layers):
                l.copy_to_alquimia(data)
                #print(l.rateconstants)
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
                        SOMtranR[n,step]=Rc['SOM'][n]
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

                ##add temperature dependance on rateconstants
                l.rateconstants=rateconstants.copy()

                Qtn=2.0**((Tw[step]-20)/10)     ##ref from Ma et al., 2017, data-constrained projections eqn(9).; and Gadfar et al., 2019, Arrhenius plot analysis, temperature coefficient and Q10 value
                for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                    l.rateconstants[reactname]=rateconstants[reactname]*Qtn
                num_cuts,Ebrate=l.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=l.rateconstants,TRc=Rc_lyr,Ebc=Eb_tmp)

                for spec in D0.keys():
                    Ridx=list(D0.keys()).index(spec)
                    Ebull[Ridx,n,step]=Ebrate[spec]*l.volume*l.porosity*l.saturation*1000      ##mol/L to mol/m3 to mol

                l.copy_from_alquimia(data)
                # Write output
                l.write_output(step+1,dt,num_cuts)

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


#cjw output
    output=convert_to_xarray(layers,t0=tstart)
    if isinstance(restart_state,xarray.Dataset):
        output=xarray.concat([restart_state,output.isel(time=slice(1,None))],dim='time')

    import datetime
    if Ebflag==1:
        ebp='EBP'
    else:
        ebp=''
    today=datetime.datetime.today()
    fname='./Outputs/wnd_temp_Methane_output_{year:04d}-{month:02d}-{day:02d}.nc'.format(year=today.year,month=today.month,day=today.day)
    fname=fname[:-3]+ebp+stnm[4:8]+'.nc'
    Lname='./Outputs/wnd_temp_LtranR_{year:04d}-{month:02d}-{day:02d}.npy'.format(year=today.year,month=today.month,day=today.day)
    Lname=Lname[:-4]+ebp+stnm[4:8]+'.npy'
    Jname='./Outputs/wnd_temp_Jswi_{year:04d}-{month:02d}-{day:02d}.npy'.format(year=today.year,month=today.month,day=today.day)
    Jname=Jname[:-4]+ebp+stnm[4:8]+'.npy'
    Ename='./Outputs/wnd_temp_Ebull_{year:04d}-{month:02d}-{day:02d}.npy'.format(year=today.year,month=today.month,day=today.day)
    Ename=Ename[:-4]+ebp+stnm[4:8]+'.npy'
    output.to_netcdf(path=fname)
    numpy.save(Lname,LtranR)
    numpy.save(Jname,Jswi)
    numpy.save(Ename,Ebull)

print('\n\n\n Simulation finished. Total time: %1.1f minutes\n'%((time.time()-starting_time)/60))
