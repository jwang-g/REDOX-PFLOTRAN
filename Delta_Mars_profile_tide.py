#Created by Ben Sulman for the Manganese soil profile, and revised by Jiaze Wang to simulate delta marsh soil profile,02/18/2021

from run_alquimia import get_alquimiavector,ffi,lib,check_status,init_alquimia,convert_condition_to_alquimia
import decomp_network
import delmar_network_tidev2 as Mar
from matplotlib import pyplot
import numpy
import xarray


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
#        self.output['CEC H+'][step]=self.total_immobile['H+']-(self.surface_site_density['>Carboxylate-']*self.mineral_volume_fraction['Rock(s)']-self.aux_doubles[len(self.total_mobile)*2+len(self.secondary_free_ion_concentration)+1+self.surface_site_names.index('>Carboxylate-')])  ##cjw--remove CEC
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

#        output_DF['CEC H+']=pandas.DataFrame(self.output['CEC H+'],index=self.output['time'])
#        output_units['CEC H+']='mol/m^3'

        self.output_DF=output_DF.reset_index(drop=True).set_index(output_DF.index/(24*3600))
        self.output_units=output_units


    def run_onestep(self,chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=0.15,truncate_concentration=0.0,rateconstants={},min_cuts=0):
        self.copy_to_alquimia(data)
        converged=False
        porosity=data.state.porosity

        max_cuts=num_cuts
        actual_dt=dt/2**num_cuts

        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]

        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            if data.state.total_mobile.data[pos] < truncate_concentration and data.state.total_mobile.data[pos] != 0.0:
                # print('Truncating concentration of {spec:s} from {conc:1.1g}'.format(spec=spec,conc=data.state.total_mobile.data[pos]))
                for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
                    if '->' in reactname and spec in reactname.split('->')[0].split():  # reactants

                        # print('Setting rate of reaction %s to zero'%reactname)
                        data.properties.aqueous_kinetic_rate_cnst.data[num]=0.0

        cut_for_flux=False
        flux=numpy.zeros(data.meta_data.primary_names.size)
        for spec in diffquo.keys():
            pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
            bc_conc=bc.total_mobile.data[pos]
            local_conc=data.state.total_mobile.data[pos]


            # flux[pos] = (bc_conc-local_conc)*diffquo[spec]*actual_dt
            # Run half step of flux first, then other half after chemistry
            flux[pos] = (bc_conc-local_conc)*(1-numpy.exp(-diffquo[spec]*actual_dt/2))

            data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + flux[pos]

            if bc_conc<0:
                raise ValueError('Boundary condition concentration of %s < 0'%spec)
            if local_conc<0:
                raise RuntimeError('Initial concentration of %s < 0'%spec)

            # Here we test if the flux was fast enough to change concentration significantly relative to bc value
            # If the change is more than flux_tol*bc, we cut the time step to resolve changes better
            # if flux[pos] != 0:
            #     print(data.state.total_mobile.data[pos],local_conc,data.state.total_mobile.data[pos] - local_conc,bc_conc)
            if abs(data.state.total_mobile.data[pos] - local_conc) > bc_conc*flux_tol:
                # print('Cutting time step to resolve flux better')
                cut_for_flux=True


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

            # Second half of flux
            flux=numpy.zeros(data.meta_data.primary_names.size)
            for spec in diffquo.keys():
                pos=get_alquimiavector(data.meta_data.primary_names).index(spec)
                bc_conc=bc.total_mobile.data[pos]
                local_conc=data.state.total_mobile.data[pos]


                # flux[pos] = (bc_conc-local_conc)*diffquo[spec]*actual_dt
                # Run half step of flux first, then other half after chemistry
                flux[pos] = (bc_conc-local_conc)*(1-numpy.exp(-diffquo[spec]*actual_dt/2))

                data.state.total_mobile.data[pos] = data.state.total_mobile.data[pos] + flux[pos]

                if bc_conc<0:
                    raise ValueError('Boundary condition concentration of %s < 0'%spec)
                if local_conc<0:
                    raise RuntimeError('Initial concentration of %s < 0'%spec)

                # Here we test if the flux was fast enough to change concentration significantly relative to bc value
                # If the change is more than flux_tol*bc, we cut the time step to resolve changes better
                if abs(data.state.total_mobile.data[pos] - local_conc) > bc_conc*flux_tol:
                    # print('Cutting time step to resolve flux better')
                    cut_for_flux=True
                    converged=False

        if converged:
            self.copy_from_alquimia(data)
            return max_cuts
        else:

            if actual_dt/2<min_dt:
                if cut_for_flux:
                    raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
                else:
                    raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))

            # data will be reset to layer contents at beginning of next call
            # Run it twice, because we cut the time step in half
            ncuts=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants)
            # If this completes, it means that run_onestep was successful
            self.copy_from_alquimia(data)
            if ncuts>max_cuts:
                max_cuts=ncuts

            # This starts from ncuts so it doesn't have to try all the ones that failed again
            for n in range(2**(ncuts-(num_cuts+1))):
                ncuts2=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants)
                self.copy_from_alquimia(data)
                if ncuts2>max_cuts:
                    max_cuts=ncuts2

            return max_cuts


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

rate_scale=2e-10
truncate_conc=1e-30
thresh=truncate_conc*1.01

reaction_network=Mar.make_network() # We will add the Mn along with leaf litter manually instead of generating it through decomposition

rateconstants={
'1.00e+00 DOM1  -> 3.33e-01 Acetate-  + 3.33e-01 HCO3-  + 6.67e-01 H+  + 1.33e+00 H2(aq)  + 3.33e-01 Tracer':rate_scale,
'1.00e+00 DOM1  + 1.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 Tracer':rate_scale,
'4.00e+00 H2(aq)  + 1.00e+00 HCO3-  + 1.00e+00 H+  -> 1.00e+00 CH4(aq)  + 3.00e+00 H2O':rate_scale,
'1.00e+00 Acetate-  + 2.00e+00 NO3-  + 1.00e+00 H+  -> 2.00e+00 HCO3-  + 1.00e+00 N2(aq)  + 0.00e+00 N2O(aq)  + 2.00e+00 H2O  + 2.00e+00 Tracer':rate_scale,
'1.00e+00 Acetate-  + 2.00e+00 O2(aq)  -> 2.00e+00 HCO3-  + 2.00e+00 H+  + 2.00e+00 Tracer':rate_scale,
'1.00e+00 CH4(aq)  + 1.00e+00 O2(aq)  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 H2O  + 1.00e+00 Tracer':rate_scale,
'2.00e+00 NH4+  + 4.00e+00 O2(aq)  -> 2.00e+00 NO3-  + 2.00e+00 H2O  + 4.00e+00 H+':rate_scale,
'1.00e+00 Acetate-  + 8.00e+00 Fe+++  -> 2.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 2.00e+00 Tracer':rate_scale,
# '1.00e+00 Fe++ + 2.50e-01 O2(aq) + 1.00e+00 H+ -> 1.00e+00 Fe+++ + 5.00e-01 H2O':rate_scale,
'1.00e+00 CH4(aq)  + 8.00e+00 Fe+++  + 3.00e+00 H2O  -> 1.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 1.00e+00 Tracer':rate_scale,
'1.00e+00 CH4(aq)  + 1.00e+00 NO3-  -> 1.00e+00 HCO3-  + 1.00e+00 NH4+  + 1.00e+00 Tracer':rate_scale,
'1.00e+00 Acetate-  + 1.00e+00 SO4--  -> 2.00e+00 HCO3-  + 2.00e+00 HS-  + 2.00e+00 Tracer':rate_scale,
'1.00e+00 CH4(aq)  + 1.00e+00 SO4--  -> 1.00e+00 HCO3-  + 1.00e+00 HS-  + 1.00e+00 H2O  + 1.00e+00 Tracer':rate_scale,
'1.00e+00 Acetate-  -> 1.00e+00 CH4(aq)  + 1.00e+00 HCO3-  + 1.00e+00 Tracer':rate_scale,
'SOM decay to CO2 (SOMDEC sandbox)': 1e-6,
'SOM decay to DOM1 (SOMDEC sandbox)':1e-7,
}

##cjw diffusion coefficient
##diffucoef={'O2(aq)':0.001**((l+1)*0.75),'CH4(aq)':0.001**((l+1)*0.75)}
#cjw
input_file='deltamarsh.in'

decomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt',input_file,log_formulation=False,truncate_concentration=truncate_conc)
#cjwdecomp_network.PF_network_writer(reaction_network).write_into_input_deck('SOMdecomp_template.txt','deltamarsh.in',length_days=30,log_formulation=False)

#cjw incubation_length=5 # Years of litter decomp

#cjw molar_mass={'N':14.007}
#cjw  molar_mass={'Mg++':24.305,'Al+++':26.982,'K+':39.098,'Ca++':40.078,'Mn':54.938,'Na+':22.99,'N':14.007}
#cjw  init_exch_cations={'Mg++':1.5,'Al+++':7.0,'K+':1.3,'Ca++':5.0,'Na+':0.2,'Mn++':0.3} # mmol/kg. From Jin et al 2010 Table 3

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
for pH in numpy.arange(6.3,7.3):

    chem,data,sizes,status=init_alquimia(input_file,hands_off=False)
        #cjw below is the codes for layer setting up.
        # Set up layers
        # Top (organic) layer should be thinner and have lower bulk density though
        # Low bulk density causes simulation to slow or crash though. Actually CEC being too low (<100 combined with BD<1) is the problem
#cjw        layers=[layer(0.05,rateconstants=rateconstants_warmed,BD=0.05,porosity=0.5)]+[layer(0.1,BD=0.25,rateconstants=rateconstants_warmed) for num in range(3)]
    layers=[layer(0.01,rateconstants=rateconstants_warmed,BD=0.08,porosity=0.95)]+[layer(0.05,BD=0.1,rateconstants=rateconstants_warmed,porosity=0.95)]+[layer(0.05,BD=0.15,rateconstants=rateconstants_warmed,porosity=0.92)]+[layer(0.05,BD=0.2,porosity=0.9,rateconstants=rateconstants_warmed)]+[layer(0.05,BD=0.25,porosity=0.85,rateconstants=rateconstants_warmed)]+[layer(0.1,BD=0.35,rateconstants=rateconstants_warmed,porosity=0.75)]
    suldpth=[0.01,0.06,0.11,0.16,0.21,0.31]   #depth in meters for sulfur inicond profile
    sulidx=0
    for l in layers:
        l.secondary_names=secondary_names
    for l in layers:
#           l.initcond=Mar.pools.copy()
        l.initcond=decomp_network.change_constraints(Mar.pools,{'H+':'%1.1f P'%pH,
                                                        'Ca++':'1e-15','Na+':'1e-15','Cl-':'1e-15',
                                                        'SO4--':'1e-10'})     #'%1.8f'%((0.75-suldpth[sulidx])*1e-15/30)})
                                                        #'O2(aq)':'0.2 G O2(g)'})               #cjw: set the sulfate profile so4 decrease with depth (0.75-suldpth[sulidx])/30)
        sulidx=sulidx+1
#cjw    for l in layers[0:3]:
#cjw        l.initcond=decomp_network.change_constraints(Mar.pools,{'Fe++':'%1.8f'%(l.volume/100)})      #cjw Fe++ profile
        # initcond=decomp_network.change_site_density(initcond, '>DOM1', 1e4)
#cjw        sed_air_interface=layers[0].initcond
#cjw        sed_water_interface=layers[0].initcond
    pools_tide=Mar.pools.copy()
    for n,p in enumerate(pools_tide):
        if p['name']=='O2(aq)':
            pools_tide[n]=pools_tide[n].copy()
            pools_tide[n].update(constraints={'initial':'0.01 G O2(g)'})
        elif p['name']=='SO4--':
            pools_tide[n]=pools_tide[n].copy()
            pools_tide[n].update(constraints={'initial':'0.0001'})
        elif p['name']=='Ca++':
            pools_tide[n]=pools_tide[n].copy()
            pools_tide[n].update(constraints={'initial':'0.0001'})
    bc=pools_tide
#cjw    bc=None
    dt=3600
    nyears=1
    nsteps=365*24//(dt//3600)*nyears
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
            bc_state.temperature=layers[0].temperature           #replace layers[0]
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
            * Starting simulation with pH = %1.1f, Ndep = %03d, Warming = %d *
            *************************************


        '''%(pH,0,20))        #cjw remove warming



    initial_HCO3 = l.total_mobile['HCO3-']
    initial_O2 = l.total_mobile['O2(aq)']

        # Flow rate cm/s = 10 L/m2/s, positive is downward
#    flow_rate=numpy.linspace(1e-6,1e-7,len(layers)) # Rate declines linearly with depth, assumes removal or accumulation in lower layers
    flow_rate=numpy.zeros(len(layers))*1e-7
#    flow_rate=numpy.linspace(1e-7,1e-8,len(layers))
    min_dt=0.1
    truncate_concentration=1e-20

#cjw update wet-dry condition for dynamic boundary conditions
##if water level > marsh elevation, then bc=sed_water_interface, else bc=sed_air_interface
###cjw compute tide water level below to
    import scipy as sp
    import mpmath as mp
    import math
    abmsl=2.59 #units in meter
    alpha=[163,154.6,176.1,153.8,37.4,30.8,37.2,19.2]    #phase angle in degrees M2,S2,N2,K2,K1,O1,P1,Q1
    aM2=[0.013,0.007,0.005,0.002,0.114,0.114,0.036,0.025]       #amplitude in meters M2,S2,N2,K2,K1,O1,P1,Q1
    w=[28.984104,30.0,28.43973,30.082138,15.041069,13.943035,14.958931,13.398661] #angual speed degree/hour
    tide_t=range(0,nsteps)        #t in hours
    Zs=2.59+0.1   #marsh elevation above mean sea level in meters
    Zt=[]
    sf=[]       #salinity in tide
    for i in tide_t:
        sz=abmsl
        for j in range(0,len(w)):
            z=aM2[j]*math.cos(w[j]*math.pi/(180*dt)*i*dt+alpha[j])
            sz=sz+z
        Zt+=[sz,]
##create salinity time series in tide, this is just an arbitrary pattern for testing
    for i in tide_t:
        sal = 36*math.cos(w[1]*i+alpha[1])
        sf+=[abs(sal),]
#cjw create a tidal forcing concentration of compounds in tide for water-sediment interface boundary
##tconc_frac is the fraction of each primary species relative to salinity in water.
#    rstcl=0.14   #unit is mg/ml:mg/ml
#    salfc=0.00180665
#ppt to sulfate concetration
#    sulfc=salfc/0.14 ## sulfate concetration is mg/L, and 1e-3*mg/L/molar_mass=M
    tconc_fc={#'SOM':0,
    #'HRimm':0,
    #'DOM1':0,
    #'Acetate-':0,
#    'HCO3-':0,
    'O2(aq)':0,
#    'NO3-':0,
#    'NH4+':0,
    'SO4--':1.29046e-2,
    'Cl-':1.80665e-3,
    'Na+':1.80665e-3*0.5769,     ##include potassium 1.13%,sodium 30.6%
    'Ca++':1.80665e-3*0.024,
#    'HS-':0,
#    'Fe+++':0,
#    'Fe++':0,
#    'CH4(aq)':0,
#    'H2(aq)':0,
#    'N2(aq)':0,
#    'N2O(aq)':0,
#    'H+':0,         #this could be related to the alklinity in porewater and tide
#    'Tracer':0,
    }

    molar_mass={#'SOM':180.156,
    #'HRimm':59.052,
    #'DOM1':180.156,
    #'Acetate-':59.052,
#    'HCO3-':61.0168,
    'O2(aq)':31.998,
#    'NO3-':62.0049,
#    'NH4+':18.039,
    'SO4--':96.06,
    'Cl-':35.454,
    'Na+':22.98977,
    'Ca++':40.0780,
#    'HS-':33.1,
#    'Fe+++':55.845,
#    'Fe++':55.845,
#    'CH4(aq)':16.04,
#    'H2(aq)':2,
#    'N2(aq)':28.0134,
#    'N2O(aq)':44.013,
#    'H+':1,         #this could be related to the alklinity in porewater and tide
#    'Tracer':1,
    }   ##molar mass is g/mol
    tide_conc={'salt':sf}
#    print('primary name %1.2f'%(tide_conc['salt'][0]))
#    for spec in layers[0].primary_names:
    for spec in tconc_fc.keys():
#        print('primary name %s'%(spec))
        if spec=='O2(aq)':
            tmp=[(i*tconc_fc[spec]+9*1e-3)/molar_mass[spec] for i in sf]
        else:
            tmp=[i*tconc_fc[spec]/molar_mass[spec] for i in sf]
#            print('primary name %s,%1.20f,%1.8f'%(spec,tmp[10],tconc_fc[spec]))
        tide_conc[spec]=tmp
#        print('primary name %s,%1.20f'%(spec,tide_conc[spec][1]))
#    print(tide_conc.keys())
### cjw Zt is tide water level
##update bc state by different layer property
#cjw    idxbc=0

#cjw        root_efolding=0.15 # e-folding depth of root biomass in m
#cjw        root_biomass_top=0.3 # gC/cm3
#cjw        root_Cfrac=0.4

    z=numpy.array([0]+[l.volume for l in layers]).cumsum()
    z_mid=(z[:-1]+z[1:])/2
    for l in range(len(layers)):
#cjw            layers[l].total_immobile['Root_biomass']=root_biomass_top*numpy.exp(-z_mid[l]/root_efolding)/(root_Cfrac*12)*100**3
#cjw update following diffusion coefficient

        layers[l].diffquo={'O2(aq)':0.001**((l+1)*0.75),
        'SO4--':0.00001**((l+1)*0.75)}
#            'H2(aq)':0.001**((l+1)*0.75),
#            'N2(aq)':0.001**((l+1)*0.75),
#            'N2O(aq)':0.001**((l+1)*0.75)
#            }
        # Treat top layer as O horizon with less root biomass and mineral Mn
#cjw        layers[0].total_immobile['Root_biomass']=root_biomass_top/(root_Cfrac*12)*100**3*1e-3
        # layers[0].mineral_volume_fraction['Mn(OH)2(am)']=1e-7
#cjw        layers[0].mineral_volume_fraction['Birnessite2']=1e-7
    layers[0].surface_site_density['>DOM1']=1e2



    for l in layers:
        l.write_output(0,dt)

    flow_in=numpy.zeros(len(layers),dtype=float)
    flow_out=numpy.zeros(len(layers),dtype=float)

    t0=time.time()
    tprev=t0

#cjw        leaf_Mn_concs=numpy.ma.masked_all(nyears)

    tstart=0

        # Restart from existing state?
#cjw        restart_state='Mn_saved.nc'
    restart_state=None
    if isinstance(restart_state,str):
        restart_state=xarray.open_dataset(restart_state)
    if isinstance(restart_state,xarray.Dataset):
        copy_to_layers(restart_state, layers)
        tstart=restart_state.dropna(dim='time')['time'].isel(time=-1).item()
    elif isinstance(restart_state,list) and isinstance(restart_state[0],layer):
        layers=restart_state



    success=True
    for step in range(nsteps):
        soillyr=[0.01,0.05,0.05,0.05,0.05,0.1]
        Ds=1e-9               #m2/s
        if Zt[step] > Zs:      #when marsh is inundated, using fick's law to have a water-sediment interface, change the code below
            #Ds=0.000761    #diffusion coef for everything in water, from Morris, 1995 (m2/s)
            Ds=2.014e-9       
            lyrdpth=[(Zt[step]-Zs)/2,0.01/2,0.06/2,0.11/2,0.16/2,0.21/2,0.31/2]    # depth of the center for each soil layer
#            soillyr=[0.01,0.05,0.05,0.05,0.05,0.1]
#            print('water level on timestep %d,level %1.5f,water conc %1.20f'%(step,Zt[step],tide_conc['salt'][step]))
#            for spec in layers[0].primary_names:
            for tspec in tconc_fc.keys():
                    #calculate first order fick's law flux (mol/m2/s) at water-sediment interface can be positive/negative
#                    print('water level on timestep %d, primary name %s,water conc %1.20f'%(step,spec,tide_conc[spec][step]))
                    #if tspec=='O2(aq)':
                    Jflx=layers[0].porosity**2*Ds*(tide_conc[tspec][step]-layers[0].total_mobile[tspec])*1e3/(soillyr[0])   #(lyrdpth[0]+lyrdpth[1]) mol/m2/s
                    spec_mol_layer=layers[0].total_mobile[tspec]*layers[0].volume*1000*layers[0].porosity*layers[0].saturation
                    spec_mol_layer=spec_mol_layer+Jflx*dt*1     ##suppose the flux happens at one unit m2 surface
                    layers[0].total_mobile[tspec]=spec_mol_layer/(layers[0].volume*1000*layers[0].porosity*layers[0].saturation)


#                if spec == 'SO4--':
#                    flow_in[:]=0.0
#                    flow_out[:]=0.0
                    #print('water level on timestep %d, water %1.2f'%(step,Zt[step]))
            for spec in layers[0].primary_names:
                if spec=='O2(aq)':
                    Ds=1e-10
                for num in range(len(layers)-1):
                    #calculate flux between soil layers
                    Fflx=layers[num+1].porosity*layers[num+1].porosity*Ds*(layers[num].total_mobile[spec]-layers[num+1].total_mobile[spec])*1e3/soillyr[num+1] #(soillyr[num+1]/2+soillyr[num]/2)
                    spec_mol_layer=layers[num+1].total_mobile[spec]*layers[num+1].volume*1000*layers[num+1].porosity*layers[num+1].saturation
                    spec_mol_layer=spec_mol_layer+Fflx*dt*1
                    layers[num+1].total_mobile[spec]=spec_mol_layer/(layers[num+1].volume*1000*layers[num+1].porosity*layers[num+1].saturation)
                        # First calculate net flow for each layer
                        #flow_out[num]=flow_rate[num]*layers[num].total_mobile[spec]*10 # mol/L*cm/s*(10L/m2)/cm -> mol/m2/s
                        #flow_in[num+1]=flow_rate[num]*layers[num].total_mobile[spec]*10
                    # Leaching out of bottom layer (to groundwater/loss)
                    #flow_out[-1]=flow_rate[-1]*layers[-1].total_mobile[spec]*10
                    #for num in range(len(layers)):
                        # Then calculate change in concentration using amount of water in the layer
                        # Units of mobile species are M (mol/L water)
                        #spec_mol_layer=layers[num].total_mobile[spec]*layers[num].volume*1000*layers[num].porosity*layers[num].saturation # mol of the species in the layer initially
                        #spec_mol_layer=spec_mol_layer+(flow_in[num]-flow_out[num])*dt
                        #layers[num].total_mobile[spec]=spec_mol_layer/(layers[num].volume*1000*layers[num].porosity*layers[num].saturation)
                        #layers[num].flow_in[spec]=flow_in[num]
                        #layers[num].flow_out[spec]=flow_out[num]
        elif Zt[step] <= Zs:
            for spec in layers[0].primary_names:
                for num in range(len(layers)-1):
                    #calculate flux between soil layers
                    Fflx=layers[num+1].porosity*layers[0].porosity*Ds*(layers[num].total_mobile[spec]-layers[num+1].total_mobile[spec])*1e3/soillyr[num+1]  ##(soillyr[num+1]/2+soillyr[num]/2)
                    spec_mol_layer=layers[num+1].total_mobile[spec]*layers[num+1].volume*1000*layers[num+1].porosity*layers[num+1].saturation
                    spec_mol_layer=spec_mol_layer+Fflx*dt*1
                    layers[num+1].total_mobile[spec]=spec_mol_layer/(layers[num+1].volume*1000*layers[num+1].porosity*layers[num+1].saturation)

        dCO2=layers[0].total_mobile['HCO3-']-initial_HCO3
        layers[0].total_mobile['HCO3-']=layers[0].total_mobile['HCO3-']-dCO2
        layers[0].total_immobile['H+']=layers[0].total_immobile['H+']-dCO2*1000*layers[0].porosity*layers[0].saturation

##cjw with changing boundary conditions every time step, first try is to move the bc setup chunk below where is right before chemistry.
#cjw        if Zt[step] > Zs:
#cjw            layers[0].diffquo={'O2(aq)':0.000001,
#cjw            'SO4--':0.2013888}

#                print('water level on timestep %d, water %1.2f'%(step,Zt[step]))
##cjw below bc could equal to either sed_air_interface or sed_water_interfece which depends on water level.


##cjw above is the first try on migrating boundary condition setup into time step loop.

            # Next do chemistry.
        try:
#cjw                for n,l in enumerate(layers+[incubation_layer]):
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
#cjw                                print('diffusion coefficient on timestep %d, dq %1.2f'%(step,dq[spec]))
                        else:
                            dq[spec]=0.0

                                # print('O2 before: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
                num_cuts=l.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=l.rateconstants)
                l.copy_from_alquimia(data)
                    # Write output
                l.write_output(step+1,dt,num_cuts)
                    # print('O2 after: %1.1g'%data.state.total_mobile.data[get_alquimiavector(data.meta_data.primary_names).index('O2(aq)')])
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
        axs[num,0].plot(t,out['Total DOM1']*1e3,c='C2',label='DOM',**kwargs)
        axs[num,1].plot(t,out['Total O2(aq)']*1e3,label='O$_2$',c='C0',**kwargs)
        axs[num,2].plot(t,out['Total Cl-']*1e3,label='Cl-',c='C1',**kwargs)
        axs[num,3].plot(t,out['Total Na+']*1e3,label='Na+',c='C5',**kwargs)
        # axs[num,0].plot(t,out['Total Sorbed DOM1']*12/100**3,label='Sorbed DOM')
        axs[num,4].plot(t,out['Total SO4--']*1e3,c='C3',label='SO$_4$',**kwargs)
        axs[num,5].plot(t,out['Total CH4(aq)']*1e3,c='C4',label='CH$_4$',**kwargs)
        axs[num,0].set_ylabel('mM')
        #axs[num,1].set_ylabel('mM')
        #axs[num,2].set_ylabel('mM')
        #axs[num,3].set_ylabel('mM')
        #axs[num,4].set_ylabel('mM')
        #if num < 2:
        #    axs[num,0].set_ylim([350,700])
        #elif num >= 2 :
        #    axs[num,0].set_ylim([350,1e4])
        axs[num,0].set_ylim([0,6])
        axs[num,1].set_ylim([0,0.02])
        #axs[num,2].set_ylim([0,0.003])
        axs[num,4].set_ylim([0,0.5])
        axs[num,5].set_ylim([0,0.6])
        #if num < 4:
        #    axs[num,4].set_ylim([0,100])
        #elif num >=4 :
        #    axs[num,4].set_ylim([0,1200])

#        axs[num,0].set_ylabel('C density\n(g C cm$^{-3}$)')

#        axs[num,1].plot(t,out['Total Mn++']*1e6/1000*porosity*saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+}$',c='C0',**kwargs)
#        axs[num,1].plot(t,out['Total Mn+++']*1e6/1000*porosity*saturation*Mn_molarmass/BD,label='Mn$^{+\!\!+\!\!+}$',c='C1',**kwargs)
        # axs[num,1].plot(t,layers[num].output_DF['Manganite VF']/molar_volume_manganite*1e6*Mn_molarmass/BD,label='Manganite')
#        axs[num,2].plot(t,out['Birnessite2 VF']*7/molar_volume_birnessite*1e6*Mn_molarmass/BD,label='Birnessite',c='C0',**kwargs)

#        annual_root_uptake=out['Total Tracer2'].groupby(numpy.floor(t)).max()
#        axs[num,1].plot(annual_root_uptake.time,annual_root_uptake*1e6/1000*porosity.mean()*saturation*Mn_molarmass/BD,label='Annual Mn$^{+\!\!+}$ root uptake',c='C2',**kwargs)

#        axs[num,3].plot(t,-numpy.log10(out['Free H+']),label='pH',c='C0',**kwargs)

        # ax=axs[num,1].twinx()
        # axs[num,2].plot(t,layers[num].output_DF['Mn(OH)2(am) VF']/molar_volume_MnOH2am*1e6*Mn_molarmass/BD,label='Mn(OH)$_2$(am)',ls='-')
#        axs[num,2].plot(t,(out['Total Mn++']+out['Total Mn+++'])*1e6/1000*porosity*saturation*Mn_molarmass/BD +
#                            (out['Birnessite2 VF']*7/molar_volume_birnessite)*1e6*Mn_molarmass/BD,c='k',label='Total Mn',**kwargs)

#        axs[num,1].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
#        axs[num,2].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
#        axs[num,3].set_ylabel('pH')
        # axs[num,1].set_ylim(*axs[0,1].get_ylim())

#        for n,cation in enumerate(['Al+++','Ca++','K+','Mg++','Na+','Mn++','Mn+++']):
#            axs[num,4].plot(t,output['Total Sorbed '+cation].isel(depth=num,time=slice(None,None,subsample)).dropna(dim='time')/(output.BD.isel(depth=2)*1e-3*100**3)*1000,c='C'+str(n),label=cation,**kwargs)
#        axs[num,4].plot(t,output['CEC H+'].isel(depth=num,time=slice(None,None,subsample)).dropna(dim='time')/(output.BD.isel(depth=2)*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),**kwargs)

#        axs[num,4].set_ylabel('Exch conc\n(mmol/kg)')



    axs[-1,0].set_xlabel('Time (days)')
    if do_legend:
        axs[0,0].legend()
        axs[0,1].legend()
        axs[0,2].legend()
        axs[0,3].legend()
        axs[0,4].legend()
        axs[0,5].legend()

    axs[-1,1].set_xlabel('Time (days)')
    axs[-1,2].set_xlabel('Time (days)')
    axs[-1,3].set_xlabel('Time (days)')
    axs[-1,4].set_xlabel('Time (days)')
    axs[-1,5].set_xlabel('Time (days)')
    axs[0,0].set_title('Total DOM1')
    axs[0,1].set_title('Total Oxygen')
    axs[0,2].set_title('Total Cl-')
    axs[0,3].set_title('Total Na+')
    axs[0,4].set_title('Total SO4')
    axs[0,5].set_title('CH4')

f,axs=pyplot.subplots(ncols=6,nrows=len(output.depth),sharex=True,clear=True,num='Simulation results',figsize=(12,6.5))
plot_output(output,axs)

pyplot.show()
#cjw networkfig=pyplot.figure('Reaction network',clear=True)
#cjw drawn=decomp_network.draw_network_with_reactions(reaction_network,omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid'],
#c        font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
#c        namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'HCO$_3^-$','DOM2':'Exposed lignin','sorbed_DOM1':'Sorbed DOM',
#c                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',})
