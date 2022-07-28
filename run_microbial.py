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
        
        
    def run_onestep(self,chem,data,dt,status,min_dt=0.1,num_cuts=0,diffquo={},bc=None,flux_tol=0.15,truncate_concentration=0.0,rateconstants={},min_cuts=0):
        self.copy_to_alquimia(data)
        converged=False
        porosity=data.state.porosity
        
        max_cuts=num_cuts
        actual_dt=dt/2**num_cuts

        #rate_nm=list(rateconstants.keys())
        
        for num,reactname in enumerate(get_alquimiavector(data.meta_data.aqueous_kinetic_names)):
            data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[reactname]
            
            #data.properties.aqueous_kinetic_rate_cnst.data[num]=rateconstants[rate_nm[num]]
            #print(data.properties.aqueous_kinetic_rate_cnst.data[num],reactname)
            #print(data.properties.aqueous_kinetic_rate_cnst.data[num],reactname,rateconstants[reactname])
        
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
                print(self)
                print('diffquo:')
                for spec in diffquo:
                    print(spec.ljust(14)+'%1.2e'%diffquo[spec])
                if cut_for_flux:
                    raise RuntimeError('Pflotran failed to converge (because of boundary fluxes) after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
                else:
                    raise RuntimeError('Pflotran failed to converge after %d cuts to dt = %1.2f s'%(num_cuts,actual_dt))
            
            # data will be reset to layer contents at beginning of next call
            # Run it twice, because we cut the time step in half
            ncuts=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=num_cuts+1,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,flux_tol=flux_tol)
            # If this completes, it means that run_onestep was successful
            self.copy_from_alquimia(data)
            if ncuts>max_cuts:
                max_cuts=ncuts
                
            # This starts from ncuts so it doesn't have to try all the ones that failed again
            for n in range(2**(ncuts-(num_cuts+1))):
                ncuts2=self.run_onestep(chem,data,dt,status,min_dt,num_cuts=ncuts,diffquo=diffquo,bc=bc,truncate_concentration=truncate_concentration,rateconstants=rateconstants,flux_tol=flux_tol)
                self.copy_from_alquimia(data)
                if ncuts2>max_cuts:
                    max_cuts=ncuts2

            return max_cuts


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
for pH in numpy.arange(6.3,7.3):

    chem,data,sizes,status=init_alquimia(input_file,hands_off=True)
        #cjw below is the codes for layer setting up.
        # Set up layers
        # Top (organic) layer should be thinner and have lower bulk density though
        # Low bulk density causes simulation to slow or crash though. Actually CEC being too low (<100 combined with BD<1) is the problem
#cjw        layers=[layer(0.05,rateconstants=rateconstants_warmed,BD=0.05,porosity=0.5)]+[layer(0.1,BD=0.25,rateconstants=rateconstants_warmed) for num in range(3)]
    layers=[layer(0.1,rateconstants=rateconstants_warmed,BD=0.08,porosity=0.95,saturation=1)]
    for l in layers:
        l.secondary_names=secondary_names
    for l in layers:
        l.initcond=Mar.pools.copy()

    for l in layers:
        l.initcond=decomp_network.change_constraints(Mar.pools,{'O2(aq)':0.002})
    bc=layers[0].initcond
    #print(bc,layers[0].initcond)

#cjw    bc=None
    dt=1
    nyears=1
    #nsteps=365*24//(dt//3600)*nyears
    nsteps=3600
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
            * Starting simulation with pH = %1.1f, Ndep = %03d, Warming = %d *
            *************************************


        '''%(pH,0,20))        #cjw remove warming



    #initial_HCO3 = l.total_mobile['HCO3-']
    #initial_O2 = l.total_mobile['O2(aq)']
        # Flow rate cm/s = 10 L/m2/s, positive is downward
    #flow_rate=numpy.linspace(1e-6,1e-7,len(layers)) # Rate declines linearly with depth, assumes removal or accumulation in lower layers
    flow_rate=numpy.zeros(len(layers))*0
#    flow_rate=numpy.linspace(1e-7,1e-8,len(layers))
    min_dt=0.1
    truncate_concentration=1e-20
    z=numpy.array([0]+[l.volume for l in layers]).cumsum()
    z_mid=(z[:-1]+z[1:])/2
    for l in range(len(layers)):
#cjw            layers[l].total_immobile['Root_biomass']=root_biomass_top*numpy.exp(-z_mid[l]/root_efolding)/(root_Cfrac*12)*100**3
#cjw update following diffusion coefficient
        layers[l].diffquo={'O2(aq)':0.01**((l+1)*0.75)}
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

#cjw 
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
                    B_new=Bgrow[micdx,step]*dt+B_old
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
# Next advect mobile species. This assumes we can separate advective transport and reactions at this time step length
            # Alternate strategy would be to pass these rates to run_onestep and spread them over the variable timesteps
        for spec in layers[0].primary_names:
            #if spec in immobile_specs:
            #    continue # Skip this one because it's not really mobile, more a workaround for Lignin decomp
            flow_in[:]=0.0
            flow_out[:]=0.0
            for num in range(len(layers)-1):
                    # First calculate net flow for each layer
                flow_out[num]=flow_rate[num]*layers[num].total_mobile[spec]*10 # mol/L*cm/s*(10L/m2)/cm -> mol/m2/s
                flow_in[num+1]=flow_rate[num]*layers[num].total_mobile[spec]*10
                # Leaching out of bottom layer (to groundwater/loss)
            flow_out[-1]=flow_rate[-1]*layers[-1].total_mobile[spec]*10
            for num in range(len(layers)):
                    # Then calculate change in concentration using amount of water in the layer
                    # Units of mobile species are M (mol/L water)
                spec_mol_layer=layers[num].total_mobile[spec]*layers[num].volume*1000*layers[num].porosity*layers[num].saturation # mol of the species in the layer initially
                spec_mol_layer=spec_mol_layer+(flow_in[num]-flow_out[num])*dt
                layers[num].total_mobile[spec]=spec_mol_layer/(layers[num].volume*1000*layers[num].porosity*layers[num].saturation)
                layers[num].flow_in[spec]=flow_in[num]
                layers[num].flow_out[spec]=flow_out[num]        
        #dCO2=layers[0].total_mobile['HCO3-']-initial_HCO3
        #layers[0].total_mobile['HCO3-']=layers[0].total_mobile['HCO3-']-dCO2
        #layers[0].total_immobile['H+']=layers[0].total_immobile['H+']-dCO2*1000*layers[0].porosity*layers[0].saturation
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
                        else:
                            dq[spec]=0.0
##cjw: update rateconstant each timestep to represent microbial activities.
                #rateconstants_microbe=rateconst.copy()
                rateconstants_microbe=rateconstants_bio.copy()   
                #print('rateconst before alquimia',rateconstants_bio)
                #rateconstants_microbe=convert_rateconstants(rateconstants_bio,reaction_network,precision=precision)
                #print('rateconst in alquimia',rateconstants_microbe)
                l.rateconstants=rateconstants_microbe
                #primarynames=l.primary_names
                #print(primarynames)
##cjw: end of microbial activities.
                num_cuts=l.run_onestep(chem,data,dt,status,min_dt=min_dt,diffquo=dq,bc=bc_state,truncate_concentration=truncate_concentration,rateconstants=l.rateconstants)
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
ax1.plot(out['Total Sorbed microbe1'],label='microbe$_1$',linewidth=2)
ax2.plot(out['Total Sorbed microbe2'],label='microbe$_2$',linewidth=2)
ax3.plot(out['Total Sorbed microbe3'],label='microbe$_3$',linewidth=2)
ax4.plot(out['Total Sorbed SOMC'],label='SOMC',linewidth=2)
ax5.plot(out['Total Sorbed SOMN'],label='SOMN',linewidth=2)
pyplot.savefig('microbes.pdf')

geochem,(ax1,ax2,ax3,ax4)=pyplot.subplots(4, 1,figsize=(8,6))
ax1.plot(out['Total Fe+++'],label='Fe+++',linewidth=2)
ax2.plot(out['Total Acetate-'],label='Acetate-',linewidth=2)
ax3.plot(out['Total DOM1'],label='DOM1',linewidth=2)
ax4.plot(out['Total NH4+'],label='NH4+',linewidth=2)
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
