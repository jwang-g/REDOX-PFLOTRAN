from pylab import *
import Arctic_Fe_CH4 as arctic
import decomp_network
import xarray
import sys
import pandas


def plot_result(result,SOM_ax=None,pH_ax=None,Fe_ax=None,CO2flux_ax=None,CH4flux_ax=None,porewater_ax=None,cation_ax=None,
                SOM_args={},pH_args={},Fe_args={},CO2flux_args={},CH4flux_args={},porewater_args={},cation_args={},
                do_legend=False,gdrywt=False,BD=None,SOC_pct=None,cellulose_SOC_frac=1.0,obs=None,obs_marker=None):
    t=result['time']
    if obs is not None:
        obs_mod=pandas.DataFrame(index=obs.Incubation_Time)
        obs_mod_pts=array([abs(t-obst).argmin() for obst in obs.Incubation_Time])
        obs_mod_pts[obs_mod_pts==0]=1
    else:
        obs_mod=None
    if gdrywt:
        if SOC_pct is None and BD is None:
            raise TypeError('SOC_pct or BD must be a number if gdrywt is True')
        SOC_mol_m3=result['Total Sorbed cellulose'].isel(time=0)/cellulose_SOC_frac # mol SOC/m3
        SOC_gC_cm3=SOC_mol_m3*12/100**3
        SOC_gC_gdwt=SOC_pct/100
        cm3_to_dwt=SOC_gC_gdwt/SOC_gC_cm3 # Conversion from /cm3 to /gdwt. This is 1/ bulk density in g/cm3
        if BD is not None:
            cm3_to_dwt=1.0/BD
        # print('cm3_to_dwt = %1.3f'%cm3_to_dwt)
    if SOM_ax is not None:
        if gdrywt:
            # Plot in %, i.e. gC/gdrywt *100.
            l=SOM_ax.plot(t,result['Total Sorbed cellulose']/SOC_mol_m3*SOC_pct+SOC_pct*(1-cellulose_SOC_frac),label='SOM',**SOM_args)[0]
            SOM_ax.set_ylabel('Concentration\n(SOC %)')
        else:
            l=SOM_ax.plot(t,result['Total Sorbed cellulose']*1e-3,label='SOM',**SOM_args)[0]
            SOM_ax.set_ylabel('Concentration\n(mmol C/cm$^{-3}$)')

        SOM_ax.set_title('SOM remaining')
        SOM_ax.set_xlabel('Time (days)')

    if obs is not None:
        obs_mod['pH_obs']=obs['pH'].set_axis(obs.Incubation_Time)
        obs_mod['pH_mod']=-log10(result['Free H+'])[obs_mod_pts]
    if pH_ax is not None:
        pH_ax.plot(t,-log10(result['Free H+']))

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
    
    if obs is not None:
        obs_mod['Fe_II_obs']=obs['Fe_II'].set_axis(obs.Incubation_Time)
        obs_mod['Fe_II_mod']=result['Total Fe++'][obs_mod_pts]*result['Porosity']*1e3*cm3_to_dwt
    if Fe_ax is not None:
        molar_volume=34.3600 # From database. cm3/mol
        molar_weight = 106.8690
        # Add sorbed Fe to this
        if gdrywt:
            # Assume we can use SOC % to convert from volume to dry weight
            l=Fe_ax.plot(t,result['Fe(OH)3 VF']/molar_volume*1e6*cm3_to_dwt   ,label='Iron oxide',**Fe_args)[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol g dwt$^{-1}$)')
            l=Fe_ax.plot(t,result['Total Fe+++']*result['Porosity']*1e3*cm3_to_dwt   ,label='Fe$^{3\!\!+}$',ls='--',**Fe_args)[0]
            
            l=Fe_ax.plot(t,result['Total Fe++']*result['Porosity']*1e3*cm3_to_dwt ,ls=':'  ,label='Fe$^{2\!\!+}$',**Fe_args)[0]
            Fe_ax.plot(t,(result['Total Sorbed Fe++']+result['Total Sorbed Fe+++'])*1e6/100**3*cm3_to_dwt,**Fe_args,label='Sorbed Fe')
            if obs is not None:
                if obs_marker is not None:
                    markers=['o','s','^','+','<','x']
                    obs_cats=unique(obs[obs_marker])
                    print(obs_cats)
                    for x in range(len(obs_cats)):
                        xx=obs[obs_marker]==obs_cats[x]
                        Fe_ax.plot(t,obs.Incubation_Time[xx],obs['Fe_II'][xx],marker=markers[x],linestyle='None',label='Measured Fe$^{2\!\!+}$',color='C2')
                else:
                    Fe_ax.plot(obs.Incubation_Time,obs['Fe_II'],marker='o',linestyle='None',label='Measured Fe$^{2\!\!+}$',color='C2')
        else:
            l=Fe_ax.plot(t,result['Fe(OH)3 VF']/molar_volume*1e6   ,label='Iron oxide (solid)')[0]
            Fe_ax.set_ylabel('Concentration\n($\mu$mol cm$^{-3}$)')
        
            # M/L to umol/cm3: 1e6/1e3=1e3
            l=Fe_ax.plot(t,result['Total Fe+++']*result['Porosity']*1e3   ,label='Fe+++',ls='--')[0]
            
            l=Fe_ax.plot(t,result['Total Fe++']*result['Porosity']*1e3 ,ls=':'  ,label='Fe++')[0]
        
        Fe_ax.set_title('Fe species')
        
        Fe_ax.set_xlabel('Time (days)')
        if do_legend:
            Fe_ax.legend(fontsize='small')
    
    if obs is not None:
        CO2mean=obs[['CO2_1','CO2_2','CO2_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff()
        obs_mod['CO2flux_obs']=CO2mean.set_axis(obs.Incubation_Time)
        obs_mod['CO2flux_mod']=(diff(result['Total Tracer']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt)[obs_mod_pts]
    if CO2flux_ax is not None:
        # gasflux_ax.set_yscale('log')
        if gdrywt:
            
            l=CO2flux_ax.plot(t[:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt,label='CO2',**CO2flux_args)[0]
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
            
            l=CO2flux_ax.plot(result['time'][:-1],diff(result['Total Tracer']*result['Porosity'])/diff(result['time'])*1e3,label='CO2',**CO2flux_args)[0]
            CO2flux_ax.set_ylabel('CO$_2$ flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        CO2flux_ax.set_title('CO$_2$ flux')
        
        CO2flux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     CO2flux_ax.legend(fontsize='small')

    if obs is not None:
        CH4mean=obs[['CH4_1','CH4_2','CH4_3']].astype(float).diff().mean(axis=1)/obs['Incubation_Time'].diff()
        obs_mod['CH4flux_obs']=CH4mean.set_axis(obs.Incubation_Time)
        obs_mod['CH4flux_mod']=(diff(result['Total CH4(aq)']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt)[obs_mod_pts]
    if CH4flux_ax is not None:
        if gdrywt:
            l=CH4flux_ax.plot(t[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result['time'])*1e3*cm3_to_dwt,label='CH4',**CH4flux_args)[0]
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
            l=CH4flux_ax.plot(t[:-1],diff(result['Total CH4(aq)']*result['Porosity'])/diff(result['time'])*1e3,label='CH4',**CH4flux_args)[0]
            
            CH4flux_ax.set_ylabel('CH$_4$ flux rate\n($\mu$mol cm$^{-3}$ day$^{-1}$)')

        CH4flux_ax.set_title('CH$_4$ flux')
        
        CH4flux_ax.set_xlabel('Time (days)')
        # if do_legend:
        #     CH4flux_ax.legend(fontsize='small')
        
    if obs is not None:
        obs_mod['DOM_obs']=obs['WEOC'].set_axis(obs.Incubation_Time)*1e-6/BD*1000/result['Porosity'][0].item()
        obs_mod['DOM_mod']=result['Total DOM1'][obs_mod_pts]
        obs_mod['Acetate_obs']=obs['TOAC'].set_axis(obs.Incubation_Time)*1e-6/BD*1000/result['Porosity'][0].item()
        obs_mod['Acetate_mod']=result['Total Acetate-'][obs_mod_pts]
    if porewater_ax is not None:
        porewater_ax.set_yscale('log')
        porewater_ax.plot(t,result['Total DOM1'],label='DOM')
        porewater_ax.plot(t,result['Total Acetate-'],label='Acetate',c='C3')
        # porewater_ax.plot(t,result['Total O2(aq)'],'--',label='O2',c='C4')
        # porewater_ax.plot(t,result['Total Fe+++'],'--',label='Fe+++',c='C1')
        # porewater_ax.plot(t,result['Free Fe+++'],':',label='Fe+++',c='C1')
        # porewater_ax.plot(t,result['Total Fe++'],':',label='Fe++',c='C2')
        
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
                    porewater_ax.plot(obs.Incubation_Time[xx],obs['TOAC'][xx]*1e-6/BD*1000/result['Porosity'][0].item(),marker=markers[x],linestyle='None',label='Measured TOAC',color='C3',**porewater_args)
                    porewater_ax.plot(obs.Incubation_Time[xx],obs['WEOC'][xx]*1e-6/BD*1000/result['Porosity'][0].item(),marker=markers[x],linestyle='None',label='Measured WEOC',color='C0',**porewater_args)
            else:
                porewater_ax.plot(obs.Incubation_Time,obs['TOAC']*1e-6/BD*1000/result['Porosity'][0].item(),marker='o',linestyle='None',label='Measured TOAC',color='C3')
                porewater_ax.plot(obs.Incubation_Time,obs['WEOC']*1e-6/BD*1000/result['Porosity'][0].item(),marker='o',linestyle='None',label='Measured WEOC',color='C0')


        
        if do_legend:
            porewater_ax.legend(fontsize='small')
            
    if cation_ax is not None:
        for n,cation in enumerate(decomp_network.pools_list_to_dict(arctic.reactions)['Cation exchange']['cations'].keys() ):
            if cation == 'H+':
                continue
            cation_ax.plot(t,result['Total Sorbed '+cation]/(BD*1e-3*100**3)*1000,c='C'+str(n),label=cation,**cation_args)
        cation_ax.plot(t,result['CEC H+']/(BD*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),**cation_args)
        
        cation_ax.set_ylabel('Exch conc\n(mmol/kg)')
        
        if do_legend:
            cation_ax.legend(fontsize='small')

    return obs_mod
        

def plot_timeseries(fname,group_suffix=''):
    if isinstance(fname,str):
        fnames=[fname]
    else:
        fnames=fname

    result_organic_trough=None
    result_organic_nottrough=None
    result_mineral_trough=None
    result_mineral_nottrough=None
    result_highO2_organic=None

    for f in fnames:
        if 'organic_trough'+group_suffix in get_groups(f):
            result_organic_trough=xarray.open_dataset(f,group='organic_trough'+group_suffix)
        if 'organic_nottrough'+group_suffix in get_groups(f):
            result_organic_nottrough=xarray.open_dataset(f,group='organic_nottrough'+group_suffix)
        if 'mineral_trough'+group_suffix in get_groups(f):
            result_mineral_trough=xarray.open_dataset(f,group='mineral_trough'+group_suffix)
        if 'mineral_nottrough'+group_suffix in get_groups(f):
            result_mineral_nottrough=xarray.open_dataset(f,group='mineral_nottrough'+group_suffix)
        if 'highO2_organic'+group_suffix in get_groups(f):
            result_highO2_organic=xarray.open_dataset(f,group='highO2_organic'+group_suffix)

    if result_organic_trough is None:
        raise ValueError('result_organic_trough'+group_suffix+' not found')
    if result_organic_nottrough is None:
        raise ValueError('result_organic_nottrough'+group_suffix+' not found')
    if result_mineral_trough is None:
        raise ValueError('result_mineral_trough'+group_suffix+' not found')
    if result_mineral_nottrough is None:
        raise ValueError('result_mineral_nottrough'+group_suffix+' not found')
    if result_highO2_organic is None:
        raise ValueError('result_highO2_organic'+group_suffix+' not found')

    import os
    # fig,axes=subplots(3,1,num='Organic horizon Anoxic',figsize=(6,8.4),clear=True)
    f,cation_axes=subplots(ncols=5,num='Cations: %s %s'%(os.path.basename(fname[0]),group_suffix),clear=True,figsize=(13,4))
    fig,axs=subplots(nrows=5,ncols=5,num='Time series: %s %s'%(os.path.basename(fname[0]),group_suffix),clear=True,figsize=(16,8),sharex=False)
    from string import ascii_lowercase
    for x in range(5):
        for y in range(5):
            axs[x,y].set_title('('+ascii_lowercase[x*5+y]+')',loc='left')

    axes=axs[:,1]
    obs=arctic.get_layer(9,'Organic')
    plot_result(result_organic_trough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[1],do_legend=False,BD=arctic.BD_organic_trough,
                gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,#BD=BD_layerest2_trough['Organic',True],
                # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')
                        # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')])
                        obs=obs)
    # axes[0].set_ylim(bottom=-0.01)

    # axes[0].set_ylim(bottom=0.9e-11,top=15)
    axes[0].set_yscale('linear')
    axes[0].text(0.5,1.3,'%s Organic horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')


    # xmax=axes[0].get_xlim()[1]
    xmax=max(arctic.simlength,90)
    for ax in axes:
        ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
        ax.set_xlim(left=-5,right=xmax)
        
    axes=axs[:,2]
    obs=arctic.get_layer(5,'Organic')
    plot_result(result_organic_nottrough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[2],do_legend=False,
                gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_organic_nottrough,
                # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')
                        # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()!='trough')])
                        obs=obs)
    # axes[0].set_ylim(bottom=-0.01)

    # axes[0].set_ylim(bottom=0.9e-11,top=15)
    axes[0].set_yscale('linear')
    axes[0].text(0.5,1.3,'%s Organic horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')


    # xmax=axes[0].get_xlim()[1]
    # xmax=90
    for ax in axes:
        ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
        ax.set_xlim(left=-5,right=xmax)


    # fig,axes=subplots(3,1,num='Mineral horizon Anoxic',figsize=(6,8.4),clear=True)
    axes=axs[:,3]
    obs=arctic.get_layer(9,'Mineral')
    plot_result(result_mineral_trough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[3],do_legend=False,
        gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_mineral_trough,
        # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')
            # &(Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()=='trough')])
            obs=obs)
    # axes[0].set_ylim(bottom=-0.01)

    # axes[0].set_ylim(bottom=0.9e-11,top=3.0)
    axes[0].set_yscale('linear')
    axes[0].text(0.5,1.3,'%s Mineral horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

    # xmax=axes[0].get_xlim()[1]
    # xmax=90
    for ax in axes:
        ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
        ax.set_xlim(left=-5,right=xmax)

    axes=axs[:,4]
    obs=arctic.get_layer(5,'Mineral')
    plot_result(result_mineral_nottrough,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[4],do_legend=False,
        gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_mineral_nottrough,
        # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Anoxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Mineral')&
            # (Barrow_synthesis_data['Incubation_Temperature']>4)&(Barrow_synthesis_data['Microtopography'].str.lower()!='trough')])
            obs=obs)
    # axes[0].set_ylim(bottom=-0.01)

    # axes[0].set_ylim(bottom=0.9e-11,top=3.0)
    axes[0].set_yscale('linear')
    axes[0].text(0.5,1.3,'%s Mineral horizon (Anoxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

    # xmax=axes[0].get_xlim()[1]
    # xmax=90
    for ax in axes:
        ax.axvspan(arctic.initfrac*arctic.simlength,xmax,color='b',alpha=0.1)
        ax.set_xlim(left=-5,right=xmax)


    # Oxic
    # fig,axes=subplots(3,1,num='Organic horizon Oxic',figsize=(6,8.4),clear=True)
    axes=axs[:,0]
    obs=arctic.get_layer(3,'Organic')
    plot_result(result_highO2_organic,Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],porewater_ax=axes[4],cation_ax=cation_axes[0],do_legend=False,
            gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_atmoO2_organic,
            # obs=Barrow_synthesis_data[(Barrow_synthesis_data.Headspace=='Oxic')&(Barrow_synthesis_data['Soil_layer'].str.capitalize()=='Organic')&(Barrow_synthesis_data['Incubation_Temperature']>4)])
            obs=obs)
    # axes[0].set_ylim(bottom=-0.01,top=30)
    axes[1].set_ylim(bottom=-0.01,top=0.2)

    axes[0].set_yscale('linear')
    # axes[0].set_ylim(-5e-2,155e-2) # This leaves out Core NGADG0003 which has 10x higher flux for some reason
    # axes[1].set_ylim(bottom=0.9e-11)
    axes[0].text(0.5,1.3,'%s Organic horizon (Oxic)'%(obs['Microtopography'].iloc[0]),fontsize='large',transform=axes[0].transAxes,ha='center')

    # axes[0].set_xlim(right=71)
    for ax in fig.axes:
        ax.set_xlim(left=-5,right=xmax)
    for ax in axs[3,:]:
        ax.set_ylim(3.8,6.1)
    for ax in axs[4,:]:
        ax.set_ylim(bottom=1e-5)

    leg=axs[2,4].legend(loc=(0.03,0.45),ncol=2,edgecolor='k')
    cation_axes[0].legend()
    axes[4].legend()

    return 


def get_groups(filename):
    import netCDF4
    with netCDF4.Dataset(filename) as d:
        groups=list(d.groups.keys())
    return groups


def mod_obs_comp(fname,group_suffix=''):
    if isinstance(fname,str):
        fnames=[fname]
    else:
        fnames=fname

    result_organic_trough=None
    result_organic_nottrough=None
    result_mineral_trough=None
    result_mineral_nottrough=None
    result_highO2_organic=None

    for f in fnames:
        if 'organic_trough'+group_suffix in get_groups(f):
            result_organic_trough=xarray.open_dataset(f,group='organic_trough'+group_suffix)
        if 'organic_nottrough'+group_suffix in get_groups(f):
            result_organic_nottrough=xarray.open_dataset(f,group='organic_nottrough'+group_suffix)
        if 'mineral_trough'+group_suffix in get_groups(f):
            result_mineral_trough=xarray.open_dataset(f,group='mineral_trough'+group_suffix)
        if 'mineral_nottrough'+group_suffix in get_groups(f):
            result_mineral_nottrough=xarray.open_dataset(f,group='mineral_nottrough'+group_suffix)
        if 'highO2_organic'+group_suffix in get_groups(f):
            result_highO2_organic=xarray.open_dataset(f,group='highO2_organic'+group_suffix)

    obs_mods={}
    if result_organic_trough is not None:
            obs_mods['organic_trough']=plot_result(result_organic_trough,BD=arctic.BD_organic_trough,
                gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,obs=arctic.get_layer(9,'Organic'))
    
    if result_organic_nottrough is not None:
            obs_mods['organic_nottrough']=plot_result(result_organic_nottrough,gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_organic_nottrough,
                        obs=arctic.get_layer(5,'Organic'))
    if result_mineral_trough is not None:
            obs_mods['mineral_trough']=plot_result(result_mineral_trough,gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,
                                    BD=arctic.BD_mineral_trough,obs=arctic.get_layer(9,'Mineral'))
    if result_mineral_nottrough is not None:
            obs_mods['mineral_nottrough']=plot_result(result_mineral_nottrough,gdrywt=True,SOC_pct=arctic.SOC_layermean['Mineral'],cellulose_SOC_frac=arctic.cellulosefrac,
                         BD=arctic.BD_mineral_nottrough,obs=arctic.get_layer(5,'Mineral'))
    if result_highO2_organic is not None:
            obs_mods['highO2_organic']=plot_result(result_highO2_organic,gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,
                        BD=arctic.BD_atmoO2_organic,obs=arctic.get_layer(3,'Organic'))

    return obs_mods
 
def getstats(obs_mod,logged=['DOM','Acetate']):
    specs=['pH','Fe_II','CO2flux','CH4flux','DOM','Acetate']
    out=pandas.DataFrame(index=specs,columns=['RMSE','nRMSE','MAE','R','R2'])
    for spec in specs:
        mod=obs_mod[spec+'_mod']
        obs=obs_mod[spec+'_obs']
        if spec in logged:
            mod=log10(mod)
            obs=log10(obs)
        out['RMSE'][spec] = sqrt(((obs-mod)**2).mean())
        out['nRMSE'][spec] = out['RMSE'][spec]/obs.mean()
        out['MAE'][spec]  = abs(obs-mod).mean()
        if spec in logged:
            out['R'][spec]    = log10(obs_mod[[spec+'_obs',spec+'_mod']]).corr().iloc[0,1]
        else:
            out['R'][spec] = obs_mod[[spec+'_obs',spec+'_mod']].corr().iloc[0,1]
        out['R2'][spec]   = out['R'][spec]**2

    return out.rename({spec:'log(%s)'%spec for spec in logged})

def save_all_figs(dirname,format='png',**kwargs):
    for fname in get_figlabels():
        fname_fixed=fname.replace('/','-')
        print(fname_fixed)
        figure(fname_fixed).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed),**kwargs)


if __name__ == '__main__':

    networkfig2=figure('Reaction network (with reactions)',clear=True,figsize=(9.4,6.6))
    drawn=decomp_network.draw_network_with_reactions(arctic.reaction_network,node_size=1200,arrowstyle='->',arrowsize=8.5,#edge_color='gray',
                omit=['NH4+','Rock(s)','gas','surf_complex','secondary','H+','implicit','Fe(II) abiotic oxidation','Fe(OH)2','Calcite'],
                namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$','CH4(aq)':'CH$_4$','HCO3-':'CO$_2$',
                            'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Fe++':r'Fe$^\mathrm{+\!\!+}$','Fe+++':r'Fe$^\mathrm{+\!\!+\!\!\!+}$',
                            'Acetate-':'Acetate','fermentation':'Fermentation','Fe(II) microbial oxidation':r'Fe$^\mathrm{+\!\!+}$'+'\noxidation',
                            'DOM aerobic respiration':'DOM aerobic\nrespiration','Hydrogen oxidation':'Hydrogen\noxidation','Aerobic decomposition':'Aerobic\ndecomposition',
                            'Acetaclastic methanogenesis':'Acetaclastic\nmethanogenesis','Acetate aerobic respiration':'Acetate\naerobic respiration',
                            'Hydrogenotrophic methanogenesis':'Hydrogenotrophic\nmethanogenesis','Gas':'Dissolved gas','Primary aqueous':'Dissolved ion',
                            'Unknown':'Cation exchange','Fe(III) reduction':r'Fe$^\mathrm{+\!\!+\!\!\!+}$'+'\nreduction','H2(aq)':'H$_2$'
                            },markers={'mineral':'8'})

    colors={'Anaerobic':'C1','Periodic':'C2','Low Fe':'C3','Aerobic':'C0'}

    if len(sys.argv)<2:
        fname='Arctic_Fe_output/results_2021-03-31.nc'
    else:
        fname=sys.argv[1:]

    data_list={}

    import netCDF4
    for filename in fname:
        with netCDF4.Dataset(filename) as d:
            groups=list(d.groups.keys())

        for g in groups:
            print(g,flush=True)
            
            d=xarray.open_dataset(filename,group=g,decode_times=False)
            if 'nperiods' in g:
                gs=g.split('_')
                d['ndryperiods']=int(gs[-3])
                d['pH']=float(gs[-1])
            
                newdims=['pH','ndryperiods']
                d=d.expand_dims(newdims).set_coords(newdims)
            
            data_list[g]=d

    results_periodic=xarray.combine_by_coords(data_list[xx] for xx in data_list if 'pH' in xx and xx.startswith('organic'))
    results_periodic_mineral=xarray.combine_by_coords(data_list[xx] for xx in data_list if 'pH' in xx and xx.startswith('mineral'))


    plot_timeseries(fname)

    # results_periodic=xarray.open_dataset(fname,group='periodic_sims')

    fig,axes=subplots(4,1,num='Periodic inundation',figsize=(6,8.4),clear=True)
    # CH4ax=axes[0].twinx()
    plot_result(results_periodic.isel(ndryperiods=1,pH=2),Fe_ax=axes[2],CO2flux_ax=axes[0],CH4flux_ax=axes[1],pH_ax=axes[3],#porewater_ax=axes[3],
            gdrywt=True,SOC_pct=arctic.SOC_layermean['Organic'],cellulose_SOC_frac=arctic.cellulosefrac,BD=arctic.BD_layerest2_trough['Organic',False])
    # axes[0].set_ylim(bottom=-0.01)
    axes[2].legend(labels=['Fe(OH)$_3$','Fe$^{3\!\!+}$','Fe$^{2\!\!+}$','Sorbed Fe'],loc=(0.83,0.2),fontsize='small')

    axes[0].set_ylim(bottom=0.9e-11)
    # axes[0].legend()
    # CH4ax.legend()
    # axes[3].legend(loc='lower left')
    t=fig.axes[0].set_title('(a)',loc='left')  
    t=fig.axes[1].set_title('(b)',loc='left')  
    t=fig.axes[2].set_title('(c)',loc='left')  
    t=fig.axes[3].set_title('(d)',loc='left')  


    xmax=axes[0].get_xlim()[1]
    O2_periodic=zeros(arctic.simlength*24)
    O2_periodic[:int(arctic.simlength*24*arctic.oxicfrac)]=arctic.dq
    for ax in axes:
        for num in range(3):
            ax.axvspan(num*arctic.simlength+nonzero(diff(O2_periodic))[0]/24,(num+1)*arctic.simlength,color='b',alpha=0.1)
        ax.set_xlim(right=xmax,left=0)







    cm=get_cmap('plasma')
    norm=matplotlib.colors.Normalize(min(results_periodic['pH']),max(results_periodic['pH']))
    porosity=results_periodic['Porosity'].isel(time=0).mean().item()
    fig,axs=subplots(ncols=2,nrows=1,num='Periodic by dry periods',clear=True,squeeze=False,figsize=(7,4))
    for pH in results_periodic['pH'].values:
        axs[0,0].plot(results_periodic['ndryperiods'].values[:],array([results_periodic.sel(ndryperiods=n,pH=pH)['Total Tracer'].isel(time=arctic.simlength*24-1) for n in results_periodic['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'o-',label=pH,c=cm(norm(pH)))
        axs[0,1].plot(results_periodic['ndryperiods'].values[:],array([results_periodic.sel(ndryperiods=n,pH=pH)['Total CH4(aq)'].isel(time=arctic.simlength*24-1) for n in results_periodic['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'o-',label=pH,c=cm(norm(pH)))
        axs[0,0].plot(results_periodic_mineral['ndryperiods'].values[:],array([results_periodic_mineral.sel(ndryperiods=n,pH=pH)['Total Tracer'].isel(time=arctic.simlength*24-1) for n in results_periodic_mineral['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'+--',c=cm(norm(pH)))
        axs[0,1].plot(results_periodic_mineral['ndryperiods'].values[:],array([results_periodic_mineral.sel(ndryperiods=n,pH=pH)['Total CH4(aq)'].isel(time=arctic.simlength*24-1) for n in results_periodic_mineral['ndryperiods'].values[:]])*porosity/arctic.BD_layerest2_trough['Organic',False],'+--',c=cm(norm(pH)))


    axs[0,0].set_xlabel('Number of oxic-anoxic cycles')
    axs[0,1].set_xlabel('Number of oxic-anoxic cycles')
    axs[0,0].set_ylabel('Cumulative CO$_2$ flux\n(mmol g dwt$^{-1}$)')
    axs[0,1].set_ylabel('Cumulative CH$_4$ flux\n(mmol g dwt$^{-1}$)')
    axs[0,0].set_title('Cumulative CO$_2$ flux')
    axs[0,1].set_title('Cumulative CH$_4$ flux')
    axs[0,0].set_ylim(bottom=0)
    axs[0,1].set_ylim(bottom=0)
    axs[0,1].legend(title='Initial pH')
    axs[0,0].set_xticks(results_periodic['ndryperiods'].values[:])
    axs[0,1].set_xticks(results_periodic['ndryperiods'].values[:])
    axs[0,0].legend(handles=[axs[0,0].lines[0],axs[0,0].lines[1]],labels=['Organic','Mineral'])

    axs[0,0].set_title('(a)',loc='left')
    axs[0,1].set_title('(b)',loc='left')


    ndry_plotted=[3]
    fig,axs=subplots(ncols=len(ndry_plotted),nrows=6,num='Periodic comparison',clear=True,figsize=(8,12),squeeze=False)

    from string import ascii_lowercase
    for nn,ndry in enumerate(ndry_plotted):
        c='k'
        for nnn,pH in enumerate(results_periodic['pH'][:3]):
            ls=[':','-','--'][nnn]
            plot_result(results_periodic.sel(ndryperiods=ndry,pH=pH),CH4flux_ax=axs[0,nn],gdrywt=True,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={'color':c,'ls':ls},CO2flux_args={'color':c,'ls':ls})
            axs[2,nn].plot(results_periodic['time'][:-1],(results_periodic.sel(ndryperiods=ndry,pH=pH)['Total Fe++'].diff(dim='time')*24*1e3*porosity+results_periodic.sel(ndryperiods=ndry,pH=pH)['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Organic',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(results_periodic['time'],results_periodic.sel(ndryperiods=ndry,pH=pH)['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c=c)
            axs[4,nn].plot(results_periodic['time'],-log10(results_periodic.sel(ndryperiods=ndry,pH=pH)['Free H+']),ls,c=c)
            axs[5,nn].plot(results_periodic['time'],results_periodic['Total Acetate-'].sel(ndryperiods=ndry,pH=pH)*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic['time'][:-1],diff(results_periodic['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic['time'])*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c='b')
            CO2flux=diff(results_periodic['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            CH4flux=diff(results_periodic['Total CH4(aq)'].sel(ndryperiods=ndry,pH=pH)*results_periodic['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            axs[1,nn].plot(results_periodic['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)

        axs[0,0].legend(labels=results_periodic['pH'][:3].values,title='Initial pH')
        axs[0,nn].set_xlim(0,arctic.simlength)
        axs[1,nn].set_xlim(0,arctic.simlength)
        axs[1,nn].set_title('CO$_2$:CH$_4$ flux ratio')
        axs[1,nn].set_ylabel('CO$_2$:CH$_4$ ratio')
        axs[1,nn].set_ylim(-2,300)
        axs[1,nn].set_xlabel('Time (days)')
        axs[0,nn].set_title('CH$_4$ flux rate')

        axs[2,nn].set_ylim(bottom=0)
        axs[2,nn].set_xlabel('Time (days)')
        axs[2,nn].set_ylabel('Fe(II) production rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
        axs[2,nn].set_title('Fe(II) production rate')

        axs[3,nn].set_xlabel('Time (days)')
        axs[3,nn].set_ylabel('Fe oxides \n(mmol g dwt$^{-1}$)')
        axs[3,nn].set_title('Fe oxide minerals ')

        axs[4,nn].set_xlabel('Time (days)')
        axs[4,nn].set_ylabel('pH')
        axs[4,nn].set_title('pH')

        axs[5,nn].set_ylim(bottom=1e-5)
        axs[5,nn].set_xlabel('Time (days)')
        axs[5,nn].set_ylabel('Concentration (mM)')
        axs[5,nn].set_title('Acetate concentration')

        for x in range(6):
            for num in range(ndry):
                axs[x,nn].axvspan(arctic.simlength*arctic.oxicfrac/ndry+num*arctic.simlength/ndry,(num+1)*arctic.simlength/ndry,color='b',alpha=0.1,zorder=-1)
            axs[x,nn].set_xlim(0,arctic.simlength)
            axs[x,nn].set_title('('+ascii_lowercase[x*len(ndry_plotted)+nn]+')',loc='left')


    fig,axs=subplots(ncols=len(ndry_plotted),nrows=6,num='Periodic comparison mineral',clear=True,figsize=(8,12),squeeze=False)

    from string import ascii_lowercase
    for nn,ndry in enumerate(ndry_plotted):
        c='k'
        for nnn,pH in enumerate(results_periodic_mineral['pH'][:3]):
            ls=[':','-','--'][nnn]
            plot_result(results_periodic_mineral.sel(ndryperiods=ndry,pH=pH),CH4flux_ax=axs[0,nn],gdrywt=True,
                    BD=arctic.BD_layerest2_trough['Organic',False],SOC_pct=arctic.SOC_layermean['Organic'],CH4flux_args={'color':c,'ls':ls},CO2flux_args={'color':c,'ls':ls})
            axs[2,nn].plot(results_periodic_mineral['time'][:-1],(results_periodic_mineral.sel(ndryperiods=ndry,pH=pH)['Total Fe++'].diff(dim='time')*24*1e3*porosity+results_periodic_mineral.sel(ndryperiods=ndry,pH=pH)['Total Sorbed Fe++'].diff(dim='time')*24*1e6/100**3)/arctic.BD_layerest2_trough['Organic',False],ls,c=c,label=str(ndry))

            axs[3,nn].plot(results_periodic_mineral['time'],results_periodic_mineral.sel(ndryperiods=ndry,pH=pH)['Fe(OH)3 VF']/arctic.molar_volume_FeOH3*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c=c)
            axs[4,nn].plot(results_periodic_mineral['time'],-log10(results_periodic_mineral.sel(ndryperiods=ndry,pH=pH)['Free H+']),ls,c=c)
            axs[5,nn].plot(results_periodic_mineral['time'],results_periodic_mineral['Total Acetate-'].sel(ndryperiods=ndry,pH=pH)*1e3,c=c,ls=ls)
            # axs[1,nn].plot(results_periodic_mineral['time'][:-1],diff(results_periodic_mineral['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_mineral['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_mineral['time'])*1e3/arctic.BD_layerest2_trough['Organic',False],ls,c='b')
            CO2flux=diff(results_periodic_mineral['Total Tracer'].sel(ndryperiods=ndry,pH=pH)*results_periodic_mineral['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_mineral['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            CH4flux=diff(results_periodic_mineral['Total CH4(aq)'].sel(ndryperiods=ndry,pH=pH)*results_periodic_mineral['Porosity'].sel(ndryperiods=ndry,pH=pH,time=0))/diff(results_periodic_mineral['time'])*1e3/arctic.BD_layerest2_trough['Organic',False]
            axs[1,nn].plot(results_periodic_mineral['time'][:-1],CO2flux/CH4flux,ls=ls,c=c)

        axs[0,0].legend(labels=results_periodic_mineral['pH'][:3].values,title='Initial pH')
        axs[0,nn].set_xlim(0,arctic.simlength)
        axs[1,nn].set_xlim(0,arctic.simlength)
        axs[1,nn].set_title('CO$_2$:CH$_4$ flux ratio')
        axs[1,nn].set_ylabel('CO$_2$:CH$_4$ ratio')
        axs[1,nn].set_ylim(-2,300)
        axs[1,nn].set_xlabel('Time (days)')
        axs[0,nn].set_title('CH$_4$ flux rate')

        axs[2,nn].set_ylim(bottom=0)
        axs[2,nn].set_xlabel('Time (days)')
        axs[2,nn].set_ylabel('Fe(II) production rate\n($\mu$mol g dwt$^{-1}$ day$^{-1}$)')
        axs[2,nn].set_title('Fe(II) production rate')

        axs[3,nn].set_xlabel('Time (days)')
        axs[3,nn].set_ylabel('Fe oxides \n(mmol g dwt$^{-1}$)')
        axs[3,nn].set_title('Fe oxide minerals ')

        axs[4,nn].set_xlabel('Time (days)')
        axs[4,nn].set_ylabel('pH')
        axs[4,nn].set_title('pH')

        axs[5,nn].set_ylim(bottom=1e-5)
        axs[5,nn].set_xlabel('Time (days)')
        axs[5,nn].set_ylabel('Concentration (mM)')
        axs[5,nn].set_title('Acetate concentration')

        for x in range(6):
            for num in range(ndry):
                axs[x,nn].axvspan(arctic.simlength*arctic.oxicfrac/ndry+num*arctic.simlength/ndry,(num+1)*arctic.simlength/ndry,color='b',alpha=0.1,zorder=-1)
            axs[x,nn].set_xlim(0,arctic.simlength)
            axs[x,nn].set_title('('+ascii_lowercase[x*len(ndry_plotted)+nn]+')',loc='left')


    save_all_figs('Arctic_Fe_output/figs')

    show()