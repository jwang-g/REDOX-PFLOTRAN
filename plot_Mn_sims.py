import matplotlib
matplotlib.use('Agg')
from pylab import *
import xarray
import sys
import netCDF4
import pandas
import numpy

filenames=sys.argv[1:-1]

if sys.argv[-1].endswith('.nc'):
    filenames.append(sys.argv[-1])
    figdir='Mn_output/Figures'
else:
    figdir=sys.argv[-1]


# result_files = [
# 'Mn_ph35_saved_2020-08-31.nc',
# 'Mn_ph40_saved_2020-08-31.nc',
#   'Mn_ph45_saved_2020-08-31.nc',
#   'Mn_ph50_saved_2020-08-31.nc',
#   'Mn_ph55_saved_2020-08-31.nc',
# ]

leafC_mass_conv=1.0 # For simulations after I started applying the C to dry mass conversion in the sims. Otherwise it should be 0.4
molar_volume_birnessite = 251.1700 # Divide by 7 below because birnessite is defined in database as 7 mols of Mn
Mn_molarmass=54.94        #g/mol


def getval_from_str(s,subs,valtype=float):
    return valtype(s[s.index(subs)+len(subs):].partition('_')[0])

def load_data(filenames):
    data_list=[]
    incubation_list=[]

    allfilenames=[]
    groupnames=[]

    for filename in sorted(filenames):
        print('Reading file %s'%filename)
        with netCDF4.Dataset(filename) as d:
            groups=list(d.groups.keys())


        for g in groups:
            print(g,flush=True)

            d=xarray.open_dataset(filename,group=g,chunks={'time':3650},decode_times=False)
            if 'soil_pH' not in d:
                d['soil_pH']=getval_from_str(g,'pH')
                d['Ndep']=getval_from_str(g,'Ndep',int)
                d['warming']=getval_from_str(g,'warming')
                d['redox_cycles']=getval_from_str(g,'anox_freq',int)
                d['anox_lenscales']=getval_from_str(g,'anox_len')
                newdims=['soil_pH','Ndep','warming','redox_cycles']
                d=d.expand_dims(newdims).set_coords(newdims)
            if 'incubation' in g:
                incubation_list.append(d)
            else:
                data_list.append(d)
                allfilenames.append(filename)
                groupnames.append(g)

    pH_sims=array([d['soil_pH'].item() for d in data_list])
    Ndep_sims=array([d['Ndep'].item() for d in data_list])
    warming_sims=array([d['warming'].item() for d in data_list])
    anox_freq_sims=array([d['redox_cycles'].item() for d in data_list])
    anox_len_sims=array([d['anox_lenscales'].item() for d in data_list])
    allsims=pandas.DataFrame({'Ndep':Ndep_sims,'warming':warming_sims,'anox_freq':anox_freq_sims,'anox_len':anox_len_sims,'pH':pH_sims})

    from Manganese_profile import setup_sims
    expected_sims=setup_sims()
    expected_sims['simnum']=expected_sims.index
    comp=allsims.merge(expected_sims,how='outer',indicator=True)
    missing=comp[comp['_merge']!='both']
    if len(missing)>0:
        print('Warning: Not all expected sims were found:')
        print(missing)
        print('Generating all nan data for missing simulations.')
        for num in missing.index:
            d=(data_list[0]*nan).assign_coords(
                soil_pH=atleast_1d(missing.loc[num]['pH']),
                Ndep=atleast_1d(missing.loc[num]['Ndep']),
                warming=atleast_1d(missing.loc[num]['warming']),
                redox_cycles=atleast_1d(missing.loc[num]['anox_freq']),
                anox_lenscales=atleast_1d(missing.loc[num]['anox_len']))
            data_list.append(d)

        pH_sims=array([d['soil_pH'].item() for d in data_list])
        Ndep_sims=array([d['Ndep'].item() for d in data_list])
        warming_sims=array([d['warming'].item() for d in data_list])
        anox_freq_sims=array([d['redox_cycles'].item() for d in data_list])
        anox_len_sims=array([d['anox_lenscales'].item() for d in data_list])


    x=(anox_freq_sims>0)&(anox_len_sims>0)
    if len(pH_sims[x])%len(unique(pH_sims)) != 0:
        raise ValueError("Simulations don't appear to divide evenly by pH and anox_freq!")

    alldata=xarray.combine_by_coords(data_list[xx] for xx in nonzero(x)[0])
    

    print('Finished combining data',flush=True)
    return alldata,incubation_list

if len(filenames)>0:
    alldata,incubation_list=load_data(filenames)
controldata=alldata.isel(Ndep=0,warming=0)

def getdata(soil_pH,Ndep,warming,redox_cycles=50,anox_lenscale=0.5):
    # return data_list[flatnonzero((pH_sims==soil_pH)&(Ndep_sims==Ndep)&(warming_sims==warming)&(anox_freq_sims==redox_cycles)&(anox_len_sims==anox_lenscale))[0]].squeeze()
    return alldata.sel(soil_pH=soil_pH,Ndep=Ndep,warming=warming,anox_lenscales=anox_lenscale,redox_cycles=redox_cycles)

def save_one_fig(f,dirname=figdir,format='png',**kwargs):
    fname_fixed=f.get_label().replace('/','-')
    savename='{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed)
    print(savename)
    f.savefig(savename,**kwargs)

def letter_label(ax,label=None,x=0.03,y=1.03,title=True,**kwargs):
    from string import ascii_lowercase
    if isinstance(label,int):
        label='('+ascii_lowercase[label]+')'
    if title:
        return ax.set_title(label,loc='left')
    else:
        return ax.text(x,y,label,transform=ax.transAxes,**kwargs)

dt=(controldata['time'][1]-controldata['time'][0]).item()
oneyr=int(365/dt)

anox_baseline_i=2
anox_baseline=controldata['anox_lenscales'][anox_baseline_i].load().item()

f,axs=subplots(num='Figure S1 Saturated time fraction',clear=True)
O2=controldata['Total O2(aq)'].sel(soil_pH=5.0).squeeze()
satfrac=(O2<O2.max()*0.9).mean(dim='time').T
h=axs.pcolormesh(satfrac['anox_lenscales'],-satfrac['depth'],satfrac*100,cmap='RdBu',vmin=0,vmax=100,shading='auto')
cb=f.colorbar(h,ax=axs)
cb.set_label('Saturated fraction of time (%)')
cb.set_ticks([0,25,50,75,100])
axs.set(title='Water-saturated fraction of time',xlabel='Drainage time scale (days)',ylabel='Depth (cm)')

save_one_fig(f)

total_litter_all=((alldata['Total Sorbed Cellulose']+alldata['Total Sorbed Lignin']).coarsen(time=oneyr,boundary='trim').mean()*12e-3*(alldata.z_bottom-alldata.z_top)/100).sum(dim='depth')

results_byph=controldata.sel(anox_lenscales=0.5)
results_byredox=controldata.squeeze()
pHs=controldata['soil_pH'].to_masked_array()


total_MAOM=((results_byredox['Total DOM3'].coarsen(time=oneyr,boundary='trim').mean()*results_byredox.saturation*results_byredox.Porosity.mean())*12e-3/1000*100**3*(results_byredox.z_bottom-results_byredox.z_top)/100).sum(dim='depth')
total_SOC=(((results_byredox['Total DOM3']+results_byredox['Total DOM1']+results_byredox['Total DOM2']).\
            coarsen(time=oneyr,boundary='trim').mean()*results_byredox.saturation*results_byredox.Porosity.mean())\
                *12e-3/1000*100**3*(results_byredox.z_bottom-results_byredox.z_top)/100).sum(dim='depth')



total_litter=((results_byredox['Total Sorbed Cellulose']+results_byredox['Total Sorbed Lignin']).load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(results_byredox.z_bottom-results_byredox.z_top)/100).sum(dim='depth').load()
Figure3,axs=subplots(nrows=3,num='Figure 3',clear=True,figsize=(3,8))
for t in [9,19,29,39]:
    axs[0].plot(results_byredox['soil_pH'],total_litter.isel(time=t,anox_lenscales=anox_baseline_i),'o-',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))
    axs[1].plot(results_byredox['soil_pH'],total_MAOM.isel(time=t,anox_lenscales=anox_baseline_i),'o-',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))
    axs[2].plot(results_byredox['soil_pH'],(total_SOC+total_litter).isel(time=t,anox_lenscales=anox_baseline_i),'o-',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))

print('litter stock max = %1.2f, min = %1.2f'%(total_litter.isel(time=t,anox_lenscales=anox_baseline_i).max(),total_litter.isel(time=t,anox_lenscales=anox_baseline_i).min()))

axs[0].set(title='Litter C',xlabel='pH',ylabel='C stock (kg C m$^{-2}$)',xticks=linspace(4,6,5))
axs[1].set(title='Mineral-associated Organic C',xlabel='pH',ylabel='C stock (kg C m$^{-2}$)',xticks=linspace(4,6,5))
axs[2].set(title='Total C (litter + soil)',xlabel='pH',ylabel='C stock (kg C m$^{-2}$)',xticks=linspace(4,6,5))
axs[0].legend()

for num in range(len(axs)):
    letter_label(axs[num],num)

save_one_fig(Figure3)

Figure3_expanded,axs=subplots(nrows=3,ncols=4,num='Figure S5 C stocks in all scenarios',clear=True,figsize=(12,10),sharey='row')
d=alldata.isel(warming=0,Ndep=0)
# MAOM_control=((d['Total DOM3'].coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())*12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
# totC_control=(d['Total Sorbed Cellulose'].coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth')+\
#             (d['Total Sorbed Lignin'].coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth') +\
#                 (((d['Total DOM3']+d['Total DOM1']+d['Total DOM2']).\
#             coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())\
#                 *12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
for num in range(4):
    w=[1,2,0,0][num]
    n=[0,0,1,2][num]
    d=alldata.isel(warming=w,Ndep=n)
    total_cellulose=(d['Total Sorbed Cellulose'].load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
    total_lignin=(d['Total Sorbed Lignin'].load().coarsen(time=oneyr,boundary='trim').mean()*12e-3*(d.z_bottom-d.z_top)/100).sum(dim='depth')

    total_MAOM=((d['Total DOM3'].load().coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())*12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
    total_SOC=(((d['Total DOM3']+d['Total DOM1']+d['Total DOM2']).load().\
            coarsen(time=oneyr,boundary='trim').mean()*d.saturation*d.Porosity.mean())\
                *12e-3/1000*100**3*(d.z_bottom-d.z_top)/100).sum(dim='depth')
    MAOM_control=total_MAOM*0.0
    totC_control=total_SOC*0.0

    for t in [9,19,29,39]:
        print(num,t)
        axs[0,num].plot(d['soil_pH'],(total_cellulose+total_lignin).isel(time=t,anox_lenscales=anox_baseline_i),'o-',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))
        axs[1,num].plot(d['soil_pH'],total_MAOM.isel(time=t,anox_lenscales=anox_baseline_i)-MAOM_control.isel(time=t,anox_lenscales=anox_baseline_i,soil_pH=-1),'o-',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))
        axs[2,num].plot(d['soil_pH'],(total_SOC+total_cellulose+total_lignin).isel(time=t,anox_lenscales=anox_baseline_i)-totC_control.isel(time=t,anox_lenscales=anox_baseline_i,soil_pH=-1),'o-',label='Year %d'%(t+1),c='%1.1f'%(1-t/39))

    axs[0,num].set(title='Warming = %d, Ndep = %d\nLitter C stock'%(d['warming'].item(),d['Ndep'].item()),xlabel='pH',ylabel='C stock (kg C m$^{-2}$)')
    axs[1,num].set(title='MAOM C stock',xlabel='pH',ylabel='C stock (kg C m$^{-2}$)')
    axs[2,num].set(title='Total soil and litter C stock',xlabel='pH',ylabel='C stock (kg C m$^{-2}$)')
    axs[0,num].tick_params(labelleft=True)
    axs[1,num].tick_params(labelleft=True)
    axs[2,num].tick_params(labelleft=True)
    axs[0,num].legend()

for num in range(axs.size):
    letter_label(axs.ravel()[num],num)

save_one_fig(Figure3_expanded)



totalMn=((results_byredox['Total Mn++']+results_byredox['Total Mn+++'])/1000*results_byredox.Porosity.mean()*results_byredox.saturation + \
            (results_byredox['Birnessite2 VF']*7/molar_volume_birnessite) + \
            results_byredox['Total Sorbed Mn++']/100**3).load().coarsen(time=oneyr,boundary='trim').mean()*100**2*(results_byredox.z_bottom-results_byredox.z_top)
birnessite=results_byredox['Birnessite2 VF'].load().coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*100**2*(results_byredox.z_bottom-results_byredox.z_top)

f,axs=subplots(nrows=len(controldata['soil_pH']),ncols=len(controldata['anox_lenscales']),clear=True,num='Figure S4 Mn bioavailability',figsize=(12,10),sharex=True,sharey=True)
maxval=1.0
minval=1e-5
for ph in range(len(controldata['soil_pH'])):
    for anox in range(len(controldata['anox_lenscales'])):
        # Total Mn in mol/m2
        row=4-anox
        h=axs[row,ph].pcolormesh(totalMn['time']/365,
                -results_byredox['z_middle'].isel(soil_pH=0,anox_lenscales=0),1-birnessite.isel(soil_pH=ph,anox_lenscales=anox)/totalMn.isel(soil_pH=ph,anox_lenscales=anox),
                shading='auto',cmap=get_cmap('viridis'),norm=matplotlib.colors.LogNorm(vmin=minval,vmax=maxval))
        axs[row,ph].set_xlabel('Time (years)')
        if (row)==0:
            axs[row,ph].set_title('pH = %1.1f'%controldata['soil_pH'][ph].load().item(),pad=10,fontsize='large',fontweight='bold')
        if ph==0:
            axs[row,ph].text(-0.6,0.5,f'Drainage time\nscale = {controldata["anox_lenscales"][anox].load().item()} days',
            rotation=90,va='center',transform=axs[row,ph].transAxes,fontsize='large',fontweight='bold')

        axs[row,ph].set_ylabel('Depth (cm)')
        axs[row,ph].tick_params(labelleft=True,labelbottom=True)
            
cb=f.colorbar(h,ax=axs)
cb.set_label('Bioavailable fraction of Mn',fontsize='large')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=-0.1)

save_one_fig(f)

f,axs=subplots(nrows=len(controldata['soil_pH']),ncols=len(controldata['anox_lenscales']),clear=True,num='Figure S3 Mn redistribution',figsize=(12,10),sharex=True,sharey=True)
maxval=(totalMn-totalMn.isel(time=0)).max().compute()
minval=(totalMn-totalMn.isel(time=0)).min().compute()
for ph in range(len(controldata['soil_pH'])):
    for anox in range(len(controldata['anox_lenscales'])):
        row=4-anox
        # Total Mn in mol/m2
        h=axs[row,ph].pcolormesh(totalMn['time']/365,
                -results_byredox['z_middle'].isel(soil_pH=0,anox_lenscales=0),totalMn.isel(soil_pH=ph,anox_lenscales=anox)-totalMn.isel(soil_pH=ph,anox_lenscales=anox,time=0),
                shading='auto',vmin=min(minval,-maxval),vmax=max(maxval,-minval),cmap=get_cmap('RdBu_r'))
        axs[row,ph].set_xlabel('Time (years)')
        if row==0:
            axs[row,ph].set_title('pH = %1.1f'%controldata['soil_pH'][ph].load().item(),pad=10,fontsize='large',fontweight='bold')
        if ph==0:
            axs[row,ph].text(-0.6,0.5,f'Drainage time\nscale = {controldata["anox_lenscales"][anox].load().item()} days',
                rotation=90,va='center',transform=axs[row,ph].transAxes,fontsize='large',fontweight='bold')

        axs[row,ph].set_ylabel('Depth (cm)')
        axs[row,ph].tick_params(labelleft=True,labelbottom=True)
            
cb=f.colorbar(h,ax=axs)
cb.set_label('Change in Mn stock (mol m$^{-2}$)',fontsize='large')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=-0.1)

save_one_fig(f)


molar_volume_manganite = 24.45 # cm3/mol
molar_volume_MnOH2am = 22.3600
molar_volume_birnessite = 251.1700
Mn_molarmass=54.94   


def plot_output(output,axs,subsample=24,do_legend=False,**kwargs):
    def fixdata(data):
        return data.coarsen(time=subsample,boundary='trim').mean().dropna(dim='time')
    outdata=output[['Total Sorbed Cellulose','Total Sorbed Lignin','Total DOM1','Total DOM3','Total Mn++',
                    'Total Mn+++','Birnessite2 VF','Total Tracer2','Free H+','Porosity','BD','saturation','CEC H+',
                    'Total Sorbed Al+++','Total Sorbed Ca++','Total Sorbed K+','Total Sorbed Mg++','Total Sorbed Na+','Total Sorbed Mn++','Total Sorbed Mn+++']].load()
    for num in range(len(output.depth)):
        out=outdata.isel(depth=num)
        t=fixdata(out.time/(365))
        porosity=out['Porosity'].mean()
        saturation=out['saturation']
        BD=out['BD']
        axs[num,0].plot(t,fixdata(out['Total Sorbed Cellulose'])*12/100**3,label='Cellulose',c='C0',**kwargs)
        axs[num,0].plot(t,fixdata(out['Total Sorbed Lignin'])*12/100**3,label='Lignin',c='C1',**kwargs)
        axs[num,0].plot(t,fixdata(out['Total DOM1']*12/1000*porosity*saturation),c='C2',label='DOM',**kwargs)
        # axs[num,0].plot(t,out['Total Sorbed DOM1']*12/100**3,label='Sorbed DOM')
        axs[num,0].plot(t,fixdata(out['Total DOM3']*12/1000*porosity*saturation),c='C3',label='Sorbed DOM',**kwargs)
        axs[num,0].set_ylabel('C density\n(g C cm$^{-3}$)')
        
        axs[num,1].plot(t,fixdata(out['Total Mn++']*1e6/1000*porosity*saturation*Mn_molarmass/BD),label='Mn$^{+\!\!+}$',c='C0',**kwargs)
        axs[num,1].plot(t,fixdata(out['Total Mn+++']*1e6/1000*porosity*saturation*Mn_molarmass/BD),label='Mn$^{+\!\!+\!\!+}$',c='C1',**kwargs)
        # axs[num,1].plot(t,layers[num].output_DF['Manganite VF']/molar_volume_manganite*1e6*Mn_molarmass/BD,label='Manganite')
        axs[num,2].plot(t,fixdata(out['Birnessite2 VF']*7/molar_volume_birnessite*1e6*Mn_molarmass/BD),label='Birnessite',c='C0',**kwargs)
        
        # annual_root_uptake=fixdata(out['Total Tracer2']).groupby(floor(t)).max() 
        annual_root_uptake=(out['Total Tracer2']).coarsen(time=oneyr,boundary='trim').max()
        axs[num,1].plot(annual_root_uptake.time/365,annual_root_uptake*1e6/1000*porosity.mean()*saturation*Mn_molarmass/BD,label='Annual Mn$^{+\!\!+}$ root uptake',c='C2',**kwargs)
        
        axs[num,3].plot(t,-log10(fixdata(out['Free H+'])),label='pH',c='C0',**kwargs)
        
        # ax=axs[num,1].twinx()
        # axs[num,2].plot(t,layers[num].output_DF['Mn(OH)2(am) VF']/molar_volume_MnOH2am*1e6*Mn_molarmass/BD,label='Mn(OH)$_2$(am)',ls='-')
        axs[num,2].plot(t,fixdata(out['Total Mn++']+out['Total Mn+++'])*1e6/1000*porosity*saturation*Mn_molarmass/BD + 
                            fixdata(out['Birnessite2 VF']*7/molar_volume_birnessite)*1e6*Mn_molarmass/BD,c='k',label='Total Mn',**kwargs)
        
        axs[num,1].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
        axs[num,2].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
        axs[num,3].set_ylabel('pH')
        # axs[num,1].set_ylim(*axs[0,1].get_ylim())

        for n,cation in enumerate(['Al+++','Ca++','K+','Mg++','Na+','Mn++','Mn+++']):
            axs[num,4].plot(t,fixdata(out['Total Sorbed '+cation])/(BD*1e-3*100**3)*1000,c='C'+str(n),label=cation,**kwargs)
        axs[num,4].plot(t,fixdata(out['CEC H+'])/(BD*1e-3*100**3)*1000,label='H+',c='C'+str(n+1),**kwargs)
        
        axs[num,4].set_ylabel('Exch conc\n(mmol/kg)')
        

        
    axs[-1,0].set_xlabel('Time (years)')
    if do_legend:
        axs[0,0].legend()
        axs[0,1].legend()
        axs[0,2].legend()
        axs[0,4].legend()

    axs[-1,1].set_xlabel('Time (years)')
    axs[-1,2].set_xlabel('Time (years)')
    axs[-1,3].set_xlabel('Time (years)')
    axs[-1,4].set_xlabel('Time (years)')
    axs[0,0].set_title('Organic Carbon')
    axs[0,1].set_title('Manganese ions')
    axs[0,2].set_title('Total Mn')
    axs[0,3].set_title('pH')
    axs[0,4].set_title('Exchangeable cations')


data_allNdeps=alldata.sel(warming=0)



data_warmings=alldata.sel(Ndep=0).squeeze()

warmings=data_warmings['warming'].to_masked_array()


f,axs=subplots(nrows=1,ncols=3,num='Figure 4 Warming time series',clear=True,sharey='row',figsize=(12,4))
litter=total_litter_all.sel(Ndep=0).squeeze()

h=axs[0].plot(litter.time/365,litter.isel(warming=0,soil_pH=0,anox_lenscales=anox_baseline_i),'b-',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][0],warmings[0]))
h=axs[0].plot(litter.time/365,litter.isel(warming=0,soil_pH=-1,anox_lenscales=anox_baseline_i),'b--',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][-1],warmings[0]))
h=axs[0].plot(litter.time/365,litter.isel(warming=1,soil_pH=0,anox_lenscales=anox_baseline_i),'y-',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][0],warmings[1]))
h=axs[0].plot(litter.time/365,litter.isel(warming=1,soil_pH=-1,anox_lenscales=anox_baseline_i),'y--',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][-1],warmings[1]))
h=axs[0].plot(litter.time/365,litter.isel(warming=2,soil_pH=0,anox_lenscales=anox_baseline_i),'r-',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][0],warmings[2]))
h=axs[0].plot(litter.time/365,litter.isel(warming=2,soil_pH=-1,anox_lenscales=anox_baseline_i),'r--',label='pH = %1.1f, Warming = %d\u00B0C'%(data_warmings['soil_pH'][-1],warmings[2]))

litter=total_litter_all.isel(Ndep=1).squeeze()

h=axs[1].plot(litter.time/365,litter.isel(warming=0,soil_pH=0,anox_lenscales=anox_baseline_i),'b-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[0]))
h=axs[1].plot(litter.time/365,litter.isel(warming=0,soil_pH=-1,anox_lenscales=anox_baseline_i),'b--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[0]))
h=axs[1].plot(litter.time/365,litter.isel(warming=1,soil_pH=0,anox_lenscales=anox_baseline_i),'y-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[1]))
h=axs[1].plot(litter.time/365,litter.isel(warming=1,soil_pH=-1,anox_lenscales=anox_baseline_i),'y--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[1]))
h=axs[1].plot(litter.time/365,litter.isel(warming=2,soil_pH=0,anox_lenscales=anox_baseline_i),'r-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[2]))
h=axs[1].plot(litter.time/365,litter.isel(warming=2,soil_pH=-1,anox_lenscales=anox_baseline_i),'r--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[2]))


litter=total_litter_all.isel(Ndep=2).squeeze()

h=axs[2].plot(litter.time/365,litter.isel(warming=0,soil_pH=0,anox_lenscales=anox_baseline_i),'b-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[0]))
h=axs[2].plot(litter.time/365,litter.isel(warming=0,soil_pH=-1,anox_lenscales=anox_baseline_i),'b--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[0]))
h=axs[2].plot(litter.time/365,litter.isel(warming=1,soil_pH=0,anox_lenscales=anox_baseline_i),'y-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[1]))
h=axs[2].plot(litter.time/365,litter.isel(warming=1,soil_pH=-1,anox_lenscales=anox_baseline_i),'y--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[1]))
h=axs[2].plot(litter.time/365,litter.isel(warming=2,soil_pH=0,anox_lenscales=anox_baseline_i),'r-',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][0],warmings[2]))
h=axs[2].plot(litter.time/365,litter.isel(warming=2,soil_pH=-1,anox_lenscales=anox_baseline_i),'r--',label='pH = %1.1f, Warming = %d'%(data_warmings['soil_pH'][-1],warmings[2]))

axs[2].set_ylabel('Total litter C (kg C m$^{-2}$)')
axs[2].set_xlabel('Year')
axs[2].set_title('N dep = 150 kg N ha$^{-1}$ year$^{-1}$')
axs[2].tick_params(labelleft=True)

axs[1].set_ylabel('Total litter C (kg C m$^{-2}$)')
axs[1].set_xlabel('Year')
axs[1].set_title('N dep = 50 kg N ha$^{-1}$ year$^{-1}$')
axs[1].tick_params(labelleft=True)

axs[0].set_ylabel('Total litter C (kg C m$^{-2}$)')
axs[0].set_xlabel('Year')
axs[0].set_title('N dep = 0 kg N ha$^{-1}$ year$^{-1}$')
axs[0].legend(ncol=2)
axs[0].tick_params(labelleft=True)

for num in range(axs.size):
    letter_label(axs.ravel()[num],num)

save_one_fig(f)

Figure2,axs=subplots(nrows=2,ncols=3,num='Figure 2 (warming)',clear=True,figsize=(12,5))
total_litter=total_litter_all.sel(Ndep=0).squeeze()
h=axs[0,0].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),total_litter.isel(warming=0,time=39).T,levels=linspace(0.3,1.2,21))
print('Litter min = %1.2f, max = %1.2f'%(total_litter.isel(warming=0,time=39).min(),total_litter.isel(warming=0,time=39).max()))
h1=axs[0,1].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=1,time=39).T-total_litter.isel(warming=0,time=39).T)/total_litter.isel(warming=0,time=39).T*100,
        levels=linspace(-50,50,21),cmap='RdBu_r')
h2=axs[0,2].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=2,time=39).T-total_litter.isel(warming=0,time=39).T)/total_litter.isel(warming=0,time=39).T*100,
        levels=linspace(-50,50,21),cmap='RdBu_r')

cb=Figure2.colorbar(h,ax=axs[0,0])
cb.set_label('Total litter C (kg C m$^{-2}$)')
cb.set_ticks([0.3,0.5,0.7,0.9,1.1])
cb=Figure2.colorbar(h2,ax=axs[0,2])
cb.set_label('Total litter C (% diff.)')

h=axs[1,0].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),data_warmings['litter_Mn'].isel(warming=0,litter_year=39).T,levels=linspace(0,120.0,11))
h1=axs[1,1].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
           ( data_warmings['litter_Mn'].isel(warming=1,litter_year=39).T-data_warmings['litter_Mn'].isel(warming=0,litter_year=39).T)/data_warmings['litter_Mn'].isel(warming=0,litter_year=39).T*100,
            levels=linspace(-95,95,31),cmap='RdBu_r')
h2=axs[1,2].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
           ( data_warmings['litter_Mn'].isel(warming=2,litter_year=39).T-data_warmings['litter_Mn'].isel(warming=0,litter_year=39).T)/data_warmings['litter_Mn'].isel(warming=0,litter_year=39).T*100,
            levels=linspace(-95,95,31),cmap='RdBu_r')

cb=Figure2.colorbar(h,ax=axs[1,0])
cb.set_label('Leaf Mn (mmol kg$^{-1}$)')
cb=Figure2.colorbar(h2,ax=axs[1,2])
cb.set_label('Leaf Mn (% diff.)')
cb.set_ticks([-75,-50,-25,0,25,50,75])

axs[0,0].set(title='Litter C stock (+0\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction')
axs[0,1].set(title='Change in litter C stock  (+2\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction')
axs[0,2].set(title='Change in litter C stock  (+5\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction')

axs[1,0].set(title='Leaf Mn concentration (+0\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction')
axs[1,1].set(title='Change in leaf Mn  (+2\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction')
axs[1,2].set(title='Change in leaf Mn  (+5\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num,x=0,title=True)

axs[0,0].text(0,1.2,'No Warming',fontsize='large',fontweight='bold',transform=axs[0,0].transAxes)
axs[0,0].text(-0.25,0.5,'Litter C Stock',fontsize='large',fontweight='bold',transform=axs[0,0].transAxes,rotation=90,va='center')
axs[0,1].text(0,1.2,'Warming = (+2\u00B0C)',fontsize='large',fontweight='bold',transform=axs[0,1].transAxes)
axs[0,2].text(0,1.2,'Warming = (+5\u00B0C)',fontsize='large',fontweight='bold',transform=axs[0,2].transAxes)
axs[1,0].text(-0.25,0.5,'Leaf Mn',fontsize='large',fontweight='bold',transform=axs[1,0].transAxes,rotation=90,va='center')

save_one_fig(Figure2)

total_litter=total_litter_all.squeeze()


f,axs=subplots(nrows=2,ncols=3,num='Figure 5 Ndep contours',clear=True,figsize=(12,5))
levs=linspace(0,1,21)
cm=get_cmap('Reds')
# h=axs[0,0].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),total_litter.isel(warming=0,Ndep=1,time=39).T,levels=linspace(0,2.0,21))
h=axs[0,0].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=0,Ndep=1,time=39).T-total_litter.isel(warming=0,Ndep=0,time=39).T),
        levels=levs,cmap=cm)
h1=axs[0,1].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=1,Ndep=1,time=39).T-total_litter.isel(warming=0,Ndep=0,time=39).T),
        levels=levs,cmap=cm)
h2=axs[0,2].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=2,Ndep=1,time=39).T-total_litter.isel(warming=0,Ndep=0,time=39).T),
        levels=levs,cmap=cm)

cb=f.colorbar(h2,ax=axs)
# cb.set_label('Total litter C diff. (%)')
cb.set_label('Total litter C diff. (kg C m$^{-2}$)')

# h=axs[1,0].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),total_litter.isel(warming=0,Ndep=2,time=39).T,levels=linspace(0,2.0,21))
h=axs[1,0].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=0,Ndep=2,time=39).T-total_litter.isel(warming=0,Ndep=0,time=39).T),
        levels=levs,cmap=cm)
h1=axs[1,1].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=1,Ndep=2,time=39).T-total_litter.isel(warming=0,Ndep=0,time=39).T),
        levels=levs,cmap=cm)
h2=axs[1,2].contourf(total_litter['soil_pH'],satfrac.mean(dim='depth'),
        (total_litter.isel(warming=2,Ndep=2,time=39).T-total_litter.isel(warming=0,Ndep=0,time=39).T),
        levels=levs,cmap=cm)


axs[0,0].set(title='Change in litter C stock  (+0\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction',xticks=linspace(4,6,5))
axs[0,1].set(title='Change in litter C stock  (+2\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction',xticks=linspace(4,6,5))
axs[0,2].set(title='Change in litter C stock  (+5\u00B0C)',xlabel='Soil pH',ylabel='Mean saturated soil fraction',xticks=linspace(4,6,5))
axs[1,0].set(xlabel='Soil pH',ylabel='Mean saturated soil fraction',xticks=linspace(4,6,5))
axs[1,1].set(xlabel='Soil pH',ylabel='Mean saturated soil fraction',xticks=linspace(4,6,5))
axs[1,2].set(xlabel='Soil pH',ylabel='Mean saturated soil fraction',xticks=linspace(4,6,5))

axs[0,0].text(0,1.2,'No Warming',fontsize='large',fontweight='bold',transform=axs[0,0].transAxes)
axs[0,0].text(-0.25,0.5,'Low N deposition',fontsize='large',fontweight='bold',transform=axs[0,0].transAxes,rotation=90,va='center')
axs[0,1].text(0,1.2,'Warming = (+2 \u00B0C)',fontsize='large',fontweight='bold',transform=axs[0,1].transAxes)
axs[0,2].text(0,1.2,'Warming = (+5 \u00B0C)',fontsize='large',fontweight='bold',transform=axs[0,2].transAxes)
axs[1,0].text(-0.25,0.5,'High N deposition',fontsize='large',fontweight='bold',transform=axs[1,0].transAxes,rotation=90,va='center')

for num in range(axs.size):
    letter_label(axs.ravel()[num],num)

save_one_fig(f)



import Manganese_network as Mn
import decomp_network
reaction_network=Mn.make_network(leaf_Mn_mgkg=0.0,Mn2_scale=1e-4,Mn_peroxidase_Mn3_leakage=1e-3,Mn3_scale=1e-12) 

networkfig=figure('Reaction network',clear=True)
pos={'Birnessite2': (644.0, 594.0),
 'Mn+++': (635.0, 522.0),
 'O2(aq)': (865.0, 450.0),
 'Cellulose': (72.0, 378.0),
 'DOM1': (346.0, 234.0),
 'DOM2': (130.0, 378.0),
 'HCO3-': (667.0, 90.0),
 'Mn++': (683.0, 378.0),
 'Tracer2': (552.0, 234.0),
 'Mg++': (204.0, 505.0),
 'Ca++': (284.0, 505.0),
 'Na+': (359.0, 505.0),
 'K+': (431.0, 505.0),
 'Al+++': (120.0, 505.0),
 'DOM3': (203.0, 90.0),
 'Hydrolysis': (147.0, 306.0),
 'DOM aerobic respiration': (884.0, 162.0),
 'DOM1 Mn+++ reduction': (667.0, 162.0),
 'DOM1 Mn+++ abiotic reduction': (422.0, 162.0),
 'DOM sorption': (216.0, 162.0),
 'Mn Peroxidase': (346.0, 306.0),
 'Bacterial Mn++ oxidation': (873.0, 306.0),
 'Cation exchange': (335.0, 420.0),
 'Root uptake of Mn++': (552.0, 306.0),
 'DOM desorption': (164.0, 18.0)}
drawn,pos=decomp_network.draw_network_with_reactions(reaction_network,
        omit=['NH4+','Rock(s)','gas','secondary','H+','>Carboxylate-','Carboxylic_acid','H2O','Root_biomass',
            'Sorption_capacity','Lignin depolymerization','Lignin exposure','Lignin','DOM2 aerobic respiration'],
        font_size='medium',node_size=1500,font_color='k',arrowstyle='->',arrowsize=10.0,edge_color='gray',node_alpha=1.0,
        namechanges={'cellulose':'Cellulose','DOM1':'DOM','O2(aq)':'O$_2$(aq)','CH4(aq)':'CH$_4$(aq)','HCO3-':'CO$_2$','DOM2':'Lignin','sorbed_DOM1':'Sorbed DOM',
                     'Fe(OH)2':'Fe(OH)$_2$','Fe(OH)3':'Fe(OH)$_3$','Mn++':r'Mn$^\mathrm{+\!\!+}$','Mn+++':r'Mn$^\mathrm{+\!\!+\!\!\!+}$','Acetate-':'Acetate',
                     'DOM3':'MAOM','Tracer2':'Root uptake','Birnessite2':'Birnessite'},connectionstyle='arc3, rad=0.2',pos=pos)

save_one_fig(networkfig)

######
# Plot annual decomp for comparison with litter decomposition studies
#####

if len(incubation_list)>0:
    incubation_data=xarray.combine_by_coords(incubation_list).squeeze()
    import pandas
    datadir='/home/b0u/Mn data/'
    Berg_massloss1=pandas.read_csv(datadir+'Berg et al 2015/Fig 4A.csv')
    Berg_massloss2=pandas.read_csv(datadir+'Berg et al 2015/Fig 4B.csv')

    f,axs=subplots(1,2,num='Figure S2 Litter annual mass loss',clear=True,figsize=(8,4.6))
    incubation_length=5
    litter_massloss=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=slice(oneyr+oneyr,None,oneyr*incubation_length))/(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1+oneyr)).to_masked_array().ravel()*100
    lignin_massloss=(1-(incubation_data['Total Sorbed Lignin']).isel(time=slice(oneyr+oneyr,None,oneyr*incubation_length),)/(incubation_data['Total Sorbed Lignin']).isel(time=1+oneyr)).to_masked_array().ravel()*100
    Mn_conc=incubation_data['litter_Mn'].to_masked_array().ravel().compressed()[::incubation_length]*leafC_mass_conv
    x=argsort(Mn_conc)
    # axs[0].plot(Mn_conc[x],litter_massloss[x],'-o',label='Model mass loss')
    axs[0].plot(Mn_conc[x],lignin_massloss[x],'-o',label='Model lignin mass loss')
    # axs[0].plot(Berg_massloss1['x'],Berg_massloss1['Massloss'],'+',label='Berg multiple species')
    axs[0].plot(Berg_massloss2['x']/(1e-3*Mn_molarmass),Berg_massloss2['Massloss'],'+',label='Berg et al. (2013) measurements')
    # axs[0].plot(Davey_data['Mn mg-g-1 DM'],Davey_data['Limit value']*(1-exp(-Davey_data['k value']*oneyr*1*100/Davey_data['Limit value'])),'x',label='Davey Oak')  
    axs[0].set_xlabel('Litter Mn concentration (mmol kg$^{-1}$)')
    axs[0].set_ylabel('Mass loss (%)')
    axs[0].legend()
    axs[0].set_title('One year mass loss')


    Davey_data=pandas.read_excel(datadir+'Davey et al 2007/Davey et al 2007 table data.xlsx',sheet_name=0)

    axs[1].set_title('Mass loss over time')
    t=linspace(0,1000,20)
    cmap=get_cmap('viridis')
    for site in range(len(Davey_data)):
        m=Davey_data['Limit value'][site]
        k=Davey_data['k value'][site]
        davey_line=axs[1].plot(t,m*(1-exp(-k*t/(m/100))),c=cmap(Davey_data['Mn mg-g-1 DM'][site]/(1e-3*Mn_molarmass)/Mn_conc.max()),ls=':')
        
    initial_mass=(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1)   
    # t=arange(0,oneyr*incubation_length-dt,dt)
    for sim in range(len(incubation_data.soil_pH)):
        for rep in range(len(incubation_data['litter_year'])//incubation_length):
            massloss=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(soil_pH=sim)[1+oneyr*rep*incubation_length:oneyr*(rep+1)*incubation_length]/initial_mass.isel(soil_pH=sim))*100
            mod_line=axs[1].plot(massloss['time']-massloss['time'][0],massloss,c=cmap(incubation_data['litter_Mn'].isel(soil_pH=sim,litter_year=rep*incubation_length)*leafC_mass_conv/Mn_conc.max()))

    axs[1].set_xlabel('Time (days)')
    axs[1].set_ylabel('Mass loss (%)')
    cb=f.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap,norm=Normalize(vmin=0,vmax=Mn_conc.max())),ax=axs[1])
    cb.set_label('Litter Mn concentration (mmol kg$^{-1}$)')
    axs[1].legend(handles=davey_line+mod_line,labels=['Davey et al. (2007) measurements','Model'])

    for num in range(axs.size):
        letter_label(axs.ravel()[num],num)

    save_one_fig(f)


def save_all_figs(dirname,format='png',**kwargs):
    for fname in get_figlabels():
        fname_fixed=fname.replace('/','-')
        print(fname_fixed)
        figure(fname_fixed).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed),**kwargs)




