from pylab import *
import xarray
import sys
import netCDF4

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

filename=sys.argv[1]

with netCDF4.Dataset(filename) as d:
    groups=list(d.groups.keys())



data_list=[]
incubation_list=[]
for g in groups:
    ph,Ndep,warming,redox=g.split('_')[:4]

    d=xarray.open_dataset(filename,group=g,chunks={'time':25},decode_times=False)
    d['soil_pH']=float(ph[2:])
    d['Ndep']=int(Ndep[4:])
    d['warming']=float(warming[7:])
    d['redox_cycles']=int(redox[5:])
    newdims=['soil_pH','Ndep','warming','redox_cycles']
    if 'incubation' in g:
        incubation_list.append(d.expand_dims(newdims).set_coords(newdims))
    else:
        data_list.append(d.expand_dims(newdims).set_coords(newdims))

alldata=xarray.combine_by_coords(data_list)
allincubations=xarray.combine_by_coords(incubation_list)
print('Done combining datasets')

dt=(alldata['time'][1]-alldata['time'][0]).item()
oneyr=int(365/dt)

def get_data(filename,ph=4.5,Ndep=0,warming=0,redox=0,incubation=False,variable=None,groupname=None,printgroup=False):
    if groupname is None:
        groupname='pH%1.1f_Ndep%03d_warming%d_redox%d'%(ph,Ndep,warming,redox)
    if incubation:
        groupname=groupname+'_incubation'
    if printgroup:
        print(groupname)
    if variable is not None:
        with xarray.open_dataset(filename,group=groupname,decode_times=False) as data:
            return data[variable]
    else:
        return xarray.open_dataset(filename,group=groupname,decode_times=False)

results_byph=alldata.isel(Ndep=0,warming=0,redox_cycles=0)
results_byredox=alldata.isel(Ndep=0,warming=0,soil_pH=0)
pH_sims=alldata['soil_pH'].to_masked_array()
f,axs=subplots(3,1,clear=True,num='Mn results by soil pH',squeeze=False,figsize=(4.5,8.75))
# These are annual means
total_cellulose=(results_byph['Total Sorbed Cellulose'].coarsen(time=oneyr,boundary='trim').mean()*12e-3*(results_byph.z_bottom-results_byph.z_top)/100).sum(dim='depth') 
total_lignin=(results_byph['Total Sorbed Lignin'].coarsen(time=oneyr,boundary='trim').mean()*12e-3*(results_byph.z_bottom-results_byph.z_top)/100).sum(dim='depth') 

ax=axs[1,0]
# ax.plot(pH_sims,total_cellulose.isel(time=39),'-o',c='C0',label='Cellulose')
ax.plot(pH_sims,total_lignin.isel(time=39),'-o',c='k',label='Lignin')
# axs[0,0].plot(pH_sims,(total_cellulose+total_lignin).isel(time=39),'-o',label='POM',c='C2')

# ax.plot(pH_sims,total_cellulose.isel(time=19),'--o',c='C0')
ax.plot(pH_sims,total_lignin.isel(time=19),'--o',c='k')
# axs[0,0].plot(pH_sims,(total_cellulose+total_lignin).isel(time=19),'--o',c='C2')

# ax.plot(pH_sims,total_cellulose.isel(time=9),':o',c='C0')
ax.plot(pH_sims,total_lignin.isel(time=9),':o',c='k')
# axs[0,0].plot(pH_sims,(total_cellulose+total_lignin).isel(time=9),':o',c='C2')

# ax.legend()
ax.set_xlabel('Soil pH')
ax.set_ylabel('OM stock (kg C m$^{-2}$)')
ax.set_title('Litter C stock')
ax.set_ylim(bottom=0)

ax=axs[2,0]
total_MAOM=((results_byph['Total DOM3'].coarsen(time=oneyr,boundary='trim').mean()*results_byph.saturation*results_byph.Porosity.mean())*12e-3/1000*100**3*(results_byph.z_bottom-results_byph.z_top)/100).sum(dim='depth')
ax.plot(pH_sims,total_MAOM.isel(time=39),'-o',c='k',label='MAOM')
ax.plot(pH_sims,total_MAOM.isel(time=19),'--o',c='k',label='MAOM')
ax.plot(pH_sims,total_MAOM.isel(time=9),':o',c='k',label='MAOM')
ax.set_xlabel('Soil pH')
ax.set_ylabel('OM stock (kg C m$^{-2}$)')
ax.set_title('Mineral-associated organic matter')
# ax.set_ylim(bottom=0)

ax=axs[0,0]
# Converting these from per kgC to kg dry mass by multiplying by 0.4
ax.plot(pH_sims,results_byph['litter_Mn'].isel(litter_year=39)*leafC_mass_conv,'o-',c='k',label='Year 40')
ax.plot(pH_sims,results_byph['litter_Mn'].isel(litter_year=19)*leafC_mass_conv,'o--',c='k',label='Year 20')
ax.plot(pH_sims,results_byph['litter_Mn'].isel(litter_year=9)*leafC_mass_conv,'o:',c='k',label='Year 10')
ax.set_xlabel('Soil pH')
ax.set_ylabel('Litter Mn conc. (mmol kg$^{-1}$)')
ax.set_title('Leaf litter Mn concentration')
ax.legend()
# 
# ax=axs[1,0]

birnessite=results_byph['Birnessite2 VF'].coarsen(time=oneyr,boundary='trim').mean()*7/molar_volume_birnessite*1e6*Mn_molarmass/results_byph.BD
totalMn=((results_byph['Total Mn++']+results_byph['Total Mn+++'])/1000*results_byph.Porosity.mean()*results_byph.saturation + \
                    (results_byph['Birnessite2 VF']*7/molar_volume_birnessite) + \
                    results_byph['Total Sorbed Mn++']/100**3).coarsen(time=oneyr,boundary='trim').mean()*Mn_molarmass/results_byph.BD*1e6

# ax.plot(pH_sims,totalMn.isel(depth=0,time=39),'-o',c='k',label='Total Mn')
# ax.plot(pH_sims,totalMn.isel(depth=0,time=19),'--o',c='k')
# ax.plot(pH_sims,totalMn.isel(depth=0,time=9),':o',c='k')
# 
# ax.plot(pH_sims,birnessite.isel(depth=0,time=39),'-^',c='C1',label='Birnessite')
# ax.plot(pH_sims,birnessite.isel(depth=0,time=19),'--^',c='C1')
# ax.plot(pH_sims,birnessite.isel(depth=0,time=9),':^',c='C1')
# 
# # ax.plot(pH_sims,(totalMn-birnessite).isel(depth=0,time=39),'-o',c='C2',label='Exchangeable Mn')
# # ax.plot(pH_sims,(totalMn-birnessite).isel(depth=0,time=19),'--o',c='C2')
# # ax.plot(pH_sims,(totalMn-birnessite).isel(depth=0,time=9),':o',c='C2')
# 
# ax.set_xlabel('Soil pH')
# ax.set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
# ax.set_title('Organic layer manganese')
# ax.legend()
# ax=axs[3,0]
# total_CO2=((result_data['Total Tracer'].coarsen(time=oneyr,boundary='trim').mean(dim='time')*result_data.saturation*result_data.Porosity.mean())*12e-3/1000*100**3*(result_data.z_bottom-result_data.z_top)/100).sum(dim='depth')
# ax.plot(pH_sims,total_CO2.isel(time=39),'-o',c='C0')
# ax.plot(pH_sims,total_CO2.isel(time=19),'--o',c='C0')
# ax.plot(pH_sims,total_CO2.isel(time=9),':o',c='C0')
# ax.set_ylabel('Cumulative CO$_2$ (kg C m$^{-2}$)')
# ax.set_title('Cumulative CO$_2$ production')
# ax.set_xlabel('Soil pH')

from string import ascii_lowercase
for num in range(len(axs)):
    axs[num,0].text(0.01,1.08,'('+ascii_lowercase[num]+')',transform=axs[num,0].transAxes)

result_data=alldata.isel(Ndep=0,warming=0,redox_cycles=0)
f,axs=subplots(1,len(result_data['soil_pH']),clear=True,num='Mn profiles',squeeze=False,figsize=(14.,3.))
cm=get_cmap('cool')
norm=matplotlib.colors.Normalize(pH_sims.min(),pH_sims.max())


for n in range(len(result_data['soil_pH'])):
    axs[0,n].plot(totalMn.isel(time=9,soil_pH=n),result_data['z_middle'].isel(soil_pH=n),':',c='C0',label='Total Mn: Year 10')
    axs[0,n].plot(totalMn.isel(time=19,soil_pH=n),result_data['z_middle'].isel(soil_pH=n),'--',c='C0',label='Year 20')
    axs[0,n].plot(totalMn.isel(time=39,soil_pH=n),result_data['z_middle'].isel(soil_pH=n),'-',c='C0',label='Year 40')
    # l=axs[0,n].plot(birnessite.isel(time=9,soil_pH=sim),result_data['z_middle'].isel(soil_pH=sim),':',c='C1')
    # axs[0,n].plot(birnessite.isel(time=19,soil_pH=sim),result_data['z_middle'].isel(soil_pH=sim),'--',c='C1')
    # axs[0,n].plot(birnessite.isel(time=30,soil_pH=sim),result_data['z_middle'].isel(soil_pH=sim),'-',label='Birnessite',c='C1')
    axs[0,n].set_ylim(result_data.depth.max(),result_data.depth.min())
    axs[0,n].set_title('pH = %1.1f'%pH_sims[n])
    axs[0,n].set_ylabel('Depth (cm)')
    axs[0,n].set_xlabel('Mn concentration ($\mu$g g$^{-1}$)')
    axs[0,n].set_ylim(top=0)

axs[0,0].legend()

f,axs=subplots(len(result_data['depth']),1,clear=True,num='Mn conc by pH',squeeze=False,figsize=(4.5,8.75))
cm=get_cmap('inferno')
for z in range(len(result_data['depth'])):
    axs[z,0].plot(pH_sims,birnessite.isel(depth=z,time=39),'o-',label='Year 40',c='C0')
    axs[z,0].plot(pH_sims,birnessite.isel(depth=z,time=19),'o--',label='Year 20',c='C0')
    axs[z,0].plot(pH_sims,birnessite.isel(depth=z,time=9),'o:',label='Year 10',c='C0')
    axs[z,0].plot(pH_sims,totalMn.isel(depth=z,time=39),'o-',c='C1')
    axs[z,0].plot(pH_sims,totalMn.isel(depth=z,time=19),'o--',c='C1')
    axs[z,0].plot(pH_sims,totalMn.isel(depth=z,time=9),'o:',c='C1')
    
    
    axs[z,0].set_ylabel('Mn concentration\n($\mu$g g$^{-1}$)')
axs[z,0].set_xlabel('Soil pH')
axs[0,0].legend()
axs[1,0].legend(handles=(axs[1,0].lines[0],axs[1,0].lines[3]),labels=('Birnessite','Total Mn'))
axs[0,0].set_title('%d-%d cm (Organic layer)'%(result_data['z_top'].isel(soil_pH=0,depth=0),result_data['z_bottom'].isel(soil_pH=0,depth=0)))
for z in range(1,len(result_data['depth'])):
    axs[z,0].set_title('%d-%d cm'%(result_data['z_top'].isel(soil_pH=0,depth=z),result_data['z_bottom'].isel(soil_pH=0,depth=z)))

######
# Plot annual decomp for comparison with litter decomposition studies
#####


import pandas
datadir='/home/b0u/Mn data/'
Berg_massloss1=pandas.read_csv(datadir+'Berg et al 2015/Fig 4A.csv')
Berg_massloss2=pandas.read_csv(datadir+'Berg et al 2015/Fig 4B.csv')


incubation_data=allincubations.isel(Ndep=0,warming=0,redox_cycles=0)
f,axs=subplots(1,2,num='Litter annual mass loss',clear=True,figsize=(8,4.6))
incubation_length=5
litter_massloss=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=slice(oneyr+oneyr,None,oneyr*incubation_length),depth=0)/(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1+oneyr,depth=0)).to_masked_array().ravel()*100
lignin_massloss=(1-(incubation_data['Total Sorbed Lignin']).isel(time=slice(oneyr+oneyr,None,oneyr*incubation_length),depth=0)/(incubation_data['Total Sorbed Lignin']).isel(time=1+oneyr,depth=0)).to_masked_array().ravel()*100
Mn_conc=result_data['litter_Mn'].to_masked_array().ravel()[::incubation_length]*leafC_mass_conv
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


# axs[1].set_title('%d year mass loss'%incubation_length)
# Berg_limitvalue=pandas.read_csv(datadir+'Berg et al 2015/Fig 7.csv')
# axs[1].plot(Berg_limitvalue['x'],Berg_limitvalue['Limit Value'],'+',label='Berg limit value')
# litter_limitvalue=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=slice(oneyr*incubation_length-1,None,oneyr*incubation_length),depth=0)/(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1,depth=0)).to_masked_array().ravel()*100
# lignin_limitvalue=(1-(incubation_data['Total Sorbed Lignin']).isel(time=slice(oneyr*incubation_length-1,None,oneyr*incubation_length),depth=0)/(incubation_data['Total Sorbed Lignin']).isel(time=1,depth=0)).to_masked_array().ravel()*100
# 
# axs[1].plot(Davey_data['Mn mg-g-1 DM'],Davey_data['Limit value']*(1-exp(-Davey_data['k value']*oneyr*5*100/Davey_data['Limit value'])),'x',label='Davey Oak')  
# 
# axs[1].plot(Mn_conc[x],litter_limitvalue[x],'-o',label='Model mass loss')
# axs[1].plot(Mn_conc[x],lignin_limitvalue[x],'-o',label='Model lignin mass loss')
# axs[1].set_xlabel('Litter Mn concentration (mg g$^{-1}$)')
# axs[1].set_ylabel('Mass loss (%)')
# axs[1].legend()

Davey_data=pandas.read_excel(datadir+'Davey et al 2007/Davey et al 2007 table data.xlsx',sheet_name=0)

axs[1].set_title('Mass loss over time')
t=linspace(0,1000,20)
cmap=get_cmap('viridis')
for site in range(len(Davey_data)):
    m=Davey_data['Limit value'][site]
    k=Davey_data['k value'][site]
    davey_line=axs[1].plot(t,m*(1-exp(-k*t/(m/100))),c=cmap(Davey_data['Mn mg-g-1 DM'][site]/(1e-3*Mn_molarmass)/Mn_conc.max()),ls=':')
    
initial_mass=(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(time=1,depth=0)   
t=arange(0,oneyr*incubation_length-dt,dt)
for sim in range(len(incubation_data.soil_pH)):
    for rep in range(len(result_data['litter_year'])//incubation_length):
        massloss=(1-(incubation_data['Total Sorbed Lignin']+incubation_data['Total Sorbed Cellulose']).isel(depth=0,soil_pH=sim)[1+oneyr*rep*incubation_length:oneyr*(rep+1)*incubation_length]/initial_mass.isel(soil_pH=sim))*100
        mod_line=axs[1].plot(t,massloss,c=cmap(result_data['litter_Mn'].isel(soil_pH=sim,litter_year=rep*incubation_length)*leafC_mass_conv/Mn_conc.max()))

axs[1].set_xlabel('Time (days)')
axs[1].set_ylabel('Mass loss (%)')
cb=f.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap,norm=Normalize(vmin=0,vmax=Mn_conc.max())),ax=axs[1])
cb.set_label('Litter Mn concentration (mmol kg$^{-1}$)')
axs[1].legend(handles=davey_line+mod_line,labels=['Davey et al. (2007) measurements','Model'])

axs[0].text(0.01,1.02,'(a)',transform=axs[0].transAxes)
axs[1].text(0.01,1.02,'(b)',transform=axs[1].transAxes)

f,axs=subplots(2,1,num='Incubations',clear=True)
for sim in range(len(incubation_data.soil_pH)):
    t=incubation_data['time']/365
    axs[0].plot(t,incubation_data['Total Mn++'].isel(depth=0,soil_pH=sim),c=get_cmap('coolwarm')(sim/len(incubation_data.soil_pH)))
    axs[1].plot(t,incubation_data['Total Sorbed Lignin'].isel(depth=0,soil_pH=sim),c=get_cmap('coolwarm')(sim/len(incubation_data.soil_pH)))
    
axs[0].set(title='Mn$^{2+}$ Concentration',xlabel='Time (years)',ylabel='Mn$^{2+}$ concentration (M)')
axs[1].set(title='Lignin Concentration',xlabel='Time (years)',ylabel='Lignin concentration (mol/m3)')

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

result_data=alldata.sel(soil_pH=4.5,Ndep=0,warming=0,redox_cycles=0)
f,axs=subplots(ncols=5,nrows=len(result_data.depth),sharex=False,clear=True,num='Simulation results',figsize=(12,6.5))
plot_output(result_data,axs,do_legend=True)
plot_output(alldata.sel(soil_pH=6.0,Ndep=0,warming=0,redox_cycles=0),axs,do_legend=False,ls='--')
plot_output(alldata.sel(soil_pH=4.5,Ndep=0,warming=0,redox_cycles=4),axs,do_legend=False,ls=':')


data_allNdeps=alldata.isel(redox_cycles=0,warming=0)

total_cellulose=(data_allNdeps['Total Sorbed Cellulose'].isel(time=slice(oneyr*35,None,None)).mean(dim='time') *12e-3*(data_allNdeps.z_bottom-data_allNdeps.z_top)/100).sum(dim='depth') 
total_lignin=(data_allNdeps['Total Sorbed Lignin'].isel(time=slice(oneyr*35,None,None)).mean(dim='time') *12e-3*(data_allNdeps.z_bottom-data_allNdeps.z_top)/100).sum(dim='depth') 

f,ax=subplots(ncols=2,num='pH by Ndep',clear=True,figsize=(12,5))
phs=data_allNdeps['soil_pH'].to_masked_array()
Ndeps=data_allNdeps['Ndep'].to_masked_array()
h=ax[0].pcolormesh(linspace(phs[0]-.25,phs[-1]+0.25,len(phs)+1),concatenate([Ndeps,[175]]),(total_cellulose+total_lignin).T,edgecolors='w')
cb=f.colorbar(h,ax=ax[0])
cb.set_label('Total litter C (kg C m$^{-2}$)')
ax[0].set_xlabel('pH')
ax[0].set_ylabel('N deposition rate (kg N ha$^{-1}$ year$^{-1}$)')
ax[0].set_title('Interaction of pH and N deposition')
ax[0].text(0.01,1.03,'(a)',transform=ax[0].transAxes)

markers=['o','s','^','h']
for Ndep in range(len(Ndeps)):
    x=data_allNdeps.litter_Mn.isel(Ndep=Ndep).mean(dim='litter_year')
    ax[1].plot(x,(total_cellulose+total_lignin).isel(Ndep=Ndep),label='Ndep = %d kg N ha$^{-1}$ year$^{-1}$'%Ndeps[Ndep],marker=markers[Ndep],lw=1.5,ms=5.0)

# ax[1].set_xlabel('pH')
ax[1].set_xlabel('Leaf Mn concentration (mmol kg$^{-1}$)')
ax[1].set_ylabel('Total litter C (kg C m$^{-2}$)')
ax[1].legend()
ax[1].set_ylim(bottom=0)
ax[1].set_title('Interaction of leaf Mn concentration and N deposition')
ax[1].text(0.01,1.03,'(b)',transform=ax[1].transAxes)

f,ax=subplots(num='pH by Ndep 2',clear=True)
for Ndep in range(len(Ndeps)):
    for ph in range(len(phs)):
        bar(Ndep+(1-ph/4)*0.8,(0.163*dt)/(total_lignin).isel(Ndep=Ndep,soil_pH=ph),width=0.8/4,color=((ph+1)/len(phs),(ph+1)/len(phs),(ph+1)/len(phs)),edgecolor='k')
xticks([0.4,1.4,2.4],['N0','N50','N150'])
legend(['%1.0f mmol kg$^{-1}$'%data_allNdeps['litter_Mn'].isel(Ndep=0,soil_pH=n,litter_year=39) for n in range(len(phs))],title='Litter Mn')
ax.set_ylabel('Lignin turnover time (year$^{-1}$)')


data_warmings=alldata.isel(redox_cycles=0,Ndep=0)
warmings=data_warmings['warming'].to_masked_array()
total_cellulose=(data_warmings['Total Sorbed Cellulose'].isel(time=slice(oneyr*35,None,None)).mean(dim='time') *12e-3*(data_warmings.z_bottom-data_warmings.z_top)/100).sum(dim='depth') 
total_lignin=(data_warmings['Total Sorbed Lignin'].isel(time=slice(oneyr*35,None,None)).mean(dim='time') *12e-3*(data_warmings.z_bottom-data_warmings.z_top)/100).sum(dim='depth') 

total_litter=((data_warmings['Total Sorbed Cellulose']+data_warmings['Total Sorbed Lignin']).coarsen(time=oneyr,boundary='trim').mean()*12e-3*(data_warmings.z_bottom-data_warmings.z_top)/100).sum(dim='depth')

f,ax=subplots(ncols=2,nrows=2,num='pH by Warming',clear=True,figsize=(11,6.8))
h=ax[0,1].pcolormesh(linspace(phs[0]-.25,phs[-1]+0.25,len(phs)+1),concatenate([warmings,[warmings[-1]+1]]),total_litter.isel(time=39).T,edgecolors='w')
cb=f.colorbar(h,ax=ax[0,1])
cb.set_label('Total litter C (kg C m$^{-2}$)')
ax[0,1].set_xlabel('pH')
ax[0,1].set_ylabel('Warming ($^\circ$C)')
ax[0,1].set_title('Interaction of pH and warming (40 years)')
ax[0,1].text(0.01,1.03,'(b)',transform=ax[0,1].transAxes)

h=ax[0,0].pcolormesh(linspace(phs[0]-.25,phs[-1]+0.25,len(phs)+1),concatenate([warmings,[warmings[-1]+1]]),total_litter.isel(time=15).T,edgecolors='w')
cb=f.colorbar(h,ax=ax[0,0])
cb.set_label('Total litter C (kg C m$^{-2}$)')
ax[0,0].set_xlabel('pH')
ax[0,0].set_ylabel('Warming ($^\circ$C)')
ax[0,0].set_title('Interaction of pH and warming (15 years)')
ax[0,0].text(0.01,1.03,'(a)',transform=ax[0,0].transAxes)


markers=['o','s','^','h']
for warming in range(len(warmings)):
    # x=data_warmings.litter_Mn.isel(warming_sim=warming).mean(dim='litter_year')
    # ax[1].plot(x,(total_cellulose+total_lignin).isel(warming_sim=warming),label='Warming = %d $^\circ$C'%warmings[warming],marker=markers[warming],lw=1.5,ms=5.0,color=get_cmap('inferno')((warming+0.5)/len(warmings)))
    ax[1,0].plot(data_warmings.litter_year,data_warmings.litter_Mn.isel(warming=warming,soil_pH=0),color=get_cmap('inferno')((warming+0.5)/len(warmings)),label='Warming: %d$^\circ$C'%warmings[warming])
    ax[1,0].plot(data_warmings.litter_year,data_warmings.litter_Mn.isel(warming=warming,soil_pH=1),color=get_cmap('inferno')((warming+0.5)/len(warmings)),ls='--')
    # ax[1,0].plot(data_warmings.litter_year,data_warmings.litter_Mn.isel(warming=warming,soil_pH=4),color=get_cmap('inferno')((warming+0.5)/len(warmings)),ls=':')
    
    ax[1,1].plot(data_warmings.litter_year,total_litter.isel(warming=warming,soil_pH=0),color=get_cmap('inferno')((warming+0.5)/len(warmings)),label='pH = 4.0')
    ax[1,1].plot(data_warmings.litter_year,total_litter.isel(warming=warming,soil_pH=1),color=get_cmap('inferno')((warming+0.5)/len(warmings)),ls='--',label='pH = 5.0')
    # ax[1,1].plot(data_warmings.litter_year,total_litter.isel(warming=warming,soil_pH=4),color=get_cmap('inferno')((warming+0.5)/len(warmings)),ls=':',label='pH = 6.0')

ax[1,0].set_ylabel('Leaf Mn concentration (mmol kg$^{-1}$)')
ax[1,0].set_xlabel('Year')
ax[1,0].legend(handles=ax[1,0].lines[::3]+ax[1,1].lines[:3],ncol=2)
ax[1,0].set_ylim(bottom=0)
ax[1,0].set_title('Leaf Mn concentration over time')
ax[1,0].text(0.01,1.03,'(c)',transform=ax[1,0].transAxes)

ax[1,1].set_ylabel('Total litter C (kg C m$^{-2}$)')
ax[1,1].set_xlabel('Year')
# ax[1,1].legend(handles=ax[1,1].lines[:3])
ax[1,1].set_ylim(bottom=0)
ax[1,1].set_title('Litter C over time')
ax[1,1].text(0.01,1.03,'(d)',transform=ax[1,1].transAxes)


# # ax[1].set_xlabel('pH')
# ax[1].set_xlabel('Leaf Mn concentration (mmol kg$^{-1}$)')
# ax[1].set_ylabel('Total litter C (kg C m$^{-2}$)')
# ax[1].legend()
# ax[1].set_ylim(bottom=0)
# ax[1].set_title('Interaction of leaf Mn concentration and warming')
# ax[1].text(0.01,1.03,'(b)',transform=ax[1].transAxes)
# 
# for phnum in range(len(total_lignin.soil_pH)):
#     c=data_warmings.litter_Mn.isel(soil_pH=ph).mean(dim='litter_year')
#     ax[2].plot(warmings,(total_cellulose+total_lignin).isel(soil_pH=phnum),'-o',label='pH = %1.1f'%phs[phnum])
# ax[2].set_xlabel('Warming (C)')
# ax[2].set_ylabel('Total litter C (kg C m$^{-2}$)')
# ax[2].legend()
# ax[2].set_ylim(bottom=0)
# ax[2].set_title('Interaction of leaf Mn concentration and warming')
# ax[2].text(0.01,1.03,'(c)',transform=ax[2].transAxes)

f,axs=subplots(ncols=5,nrows=len(result_data.depth),sharex=False,clear=True,num='Warming results low ph',figsize=(12,6.5))
plot_output(data_warmings.isel(soil_pH=0,warming=0),axs,do_legend=True)
plot_output(data_warmings.isel(soil_pH=0,warming=1),axs,do_legend=False,ls='--')
plot_output(data_warmings.isel(soil_pH=0,warming=2),axs,do_legend=False,ls=':')

f,axs=subplots(ncols=5,nrows=len(result_data.depth),sharex=False,clear=True,num='Warming results high ph',figsize=(12,6.5))
plot_output(data_warmings.isel(soil_pH=1,warming=0),axs,do_legend=True)
plot_output(data_warmings.isel(soil_pH=1,warming=1),axs,do_legend=False,ls='--')
plot_output(data_warmings.isel(soil_pH=1,warming=2),axs,do_legend=False,ls=':')

f,axs=subplots(ncols=5,nrows=len(result_data.depth),sharex=False,clear=True,num='Warming results redox cycles',figsize=(12,6.5))
plot_output(alldata.isel(soil_pH=0,warming=0,Ndep=0,redox_cycles=2),axs,do_legend=True)
plot_output(alldata.isel(soil_pH=0,warming=1,Ndep=0,redox_cycles=2),axs,do_legend=False,ls='--')
plot_output(alldata.isel(soil_pH=0,warming=2,Ndep=0,redox_cycles=2),axs,do_legend=False,ls=':')


f,axs=subplots(ncols=5,nrows=len(result_data.depth),sharex=False,clear=True,num='Redox cycles results',figsize=(12,6.5))
plot_output(alldata.isel(soil_pH=0,warming=0,Ndep=0,redox_cycles=0),axs,do_legend=True)
plot_output(alldata.isel(soil_pH=0,warming=0,Ndep=0,redox_cycles=2),axs,do_legend=False,ls='--')
plot_output(alldata.isel(soil_pH=0,warming=0,Ndep=0,redox_cycles=3),axs,do_legend=False,ls=':')



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

