from pylab import *

import pandas

def read_tecfile(filename,convert_units_to=None,saturation=None):
    f=open(filename,'r')
    headerline=f.readline()
    f.close()
    colnames=headerline.split(',')
    header=[]
    units={}
    for col in colnames:
        if '[' in col:
            name=col[:col.find('[')].strip(' "')
            unit=col[col.find('[')+1:col.find(']')]
        else:
            name=col.strip('"').split()[0]
            unit='NA'
        header.append(name)
        units[name]=unit


    # This converter stuff is necessary because pflotran doesn't produce the right text format when exponent is <= E-100
    def convertfunc(strval):
        if '-' in strval and 'E' not in strval:
            val,exponent=strval.split('-')
            return float(val)*10**-float(exponent)  
        else:
            return float(strval)

    converters=dict(((n,convertfunc) for n in range(len(header))))
    
    data=pandas.read_table(filename,skiprows=1,names=header,header=None,delim_whitespace=True,index_col=0,converters=converters)
    
    if convert_units_to is not None:
        data,units=convert_units(data,units,convert_units_to)
    else:
        print('Note: not converting units')
        
    return data,units
    
def convert_units(data,units,convert_units_to):
    out=data.copy()
    units_out=units.copy()
    water_molar_density=55.345 # mol/L
    saturation = 1.0 # Assuming we are always doing chemistry under saturated conditions (even if it is tricking Pflotran)
    if 'Porosity' in data.columns:
        porosity=data['Porosity']
    else:
        raise ValueError('Porosity must be included in output file to calculate unit changes')
    
    if convert_units_to=='M':
        print('Converting all concentration units to M equivalent assuming saturation')
        for col in data.columns:
            if units[col] == 'mol/m^3':
                print(col)
                # Convert from mol[x]/m^3[bulk] to equivalent mol[x]/L[H2O]:  /(1000L/m^3*saturation*porosity) (m^3[H2O]/m^3[bulk])
                out[col]=data[col]/(1000*saturation/porosity)
                units_out[col]='M'
                
    elif convert_units_to=='mol/m^3':
        print('Converting all concentration units to mol/m^3')
        for col in data.columns:
            if units[col] == 'M':
                print(col)
                # Convert from mol[x]/L[H2O] to mol[x]/m^3[bulk]: Divide by porosity*saturation*1000
                out[col]=data[col]/(1000*saturation*porosity)
                units_out[col]='mol/m^3'

    else:
        raise ValueError('Units %s not supported. Must be "M" or "mol/m^3"'%convert_units_to)
        
    return out,units_out
    


def read_hdf5file(filename):
    import h5py
    data_file=h5py.File(filename)
    element_names=list(data_file)
    timepoints=element_names[2:]
    output_fields=list(data_file[timepoints[0]])
    import xarray,numpy
    time_axis=numpy.array([float(t[t.find(':')+1:-1].strip()) for t in timepoints])
    coords=list(data_file['Coordinates'])
    datashape=data_file[timepoints[0]][output_fields[0]].value.shape
    data_coords={'time':time_axis}
    for num,coord in enumerate(coords):
        data_coords[coord]=numpy.arange(datashape[num])
    data_out=xarray.Dataset(coords=data_coords)
    data_dict={}
    for field in output_fields:
        data_dict[field]=(numpy.zeros((len(time_axis),datashape[0],datashape[1],datashape[2])))
    for tnum,timepoint in enumerate(timepoints):
        print('Timepoint: ',timepoint)
        vals=data_file[timepoint]
        for field in output_fields:
            data_dict[field][tnum,:,:,:]=vals[field].value
    
    for field in output_fields:
        data_out[field]=(('time',coords[0],coords[1],coords[2]),data_dict[field])
    return data_out.sortby('time')

if __name__=='__main__':
    # data=read_tecfile('CLM-CN-obs-0.tec')*100**3
    #
    # figure(1);clf()
    # subplot(211)
    # plot(data['C'],'k--',label='CO2')
    # plot(data['SOM1'])
    # plot(data['SOM2'])
    # plot(data['SOM3'])
    # plot(data['SOM4'])
    # plot(data['Lit1C'])
    # plot(data['Lit2C'])
    # plot(data['Lit3C'])
    # title('C pools')
    # legend()
    # xlabel('Time (days)')
    # ylabel('Concentration (mol cm$^{-3}$)')
    #
    #
    # subplot(212)
    # plot(data['N'],'k--',label='Mineral N')
    # plot(data['Lit1N'])
    # plot(data['Lit2N'])
    # plot(data['Lit3N'])
    # legend()
    # title('N pools')
    # xlabel('Time (days)')
    # ylabel('Concentration (mol cm$^{-3}$)')
    #
    #
    # tight_layout()
    
    import sys
    if len(sys.argv)==1:
        filename='test3d-obs-0.tec'
    else:
        filename=sys.argv[1]

    if filename.endswith('.tec'):
        data,units=read_tecfile(filename,convert_units_to='M')
        data.loc[:,pandas.Series(units)=='M'].plot()
    elif filename.endswith('mas.dat'):
        data,units=read_tecfile(filename)
        data.loc[:,pandas.Series(units)=='mol'].plot()
    elif filename.endswith('.h5'):
        data=read_hdf5file(filename)
        for var in data.variables:
            data[var].plot()

    show()
