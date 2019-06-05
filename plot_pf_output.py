from pylab import *

import pandas

def read_tecfile(filename):
    f=open(filename,'r')
    headerline=f.readline()
    header=[x[:x.find('[')].strip(' "') for x in headerline.split(',')]
    lines=f.readlines()

    f.close()

    # This converter stuff is necessary because pflotran doesn't produce the right text format when exponent is <= E-100
    def convertfunc(strval):
        if '-' in strval and 'E' not in strval:
            return 0.0
        else:
            return float(strval)

    converters=dict(((n,convertfunc) for n in range(len(header))))
    return pandas.read_table(filename,skiprows=1,names=header,header=None,delim_whitespace=True,index_col=0,converters=converters)

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
        data=read_tecfile(filename)
        data.plot()
    elif filename.endswith('.h5'):
        data=read_hdf5file(filename).plot()
        for var in data:
            data[var].plot()

    show()
