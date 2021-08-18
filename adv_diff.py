# Advection-diffusion solver adapted from ELM SoilLittVertTranspMod.F90
# Based on S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5

import numpy
from scipy.linalg import solve_banded

def adv_diff(conc_trcr,zsoi,dtime,diffus=1e-50,adv_flux=0.0,source=0.0):
    '''Advection-diffusion solver adapted from ELM SoilLittVertTranspMod.F90
    Based on S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5
    '''

    # Set statement functions
    def aaa(pe):
        # A function from Patankar, Table 5.2, pg 95
        return  max (0., (1. - 0.1 * abs(pe))**5)  

    nlayers=len(zsoi)

    # Set the distance between the node and the one ABOVE it   
    dz=numpy.diff(zsoi)
    dz=numpy.append(dz,dz[-1])
    dz_node=numpy.zeros(nlayers)
    dz_node[0] = zsoi[0]
    for j in range(1,nlayers):
        dz_node[j]= zsoi[j] - zsoi[j-1]
    
    zisoi=numpy.zeros(nlayers+1)
    zisoi[1:-1]=(zsoi[:-1]+zsoi[1:])*0.5
    zisoi[-1]=zisoi[-2]+dz_node[-1]

    d_m1_zm1 = numpy.zeros(nlayers)
    d_p1_zp1 = numpy.zeros(nlayers)
    f_m1     = numpy.zeros(nlayers)
    f_p1 =     numpy.zeros(nlayers)
    pe_m1 =    numpy.zeros(nlayers)
    pe_p1 =    numpy.zeros(nlayers)

    a_tri = numpy.zeros(nlayers+2)
    b_tri = numpy.zeros(nlayers+2)
    c_tri = numpy.zeros(nlayers+2)
    r_tri = numpy.zeros(nlayers+2)

    if numpy.ndim(diffus) == 0:
        diffus = numpy.zeros(nlayers)+diffus
        # diffus[0]=1e-50
        # diffus[-1]=1e-50
    if numpy.ndim(adv_flux) == 0:
        adv_flux = numpy.zeros(nlayers)+adv_flux
        adv_flux[-1]=0.0
    if numpy.ndim(source) == 0:
        source = numpy.zeros(nlayers)+source


    # Set Pe (Peclet #) and D/dz throughout column

    for j in range(nlayers):

        # dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
        # dz_node_tracer is difference between cell centers 

        # Calculate the D and F terms in the Patankar algorithm
        if (j == 0) :
            d_m1_zm1[j] = 0.
            w_p1 = (zsoi[j+1] - zisoi[j]) / dz_node[j+1]
            if ( diffus[j+1] > 0. and diffus[j] > 0.) :
                d_p1 = 1. / ((1. - w_p1) / diffus[j] + w_p1 / diffus[j+1]) # Harmonic mean of diffus
            else:
                d_p1 = 0.
            
            d_p1_zp1[j] = d_p1 / dz_node[j+1]
            f_m1[j] = adv_flux[j]  # Include infiltration here
            f_p1[j] = adv_flux[j+1]
            pe_m1[j] = 0.
            pe_p1[j] = f_p1[j] / d_p1_zp1[j] # Peclet #
        elif (j == nlayers-1) :
            # At the bottom, assume no gradient in d_z (i.e., they're the same)
            w_m1 = (zisoi[j-1] - zsoi[j-1]) / dz_node[j]
            if ( diffus[j] > 0. and diffus[j-1] > 0.) :
                d_m1 = 1. / ((1. - w_m1) / diffus[j] + w_m1 / diffus[j-1]) # Harmonic mean of diffus
            else:
                d_m1 = 0.

            d_m1_zm1[j] = d_m1 / dz_node[j]
            d_p1_zp1[j] = d_m1_zm1[j] # Set to be the same
            f_m1[j] = adv_flux[j]
            #f_p1[j] = adv_flux(c,j+1)
            f_p1[j] = 0.
            pe_m1[j] = f_m1[j] / d_m1_zm1[j] # Peclet #
            pe_p1[j] = f_p1[j] / d_p1_zp1[j] # Peclet #
        else:
            # Use distance from j-1 node to interface with j divided by distance between nodes
            w_m1 = (zisoi[j-1] - zsoi[j-1]) / dz_node[j]
            if ( diffus[j-1] > 0. and diffus[j] > 0.) :
                d_m1 = 1. / ((1. - w_m1) / diffus[j] + w_m1 / diffus[j-1]) # Harmonic mean of diffus
            else:
                d_m1 = 0.
            
            w_p1 = (zsoi[j+1] - zisoi[j]) / dz_node[j+1]
            if ( diffus[j+1] > 0. and diffus[j] > 0.) :
                d_p1 = 1. / ((1. - w_p1) / diffus[j] + w_p1 / diffus[j+1]) # Harmonic mean of diffus
            else:
                d_p1 = (1. - w_m1) * diffus[j] + w_p1 * diffus[j+1] # Arithmetic mean of diffus
            
            d_m1_zm1[j] = d_m1 / dz_node[j]
            d_p1_zp1[j] = d_p1 / dz_node[j+1]
            f_m1[j] = adv_flux[j]
            f_p1[j] = adv_flux[j+1]
            pe_m1[j] = f_m1[j] / d_m1_zm1[j] # Peclet #
            pe_p1[j] = f_p1[j] / d_p1_zp1[j] # Peclet #



    # Calculate the tridiagonal coefficients
    # top layer (atmosphere)
    a_tri[0] = 0.
    b_tri[0] = 1.
    c_tri[0] = -1.
    r_tri[0] = 0.


    for j in range(nlayers+1):

        
        if (j == 0) :
            a_p_0 =  dz[j] / dtime
            a_tri[j+1] = -(d_m1_zm1[j] * aaa(pe_m1[j]) + max( f_m1[j], 0.)) # Eqn 5.47 Patankar
            c_tri[j+1] = -(d_p1_zp1[j] * aaa(pe_p1[j]) + max(-f_p1[j], 0.))
            b_tri[j+1] = -a_tri[j+1] - c_tri[j+1] + a_p_0
            r_tri[j+1] = source[j] * dz[j] /dtime + (a_p_0 - adv_flux[j]*0) * conc_trcr[j] 
        elif (j < nlayers) :
            a_p_0 =  dz[j] / dtime
            a_tri[j+1] = -(d_m1_zm1[j] * aaa(pe_m1[j]) + max( f_m1[j], 0.)) # Eqn 5.47 Patankar
            c_tri[j+1] = -(d_p1_zp1[j] * aaa(pe_p1[j]) + max(-f_p1[j], 0.))
            b_tri[j+1] = -a_tri[j+1] - c_tri[j+1] + a_p_0
            r_tri[j+1] = source[j] * dz[j] /dtime + a_p_0 * conc_trcr[j]
        else: #  0 concentration gradient at bottom
            a_tri[j+1] = -1.
            b_tri[j+1] = 1.
            c_tri[j+1] = 0. 
            r_tri[j+1] = 0.

    # Solve for concentration profile
    conc_after = solve_banded((1,1),numpy.row_stack((c_tri,b_tri,a_tri)),r_tri)

    return conc_after

def run_steps(conc,nsteps,*args,**kwargs):
    for step in range(nsteps):
        conc=adv_diff(conc,*args,**kwargs)[1:-1]
    return conc

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    z=numpy.linspace(0.05,10,20)
    # z=numpy.array([7.10063521e-03, 2.79249996e-02, 6.22585751e-02, 1.18865065e-01, # ELM layer thicknesses
    #    2.12193400e-01, 3.66065800e-01, 6.19758487e-01, 1.03802705e+00,
    #    1.72763526e+00, 2.86460710e+00, 4.73915672e+00, 7.82976627e+00,
    #    1.29253206e+01, 2.13264694e+01, 3.51776199e+01])
    nlayers=len(z)
    conc=numpy.linspace(0.1,0,nlayers)**2
    conc[3]=0.2
    
    # conc[-1]=0.0
    # conc[10]=1.0

    dz=numpy.diff(z)
    dz=numpy.append(dz,dz[-1])

    dt=3600
    adv=5e-5
    diffco=5e-5

    f,a=plt.subplots(2,3,num='Result',clear=True)
    for num in range(3):
        a[0,num].plot(z,conc,label='Original',c='k')
        a[1,num].axhline(0.0,ls='--',c='k',lw=0.5)
        a[1,num].axhline((conc*dz).sum(),ls='--',c='k',lw=0.5)
    a[0,0].plot(z,adv_diff(conc,z,dt,adv_flux=adv)[1:-1],label='Adv',c='C1')
    a[0,1].plot(z,adv_diff(conc,z,dt,diffus=diffco)[1:-1],label='Diff',c='C2')
    a[0,2].plot(z,adv_diff(conc,z,dt,diffus=diffco,adv_flux=adv)[1:-1],label='Diff and adv',c='C3')

    a[0,0].plot(z,run_steps(conc,4,z,dt,adv_flux=adv),ls='--',c='C1')
    a[0,1].plot(z,run_steps(conc,4,z,dt,diffus=diffco),ls='--',c='C2')
    a[0,2].plot(z,run_steps(conc,4,z,dt,diffus=diffco,adv_flux=adv),ls='--',c='C3')

    a[0,0].plot(z,run_steps(conc,10,z,dt,adv_flux=adv),ls=':',c='C1')
    a[0,1].plot(z,run_steps(conc,10,z,dt,diffus=diffco),ls=':',c='C2')
    a[0,2].plot(z,run_steps(conc,10,z,dt,diffus=diffco,adv_flux=adv),ls=':',c='C3')



    a[1,0].plot(z,(conc*dz).cumsum(),label='Original',c='C0')
    a[1,0].plot(z,(adv_diff(conc,z,dt,adv_flux=adv)[1:-1]*dz).cumsum(),ls='-',c='C1')
    a[1,1].plot(z,(adv_diff(conc,z,dt,diffus=diffco)[1:-1]*dz).cumsum(),ls='-',c='C2')
    a[1,2].plot(z,(adv_diff(conc,z,dt,diffus=diffco,adv_flux=adv)[1:-1]*dz).cumsum(),ls='-',c='C3')

    a[1,0].plot(z,(run_steps(conc,4,z,dt,adv_flux=adv)*dz).cumsum(),label='Adv = 5e-4',c='C1',ls='--')
    a[1,1].plot(z,(run_steps(conc,4,z,dt,diffus=diffco)*dz).cumsum(),label='Diff = 1e-5',c='C2',ls='--')
    a[1,2].plot(z,(run_steps(conc,4,z,dt,diffus=diffco,adv_flux=adv)*dz).cumsum(),label='Diff and adv',c='C3',ls='--')

    a[1,0].plot(z,(run_steps(conc,10,z,dt,adv_flux=adv)*dz).cumsum(),label='Adv = 5e-4',c='C1',ls=':')
    a[1,1].plot(z,(run_steps(conc,10,z,dt,diffus=diffco)*dz).cumsum(),label='Diff = 1e-5',c='C2',ls=':')
    a[1,2].plot(z,(run_steps(conc,10,z,dt,diffus=diffco,adv_flux=adv)*dz).cumsum(),label='Diff and adv',c='C3',ls=':')


    a[0,0].set_title('Advection')
    a[0,1].set_title('Diffusion')
    a[0,2].set_title('Advection+Diffusion')
    a[0,0].legend()
    a[0,0].set_ylabel('Concentration')
    a[1,0].set_ylabel('Cumulative concentration')
    a[1,0].set_xlabel('Depth')
    a[1,1].set_xlabel('Depth')
    a[1,2].set_xlabel('Depth')

    plt.show()