# Advection-diffusion solver adapted from ELM SoilLittVertTranspMod.F90
# Based on S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5

import numpy
from scipy.linalg import solve_banded

def adv_diff(conc_trcr,zsoi,dz,dtime,diffus=1e-50,adv_flux=0.0,source=0.0,printmatrix=False,return_matrix=False,solver='scipy'):
    '''Advection-diffusion solver adapted from ELM SoilLittVertTranspMod.F90
    Based on S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5
    '''

    # Set statement functions
    def aaa(pe):
        # A function from Patankar, Table 5.2, pg 95
        # This essentially determines the relative influence of the two adjacent points on the central point which depends on flow and diffusion rates and direction
        # When Pe is large, flow dominates in a particular direction. When Pe is small, diffusion dominates and the influence is more symmetrical
        return  max (0., (1. - 0.1 * abs(pe))**5)  

    nlayers=len(zsoi)


    dz_node=numpy.zeros(nlayers)
    dz_node[0] = zsoi[0]
    for j in range(1,nlayers):
        dz_node[j]= zsoi[j] - zsoi[j-1]
    
    zisoi=numpy.zeros(nlayers+1) # Starts at zero. zsoi starts at 1
    zisoi[1:-1]=(zsoi[:-1]+zsoi[1:])*0.5
    zisoi[-1]=zisoi[-2]+dz_node[-1]

    # d: diffusivity
    # m: layer above
    # p: layer below
    # pe: Peclet number (ratio of convection to diffusion)
    d_m1_zm1 = numpy.zeros(nlayers) # Diffusivity in layer above divided by layer above thickness
    d_p1_zp1 = numpy.zeros(nlayers) # Diffusivity in layer below divided by layer below thickness
    f_m1     = numpy.zeros(nlayers) # Flow out of above layer
    f_p1 =     numpy.zeros(nlayers) # Flow out of below layer
    pe_m1 =    numpy.zeros(nlayers) # Peclet number of above layer
    pe_p1 =    numpy.zeros(nlayers) # Peclet number of below layer

    # Coefficients of tridiagonal problem: a_i*x_(i-1) + b_i*(x_i) + c_i*x_(i+1) = r_i
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
        adv_flux[0]=0.0
    if numpy.ndim(source) == 0:
        source = numpy.zeros(nlayers)+source


    # Set Pe (Peclet #) and D/dz throughout column
    # This is setting weights and coefficients that will be used for calculating mass balance equation in the tridiagonal step
    for j in range(nlayers):

        # dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
        # dz_node_tracer is difference between cell centers 

        # Calculate the D and F terms in the Patankar algorithm
        if (j == 0) :
            d_m1_zm1[j] = 0.
            w_p1 = (zsoi[j+1] - zisoi[j+1]) / dz_node[j+1]
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
            w_m1 = (zisoi[j] - zsoi[j-1]) / dz_node[j]
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
            w_m1 = (zisoi[j] - zsoi[j-1]) / dz_node[j]
            if ( diffus[j-1] > 0. and diffus[j] > 0.) :
                d_m1 = 1. / ((1. - w_m1) / diffus[j] + w_m1 / diffus[j-1]) # Harmonic mean of diffus
            else:
                d_m1 = 0.
            
            w_p1 = (zsoi[j+1] - zisoi[j+1]) / dz_node[j+1]
            if ( diffus[j+1] > 0. and diffus[j] > 0.) :
                d_p1 = 1. / ((1. - w_p1) / diffus[j] + w_p1 / diffus[j+1]) # Harmonic mean of diffus
            else:
                d_p1 = (1. - w_p1) * diffus[j] + w_p1 * diffus[j+1] # Arithmetic mean of diffus
            
            d_m1_zm1[j] = d_m1 / dz_node[j]
            d_p1_zp1[j] = d_p1 / dz_node[j+1]
            f_m1[j] = adv_flux[j]
            f_p1[j] = adv_flux[j+1]
            pe_m1[j] = f_m1[j] / d_m1_zm1[j] # Peclet #
            pe_p1[j] = f_p1[j] / d_p1_zp1[j] # Peclet #



    # Calculate the tridiagonal coefficients
    # Coefficients of tridiagonal problem: a_i*x_(i-1) + b_i*(x_i) + c_i*x_(i+1) = r_i
    # Here, this is equivalent to Patankar equation 5.56 and 5.57 (but in one dimension):
    # a_P*phi_P = a_E*phi_E + a_W*phi_W + b [phi is concentration, = x in tridiagonal]. Converting East/West to above/below
    # -> -a_E*phi_E + a_P*phi_P - a_W+phi_W = b
    # -a_tri = a_above = D_above*A(Pe)+max(-F_above,0); D_above=diffus_above/dz
    # b_tri = a_above+a_below+(F_above-F_below)+rho*dz/dt
    # -c_tri = D_below*A(Pe)+max(F_below,0); D_below = diffus_below/dz
    # r_tri = b = source_const*dz + conc*rho*dz/dt

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
            r_tri[j+1] = source[j] * dz[j] + (a_p_0 - adv_flux[j]) * conc_trcr[j] 
        elif (j < nlayers) :
            a_p_0 =  dz[j] / dtime
            a_tri[j+1] = -(d_m1_zm1[j] * aaa(pe_m1[j]) + max( f_m1[j], 0.)) # Eqn 5.47 Patankar
            c_tri[j+1] = -(d_p1_zp1[j] * aaa(pe_p1[j]) + max(-f_p1[j], 0.))
            b_tri[j+1] = -a_tri[j+1] - c_tri[j+1] + a_p_0
            r_tri[j+1] = source[j] * dz[j] + a_p_0 * conc_trcr[j]
        else: #  0 concentration gradient at bottom
            a_tri[j+1] = -1.
            b_tri[j+1] = 1.
            c_tri[j+1] = 0. 
            r_tri[j+1] = 0.

    # Solve for concentration profile
    if solver=='scipy':
        conc_after = solve_banded((1,1),numpy.row_stack((a_tri,b_tri,c_tri)),r_tri)
    elif solver=='ELM_tridiag':
        conc_after=tridiag_ELM(a_tri,b_tri,c_tri,r_tri)
    else:
        raise ValueError('Unknown solver %s'%solver)

    if printmatrix:
        import pandas
        df=pandas.DataFrame(data={'a':a_tri,'b':b_tri,'c':c_tri,'r':r_tri,})
        df2=pandas.DataFrame(data={'ap0':dz/dtime,'pe_m':pe_m1,'pe_p':pe_p1,'f_m':f_m1,'f_p':f_p1,'d_m':d_m1_zm1*dz_node,'d_p':d_p1_zp1*numpy.append(dz_node[1:],dz_node[-1]),'diff':conc_after[1:-1]-conc_trcr},index=numpy.arange(1,nlayers+1))
        print(df.merge(df2,left_index=True,right_index=True,how='outer').to_string(float_format='%1.6e'))
        print('Conc diff = %1.4e'%((conc_after[1:-1]*dz).sum()-(conc_trcr*dz).sum()))

    if return_matrix:
        return(conc_after,df)
    return conc_after

def run_steps(conc,nsteps,*args,**kwargs):
    for step in range(nsteps):
        conc=adv_diff(conc,*args,**kwargs)[1:-1]
    return conc

def conv(s):
   return numpy.array([float(n) for n in s.split()])

def tridiag_ELM(a,b,c,r):
    jtop=0
    lbj=0
    ubj=len(a)

    gam=numpy.zeros(len(a))
    u=numpy.zeros(len(a))

    bet = b[jtop]

    for j in range(lbj, ubj):
        if (j >= jtop):
            if (j == jtop):
                u[j] = r[j] / bet
            else:
                gam[j] = c[j-1] / bet
                bet = b[j] - a[j] * gam[j]
                u[j] = (r[j] - a[j]*u[j-1]) / bet

    for j in range(ubj-2,lbj,-1):
        if (j >= jtop):
            u[j] = u[j] - gam[j+1] * u[j+1]

    return u



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    z=numpy.linspace(0.05,10,20)
    # z=numpy.array([7.10063521e-03, 2.79249996e-02, 6.22585751e-02, 1.18865065e-01, # ELM layer thicknesses
    #    2.12193400e-01, 3.66065800e-01, 6.19758487e-01, 1.03802705e+00,
    #    1.72763526e+00, 2.86460710e+00, 4.73915672e+00, 7.82976627e+00,
    #    1.29253206e+01, 2.13264694e+01, 3.51776199e+01])[:10]
    nlayers=len(z)
    conc=numpy.linspace(0.1,0,nlayers)**2
    # conc=numpy.zeros(nlayers)
    conc[3]=0.2
    conc[-4]=0.2


    # conc[-1]=0.0
    # conc[10]=1.0

    dz=numpy.diff(z)
    dz=numpy.append(dz,dz[-1])

    dt=3600/2
    adv=-5e-5 # Negative means downward velocity. Layer zero flux is negative of infiltration
    adv=numpy.linspace(-5e-5,0,nlayers)
    adv[0]=0.0
    adv[-1]=0.0
    diffco=5e-5
    diffco=numpy.zeros(nlayers)
    diffco[:8]=1e-3
    diffco[8:]=1e-5

    conc=conv('5.8473004272365748E-003   5.3766436140742467E-003   4.0963920138112264E-003   3.3715881363458491E-003   3.0643021179267644E-003   2.8091586335862299E-003   2.7519242675744958E-003   2.7933799778690465E-003   3.0975783630216246E-003   3.0975783630216246E-003')
    adv=conv('-0.0000000000000000        3.6156883648582949E-009   1.0751958543366249E-009   4.3556589863946244E-010   5.2408802310046698E-010   1.2948885408737850E-009   4.7926796534336999E-010  -4.9948235117437907E-010  -6.8871993914108431E-010  -4.2306083932048843E-012   5.7766670618938444E-011')[:-1]
    diffco=conv('1.6076247457285592E-011   4.7824915654421890E-011   5.5988558307227856E-011   5.6023462843649933E-011   5.6033643445094697E-011   5.6013057307543900E-011   5.6015241263449222E-011   5.6020423167321922E-011   5.6022624794902403E-011   5.6021196434299971E-011')
    # dz=conv('1.7512817916255204E-002   2.7578969259676251E-002   4.5470033242413201E-002   7.4967410986208557E-002  0.12360036510228053       0.20378255101043175       0.33598062644843263       0.55393840536868488       0.91329003158906108        1.5057607013992766        2.4825796969813321        4.0930819526214002        6.7483512780057175        11.126150294204420        13.851152141963599')[:11]
    z=conv('7.1006354171935350E-003   2.7925000415316870E-002   6.2258573936546040E-002  0.11886506690014327       0.21219339590896316       0.36606579710470433       0.61975849792982662        1.0380270500015696        1.7276353086671965        2.8646071131796917        4.7391567114657498        7.8297665071423559        12.925320616708550        21.326469063153791        35.177621205117390')[:10]
    dz=conv('1.7512817916255204E-002   2.7578969259676251E-002   4.5470033242413201E-002   7.4967410986208557E-002  0.12360036510228053       0.20378255101043175       0.33598062644843263       0.55393840536868488       0.91329003158906108        1.5057607013992766        2.4825796969813321        4.0930819526214002        6.7483512780057175        11.126150294204420        13.851152141963599')[:10]

    source=numpy.zeros(nlayers)
    # source[4]=1e-5
    # source[10]=1e-5
    
    f,a=plt.subplots(2,3,num='Result',clear=True)
    for num in range(3):
        a[0,num].plot(z,conc,label='Original',c='k')
        a[1,num].axhline(0.0,ls='--',c='k',lw=0.5)
        a[1,num].axhline((conc*dz).sum(),ls='--',c='k',lw=0.5)
    a[0,0].plot(z,adv_diff(conc,z,dz,dt,adv_flux=adv,source=source,printmatrix=False)[1:-1],label='Adv',c='C1')
    a[0,1].plot(z,adv_diff(conc,z,dz,dt,diffus=diffco,source=source,printmatrix=False)[1:-1],label='Diff',c='C1')
    a[0,2].plot(z,adv_diff(conc,z,dz,dt,diffus=diffco,adv_flux=adv,source=source,printmatrix=True)[1:-1],label='Diff and adv',c='C1')

    a[0,0].plot(z,run_steps(conc,4,z,dz,dt,adv_flux=adv,source=source),ls='--',c='C1')
    a[0,1].plot(z,run_steps(conc,4,z,dz,dt,diffus=diffco,source=source),ls='--',c='C1')
    a[0,2].plot(z,run_steps(conc,4,z,dz,dt,diffus=diffco,adv_flux=adv,source=source),ls='--',c='C1')

    a[0,0].plot(z,run_steps(conc,10,z,dz,dt,adv_flux=adv,source=source),ls=':',c='C1')
    a[0,1].plot(z,run_steps(conc,10,z,dz,dt,diffus=diffco,source=source),ls=':',c='C1')
    a[0,2].plot(z,run_steps(conc,10,z,dz,dt,diffus=diffco,adv_flux=adv,source=source),ls=':',c='C1')

    a[0,0].plot(z,adv_diff(conc,z,dz,dt,adv_flux=-adv,source=source,printmatrix=False,solver='ELM_tridiag')[1:-1],label='Adv (ELM)',c='C2')
    a[0,1].plot(z,adv_diff(conc,z,dz,dt,diffus=diffco,source=source,printmatrix=False,solver='ELM_tridiag')[1:-1],label='Diff (ELM)',c='C2')
    a[0,2].plot(z,adv_diff(conc,z,dz,dt,diffus=diffco,adv_flux=-adv,source=source,printmatrix=True,solver='ELM_tridiag')[1:-1],label='Diff and adv (ELM)',c='C2')

    a[0,0].plot(z,run_steps(conc,4,z,dz,dt,adv_flux=-adv,source=source,solver='ELM_tridiag'),ls='--',c='C2')
    a[0,1].plot(z,run_steps(conc,4,z,dz,dt,diffus=diffco,source=source,solver='ELM_tridiag'),ls='--',c='C2')
    a[0,2].plot(z,run_steps(conc,4,z,dz,dt,diffus=diffco,adv_flux=-adv,source=source,solver='ELM_tridiag'),ls='--',c='C2')

    a[0,0].plot(z,run_steps(conc,10,z,dz,dt,adv_flux=-adv,source=source,solver='ELM_tridiag'),ls=':',c='C2')
    a[0,1].plot(z,run_steps(conc,10,z,dz,dt,diffus=diffco,source=source,solver='ELM_tridiag'),ls=':',c='C2')
    a[0,2].plot(z,run_steps(conc,10,z,dz,dt,diffus=diffco,adv_flux=-adv,source=source,solver='ELM_tridiag'),ls=':',c='C2')



    a[1,0].plot(z,(conc*dz).cumsum(),label='Original',c='C0')
    a[1,0].plot(z,(adv_diff(conc,z,dz,dt,adv_flux=adv,source=source)[1:-1]*dz).cumsum(),ls='-',c='C1')
    a[1,1].plot(z,(adv_diff(conc,z,dz,dt,diffus=diffco,source=source)[1:-1]*dz).cumsum(),ls='-',c='C2')
    a[1,2].plot(z,(adv_diff(conc,z,dz,dt,diffus=diffco,adv_flux=adv,source=source)[1:-1]*dz).cumsum(),ls='-',c='C3')

    a[1,0].plot(z,(run_steps(conc,4,z,dz,dt,adv_flux=adv,source=source)*dz).cumsum(),label='Adv = 5e-4',c='C1',ls='--')
    a[1,1].plot(z,(run_steps(conc,4,z,dz,dt,diffus=diffco,source=source)*dz).cumsum(),label='Diff = 1e-5',c='C2',ls='--')
    a[1,2].plot(z,(run_steps(conc,4,z,dz,dt,diffus=diffco,adv_flux=adv,source=source)*dz).cumsum(),label='Diff and adv',c='C3',ls='--')

    a[1,0].plot(z,(run_steps(conc,10,z,dz,dt,adv_flux=adv,source=source)*dz).cumsum(),label='Adv = 5e-4',c='C1',ls=':')
    a[1,1].plot(z,(run_steps(conc,10,z,dz,dt,diffus=diffco,source=source)*dz).cumsum(),label='Diff = 1e-5',c='C2',ls=':')
    a[1,2].plot(z,(run_steps(conc,10,z,dz,dt,diffus=diffco,adv_flux=adv,source=source)*dz).cumsum(),label='Diff and adv',c='C3',ls=':')


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