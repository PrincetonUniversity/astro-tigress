def calc_IQU(nH,Bx,By,Bz,deltas,los='y'):
    '''
        Flat sky projection of dust polarization at 353 GHz
        Parameters assumed are identical to Kim, Choi, and Flauger (2019)
        Tdust = 18K
        sigma_d,353 = 1.2e-26 cm^2
        p0 = 0.2
    '''
    if los=='x':
        Blos = Bx
        Bperp_x = By
        Bperp_y = Bz
        ilos = 2
    elif los=='y':
        Blos = By
        Bperp_x = -Bx
        Bperp_y = Bz
        ilos = 1
    elif los=='z':
        Blos = Bz
        Bperp_x = Bx
        Bperp_y = By
        ilos = 0

    args={'Bnu':41495.876171482356, 'sigma':1.2e-26, 'p0':0.2}
    Bnu=args['Bnu']
    p0=args['p0']
    sigma=args['sigma']

    Bperp2=Bperp_x*Bperp_x+Bperp_y*Bperp_y
    B2=Bperp2+Blos*Blos
    cos2phi=(Bperp_y*Bperp_y-Bperp_x*Bperp_x)/Bperp2
    sin2phi=-Bperp_x*Bperp_y/Bperp2*2
    cosgam2=Bperp2/B2

    ds=deltas*3.085677581467192e+18
    dtau=sigma*nH*ds

    I=Bnu*(1.0-p0*(cosgam2-2./3.0))*dtau
    Q=p0*Bnu*cos2phi*cosgam2*dtau
    U=p0*Bnu*sin2phi*cosgam2*dtau

    return I.sum(axis=ilos),Q.sum(axis=ilos),U.sum(axis=ilos)
