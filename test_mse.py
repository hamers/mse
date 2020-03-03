import numpy as np

from mse import MSE,Particle,Tools

"""
Several routines for testing the code/installation. 
To use, run `python test_MSE.py i', where i is
the test number.

Adrian Hamers, June 2019
"""

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class test_mse():
    def test0(self):
        code = MSE()

    def test1(self):
        """
        test reference system of Naoz et al. (2009)
        """

        particles = Tools.create_nested_multiple(3, [1.0,1.0e-3,40.0e-3],[6.0,100.0],[0.001,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])

        code = MSE()
        code.add_particles(particles)
        inner_binary = code.particles[3]
        outer_binary = code.particles[4]
        
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        
        t = 0.0
        N = 1000
        tend = 3.0e7
        dt = tend/float(N)
        while t<=tend:
            code.evolve_model(t)
            t+=dt
        
            print( 't',t,'e',inner_binary.e,'INCL',inner_binary.INCL,outer_binary.INCL)
                        
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            rel_INCL_print.append(inner_binary.INCL_parent)
            e_print.append(inner_binary.e)
            INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        t_print = np.array(t_print)
        rel_INCL_print = np.array(rel_INCL_print)
        e_print = np.array(e_print)
        
        if HAS_MATPLOTLIB==True:
            fig=pyplot.figure()
            plot1=fig.add_subplot(2,1,1)
            plot2=fig.add_subplot(2,1,2,yscale="log")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            plot1.plot(1.0e-6*t_print,rel_INCL_print*180.0/np.pi)
            plot2.plot(1.0e-6*t_print,1.0-e_print)
            pyplot.show()

    def test2(self):
        pass
        """
        test 1PN precession in 2-body system
        """
        particles = Tools.create_nested_multiple(2,[1.0, 1.0], [1.0], [0.99], [0.01], [0.01], [0.01])
        binaries = [x for x in particles if x.is_binary == True]

        for b in binaries:
            b.include_pairwise_1PN_terms = True

        code = MSE()
        code.add_particles(particles)

        binaries = [particles[2]]

        t = 0.0
        N=100
        tend = 1.0e6
        dt=tend/float(N)

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)

            print( 't/Myr',t,'omega',binaries[0].AP)
            t_print_array.append(t)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)
            AP_print_array.append(binaries[0].AP)
        t_print_array = np.array(t_print_array)
        a_print_array = np.array(a_print_array)
        e_print_array = np.array(e_print_array)
        AP_print_array = np.array(AP_print_array)
        
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C

        a = binaries[0].a
        e = binaries[0].e
        M = binaries[0].mass
        rg = CONST_G*M/(CONST_C**2)
        P = 2.0*np.pi*np.sqrt(a**3/(CONST_G*M))
        t_1PN = (1.0/3.0)*P*(1.0-e**2)*(a/rg)
        
        
        if HAS_MATPLOTLIB == True:
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array*1.0e-6,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,N)
            AP = 0.01 +2.0*np.pi*points/(t_1PN*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g')

            plot2.plot(t_print_array*1.0e-6,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array*1.0e-6,np.fabs((a-a_print_array)/a), color='r')
            plot4.plot(t_print_array*1.0e-6,np.fabs((e-e_print_array)/e), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)

            pyplot.show()

    def test3(self):
        """
        test GW emission in 2-body system + collision detection
        """
        a0 = 1.0
        e0 = 0.999
        m1 = 1.0
        m2 = 1.0
        particles = Tools.create_nested_multiple(2,[m1, m2], [a0], [e0], [0.01], [0.01], [0.01])

        binary = particles[2]
        #stars = particles - binaries
        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        
        rg = (m1+m2)*CONST_G/(CONST_C**2)
        print( 'rg/AU',rg)
        for i in range(2):
            particles[i].radius = 100.0*rg

        binary.check_for_physical_collision_or_orbit_crossing = True
        binary.include_pairwise_25PN_terms = True



        code.add_particles(particles)
        binary = code.particles[2]

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        tend = 1.0e8
        N = 0
        t = 0.0
        dt = 1.0e6
        while (t<tend):
            t+=dt
            N+=1
            code.evolve_model(t)
            flag = code.flag

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)


            print( 'e',binary.e,'a',binary.a)

            if flag == 2:
                print( 'root found')
                break

        if HAS_MATPLOTLIB == True:
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(2,1,1)
            plot2 = fig.add_subplot(2,1,2,yscale="log")

            plot1.plot(t_print_array,e_print_array, color='r')

            plot2.plot(t_print_array,a_print_array, color='r')

            ### Peters 1964 ###
            c0 = a0*(1.0-e0**2)/( pow(e0,12.0/19.0)*pow(1.0 + (121.0/304.0)*e0**2,870.0/2299.0))
            a_an = c0*pow(e_print_array,12.0/19.0)*pow(1.0+(121.0/304.0)*e_print_array**2,870.0/2299.0)/(1.0-e_print_array**2)
            beta = (64.0/5.0)*CONST_G**3*m1*m2*(m1+m2)/(CONST_C**5)
            #T = c0**4*pow(e0,48.0/19.0)/(4.0*beta)
            T_c = a0**4/(4.0*beta)
            T = (768.0/425.0)*T_c*pow(1.0-e0**2,7.0/2.0)
            print( 'T/Myr (approx)',T*1.0e-6)
            plot2.plot(t_print_array,a_an,color='g',linestyle='dashed',linewidth=2)

            fontsize = 15
            plot1.set_ylabel("$e$",fontsize=fontsize)
            plot2.set_ylabel("$a/\mathrm{AU}$",fontsize=fontsize)
            plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    
    def test4(self):
        """
        test tidal friction in 2-body system
        """
        
        code = MSE()
        code.enable_tides = True
        code.enable_root_finding = True

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        M = 0.0009546386983890755 ### Jupiter mass
        R = 40.0*0.1027922358015816*CONST_R_SUN ### Jupiter radius = 0.1 R_SUN
        m_per = 1.0
        mu = m_per*M/(m_per+M)
        a0 = 0.1
        e0 = 0.3
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        aF = a0*(1.0-e0**2)
        nF = np.sqrt( CONST_G*(M+m_per)/(aF**3) )

        particles = Tools.create_nested_multiple(2, [m_per, M], [a0], [e0], [0.01], [0.01], [0.01])
        binary = particles[2]
        
        
        particles[0].radius = CONST_R_SUN
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0
        particles[1].spin_vec_y = 0.0
        particles[1].spin_vec_z = 4.0e-2/day

        k_L = 0.38
        k_AM = k_L/2.0
        rg = 0.25
        tau = 0.66*second

        I = rg*M*R**2
        alpha = I/(mu*a0**2)
        T = R**3/(CONST_G*M*tau)
        t_V = 3.0*(1.0 + 1.0/k_L)*T
        print( 't_V',t_V,'M',M,'R',R)

        particles[0].include_tidal_friction_terms = False
        
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = True
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-8

        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_viscous_time_scale = t_V
        particles[1].tides_gyration_radius = rg

        tD = M*aF**8/(3.0*k_L*tau*CONST_G*m_per*(M+m_per)*R**5)
        particles[2].check_for_physical_collision_or_orbit_crossing = True

        code.add_particles(particles)
        binary = code.particles[2]
        #code.parameters.relative_tolerance = 1.0e-14


        t = 0.0
        N=100
        tend = 1.0e2
        dt = tend/float(N)

        t_print_array = []
        a_print_array = []
        n_print_array = []
        e_print_array = []
        AP_print_array = []
        spin_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            print( 'flag',code.flag,t,'a/AU',binary.a,'e',binary.e)


            t_print_array.append(t)
            a_print_array.append(binary.a)
            n_print_array.append(np.sqrt(CONST_G*(M+m_per)/(binary.a**3)))
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)
            spin_print_array.append( np.sqrt( particles[1].spin_vec_x**2 + particles[1].spin_vec_y**2 + particles[1].spin_vec_z**2) )

            bodies = particles[0:2]
            for body in bodies:
                print( 'S_x',body.spin_vec_x)
                print( 'S_y',body.spin_vec_y)
                print( 'S_z',body.spin_vec_z)
            print( '='*50)

        if HAS_MATPLOTLIB == True:
            fig = pyplot.figure()
            fontsize=12
    
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)
            spin_print_array = np.array(spin_print_array)
            n_print_array = np.array(n_print_array)
            
            N_p = 4
            plot1 = fig.add_subplot(N_p,1,1)
            plot1.plot(t_print_array*1.0e-6,a_print_array, color='r')
            plot1.set_ylabel("$a/\mathrm{AU}$",fontsize=fontsize)

            plot2 = fig.add_subplot(N_p,1,2)
            plot2.plot(t_print_array*1.0e-6,e_print_array,color='k')
            plot2.set_ylabel("$e$",fontsize=fontsize)

            plot3 = fig.add_subplot(N_p,1,3,yscale="log")

            plot3.plot(t_print_array*1.0e-6,a_print_array*(1.0-e_print_array**2),color='k')
            plot3.axhline(y = a0*(1.0 - e0**2), color='k')
            plot3.set_ylabel("$a(1-e^2)/\mathrm{AU}$",fontsize=fontsize)

            plot4 = fig.add_subplot(N_p,1,4)
            plot4.plot(t_print_array*1.0e-6,spin_print_array/n_print_array)
            plot4.set_ylabel("$\Omega/n$",fontsize=fontsize)

            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    def test5(self):
        """
        test precession due to tidal bulges
        """
        code = MSE()
        code.enable_root_finding = True
        code.enable_tides = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        M = 0.0009546386983890755 ### Jupiter mass
        R = 1.0*0.1027922358015816*CONST_R_SUN ### Jupiter radius ~ 0.1 R_SUN

        m_per = 1.0
        a0 = 30.0
        e0 = 0.999
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        particles = Tools.create_nested_multiple(2, [m_per, M], [a0], [e0], [0.01], [0.01], [0.01])
        binary = particles[2]
        particles[0].radius = 1.0*CONST_R_SUN
        particles[1].radius = R

        k_L = 0.41
        k_AM = k_L/2.0

        particles[1].tides_method = 0
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = True
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_gyration_radius = 0.25
        
        code.add_particles(particles)

        t = 0.0
        dt = 1.0e6
        tend = 1.0e8

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        g_dot_TB = (15.0/8.0)*n0*(8.0+12.0*e0**2+e0**4)*(m_per/M)*k_AM*pow(R/a0,5.0)/pow(1.0-e0**2,5.0)
        t_TB = 2.0*np.pi/g_dot_TB

        N=0
        while (t<tend):
            t+=dt
            code.evolve_model(t)
            print( 'flag',code.flag,t,binary.a,binary.e)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

            N+=1
        if HAS_MATPLOTLIB == True:
            
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)
            
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,N)
            AP = 0.01 +2.0*np.pi*points/(t_TB*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',linestyle='dotted',linewidth=2)

            plot2.plot(t_print_array,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array,np.fabs((a0-a_print_array)/a0), color='r')
            plot4.plot(t_print_array,np.fabs((e0-e_print_array)/e0), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)

            pyplot.show()
        #exit(-1)
    def test6(self):
        """
        test precession due to rotation
        """

        code = MSE()
        code.enable_root_finding = True
        code.enable_tides = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)
        
        M = 0.0009546386983890755 ### Jupiter mass
        R = 1.5*0.1027922358015816*CONST_R_SUN ### Jupiter radius ~ 0.1 R_SUN
        m_per = 1.0
        a0 = 30.0
        e0 = 0.999
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        aF = a0*(1.0-e0**2)
        nF = np.sqrt( CONST_G*(M+m_per)/(aF**3) )

        particles = Tools.create_nested_multiple(2, [m_per, M], [a0], [e0], [1.0e-5], [1.0e-5], [1.0e-5])
        binary = particles[2]
        particles[0].radius = 1.0
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0
        particles[1].spin_vec_y = 0.0
        Omega_PS0 = n0*(33.0/10.0)*pow(a0/aF,3.0/2.0)
        particles[1].spin_vec_z = Omega_PS0


        k_L = 0.51
        k_AM = k_L/2.0
        rg = 0.25
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = True
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_gyration_radius = rg

        code.add_particles(particles)

        t = 0.0
        dt = 1.0e6
        tend = 1.0e8

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        Omega_vec = [particles[1].spin_vec_x,particles[1].spin_vec_y,particles[1].spin_vec_z]
        Omega = np.sqrt(Omega_vec[0]**2 + Omega_vec[1]**2 + Omega_vec[2]**2)
        print( 'Omega/n',Omega/n0)

        g_dot_rot = n0*(1.0 + m_per/M)*k_AM*pow(R/a0,5.0)*(Omega/n0)**2/((1.0-e0**2)**2)
        t_rot = 2.0*np.pi/g_dot_rot
        print( 't_rot/Myr',t_rot*1.0e-6)

        N=0
        while (t<tend):
            t+=dt
            code.evolve_model(t)
            print( 'flag',code.flag,t,binary.a,binary.e)


            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

            N+=1
        if HAS_MATPLOTLIB == True:
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)

            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,N)
            AP = 0.01 +2.0*np.pi*points/(t_rot*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',linestyle='dotted',linewidth=2)

            plot2.plot(t_print_array,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array,np.fabs((a0-a_print_array)/a0), color='r')
            plot4.plot(t_print_array,np.fabs((e0-e_print_array)/e0), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)

            pyplot.show()

    def test7(self):
        """
        test collision detection in 3-body system
        """
        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)
        

        particles = Tools.create_nested_multiple(3,[1.0, 1.2, 0.9], [1.0, 100.0], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        binaries[0].check_for_physical_collision_or_orbit_crossing = True
        for body in bodies:
            body.radius = 0.03

        code.add_particles(particles)

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6

        t_print_array = []
        a_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag

            print( 'secular_breakdown_has_occurred',binaries[0].secular_breakdown_has_occurred)
            print( 'dynamical_instability_has_occurred',binaries[0].dynamical_instability_has_occurred)
            print( 'physical_collision_or_orbit_crossing_has_occurred',binaries[0].physical_collision_or_orbit_crossing_has_occurred)
            print( 'minimum_periapse_distance_has_occurred',binaries[0].minimum_periapse_distance_has_occurred)
            print( 'RLOF_at_pericentre_has_occurred',binaries[0].RLOF_at_pericentre_has_occurred)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                print( 'root found')
                break
            print( 't_end',code.model_time)


        if HAS_MATPLOTLIB == True:
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure()
            plot = fig.add_subplot(2,1,1)
            plot.plot(t_print_array,e_print_array)

            plot = fig.add_subplot(2,1,2)
            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))
            plot.axhline(y = bodies[0].radius + bodies[1].radius,color='k')

            pyplot.show()

    def test8(self):
        """
        test minimum periapse occurrence
        """

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        particles = Tools.create_nested_multiple(3,[1.0, 1.2, 0.9], [1.0, 100.0], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        binaries[0].check_for_minimum_periapse_distance = True
        rp_min = 0.1
        binaries[0].check_for_minimum_periapse_distance_value = rp_min
        
        code.add_particles(particles)

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6

        t_print_array = []
        a_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag

            print( 'secular_breakdown_has_occurred',binaries[0].secular_breakdown_has_occurred)
            print( 'dynamical_instability_has_occurred',binaries[0].dynamical_instability_has_occurred)
            print( 'physical_collision_or_orbit_crossing_has_occurred',binaries[0].physical_collision_or_orbit_crossing_has_occurred)
            print( 'minimum_periapse_distance_has_occurred',binaries[0].minimum_periapse_distance_has_occurred)
            print( 'RLOF_at_pericentre_has_occurred',binaries[0].RLOF_at_pericentre_has_occurred)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                print( 'root found')
                break
            print( 't_end',code.model_time)


        if HAS_MATPLOTLIB == True:
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure()
            plot = fig.add_subplot(2,1,1)
            plot.plot(t_print_array,e_print_array)

            plot = fig.add_subplot(2,1,2)
            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))
            plot.axhline(y = rp_min,color='k')

            pyplot.show()

    def test9(self):
        """
        test adiabatic mass loss
        """

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        m1i = 1.0
        m2i = 1.2
        m3i = 0.9
        a1i = 1.0
        a2i = 100.0
        particles = Tools.create_nested_multiple(3,[m1i,m2i,m3i], [a1i,a2i], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

#        for index,body in enumerate(bodies):
#            if index==0:
        m1dot = -1.0e-7
        m2dot = -1.0e-7
        m3dot = -1.0e-7
        bodies[0].mass_dot = m1dot
        bodies[1].mass_dot = m2dot
        bodies[2].mass_dot = m3dot

        code.add_particles(particles)
        #particles = code.particles
        #bodies = [x for x in particles if x.is_binary==False]
        #binaries = [x for x in particles if x.is_binary==True]

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6

        t_print_array = []
        a1_print_array = []
        a2_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag

            t_print_array.append(t*1.0e-6)
            a1_print_array.append(binaries[0].a)
            a2_print_array.append(binaries[1].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                print( 'root found')
                break
            print( 't_end',code.model_time,bodies[0].mass,bodies[1].mass,bodies[2].mass,binaries[0].a)


        if HAS_MATPLOTLIB == True:
            a1_print_array = np.array(a1_print_array)
            a2_print_array = np.array(a2_print_array)
            e_print_array = np.array(e_print_array)

            t = pow(10.0,np.linspace(4.0,6.0,100))
            m1 = m1i + m1dot*t
            m2 = m2i + m2dot*t
            m3 = m3i + m3dot*t
            
            a1 = a1i*(m1i+m2i)/(m1+m2)
            a2 = a2i*(m1i+m2i+m3i)/(m1+m2+m3)

            fig = pyplot.figure()
            plot = fig.add_subplot(1,1,1,yscale="log")
            plot.plot(t_print_array,a1_print_array,color='k',linestyle='solid',label='$a_1$')
            plot.plot(t_print_array,a2_print_array,color='k',linestyle='solid',label='$a_2$')
            
            plot.plot(t*1.0e-6,a1,color='g',linestyle='dashed',linewidth=2,label='$a_1\,(\mathrm{an})$')
            plot.plot(t*1.0e-6,a2,color='g',linestyle='dashed',linewidth=2,label='$a_2\,(\mathrm{an})$')

            handles,labels = plot.get_legend_handles_labels()
            plot.legend(handles,labels,loc="lower left",fontsize=18)

#            plot = fig.add_subplot(2,1,2)
#            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))


            pyplot.show()

    def test10(self):
        """
        test flybys module: numerically integrating over perturber's orbit
        """

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        a = 10.0
        e = 0.1
        masses = [1.0, 0.8]
        m = masses[0] + masses[1]
        M = 1.0
        E = 10.0
        Q = 100.0
        INCL = 0.5*np.pi
        AP = 0.25*np.pi
        LAN = 0.25*np.pi
        particles = Tools.create_nested_multiple(2,masses, [a], [e], [INCL], [AP], [LAN])
        
        t = 0.0
        #t_ref = 1.0e-2 | units.Myr
        t_ref = 1.0e5
        tend = 2.0*t_ref
        Nsteps = 100
        dt = tend/float(Nsteps)

        import np.random as randomf
        randomf.seed(0)
        N=1
        external_particle = Particle(mass = M, is_binary=True, is_external=True, external_t_ref=t_ref, e=E, external_r_p = Q, INCL = 1.0e-10, AP = 1.0e-10, LAN = 1.0e-10)

        particles.append(external_particle)

        
        
        code = MSE()
        code.add_particles(particles)
        binary = code.particles[2]
        
        #code.parameters.include_quadrupole_order_terms = True
        #code.parameters.include_octupole_order_binary_pair_terms = True
        #code.parameters.include_hexadecupole_order_binary_pair_terms = True
        #code.parameters.include_dotriacontupole_order_binary_pair_terms = True

        t_print_array = []
        a_print_array = []
        e_print_array = []
        INCL_print_array = []
        AP_print_array = []
        LAN_print_array = []
        
        N = 0
        while (t<tend):

            t+=dt
            N+=1

            code.evolve_model(t)
             
            
            t_print_array.append(t)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
#            INCL_print_array.append(binaries[0].inclination | units.none)            
#            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)            
#            LAN_print_array.append(binaries[0].longitude_of_ascending_node | units.none)            
            
        Delta_e = binary.e-e
        print( 'Delta_e',Delta_e)
            
        if HAS_MATPLOTLIB == True:
            fig = pyplot.figure(figsize=(8,6))
            plot1 = fig.add_subplot(1,1,1)
            
            t_print_array = np.array(t_print_array)
            e_print_array = np.array(e_print_array)
            
            plot1.plot(t_print_array*1.0e-6,e_print_array, color='k')
            #plot2.plot(t_print_array*1.0e-6,(180.0/np.pi)*INCL_print_array.value_in(units.none), color='k')
            #plot3.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*AP_print_array.value_in(units.none), color='k')
            #plot4.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*LAN_print_array.value_in(units.none), color='k')
            
            
            fontsize = 15
            plot1.set_ylabel("$e$",fontsize=fontsize)
            
            
            pyplot.show()
        exit(-1)
    
    def test11(self):
        """
        test flybys module: using analytic formulae
        """

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        a = 10.0
        e = 0.1
        masses = [1.0, 0.8]
        m = masses[0] + masses[1]
        M = 1.0
        E = 10.0
        Q = 100.0
        INCL = 0.5*np.pi
        AP = 0.25*np.pi
        LAN = 0.25*np.pi
        particles = Tools.create_nested_multiple(2,masses, [a], [e], [INCL], [AP], [LAN])
        
        t = 0.0
        #t_ref = 1.0e-2 | units.Myr
        t_ref = 1.0e5
        tend = 2.0*t_ref
        Nsteps = 100
        dt = tend/float(Nsteps)

        import np.random as randomf
        randomf.seed(0)
        N=1
        external_particle = Particle(mass = M, is_binary=True, is_external=True, external_t_ref=t_ref, e=E, external_r_p = Q, INCL = 1.0e-10, AP = 1.0e-10, LAN = 1.0e-10)

        particles.append(external_particle)

        
        
        code = MSE()
        code.add_particles(particles)
        binary = code.particles[2]
        
        #code.parameters.include_quadrupole_order_terms = True
        #code.parameters.include_octupole_order_binary_pair_terms = True
        #code.parameters.include_hexadecupole_order_binary_pair_terms = True
        #code.parameters.include_dotriacontupole_order_binary_pair_terms = True
        
        code.apply_external_perturbation_assuming_integrated_orbits()
            
        Delta_e = binary.e-e
        print( 'Delta_e',Delta_e)

        exit(-1)

    def test12(self):
        """
        test flybys module: instantaneous change -- binary
        """

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        a = 10.0
        e = 0.1
        masses = [1.0, 0.8]
        m = masses[0] + masses[1]
        INCL = 0.1
        AP = 0.2
        LAN = 0.3
        particles = Tools.create_nested_multiple(2,masses, [a], [e], [INCL], [AP], [LAN])
        
        binary = particles[2]
        binary.sample_orbital_phase_randomly = 0
        binary.TA = 60.5*np.pi/180.0
        
        delta_m1 = -0.5
        #delta_m1 = 0.0
        #delta_m1 = 0.0 | units.MSun
        
        
        
        km_p_s_to_AU_p_yr = 0.21094502112788768
        V_k_vec = np.array([0.0,1.0,1.0])*km_p_s_to_AU_p_yr
        #V_k_vec = [0.0,0.0,0.0]
        particles[0].instantaneous_perturbation_delta_mass = delta_m1
        particles[0].instantaneous_perturbation_delta_vx = V_k_vec[0]
        particles[0].instantaneous_perturbation_delta_vy = V_k_vec[1]
        particles[0].instantaneous_perturbation_delta_vz = V_k_vec[2]

        #particles[0].instantaneous_perturbation_delta_mass = 0 | units.MSun
        
#        import np.random as randomf
#        randomf.seed(0)

        #binaries.include_pairwise_25PN_terms = True
        
        code.add_particles(particles)
        binary = code.particles[2]
        #code.external_particles.add_particles(external_particles)
        
        #code.parameters.include_quadrupole_order_terms = True
        #code.parameters.include_octupole_order_binary_pair_terms = True
        #code.parameters.include_hexadecupole_order_binary_pair_terms = True
        #code.parameters.include_dotriacontupole_order_binary_pair_terms = True

        print( '='*50)
        print( 'pre')
        print( 'a',binary.a,'e',binary.e,'INCL',binary.INCL*180.0/np.pi,'AP',binary.AP*180.0/np.pi,'LAN',binary.LAN*180.0/np.pi,binary.mass)
#        print 'r',particles.position_x,particles.position_y,particles.position_z
#        print 'v',particles.velocity_x.value_in(units.km/units.s),particles.velocity_y.value_in(units.km/units.s),particles.velocity_z.value_in(units.km/units.s)
        
#        c1 = particles[0]
#        c2 = particles[1]
        
#        r1_vec = [c1.position_x - c2.position_x,c1.position_y - c2.position_y,c1.position_z - c2.position_z]
#        v1_vec = [c1.velocity_x - c2.velocity_x,c1.velocity_y - c2.velocity_y,c1.velocity_z - c2.velocity_z]

#        code.parameters.orbital_phases_random_seed = 1
        
        code.apply_user_specified_instantaneous_perturbation()
        
        #channel_from_code_external.copy()

        print( '='*50)
        print( 'post')
        print( 'a',binary.a,'e',binary.e,'INCL',binary.INCL*180.0/np.pi,'AP',binary.AP*180.0/np.pi,'LAN',binary.LAN*180.0/np.pi,binary.mass)
      
        exit(-1)
        print( 'r',particles.position_x,particles.position_y,particles.position_z)
        print( 'v',particles.velocity_x.value_in(units.km/units.s),particles.velocity_y.value_in(units.km/units.s),particles.velocity_z.value_in(units.km/units.s))

        #code.apply_external_perturbation_assuming_integrated_orbits()

        
            
        t += 1.0 | units.Myr ### arbitrary
        t_print_array.append(t)
        a_print_array.append(binaries[0].semimajor_axis)
        e_print_array.append(binaries[0].eccentricity | units.none)
        INCL_print_array.append(binaries[0].inclination | units.none)            
        AP_print_array.append(binaries[0].argument_of_pericenter | units.none)            
        LAN_print_array.append(binaries[0].longitude_of_ascending_node | units.none)            

        print( 'a_print_array',a_print_array.value_in(units.AU))
        print( 'e_print_array',e_print_array)
        print( 'INCL_print_array',INCL_print_array*(180.0/np.pi))
        print( 'AP_print_array',AP_print_array*(180.0/np.pi))
        print( 'LAN_print_array',LAN_print_array*(180.0/np.pi))
        
        f1 = binaries.true_anomaly
        print( 'true_anomaly',f1)
        
        m1 = masses[0]
        m2 = masses[1]
        a1 = semimajor_axes[0]
        e1 = eccentricities[0]
        
        r1 = a1*(1.0-e1**2)/(1.0 + e1*np.cos(f1))
        

        r1_dot_v1 = np.sum([x*y for x,y in zip(r1_vec,v1_vec)])
        r1_dot_V_k = np.sum([x*y for x,y in zip(r1_vec,V_k_vec)])
        v1_dot_V_k = np.sum([x*y for x,y in zip(v1_vec,V_k_vec)])
        V_k_dot_V_k = np.sum([x*y for x,y in zip(V_k_vec,V_k_vec)])
        v1c_sq =constants.G*(m1+m2)/a1
        
        
        a1_p = a1*(1.0 + delta_m1/(m1+m2))*pow( 1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq, -1.0)
        j1_p = pow( (m1+m2)/(m1+m2+delta_m1), 2.0)*(1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq)*( 1.0 - e1**2 \
            + (1.0/(constants.G*(m1+m2)*a1))*( r1**2*( 2.0*v1_dot_V_k + V_k_dot_V_k) - 2.0*r1_dot_v1*r1_dot_V_k - r1_dot_V_k**2) )
        e1_p = np.sqrt(1.0 - j1_p)
        
        print( 'analytic results Toonen+16',a1_p.value_in(units.AU),e1_p)
        
        if HAS_MATPLOTLIB == True and 1==0:
            fig = pyplot.figure(figsize=(8,6))
            plot1 = fig.add_subplot(2,2,1)
            plot2 = fig.add_subplot(2,2,2)
            plot3 = fig.add_subplot(2,2,3)
            plot4 = fig.add_subplot(2,2,4)
            
            plot1.plot(t_print_array.value_in(units.Myr),e_print_array.value_in(units.none), color='k')
            plot2.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*INCL_print_array.value_in(units.none), color='k')
            plot3.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*AP_print_array.value_in(units.none), color='k')
            plot4.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*LAN_print_array.value_in(units.none), color='k')
            
            fontsize = 15
            plot1.set_ylabel("$e$",fontsize=fontsize)
            plot2.set_ylabel("$i$",fontsize=fontsize)
            plot3.set_ylabel("$\omega$",fontsize=fontsize)
            plot4.set_ylabel("$\Omega$",fontsize=fontsize)
            
#            fontsize = 15
#            plot1.set_ylabel("$e$",fontsize=fontsize)
            
            pyplot.show()


    def test13(self):
        """
        test flybys module: instantaneous change -- triple
        """

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        a1 = 10.0
        a2 = 100.0
        e1 = 0.1
        e2 = 0.3
        masses = [1.0, 0.8,1.2]
        m = masses[0] + masses[1]
        INCLs = [0.1,0.5]
        APs = [0.1,1.0]
        LANs = [0.1,2.0]
        particles = Tools.create_nested_multiple(3,masses, [a1,a2], [e1,e2], INCLs, APs, LANs)
        
        inner_binary = particles[3]
        outer_binary = particles[4]
        inner_binary.sample_orbital_phase_randomly = 0
        outer_binary.sample_orbital_phase_randomly = 0
        inner_binary.TA = 60.5*np.pi/180.0
        outer_binary.TA = 30.5*np.pi/180.0
        

        delta_m1 = -0.3
        #delta_m1 = 0.0
        #delta_m1 = 0.0 | units.MSun
        
        km_p_s_to_AU_p_yr = 0.21094502112788768
        V_k_vec = np.array([1.0,2.0,2.0])*km_p_s_to_AU_p_yr
        #V_k_vec = [0.0,0.0,0.0]
        particles[0].instantaneous_perturbation_delta_mass = delta_m1
        particles[0].instantaneous_perturbation_delta_vx = V_k_vec[0]
        particles[0].instantaneous_perturbation_delta_vy = V_k_vec[1]
        particles[0].instantaneous_perturbation_delta_vz = V_k_vec[2]

        #particles[0].instantaneous_perturbation_delta_mass = 0 | units.MSun
        
#        import np.random as randomf
#        randomf.seed(0)

        #binaries.include_pairwise_25PN_terms = True
        
        code.add_particles(particles)
        binary = code.particles[2]
        #code.external_particles.add_particles(external_particles)
        
        #code.parameters.include_quadrupole_order_terms = True
        #code.parameters.include_octupole_order_binary_pair_terms = True
        #code.parameters.include_hexadecupole_order_binary_pair_terms = True
        #code.parameters.include_dotriacontupole_order_binary_pair_terms = True

        print( '='*50)
        print( 'pre')
        print( 'inner','a',inner_binary.a,'e',inner_binary.e,'INCL',inner_binary.INCL*180.0/np.pi,'AP',inner_binary.AP*180.0/np.pi,'LAN',inner_binary.LAN*180.0/np.pi,'m',inner_binary.mass)
        print( 'outer','a',outer_binary.a,'e',outer_binary.e,'INCL',outer_binary.INCL*180.0/np.pi,'AP',outer_binary.AP*180.0/np.pi,'LAN',outer_binary.LAN*180.0/np.pi,'m',outer_binary.mass)
#        print 'r',particles.position_x,particles.position_y,particles.position_z
#        print 'v',particles.velocity_x.value_in(units.km/units.s),particles.velocity_y.value_in(units.km/units.s),particles.velocity_z.value_in(units.km/units.s)
        
#        c1 = particles[0]
#        c2 = particles[1]
        
#        r1_vec = [c1.position_x - c2.position_x,c1.position_y - c2.position_y,c1.position_z - c2.position_z]
#        v1_vec = [c1.velocity_x - c2.velocity_x,c1.velocity_y - c2.velocity_y,c1.velocity_z - c2.velocity_z]

#        code.parameters.orbital_phases_random_seed = 1
        
        code.apply_user_specified_instantaneous_perturbation()
        
        #channel_from_code_external.copy()

        print( '='*50)
        print( 'post')
        print( 'inner','a',inner_binary.a,'e',inner_binary.e,'INCL',inner_binary.INCL*180.0/np.pi,'AP',inner_binary.AP*180.0/np.pi,'LAN',inner_binary.LAN*180.0/np.pi,'m',inner_binary.mass)
        print( 'outer','a',outer_binary.a,'e',outer_binary.e,'INCL',outer_binary.INCL*180.0/np.pi,'AP',outer_binary.AP*180.0/np.pi,'LAN',outer_binary.LAN*180.0/np.pi,'m',outer_binary.mass)

      
        exit(-1)
        print( 'r',particles.position_x,particles.position_y,particles.position_z)
        print( 'v',particles.velocity_x.value_in(units.km/units.s),particles.velocity_y.value_in(units.km/units.s),particles.velocity_z.value_in(units.km/units.s))

        #code.apply_external_perturbation_assuming_integrated_orbits()

        
            
        t += 1.0 | units.Myr ### arbitrary
        t_print_array.append(t)
        a_print_array.append(binaries[0].semimajor_axis)
        e_print_array.append(binaries[0].eccentricity | units.none)
        INCL_print_array.append(binaries[0].inclination | units.none)            
        AP_print_array.append(binaries[0].argument_of_pericenter | units.none)            
        LAN_print_array.append(binaries[0].longitude_of_ascending_node | units.none)            

        print( 'a_print_array',a_print_array.value_in(units.AU))
        print( 'e_print_array',e_print_array)
        print( 'INCL_print_array',INCL_print_array*(180.0/np.pi))
        print( 'AP_print_array',AP_print_array*(180.0/np.pi))
        print( 'LAN_print_array',LAN_print_array*(180.0/np.pi))
        
        f1 = binaries.true_anomaly
        print( 'true_anomaly',f1)
        
        m1 = masses[0]
        m2 = masses[1]
        a1 = semimajor_axes[0]
        e1 = eccentricities[0]
        
        r1 = a1*(1.0-e1**2)/(1.0 + e1*np.cos(f1))
        

        r1_dot_v1 = np.sum([x*y for x,y in zip(r1_vec,v1_vec)])
        r1_dot_V_k = np.sum([x*y for x,y in zip(r1_vec,V_k_vec)])
        v1_dot_V_k = np.sum([x*y for x,y in zip(v1_vec,V_k_vec)])
        V_k_dot_V_k = np.sum([x*y for x,y in zip(V_k_vec,V_k_vec)])
        v1c_sq =constants.G*(m1+m2)/a1
        
        
        a1_p = a1*(1.0 + delta_m1/(m1+m2))*pow( 1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq, -1.0)
        j1_p = pow( (m1+m2)/(m1+m2+delta_m1), 2.0)*(1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq)*( 1.0 - e1**2 \
            + (1.0/(constants.G*(m1+m2)*a1))*( r1**2*( 2.0*v1_dot_V_k + V_k_dot_V_k) - 2.0*r1_dot_v1*r1_dot_V_k - r1_dot_V_k**2) )
        e1_p = np.sqrt(1.0 - j1_p)
        
        print( 'analytic results Toonen+16',a1_p.value_in(units.AU),e1_p)
        
        if HAS_MATPLOTLIB == True and 1==0:
            fig = pyplot.figure(figsize=(8,6))
            plot1 = fig.add_subplot(2,2,1)
            plot2 = fig.add_subplot(2,2,2)
            plot3 = fig.add_subplot(2,2,3)
            plot4 = fig.add_subplot(2,2,4)
            
            plot1.plot(t_print_array.value_in(units.Myr),e_print_array.value_in(units.none), color='k')
            plot2.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*INCL_print_array.value_in(units.none), color='k')
            plot3.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*AP_print_array.value_in(units.none), color='k')
            plot4.plot(t_print_array.value_in(units.Myr),(180.0/np.pi)*LAN_print_array.value_in(units.none), color='k')
            
            fontsize = 15
            plot1.set_ylabel("$e$",fontsize=fontsize)
            plot2.set_ylabel("$i$",fontsize=fontsize)
            plot3.set_ylabel("$\omega$",fontsize=fontsize)
            plot4.set_ylabel("$\Omega$",fontsize=fontsize)
            
#            fontsize = 15
#            plot1.set_ylabel("$e$",fontsize=fontsize)
            
            pyplot.show()

    def development_test(self):
        """
        Tests related to development -- not for production tests
        """
        SM = MSE()
        
        #primary = Particle(is_binary=False, mass=1.0)
        #secondary = Particle(is_binary=False, mass=0.5)
        #binary = Particle(is_binary=True, child1=primary,child2=secondary,a=1.0,e=0.0,INCL=0.0,AP=0.0,LAN=0.0)

    #    primary = Particle(is_binary=False, mass=1.0)
    #    secondary = Particle(is_binary=False, mass=1.0e-3)
    #    tertiary = Particle(is_binary=False, mass=40.0e-3)
    #    inner_binary = Particle(is_binary=True, child1=primary,child2=secondary,a=6.0,e=0.001,INCL=0.0001,AP=45.0*np.pi/180.0,LAN=0.01)
    #    outer_binary = Particle(is_binary=True, child1=inner_binary,child2=tertiary,a=100.0,e=0.6,INCL=65.0*np.pi/180.0,AP=0.01,LAN=0.01)
        
        #particles = create_nested_multiple(3,[1.0|units.MSun, 1.0|units.MJupiter, 40.0|units.MJupiter], [6.0|units.AU,100.0|units.AU], [0.001,0.6], [0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.0001],[0.01,0.01])
        #print primary.is_binary
    #    print 'p',primary
    #    particles = [primary,secondary,tertiary,inner_binary,outer_binary]

        code = MSE()
        code.enable_root_finding = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        
        a = 10.0
        e = 0.1
        masses = [1.0, 0.8]
        m = masses[0] + masses[1]
        M = 1.0
        E = 10.0
        Q = 100.0
        INCL = 0.5*np.pi
        AP = 0.4*np.pi
        LAN = 0.25*np.pi
        particles = Tools.create_nested_multiple(2,masses, [a], [e], [INCL], [AP], [LAN])
        
        t = 0.0
        #t_ref = 1.0e-2 | units.Myr
        t_ref = 1.0e5
        tend = 2.0*t_ref
        Nsteps = 100
        dt = tend/float(Nsteps)

        import np.random as randomf
        randomf.seed(0)
        N=1
        external_particle = Particle(mass = M, is_binary=True, is_external=True, external_t_ref=t_ref, e=E, external_r_p = Q, INCL = 1.0e-10, AP = 1.0e-10, LAN = 1.0e-10)

        particles.append(external_particle)

        code = MSE()
        code.add_particles(particles)
        binary = code.particles[2]
        
        #code.parameters.include_quadrupole_order_terms = True
        #code.parameters.include_octupole_order_binary_pair_terms = True
        #code.parameters.include_hexadecupole_order_binary_pair_terms = True
        #code.parameters.include_dotriacontupole_order_binary_pair_terms = True

        t_print_array = []
        a_print_array = []
        e_print_array = []
        INCL_print_array = []
        AP_print_array = []
        LAN_print_array = []
        
        if 1==0:
            N = 0
            while (t<tend):

                t+=dt
                N+=1

                code.evolve_model(t)
                 
                
                t_print_array.append(t)
                a_print_array.append(binary.a)
                e_print_array.append(binary.e)
        #            INCL_print_array.append(binaries[0].inclination | units.none)            
        #            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)            
        #            LAN_print_array.append(binaries[0].longitude_of_ascending_node | units.none)            
                
                print( 't',t,'e',binary.e)
                                

            Delta_e_num = e_print_array[-1] - e
        else:
            code.apply_external_perturbation_assuming_integrated_orbits()
            Delta_e_num = binary.e - e
            
        eps_SA = np.sqrt( (M**2/(m*(m+M)))*(a/Q)**3*1.0/((1.0 + E)**3))
        i=INCL
        Delta_e_an = (-5*np.sqrt(e**2 - e**4)*eps_SA*(-3*E*np.arccos(-1.0/E)*np.cos(AP)*np.sin(i)**2*np.sin(AP) + np.sqrt(1 - E**(-2))*(-((-1 + E**2)*np.cos(i)**2*np.cos(2*LAN)*np.sin(2*AP)) - ((3*E**2 + (-1 + E**2)*np.cos(2*LAN))*np.sin(i)**2*np.sin(2*AP))/2. + ((-1 + E**2)*np.cos(i)*(-2 - 2*np.cos(2*i) + np.cos(2*(i - AP)) - 6*np.cos(2*AP) + np.cos(2*(i + AP)))*np.sin(2*LAN))/8. + (-1 + E**2)*np.cos(i)**3*np.sin(AP)**2*np.sin(2*LAN))))/(2.*E)
        
        print( 'Delta e num',Delta_e_num,'Delta_e_an',Delta_e_an)
        
        exit(0)


        particles = Tools.create_nested_multiple(3, [1.0,1.0e-3,40.0e-3],[6.0,100.0],[0.001,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(3, [1.0,1.0e-3,1.0e-3,1.0e-3],[1.0,6.0,40.0],[0.001,0.2,0.4],[0.0001,50.0*np.pi/180.0,55.0*np.pi/180.0],[45.0*np.pi/180.0,0.01,0.01*np.pi/180.0],[0.01,0.01,0.01])
        binaries = [x for x in particles if x.is_binary==True]
        #binaries[0].include_1PN_terms = True
        #print 'b',particles

        SM.add_particles(particles)
        primary = SM.particles[0]
        inner_binary = SM.particles[3]
        outer_binary = SM.particles[4]

        #SM.enable_tides = True
        #SM.enable_root_finding = True
        
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        
        t = 0.0
        N = 1000
        tend = 3.0e7
        dt = tend/float(N)
        import time
        start = time.time()
        while t<=tend:
            print( 't',t)
            SM.evolve_model(t)
            t+=dt
        
            print( 't',t,'e',inner_binary.e,'INCL',inner_binary.INCL,outer_binary.INCL,primary.spin_vec_x)
            
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            rel_INCL_print.append(inner_binary.INCL_parent)
            e_print.append(inner_binary.e)
            INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        print( 'wall time',time.time()-start)
        t_print = np.array(t_print)
        rel_INCL_print = np.array(rel_INCL_print)
        e_print = np.array(e_print)
        
        from matplotlib import pyplot
        fig=pyplot.figure()
        plot1=fig.add_subplot(2,1,1)
        plot2=fig.add_subplot(2,1,2,yscale="log")
        #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
        plot1.plot(1.0e-6*t_print,rel_INCL_print*180.0/np.pi)
        plot2.plot(1.0e-6*t_print,1.0-e_print)
        pyplot.show()
        
        
        #print 'p',primary

    #    print 'test 1',SM.particles

    #    SM.delete_particle(inner_binary)
    #    print 'test 2',SM.particles

        #print 'C1',SM.CONST_M_SUN
        #SM.CONST_G = 2.0
        #print 'C2',SM.CONST_C
    

    def test14(self):
        print('test14')
        
        N_bodies=4
        
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5],[30.0,400.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [4.0,2.8,1.5],[3000.0,400000.0],[0.1,0.3],[0.0001,89.9*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        particles = Tools.create_nested_multiple(N_bodies, [6.0,4.6,3.5,2.5],[10.0,1000.0,12000.0],[0.1,0.3,0.3],[0.0001,51.0*np.pi/180.0,123.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0,1.0],[10.0,200.0],[0.1,0.6],[0.0001,89.8*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0],[1.0],[0.1],[0.0001,],[45.0*np.pi/180.0,],[0.01])
        
        #particles = [Particle(is_binary=False,mass=8.5,metallicity=0.02)]
        #particles = [Particle(is_binary=False,mass=10.0),Particle(is_binary=False,mass=10.1)]
        orbits = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]

        N_orbits = len(orbits)
        
       # for orbit in orbits:
       #     orbit.check_for_physical_collision_or_orbit_crossing = True

       #     orbit.include_tidal_friction_terms = False

        for b in bodies:
            b.kick_distribution_sigma = 0.0
       #     b.metallicity = 0.02
       #     b.include_tidal_friction_terms = True

        
#        masses=[10.0,9.5,8.5]
#        N_bodies=len(masses)
#        for index in range(N_bodies):
#            particle = Particle(is_binary=False,mass=masses[index])
#            particle.evolve_as_star = True
#            particles.append(particle)

        
        
        code = MSE()
        code.enable_root_finding = True
        code.add_particles(particles)
        code.orbital_phases_random_seed = 1
        
        code.include_flybys = True
        #code.flybys_reference_binary = orbits[0].index
        #print 'a',orbits[0].a
        code.flybys_stellar_density = 0.1*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e5


        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        k_print = [[] for x in range(N_bodies)]
        m_print = [[] for x in range(N_bodies)]
        R_print = [[] for x in range(N_bodies)]
        Rc_print = [[] for x in range(N_bodies)]
        R_L_print = [[] for x in range(N_bodies)]
        t_V_print = [[] for x in range(N_bodies)]
        mc_print = [[] for x in range(N_bodies)]
        a_print = [[] for x in range(N_orbits)]
        e_print = [[] for x in range(N_orbits)]
        rel_INCL_print = [[] for x in range(N_orbits)]
        
        t = 0.0

        if 1==0:
            for index in range(N_bodies):
                m_print[index].append(particles[index].mass)
                k_print[index].append(particles[index].stellar_type)
                R_print[index].append(particles[index].radius)
                t_V_print[index].append(particles[index].tides_viscous_time_scale)
                Rc_print[index].append(particles[index].convective_envelope_radius)
                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        N = 10
        tend = 1e10

        N =10000
        tend = 1e9
        #tend = 1.4e10

        seed=1
        dt = tend/float(N)
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state
        
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag)

            for index in range(N_orbits):
                rel_INCL_print[index].append(orbits[index].INCL_parent)
                e_print[index].append(orbits[index].e)
                a_print[index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[index].append(particles[index].mass)
                k_print[index].append(particles[index].stellar_type)
                R_print[index].append(particles[index].radius)
                t_V_print[index].append(particles[index].tides_viscous_time_scale)
                Rc_print[index].append(particles[index].convective_envelope_radius)
                R_L_print[index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
            t_print.append(t)        
            if flag==2:
                print("Root",t,[o.RLOF_at_pericentre_has_occurred for o in orbits],[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1

        t_print = np.array(t_print)
        for index in range(N_bodies):
                m_print[index] = np.array(m_print[index])
                t_V_print[index] = np.array(t_V_print[index])
        for index in range(N_orbits):

            rel_INCL_print[index] = np.array(rel_INCL_print[index])
            e_print[index] = np.array(e_print[index])
            a_print[index] = np.array(a_print[index])

        #rel_INCL_print = np.array(rel_INCL_print)
        #e_print = np.array(e_print)
        
        print("Final properties -- ","masses/MSun",[m_print[i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            linewidth=1.0
            for index in range(N_bodies):
                plot1.plot(1.0e-6*t_print,m_print[index],color='k',linewidth=linewidth)
                plot1.plot(1.0e-6*t_print,mc_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                plot3.plot(1.0e-6*t_print,k_print[index],color='k',linestyle='solid',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,R_print[index],color='r',linestyle='solid',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,Rc_print[index],color='r',linestyle='dotted',linewidth=linewidth)
                #if index in [0,1]:
                    #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                linewidth+=0.8
            linewidth=1.0
            for index in range(N_orbits):
                plot2.plot(1.0e-6*t_print,a_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,a_print[index]*(1.0-e_print[index]),color='k',linestyle='solid',linewidth=linewidth)
                #plot3.plot(1.0e-6*t_print,(180.0/np.pi)*rel_INCL_print[index],color='k',linestyle='solid',linewidth=linewidth)
                linewidth+=0.8

            #plot3.plot(1.0e-6*t_print,a_print[0]*(m_print[0]+m_print[1]),color='k')
            #plot3.plot(1.0e-6*t_print,a_print[1]*(m_print[0]+m_print[1]+m_print[2]),color='g')
            #plot3.plot(1.0e-6*t_print,a_print[2]*(m_print[0]+m_print[1]+m_print[2]+m_print[3]),color='b')
            fontsize=18
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            plot3.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            fig.savefig("test14.pdf")
            pyplot.show()
        


    def test15(self):
        print('test15')
        
        #ms = np.linspace(1.0,20.0,20)
        ms = np.linspace(1.0,20.0,1)
        
        for index,m in enumerate(ms):
            print( 'index',index,'m',m)
            N_bodies=4
            
            #particles = Tools.create_nested_multiple(N_bodies, [m*5,m*4,m*3],[20.0,400.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
            
            particles = Tools.create_nested_multiple(N_bodies, [m*5,m*4,m*3,m*2],[35.0,1000.0,10000.0],[0.1,0.6,0.3],[0.0001,65.0*np.pi/180.0,43.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
            #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0,1.0],[10.0,200.0],[0.1,0.6],[0.0001,89.8*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
            #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0],[1.0],[0.1],[0.0001,],[45.0*np.pi/180.0,],[0.01])
            
            #particles = [Particle(is_binary=False,mass=m,metallicity=0.02)]
            #particles = [Particle(is_binary=False,mass=10.0),Particle(is_binary=False,mass=10.1)]
            orbits = [x for x in particles if x.is_binary==True]
            bodies = [x for x in particles if x.is_binary==False]

            N_orbits = len(orbits)
            
            for orbit in orbits:
                orbit.check_for_physical_collision_or_orbit_crossing = True

            for b in bodies:
                b.metallicity = 0.02

            
    #        masses=[10.0,9.5,8.5]
    #        N_bodies=len(masses)
    #        for index in range(N_bodies):
    #            particle = Particle(is_binary=False,mass=masses[index])
    #            particle.evolve_as_star = True
    #            particles.append(particle)

            
            
            code = MSE()
            code.enable_root_finding = True
            code.add_particles(particles)
            code.orbital_phases_random_seed = 1
            
            e_print = []
            INCL_print = []
            rel_INCL_print = []
            t_print = []
            k_print = [[] for x in range(N_bodies)]
            m_print = [[] for x in range(N_bodies)]
            R_print = [[] for x in range(N_bodies)]
            Rc_print = [[] for x in range(N_bodies)]
            mc_print = [[] for x in range(N_bodies)]
            a_print = [[] for x in range(N_orbits)]
            e_print = [[] for x in range(N_orbits)]
            rel_INCL_print = [[] for x in range(N_orbits)]
            t = 0.0

            if 1==0:
                for index in range(N_bodies):
                    m_print[index].append(particles[index].mass)
                    k_print[index].append(particles[index].stellar_type)
                    R_print[index].append(particles[index].radius)
                    Rc_print[index].append(particles[index].convective_envelope_radius)
                    mc_print[index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
                t_print.append(t)
            N = 10
            tend = 1e10

            N =5000
            tend = 1e9
            #tend = 1.4e10

            seed=1
            dt = tend/float(N)
            while t<tend:

                t+=dt
                code.evolve_model(t)
                
                flag = code.CVODE_flag
                state = code.state
            
                #print( 't',t,orbits[0].e,flag)

                for index in range(N_orbits):
                    rel_INCL_print[index].append(orbits[index].INCL_parent)
                    e_print[index].append(orbits[index].e)
                    a_print[index].append(orbits[index].a)
                for index in range(N_bodies):
                    m_print[index].append(particles[index].mass)
                    k_print[index].append(particles[index].stellar_type)
                    R_print[index].append(particles[index].radius)
                    Rc_print[index].append(particles[index].convective_envelope_radius)
                    mc_print[index].append(particles[index].convective_envelope_mass)
                #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
                #rel_INCL_print.append(inner_binary.INCL_parent)
                #e_print.append(inner_binary.e)
                #INCL_print.append(inner_binary.INCL)
                t_print.append(t)        
                if flag==2:
                    print("Root",t,[o.RLOF_at_pericentre_has_occurred for o in orbits],[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits])
                    break
                if flag==3:
                    print("SNe")
                    break
                if state==3:
                    print("unbound")
                    break
                code.random_seed = seed
                seed += 1

            code.reset()
            t_print = np.array(t_print)
            for index in range(N_bodies):
                    m_print[index] = np.array(m_print[index])
            for index in range(N_orbits):

                e_print[index] = np.array(e_print[index])
                a_print[index] = np.array(a_print[index])

            #rel_INCL_print = np.array(rel_INCL_print)
            #e_print = np.array(e_print)
            
            print("Final properties -- ","masses/MSun",[m_print[i][-1] for i in range(N_bodies)])
            
            if HAS_MATPLOTLIB==True:
                fig=pyplot.figure(figsize=(8,8))
                plot1=fig.add_subplot(2,1,1)
                plot2=fig.add_subplot(2,1,2,yscale="log")
                #plot3=fig.add_subplot(2,1,3,yscale="log")
                #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
                linewidth=1.0
                for index in range(N_bodies):
                    plot1.plot(1.0e-6*t_print,m_print[index],color='k')
                    plot1.plot(1.0e-6*t_print,mc_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                    #plot1.plot(1.0e-6*t_print,k_print[index],color='r',linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print,R_print[index],color='r',linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print,Rc_print[index],color='r',linestyle='dotted',linewidth=linewidth)
                    linewidth+=1.0
                linewidth=1.0
                for index in range(N_orbits):
                    plot2.plot(1.0e-6*t_print,a_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print,a_print[index]*(1.0-e_print[index]),color='k',linestyle='solid',linewidth=linewidth)
                    linewidth+=1.0

                #plot3.plot(1.0e-6*t_print,a_print[0]*(m_print[0]+m_print[1]),color='k')
                #plot3.plot(1.0e-6*t_print,a_print[1]*(m_print[0]+m_print[1]+m_print[2]),color='g')
                #plot3.plot(1.0e-6*t_print,a_print[2]*(m_print[0]+m_print[1]+m_print[2]+m_print[3]),color='b')
                fontsize=18
                plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
                plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
                plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
                plot2.set_ylim(1.0e-2,1.0e4)
                #pyplot.show()
            
        pyplot.show()


    def test16(self):
        print('test16')
        
        N_bodies=3
        
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5],[30.0,400.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        particles = Tools.create_nested_multiple(N_bodies, [20.0,10.0,4.0e6],[10.0,1.0e4],[0.1,0.3],[0.0001,80.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5,8.5],[35.0,400.0,40000.0],[0.1,0.6,0.3],[0.0001,65.0*np.pi/180.0,43.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0,1.0],[10.0,200.0],[0.1,0.6],[0.0001,89.8*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0],[1.0],[0.1],[0.0001,],[45.0*np.pi/180.0,],[0.01])
        
        #particles = [Particle(is_binary=False,mass=8.5,metallicity=0.02)]
        #particles = [Particle(is_binary=False,mass=10.0),Particle(is_binary=False,mass=10.1)]
        orbits = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]

        N_orbits = len(orbits)
        
       # for orbit in orbits:
       #     orbit.check_for_physical_collision_or_orbit_crossing = True

       #     orbit.include_tidal_friction_terms = False


       #     b.metallicity = 0.02
       #     b.include_tidal_friction_terms = True

        
#        masses=[10.0,9.5,8.5]
#        N_bodies=len(masses)
#        for index in range(N_bodies):
#            particle = Particle(is_binary=False,mass=masses[index])
#            particle.evolve_as_star = True
#            particles.append(particle)

        
        
        code = MSE()
        CONST_G=code.CONST_G
        CONST_C=code.CONST_C
        CONST_KM_PER_S=code.CONST_KM_PER_S
        code.enable_root_finding = True
        code.add_particles(particles)
        code.orbital_phases_random_seed = 1
        
        for b in bodies:
            b.evolve_as_star = False
            b.radius = CONST_G*b.mass/(CONST_C**2)
            
        code.include_flybys = True
        code.flybys_reference_binary = orbits[0].index
        #print 'a',orbits[0].a
        code.flybys_stellar_density = 1e5*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e2
        code.flybys_stellar_relative_velocity_dispersion = 1000.0*CONST_KM_PER_S

        
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        k_print = [[] for x in range(N_bodies)]
        m_print = [[] for x in range(N_bodies)]
        R_print = [[] for x in range(N_bodies)]
        Rc_print = [[] for x in range(N_bodies)]
        t_V_print = [[] for x in range(N_bodies)]
        mc_print = [[] for x in range(N_bodies)]
        a_print = [[] for x in range(N_orbits)]
        e_print = [[] for x in range(N_orbits)]
        rel_INCL_print = [[] for x in range(N_orbits)]
        
        t = 0.0

        if 1==0:
            for index in range(N_bodies):
                m_print[index].append(particles[index].mass)
                k_print[index].append(particles[index].stellar_type)
                R_print[index].append(particles[index].radius)
                t_V_print[index].append(particles[index].tides_viscous_time_scale)
                Rc_print[index].append(particles[index].convective_envelope_radius)
                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        N = 10
        tend = 1e10

        N =1000
        tend = 1e8
        #tend = 1.4e10

        seed=1
        dt = tend/float(N)
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state
        
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag)

            for index in range(N_orbits):
                rel_INCL_print[index].append(orbits[index].INCL_parent)
                e_print[index].append(orbits[index].e)
                a_print[index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[index].append(particles[index].mass)
                k_print[index].append(particles[index].stellar_type)
                R_print[index].append(particles[index].radius)
                t_V_print[index].append(particles[index].tides_viscous_time_scale)
                Rc_print[index].append(particles[index].convective_envelope_radius)
                mc_print[index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
            t_print.append(t)        
            if flag==2:
                print("Root",t,[o.RLOF_at_pericentre_has_occurred for o in orbits],[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1

        t_print = np.array(t_print)
        for index in range(N_bodies):
                m_print[index] = np.array(m_print[index])
                t_V_print[index] = np.array(t_V_print[index])
        for index in range(N_orbits):

            rel_INCL_print[index] = np.array(rel_INCL_print[index])
            e_print[index] = np.array(e_print[index])
            a_print[index] = np.array(a_print[index])

        #rel_INCL_print = np.array(rel_INCL_print)
        #e_print = np.array(e_print)
        
        print("Final properties -- ","masses/MSun",[m_print[i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True:
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            linewidth=1.0
            for index in range(N_bodies):
                plot1.plot(1.0e-6*t_print,m_print[index],color='k')
                plot1.plot(1.0e-6*t_print,mc_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                #plot1.plot(1.0e-6*t_print,k_print[index],color='r',linestyle='dotted',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,R_print[index],color='r',linestyle='solid',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,Rc_print[index],color='r',linestyle='dotted',linewidth=linewidth)
                #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                
                linewidth+=1.0
            linewidth=1.0
            for index in range(N_orbits):
                plot2.plot(1.0e-6*t_print,a_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,a_print[index]*(1.0-e_print[index]),color='k',linestyle='solid',linewidth=linewidth)
                plot3.plot(1.0e-6*t_print,(180.0/np.pi)*rel_INCL_print[index],color='k',linestyle='solid',linewidth=linewidth)
                linewidth+=1.0

            #plot3.plot(1.0e-6*t_print,a_print[0]*(m_print[0]+m_print[1]),color='k')
            #plot3.plot(1.0e-6*t_print,a_print[1]*(m_print[0]+m_print[1]+m_print[2]),color='g')
            #plot3.plot(1.0e-6*t_print,a_print[2]*(m_print[0]+m_print[1]+m_print[2]+m_print[3]),color='b')
            fontsize=18
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
            plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-2,1.0e4)
            pyplot.show()

    def test17(self):
        print('test17')
        """
        test mass transfer
        """
        
        N_bodies=3
        
        particles = Tools.create_nested_multiple(N_bodies, [10.0,5.0,10.0],[5.0,1.0e2],[0.1,0.3],[0.0001,59.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01]) 
        
        #particles = Tools.create_nested_multiple(N_bodies, [10.0,5.0,10.0],[5.0,1.0e2],[0.1,0.5],[0.0001,79.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01]) ###!!!!!
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5,8.5],[35.0,400.0,40000.0],[0.1,0.6,0.3],[0.0001,65.0*np.pi/180.0,43.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0,1.0],[10.0,200.0],[0.1,0.6],[0.0001,89.8*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0],[1.0],[0.1],[0.0001,],[45.0*np.pi/180.0,],[0.01])
        
        #particles = [Particle(is_binary=False,mass=8.5,metallicity=0.02)]
        #particles = [Particle(is_binary=False,mass=10.0),Particle(is_binary=False,mass=10.1)]
        orbits = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]

        N_orbits = len(orbits)
        
       # for orbit in orbits:
       #     orbit.check_for_physical_collision_or_orbit_crossing = True

       #     orbit.include_tidal_friction_terms = False


       #     b.metallicity = 0.02
       #     b.include_tidal_friction_terms = True

        
#        masses=[10.0,9.5,8.5]
#        N_bodies=len(masses)
#        for index in range(N_bodies):
#            particle = Particle(is_binary=False,mass=masses[index])
#            particle.evolve_as_star = True
#            particles.append(particle)

        
        
        code = MSE()
        CONST_G=code.CONST_G
        CONST_C=code.CONST_C
        CONST_KM_PER_S=code.CONST_KM_PER_S
        code.enable_root_finding = True
        code.add_particles(particles)
        code.orbital_phases_random_seed = 1
        
#        for b in bodies:
#            b.evolve_as_star = False
#            b.radius = CONST_G*b.mass/(CONST_C**2)
            
        code.include_flybys = True
        #code.flybys_reference_binary = orbits[0].index
        #print 'a',orbits[0].a
        code.flybys_stellar_density = 0.1*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e5
        code.flybys_stellar_relative_velocity_dispersion = 40.0*CONST_KM_PER_S

        
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        k_print = [[] for x in range(N_bodies)]
        m_print = [[] for x in range(N_bodies)]
        R_print = [[] for x in range(N_bodies)]
        Rc_print = [[] for x in range(N_bodies)]
        R_L_print = [[] for x in range(N_bodies)]
        t_V_print = [[] for x in range(N_bodies)]
        mc_print = [[] for x in range(N_bodies)]
        a_print = [[] for x in range(N_orbits)]
        e_print = [[] for x in range(N_orbits)]
        rel_INCL_print = [[] for x in range(N_orbits)]
        
        t = 0.0

        if 1==0:
            for index in range(N_bodies):
                m_print[index].append(particles[index].mass)
                k_print[index].append(particles[index].stellar_type)
                R_print[index].append(particles[index].radius)
                t_V_print[index].append(particles[index].tides_viscous_time_scale)
                Rc_print[index].append(particles[index].convective_envelope_radius)
                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        N = 10
        tend = 1e10

        N =10000
        tend = 1e8
        #tend = 1.4e10

        seed=1
        dt = tend/float(N)
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state
        
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag)

            for index in range(N_orbits):
                rel_INCL_print[index].append(orbits[index].INCL_parent)
                e_print[index].append(orbits[index].e)
                a_print[index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[index].append(particles[index].mass)
                k_print[index].append(particles[index].stellar_type)
                R_print[index].append(particles[index].radius)
                t_V_print[index].append(particles[index].tides_viscous_time_scale)
                Rc_print[index].append(particles[index].convective_envelope_radius)
                R_L_print[index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
            t_print.append(t)        
            if flag==2:
                print("Root",t,[o.RLOF_at_pericentre_has_occurred for o in bodies],[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1

        t_print = np.array(t_print)
        for index in range(N_bodies):
                m_print[index] = np.array(m_print[index])
                t_V_print[index] = np.array(t_V_print[index])
        for index in range(N_orbits):

            rel_INCL_print[index] = np.array(rel_INCL_print[index])
            e_print[index] = np.array(e_print[index])
            a_print[index] = np.array(a_print[index])

        #rel_INCL_print = np.array(rel_INCL_print)
        #e_print = np.array(e_print)
        
        print("Final properties -- ","masses/MSun",[m_print[i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True:
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            linewidth=0.5
            for index in range(N_bodies):
                plot1.plot(1.0e-6*t_print,m_print[index],color='k')
                plot1.plot(1.0e-6*t_print,mc_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                plot3.plot(1.0e-6*t_print,k_print[index],color='k',linestyle='solid',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,R_print[index],color='r',linestyle='solid',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,Rc_print[index],color='r',linestyle='dotted',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                
                linewidth+=1.5
            linewidth=0.5
            for index in range(N_orbits):
                plot2.plot(1.0e-6*t_print,a_print[index],color='k',linestyle='dotted',linewidth=linewidth)
                plot2.plot(1.0e-6*t_print,a_print[index]*(1.0-e_print[index]),color='k',linestyle='solid',linewidth=linewidth)
                #plot3.plot(1.0e-6*t_print,(180.0/np.pi)*rel_INCL_print[index],color='k',linestyle='solid',linewidth=linewidth)
                #plot3.plot(1.0e-6*t_print,e_print[index],color='k',linestyle='solid',linewidth=linewidth)
                linewidth+=1.5

            #plot3.plot(1.0e-6*t_print,a_print[0]*(m_print[0]+m_print[1]),color='k')
            #plot3.plot(1.0e-6*t_print,a_print[1]*(m_print[0]+m_print[1]+m_print[2]),color='g')
            #plot3.plot(1.0e-6*t_print,a_print[2]*(m_print[0]+m_print[1]+m_print[2]+m_print[3]),color='b')
            fontsize=18
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
            plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-2,1.0e4)
            fig.savefig("test.pdf")
            pyplot.show()



    def test18(self):
        """
        Test running on parallel environment using mpi4py
        """

        print('test18')



        #ms = np.linspace(1.0,20.0,20)
        ms = np.linspace(1.0,5.0,1000)
        system_particles = []
        
        for index,m in enumerate(ms):
            N_bodies=4
            particles = Tools.create_nested_multiple(N_bodies, [m*5,m*4,m*3,m*2],[35.0,1000.0,10000.0],[0.1,0.6,0.3],[0.0001,65.0*np.pi/180.0,43.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
            system_particles.append(particles)

        import math
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        N_threads = comm.Get_size()
        N_jobs = len(system_particles)
        N_jobs_per_thread = int(math.ceil( float(N_jobs)/float(N_threads) ))
        py_name = "test_mse.py"
        print(py_name,'-- rank',rank,'; N_threads = ',N_threads,'; N_jobs_per_thread',N_jobs_per_thread)

        indices = range(rank,N_jobs,N_threads)
#        print(py_name,'--rank',rank,'-- will run systems with indices ',[x.index for x in [simulation_parameters_particles[i] for i in indices]])

        for index_rank,index_sys in enumerate(indices):
            particles = system_particles[index_sys]

            run_system(index_sys,particles)

            #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0,1.0],[10.0,200.0],[0.1,0.6],[0.0001,89.8*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
            #particles = Tools.create_nested_multiple(N_bodies, [1.0,1.0],[1.0],[0.1],[0.0001,],[45.0*np.pi/180.0,],[0.01])
            
            #particles = [Particle(is_binary=False,mass=m,metallicity=0.02)]
            #particles = [Particle(is_binary=False,mass=10.0),Particle(is_binary=False,mass=10.1)]

def run_system(index_system,particles):
    
    orbits = [x for x in particles if x.is_binary==True]
    bodies = [x for x in particles if x.is_binary==False]

    print("mass",bodies[0].mass)
    N_orbits = len(orbits)
    N_bodies = len(bodies)
    
    for orbit in orbits:
        orbit.check_for_physical_collision_or_orbit_crossing = True

    for b in bodies:
        b.metallicity = 0.02
    
#        masses=[10.0,9.5,8.5]
#        N_bodies=len(masses)
#        for index in range(N_bodies):
#            particle = Particle(is_binary=False,mass=masses[index])
#            particle.evolve_as_star = True
#            particles.append(particle)
    
    code = MSE()
    code.enable_root_finding = True
    code.add_particles(particles)
    code.orbital_phases_random_seed = 1
    
    e_print = []
    INCL_print = []
    rel_INCL_print = []
    t_print = []
    k_print = [[] for x in range(N_bodies)]
    m_print = [[] for x in range(N_bodies)]
    R_print = [[] for x in range(N_bodies)]
    Rc_print = [[] for x in range(N_bodies)]
    mc_print = [[] for x in range(N_bodies)]
    a_print = [[] for x in range(N_orbits)]
    e_print = [[] for x in range(N_orbits)]
    rel_INCL_print = [[] for x in range(N_orbits)]
    t = 0.0

    if 1==0:
        for index in range(N_bodies):
            m_print[index].append(particles[index].mass)
            k_print[index].append(particles[index].stellar_type)
            R_print[index].append(particles[index].radius)
            Rc_print[index].append(particles[index].convective_envelope_radius)
            mc_print[index].append(particles[index].convective_envelope_mass)
    #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
    #rel_INCL_print.append(inner_binary.INCL_parent)
    #e_print.append(inner_binary.e)
    #INCL_print.append(inner_binary.INCL)
        t_print.append(t)
    N = 10
    tend = 1e10

    N =5000
    tend = 1e8
    #tend = 1.4e10

    seed=1
    dt = tend/float(N)
    while t<tend:

        t+=dt
        code.evolve_model(t)
        
        flag = code.CVODE_flag
        state = code.state
    
        #print( 't',t,orbits[0].e,flag)

        for index in range(N_orbits):
            rel_INCL_print[index].append(orbits[index].INCL_parent)
            e_print[index].append(orbits[index].e)
            a_print[index].append(orbits[index].a)
        for index in range(N_bodies):
            m_print[index].append(particles[index].mass)
            k_print[index].append(particles[index].stellar_type)
            R_print[index].append(particles[index].radius)
            Rc_print[index].append(particles[index].convective_envelope_radius)
            mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
        t_print.append(t)        
        if flag==2:
            print("Root",t,[o.RLOF_at_pericentre_has_occurred for o in orbits],[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits])
            break
        if flag==3:
            print("SNe")
            break
        if state==3:
            print("unbound")
            break
        code.random_seed = seed
        seed += 1

    code.reset()
    t_print = np.array(t_print)
    for index in range(N_bodies):
            m_print[index] = np.array(m_print[index])
    for index in range(N_orbits):

        e_print[index] = np.array(e_print[index])
        a_print[index] = np.array(a_print[index])

    #rel_INCL_print = np.array(rel_INCL_print)
    #e_print = np.array(e_print)
    
    print("Final properties -- ","masses/MSun",[m_print[i][-1] for i in range(N_bodies)])
    
    if HAS_MATPLOTLIB==True:
        fig=pyplot.figure(figsize=(8,8))
        plot1=fig.add_subplot(2,1,1)
        plot2=fig.add_subplot(2,1,2,yscale="log")
        #plot3=fig.add_subplot(2,1,3,yscale="log")
        #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
        linewidth=1.0
        for index in range(N_bodies):
            plot1.plot(1.0e-6*t_print,m_print[index],color='k')
            plot1.plot(1.0e-6*t_print,mc_print[index],color='k',linestyle='dotted',linewidth=linewidth)
            #plot1.plot(1.0e-6*t_print,k_print[index],color='r',linestyle='dotted',linewidth=linewidth)
            plot2.plot(1.0e-6*t_print,R_print[index],color='r',linestyle='solid',linewidth=linewidth)
            plot2.plot(1.0e-6*t_print,Rc_print[index],color='r',linestyle='dotted',linewidth=linewidth)
            linewidth+=1.0
        linewidth=1.0
        for index in range(N_orbits):
            plot2.plot(1.0e-6*t_print,a_print[index],color='k',linestyle='dotted',linewidth=linewidth)
            plot2.plot(1.0e-6*t_print,a_print[index]*(1.0-e_print[index]),color='k',linestyle='solid',linewidth=linewidth)
            linewidth+=1.0

        #plot3.plot(1.0e-6*t_print,a_print[0]*(m_print[0]+m_print[1]),color='k')
        #plot3.plot(1.0e-6*t_print,a_print[1]*(m_print[0]+m_print[1]+m_print[2]),color='g')
        #plot3.plot(1.0e-6*t_print,a_print[2]*(m_print[0]+m_print[1]+m_print[2]+m_print[3]),color='b')
        fontsize=18
        plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
        plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
        plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
        plot2.set_ylim(1.0e-2,1.0e4)
        #pyplot.show()
        
        fig.savefig("figs/system_" + str(index_system) + ".pdf")



if __name__ == '__main__':
    import sys
    t=test_mse()
    if len(sys.argv)>1:
        i = int(sys.argv[1])
        if i<0:
            t.development_test()
        else:
            print( 'Running test',i)
            function = getattr(t, 'test%s'%i)
            function()

