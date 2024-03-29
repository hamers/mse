import numpy as np
import argparse
import time
        
from mse import MSE,Particle,Tools

"""
Several routines for testing the code/installation. 
To run all tests, simply run `python test_mse.py'.
Specific tests can be run with the command line --t i, where i is the
number of the test. Use --verbose for verbose terminal output, and --plot to
make and show plots if applicable (required Matplotlib).

"""

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def add_bool_arg(parser, name, default=False,help=None):
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--' + name, dest=name, action='store_true',help="Enable %s"%help)
    group.add_argument('--no-' + name, dest=name, action='store_false',help="Disable %s"%help)
    parser.set_defaults(**{name:default})

def parse_arguments():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--t",                           type=int,     dest="test",                        default=0,              help="Test number")
    parser.add_argument("--m","--mode",                  type=int,     dest="mode",                        default=0,              help="Mode -- 0: standard; 1: detailed (more extensive tests, but takes more time to run)")
    
    
    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                         default=False,         help="Verbose terminal output")
    add_bool_arg(parser, 'plot',                            default=False,         help="Make plots")
    add_bool_arg(parser, 'fancy_plots',                     default=False,         help="Use LaTeX fonts for plots (slow)")
    
    args = parser.parse_args()

    return args

class test_mse():

    def test1(self,args):
        print("Test secular evolution")
        print("Test 1a: secular equations of motion using reference triple system of Naoz et al. (2009)")

        particles = Tools.create_fully_nested_multiple(3, [1.0,1.0e-3,40.0e-3],[6.0,100.0],[0.001,0.6],[0.0,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.0],[0.0,0.0],metallicities=[0.02,0.02,0.02],stellar_types=[1,1,1],object_types=[2,2,2])
        binaries = [x for x in particles if x.is_binary==True]
        inner_binary = binaries[0]
        outer_binary = binaries[1]
        bodies = [x for x in particles if x.is_binary==False]
        
        for b in bodies:
            b.evolve_as_star = False
            b.include_mass_transfer_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False
        
        code = MSE()
        code.add_particles(particles)

        code.relative_tolerance = 1.0e-14
        code.absolute_tolerance_eccentricity_vectors = 1.0e-14

        code.include_flybys = False
        code.enable_tides = False
        code.enable_root_finding = False
        code.verbose_flag = 0
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        
        start = time.time()
    
        t = 0.0
        N = 1000
        tend = 3.0e7

        dt = tend/float(N)
        while t<=tend:

            code.evolve_model(t)
            t+=dt
            
            particles=code.particles
            binaries = [x for x in particles if x.is_binary==True]
            inner_binary = binaries[0]
            outer_binary = binaries[1]
            bodies = [x for x in particles if x.is_binary==False]
        
            if args.verbose==True:
                print( 't',t,'es',[x.e for x in binaries],'INCL_parent',inner_binary.INCL_parent,[x.mass for x in bodies])

            rel_INCL_print.append(inner_binary.INCL_parent)
            e_print.append(inner_binary.e)
            INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        
        if args.verbose==True:
            print('wall time',time.time()-start)
            print("e_print[-1]",e_print[-1],"rel_INCL_print[-1]",rel_INCL_print[-1])
        
        t_print = np.array(t_print)
        rel_INCL_print = np.array(rel_INCL_print)
        e_print = np.array(e_print)

        assert round(e_print[-1],2) == 0.20
        assert round(rel_INCL_print[-1],2) == 1.22

        print("Test 1a passed")

        code.reset()
                
        if HAS_MATPLOTLIB==True and args.plot==True:
            fig=pyplot.figure()
            plot1=fig.add_subplot(2,1,1)
            plot2=fig.add_subplot(2,1,2,yscale="log")
            plot1.plot(1.0e-6*t_print,rel_INCL_print*180.0/np.pi)
            plot2.plot(1.0e-6*t_print,1.0-e_print)
            plot2.set_xlabel("$t/\mathrm{Myr}$")
            plot1.set_ylabel("$i_\mathrm{rel}/\mathrm{deg}$")
            plot2.set_ylabel("$1-e_\mathrm{in}$")
            pyplot.show()


        print("Test 1b: test of canonical e_max expression (quadruple order; test particle; zero initial inner eccentricity)")
        
        i_rel = 85.0*np.pi/180.0
        particles = Tools.create_fully_nested_multiple(3, [1.0,0.001,1.0],[1.0,100.0],[0.0001,0.0001],[0.0,i_rel],[45.0*np.pi/180.0,0.0],[0.0,0.0],metallicities=[0.02,0.02,0.02],stellar_types=[1,1,1],object_types=[2,2,2])
        binaries = [x for x in particles if x.is_binary==True]
        inner_binary = binaries[0]
        outer_binary = binaries[1]
        bodies = [x for x in particles if x.is_binary==False]
        
        for b in bodies:
            b.evolve_as_star = False
            b.include_mass_transfer_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False
        
        code = MSE()
        code.add_particles(particles)

        code.relative_tolerance = 1.0e-14
        code.absolute_tolerance_eccentricity_vectors = 1.0e-14

        code.include_flybys = False
        code.enable_tides = False
        code.enable_root_finding = False
        code.verbose_flag = 0
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        
        start = time.time()
    
        t = 0.0
        N = 1000
        tend = 2.0e6

        dt = tend/float(N)
        while t<=tend:

            code.evolve_model(t)
            t+=dt
            
            particles=code.particles
            binaries = [x for x in particles if x.is_binary==True]
            inner_binary = binaries[0]
            outer_binary = binaries[1]
            bodies = [x for x in particles if x.is_binary==False]
        
            if args.verbose==True:
                print( 't',t,'es',[x.e for x in binaries],'INCL_parent',inner_binary.INCL_parent,[x.mass for x in bodies])

            rel_INCL_print.append(inner_binary.INCL_parent)
            e_print.append(inner_binary.e)
            INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        
        t_print = np.array(t_print)
        rel_INCL_print = np.array(rel_INCL_print)
        e_print = np.array(e_print)

        e_max_an = np.sqrt(1.0 - (5.0/3.0)*np.cos(i_rel)**2)
        e_max_num = np.amax(e_print)

        if args.verbose==True:
            print('wall time',time.time()-start)
            print("e_max_an",e_max_an,"e_max_num",e_max_num)

        N_r = 3
        assert round(e_max_an,N_r) == round(e_max_num,N_r)

        print("Test 1b passed")

        code.reset()

        if HAS_MATPLOTLIB==True and args.plot==True:
            fig=pyplot.figure()
            plot1=fig.add_subplot(2,1,1)
            plot2=fig.add_subplot(2,1,2,yscale="log")
            plot1.plot(1.0e-6*t_print,rel_INCL_print*180.0/np.pi)
            plot2.plot(1.0e-6*t_print,1.0-e_print)
            plot2.set_xlabel("$t/\mathrm{Myr}$")
            plot1.set_ylabel("$i_\mathrm{rel}/\mathrm{deg}$")
            plot2.set_ylabel("$1-e_\mathrm{in}$")
            pyplot.show()
        

    def test2(self,args):
        print("Test 1PN precession in 2-body system")

        particles = Tools.create_fully_nested_multiple(2,[1.0, 1.0], [1.0], [0.99], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[2,2])
        binaries = [x for x in particles if x.is_binary == True]
        bodies = [x for x in particles if x.is_binary == False]

        for b in bodies:
            b.evolve_as_star = False
            b.include_mass_transfer_terms = False
            b.radius = 1.0e-10

        for b in binaries:
            b.include_pairwise_1PN_terms = True
            b.include_pairwise_25PN_terms = False
            b.exclude_1PN_precession_in_case_of_isolated_binary = False ### By default, 1PN apsidal motion is not calculated for isolated binaries; override for this test
        
        code = MSE()
        code.add_particles(particles)

        code.include_flybys = False
        code.enable_tides = False
        code.enable_root_finding = False
        code.verbose_flag=0
        code.relative_tolerance = 1.0e-14
        code.absolute_tolerance_eccentricity_vectors = 1.0e-14 ### need to set lower than default to get more accurate result and compare to analytic expression
        t = 0.0
        N=1000
        tend = 1.0e7
        dt=tend/float(N)

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            
            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]

            if args.verbose==True:
                print( 't/Myr',t,'AP',binaries[0].AP)

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

        # Theoretical prediction #
        a = binaries[0].a
        e = binaries[0].e
        M = binaries[0].mass
        rg = CONST_G*M/(CONST_C**2)
        P = 2.0*np.pi*np.sqrt(a**3/(CONST_G*M))
        t_1PN = (1.0/3.0)*P*(1.0-e**2)*(a/rg)

        AP = 0.01 +2.0*np.pi*tend/t_1PN
        AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi

        if args.verbose == True:
            print("AP num",AP_print_array[-1], "AP an",AP)
        
        N_r=4
        assert round(AP_print_array[-1],N_r) == round(AP,N_r)
        print("Test passed")

        code.reset()
                
        if HAS_MATPLOTLIB == True and args.plot==True:
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array*1.0e-6,AP_print_array, color='r',label="$\mathrm{MSE}$")
            points = np.linspace(0.0,tend*1.0e-6,N)
            AP = 0.01 +2.0*np.pi*points/(t_1PN*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',label="$\mathrm{Analytic}$")

            plot2.plot(t_print_array*1.0e-6,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array*1.0e-6,np.fabs((a-a_print_array)/a), color='r')
            plot4.plot(t_print_array*1.0e-6,np.fabs((e-e_print_array)/e), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)
            
            handles,labels = plot1.get_legend_handles_labels()
            plot1.legend(handles,labels,loc="upper left",fontsize=0.6*fontsize)

            pyplot.show()

    def test3(self,args):
        print("Test GW emission in 2-body system")

        code = MSE()
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C

        a0 = 1.0
        e0 = 0.999
        m1 = 1.0
        m2 = 1.0
        particles = Tools.create_fully_nested_multiple(2,[m1,m2], [a0], [e0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[2,2])
        
        bodies = [x for x in particles if x.is_binary == False]
        binaries = [x for x in particles if x.is_binary == True]

        for b in bodies:
            b.evolve_as_star = False
            b.include_mass_transfer_terms = False
            b.check_for_RLOF_at_pericentre = False

        code.include_flybys = False
        code.enable_tides = False
        code.enable_root_finding = True
        code.verbose_flag = 0
        
        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = True
            b.check_for_physical_collision_or_orbit_crossing = True
        
        rg = (m1+m2)*CONST_G/(CONST_C**2)
        
        for b in bodies:
            b.radius = 100.0*rg

        binary = binaries[0]

        code.add_particles(particles)

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        tend = 0.97e8
        N = 1000
        dt = tend/float(N)
        t = 0.0
        while (t<tend):
            t+=dt

            code.evolve_model(t)
            flag = code.CVODE_flag

            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]

            if len(binaries)>0:
                binary = binaries[0]

                t_print_array.append(t*1.0e-6)
                a_print_array.append(binary.a)
                e_print_array.append(binary.e)
                AP_print_array.append(binary.AP)

                if args.verbose==True:
                    print("t",t*1e-6,'a/au',binary.a,'e',binary.e)

        a_print_array = np.array(a_print_array)
        e_print_array = np.array(e_print_array)
        
        ### Peters 1964 ###
        c0 = a0*(1.0-e0**2)/( pow(e0,12.0/19.0)*pow(1.0 + (121.0/304.0)*e0**2,870.0/2299.0))
        a_an = c0*pow(e_print_array,12.0/19.0)*pow(1.0+(121.0/304.0)*e_print_array**2,870.0/2299.0)/(1.0-e_print_array**2)
        beta = (64.0/5.0)*CONST_G**3*m1*m2*(m1+m2)/(CONST_C**5)
        T_c = a0**4/(4.0*beta)
        T = (768.0/425.0)*T_c*pow(1.0-e0**2,7.0/2.0)

        N_r = 5
        if args.verbose == True:
            print("a_print_array[-1]",a_print_array[-1],"a_an[-1]",a_an[-1])
        assert(round(a_print_array[-1],N_r) == round(a_an[-1],N_r))
        
        print("Test passed")
        
        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(2,1,1)
            plot2 = fig.add_subplot(2,1,2,yscale="log")

            plot1.plot(t_print_array,e_print_array, color='r')

            plot2.plot(t_print_array,a_print_array, color='r',label='$\mathrm{MSE}$')

            plot2.plot(t_print_array,a_an,color='g',linestyle='dashed',linewidth=2,label='$\mathrm{Semi-analytic\,Peters\,(1964)}$')

            fontsize = 15
            plot1.set_ylabel("$e$",fontsize=fontsize)
            plot2.set_ylabel("$a/\mathrm{au}$",fontsize=fontsize)
            plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            handles,labels = plot2.get_legend_handles_labels()
            plot2.legend(handles,labels,loc="upper left",fontsize=0.6*fontsize)

            pyplot.show()

    def test4(self,args):
        print("Test tidal friction in 2-body system")
        
        code = MSE()
        code.enable_tides = True
        code.include_flybys = False
        code.enable_root_finding = False

        code.verbose_flag = 0
        code.relative_tolerance = 1.0e-10
        code.absolute_tolerance_eccentricity_vectors = 1.0e-14

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        M = 0.0009546386983890755 ### Jupiter mass
        R = 40.0*0.1027922358015816*CONST_R_SUN ### 40 R_J
        m_per = 1.0
        mu = m_per*M/(m_per+M)
        a0 = 1.0
        e0 = 0.3
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        omega_crit = np.sqrt(CONST_G*M/(R**3))

        aF = a0*(1.0-e0**2)
        nF = np.sqrt( CONST_G*(M+m_per)/(aF**3) )

        particles = Tools.create_fully_nested_multiple(2,[m_per, M], [a0], [e0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[2,2])
        binaries = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]
        binary = particles[2]

        for b in bodies:
            b.include_mass_transfer_terms = False
            b.check_for_RLOF_at_pericentre = False
            b.include_spin_orbit_1PN_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False

        particles[0].radius = CONST_R_SUN
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0
        particles[1].spin_vec_y = 0.0
        particles[1].spin_vec_z = 4.0e-2/day



        k_L = 0.38
        k_AM = k_L/2.0
        rg = 0.25
        tau = 1e5*0.66*second

        I = rg*M*R**2
        alpha = I/(mu*a0**2)
        T = R**3/(CONST_G*M*tau)
        t_V = 3.0*(1.0 + 2.0*k_AM)**2*T/k_AM
        
        if args.verbose==True:
            print( 't_V',t_V,'M',M,'R',R)
            print("n",n0,"omega_crit",omega_crit,"omega_init",particles[1].spin_vec_z)
            
        particles[0].include_tidal_friction_terms = False
        particles[0].include_tidal_bulges_precession_terms = False
        particles[0].include_rotation_precession_terms = False

        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = True
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-8

        particles[1].apsidal_motion_constant = k_AM
        particles[1].tides_viscous_time_scale = t_V
        particles[1].gyration_radius = rg
        particles[1].tides_viscous_time_scale_prescription = 0

        tD = M*aF**8/(3.0*k_L*tau*CONST_G*m_per*(M+m_per)*R**5)
        particles[2].check_for_physical_collision_or_orbit_crossing = True

        code.add_particles(particles)
        code.verbose_flag = 0
        t = 0.0
        N=100
        tend = 1.0e7
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

            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]
            binary = binaries[0]

            if args.verbose==True:
                print( 'flag',code.CVODE_flag,'t/yr',t,'a/au',binary.a,'e',binary.e)


            t_print_array.append(t)
            a_print_array.append(binary.a)
            n_print_array.append(np.sqrt(CONST_G*(M+m_per)/(binary.a**3)))
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)
            spin_print_array.append( np.sqrt( particles[1].spin_vec_x**2 + particles[1].spin_vec_y**2 + particles[1].spin_vec_z**2) )

            bodies = particles[0:2]
            if args.verbose==True:
                for body in bodies:
                    print( 'S_x',body.spin_vec_x)
                    print( 'S_y',body.spin_vec_y)
                    print( 'S_z',body.spin_vec_z)
                print( '='*50)

        if args.verbose == True:
            print("spin_print_array[-1]",spin_print_array[-1],"n_print_array[-1]",n_print_array[-1])
            print("a_print_array[-1]",a_print_array[-1],"a0(1-e0^2)",a0*(1.0-e0**2))
            
        N_r = 3
        assert round(spin_print_array[-1],N_r) == round(n_print_array[-1],N_r)
        assert round(a_print_array[-1],N_r) == round(aF,N_r)
        #assert len([x for x in range(len(t_print_array)) if round(aF,N_r) not in [round(a*(1.0-e**2),N_r) for a,e in zip( a_print_array,e_print_array)] ] ) == 0
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            fig = pyplot.figure(figsize=(10,10))
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
            plot1.set_ylabel("$a/\mathrm{au}$",fontsize=fontsize)

            plot2 = fig.add_subplot(N_p,1,2)
            plot2.plot(t_print_array*1.0e-6,e_print_array,color='k')
            plot2.set_ylabel("$e$",fontsize=fontsize)

            plot3 = fig.add_subplot(N_p,1,3,yscale="log")

            plot3.plot(t_print_array*1.0e-6,a_print_array*(1.0-e_print_array**2),color='k')
            plot3.axhline(y = a0*(1.0 - e0**2), color='k')
            plot3.set_ylabel("$a(1-e^2)/\mathrm{au}$",fontsize=fontsize)

            plot4 = fig.add_subplot(N_p,1,4)
            plot4.plot(t_print_array*1.0e-6,spin_print_array/n_print_array)
            plot4.set_ylabel("$\Omega/n$",fontsize=fontsize)

            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    def test5(self,args):
        print("Test apsidal motion due to tidal bulges in binary")

        code = MSE()
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

        particles = Tools.create_fully_nested_multiple(2,[m_per, M], [a0], [e0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[2,2])
        binaries = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]
        
        binary = particles[2]
        particles[0].radius = 1.0e-10*R
        particles[1].radius = R

        for b in bodies:
            b.include_mass_transfer_terms = False
            b.check_for_RLOF_at_pericentre = False
            b.include_spin_orbit_1PN_terms = False
            b.spin_vec_x = 0.0
            b.spin_vec_y = 0.0
            b.spin_vec_z = 1.0e-15
            
        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False
            b.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = False

        k_L = 0.41
        k_AM = k_L/2.0


        particles[0].include_tidal_friction_terms = False
        particles[0].include_tidal_bulges_precession_terms = False
        particles[0].include_rotation_precession_terms = False

        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = True
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].apsidal_motion_constant = k_AM
        particles[1].gyration_radius = 0.25
        
        code.add_particles(particles)

        code.relative_tolerance = 1.0e-14
        code.absolute_tolerance_eccentricity_vectors = 1.0e-14
        code.absolute_tolerance_spin_vectors = 1.0e-4
        code.absolute_tolerance_angular_momentum_vectors = 1.0e-4
        code.include_flybys = False
        code.verbose_flag = 0
        
        t = 0.0
        dt = 1.0e6
        tend = 1.0e8

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        g_dot_TB = (15.0/8.0)*n0*(8.0+12.0*e0**2+e0**4)*(m_per/M)*k_AM*pow(R/a0,5.0)/pow(1.0-e0**2,5.0)
        t_TB = 2.0*np.pi/g_dot_TB

        while (t<tend):
            t+=dt
            code.evolve_model(t)

            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]
            #print("?",[b.include_spin_orbit_1PN_terms for b in bodies])
            binary = binaries[0]

            if args.verbose==True:
                print( 'flag',code.CVODE_flag,'t',t,'a/au',binary.a,'e',binary.e,"AP",binary.AP)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

        AP = 0.01 + 2.0*np.pi*tend/(t_TB)
        AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
        
        if args.verbose == True:
            print("Predicted AP",AP,"AP_print_array[-1]",AP_print_array[-1])        

        N_r = 2
        assert round(AP,N_r) == round(AP_print_array[-1],N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)
            
            fig = pyplot.figure(figsize=(10,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,len(t_print_array))
            AP = 0.01 +2.0*np.pi*points/(t_TB*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',linestyle='dotted',linewidth=2)

            plot2.plot(t_print_array,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array,np.fabs((a0-a_print_array)/a0), color='r')
            plot4.plot(t_print_array,np.fabs((e0-e_print_array)/e0), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega/\mathrm{rad}$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)
            
            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    def test6(self,args):
        print("Test apsidal motion due to rotation in binary")

        code = MSE()
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

        particles = Tools.create_fully_nested_multiple(2,[m_per, M], [a0], [e0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[2,2])
        binaries = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]
        
        for b in bodies:
            b.include_mass_transfer_terms = False
            b.check_for_RLOF_at_pericentre = False
            b.include_spin_orbit_1PN_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False
            b.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = False
        
        particles[0].radius = CONST_R_SUN

        particles[0].include_tidal_friction_terms = False
        particles[0].include_tidal_bulges_precession_terms = False
        particles[0].include_rotation_precession_terms = False

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
        particles[1].apsidal_motion_constant = k_AM
        particles[1].gyration_radius = rg

        code.add_particles(particles)

        code.relative_tolerance = 1.0e-14
        code.absolute_tolerance_eccentricity_vectors = 1.0e-14
        code.absolute_tolerance_spin_vectors = 1.0e-4
        code.absolute_tolerance_angular_momentum_vectors = 1.0e-4
        code.include_flybys = False
        code.verbose_flag = 0
        
        t = 0.0
        dt = 1.0e6
        tend = 1.0e7

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        Omega_vec = [particles[1].spin_vec_x,particles[1].spin_vec_y,particles[1].spin_vec_z]
        Omega = np.sqrt(Omega_vec[0]**2 + Omega_vec[1]**2 + Omega_vec[2]**2)
        if args.verbose==True:
            print( 'Omega/n',Omega/n0)

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            
            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]
            binary = binaries[0]
            
            if args.verbose==True:
                print( 'flag',code.CVODE_flag,'t',t,'a',binary.a,'e',binary.e,"AP",binary.AP)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

        g_dot_rot = n0*(1.0 + m_per/M)*k_AM*pow(R/a0,5.0)*(Omega/n0)**2/((1.0-e0**2)**2)
        t_rot = 2.0*np.pi/g_dot_rot

        AP = 0.01 + 2.0*np.pi*tend/(t_rot)
        AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi

        if args.verbose == True:
            print("Predicted AP",AP,"AP_print_array[-1]",AP_print_array[-1])        

        N_r = 3
        assert round(AP,N_r) == round(AP_print_array[-1],N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)

            fig = pyplot.figure(figsize=(10,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,len(t_print_array))
            AP = 0.01 +2.0*np.pi*points/(t_rot*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',linestyle='dotted',linewidth=2)

            plot2.plot(t_print_array,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array,np.fabs((a0-a_print_array)/a0), color='r')
            plot4.plot(t_print_array,np.fabs((e0-e_print_array)/e0), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega/\mathrm{rad}$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)

            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            pyplot.show()

    def test7(self,args):
        print("Test ODE collision detection in 3-body system")

        code = MSE()

        particles = Tools.create_fully_nested_multiple(3,[1.0, 1.2, 0.9], [1.0, 100.0], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01], metallicities=[0.02,0.02,0.02],stellar_types=[1,1,1],object_types=[2,2,2])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        R = 0.03 ### radius of individual objects; collision distance is 2*R

        for b in bodies:
            b.include_mass_transfer_terms = False
            b.check_for_RLOF_at_pericentre = False
            b.include_spin_orbit_1PN_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False
            b.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = False
        

        code.add_particles(particles)
        code.stop_after_root_found = True
        code.enable_tides = False
        code.verbose_flag = 0
        
        t = 0.0
        dt = 1.0e4
        tend = 1.0e6
        t_root = 0.0
        
        t_print_array = []
        a_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.CVODE_flag

            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]

            for b in bodies:
                b.include_mass_transfer_terms = False
                b.check_for_RLOF_at_pericentre = False
                b.include_spin_orbit_1PN_terms = False
                
            binaries[0].check_for_minimum_periapse_distance = False
            binaries[0].check_for_physical_collision_or_orbit_crossing = True

            for body in bodies:
                body.radius = R ### Force stellar radii to custom value

            if args.verbose==True:
                print("="*50)
                print("t/Myr",t*1e-6,"a",binaries[0].a,"e",binaries[0].e,"rp/au",binaries[0].a*(1.0 - binaries[0].e) )
                print( 'secular_breakdown_has_occurred',binaries[0].secular_breakdown_has_occurred)
                print( 'dynamical_instability_has_occurred',binaries[0].dynamical_instability_has_occurred)
                print( 'physical_collision_or_orbit_crossing_has_occurred',binaries[0].physical_collision_or_orbit_crossing_has_occurred)
                print( 'minimum_periapse_distance_has_occurred',binaries[0].minimum_periapse_distance_has_occurred)
                print( 'RLOF_at_pericentre_has_occurred',binaries[0].RLOF_at_pericentre_has_occurred)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                t_root = code.model_time
                if args.verbose==True:
                    print( 'root found at t=',t_root)
                break
        
        if args.verbose == True:
            print("num rp ",a_print_array[-1]*(1.0 - e_print_array[-1])," 2*R ",2*R)

        N_r = 10       
        assert round(a_print_array[-1]*(1.0 - e_print_array[-1]),N_r) == round(2.0*R,N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure()
            plot = fig.add_subplot(1,1,1)
            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))
            plot.axhline(y = bodies[0].radius + bodies[1].radius,color='k')
            plot.set_ylabel("$r_\mathrm{p}/\mathrm{au}$")
            plot.set_xlabel("$t/\mathrm{Myr}$")

            pyplot.show()

    def test8(self,args):
        print("Test ODE minimum periapsis distance root finding")

        code = MSE()

        particles = Tools.create_fully_nested_multiple(3,[1.0, 1.2, 0.9], [1.0, 100.0], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01], metallicities=[0.02,0.02,0.02],stellar_types=[1,1,1],object_types=[2,2,2])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        rp_min = 0.1
        R = 1.0e-10

        for b in bodies:
            b.include_mass_transfer_terms = False
            b.check_for_RLOF_at_pericentre = False
            b.include_spin_orbit_1PN_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False
            b.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = False
        

        code.add_particles(particles)
        code.stop_after_root_found = True
        code.enable_tides = False
        code.verbose_flag = 0
        
        t = 0.0
        dt = 1.0e4
        tend = 1.0e6
        t_root = 0.0
        
        t_print_array = []
        a_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.CVODE_flag

            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]

            for b in bodies:
                b.include_mass_transfer_terms = False
                b.check_for_RLOF_at_pericentre = False
                b.include_spin_orbit_1PN_terms = False
                
            binaries[0].check_for_minimum_periapse_distance = True
            binaries[0].check_for_minimum_periapse_distance_value = rp_min
            binaries[0].check_for_physical_collision_or_orbit_crossing = False

            for body in bodies:
                body.radius = R ### Force stellar radii to custom value

            if args.verbose==True:
                print("="*50)
                print("t/Myr",t*1e-6,"a",binaries[0].a,"e",binaries[0].e,"rp/au",binaries[0].a*(1.0 - binaries[0].e) )
                print( 'secular_breakdown_has_occurred',binaries[0].secular_breakdown_has_occurred)
                print( 'dynamical_instability_has_occurred',binaries[0].dynamical_instability_has_occurred)
                print( 'physical_collision_or_orbit_crossing_has_occurred',binaries[0].physical_collision_or_orbit_crossing_has_occurred)
                print( 'minimum_periapse_distance_has_occurred',binaries[0].minimum_periapse_distance_has_occurred)
                print( 'RLOF_at_pericentre_has_occurred',binaries[0].RLOF_at_pericentre_has_occurred)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                t_root = code.model_time
                if args.verbose==True:
                    print( 'root found at t=',t_root)
                break
        
        if args.verbose == True:
            print("num rp ",a_print_array[-1]*(1.0 - e_print_array[-1]),"rp_min",rp_min)

        N_r = 10            
        assert round(a_print_array[-1]*(1.0 - e_print_array[-1]),N_r) == round(rp_min,N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure()
            plot = fig.add_subplot(1,1,1)
            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))
            plot.axhline(y = rp_min,color='k')
            plot.set_ylabel("$r_\mathrm{p}/\mathrm{au}$")
            plot.set_xlabel("$t/\mathrm{Myr}$")

            pyplot.show()

    def test9(self,args):
        print("Test adiabatic mass loss")

        code = MSE()
        
        m1i = 15.0
        m2i = 10.0
        m3i = 1.0
        a_in_i = 1.0
        a_out_i = 1000.0
        particles = Tools.create_fully_nested_multiple(3,[m1i,m2i,m3i], [a_in_i,a_out_i], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01],metallicities=[0.02,0.02,0.02],stellar_types=[1,1,1],object_types=[1,1,1])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        #m1dot = -1.0e-7
        #m2dot = -1.0e-8
        #m3dot = -1.0e-9
        #mdots = [m1dot,m2dot,m3dot]
        #for index,body in enumerate(bodies):
        #    body.mass_dot = mdots[index]

        code.add_particles(particles)
        code.enable_tides = False

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6

        t_print_array = []
        m_in_print_array = []
        m_out_print_array = []
        a_in_print_array = []
        a_out_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.CVODE_flag

            particles = code.particles
            binaries = [x for x in particles if x.is_binary == True]
            bodies = [x for x in particles if x.is_binary == False]

            for b in bodies:
                b.include_mass_transfer_terms = False
                b.check_for_RLOF_at_pericentre = False
                b.include_spin_orbit_1PN_terms = False

            if args.verbose==True:
                print( 't/Myr',t*1e-6,'m1',bodies[0].mass,'m2',bodies[1].mass,'m3',bodies[2].mass,'a/au',[b.a for b in binaries])

            t_print_array.append(t*1.0e-6)
            m_in_print_array.append(bodies[0].mass + bodies[1].mass)
            m_out_print_array.append(bodies[0].mass + bodies[1].mass + bodies[2].mass)
            a_in_print_array.append(binaries[0].a)
            a_out_print_array.append(binaries[1].a)
            e_print_array.append(binaries[0].e)

        m_in_print_array = np.array(m_in_print_array)
        m_out_print_array = np.array(m_out_print_array)
        a_in_print_array = np.array(a_in_print_array)
        a_out_print_array = np.array(a_out_print_array)

        inner_adiabat = m_in_print_array*a_in_print_array
        outer_adiabat = m_out_print_array*a_out_print_array

        if args.verbose == True:
            print("a_in_i*(m1i+m2i)",a_in_i*(m1i+m2i),"num ",inner_adiabat[1])
            print("a_out_i*(m1i+m2i+m3i)",a_out_i*(m1i+m2i+m3i),"num ",outer_adiabat[1])
            
        N_r = 4
        assert round(a_in_i*(m1i+m2i),N_r) == round(inner_adiabat[-1],N_r)
        assert round(a_out_i*(m1i+m2i+m3i),N_r) == round(outer_adiabat[-1],N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            fig = pyplot.figure(figsize=(10,8))
            plot1 = fig.add_subplot(2,1,1,yscale="log")
            plot2 = fig.add_subplot(2,1,2,yscale="log")
            plot1.plot(t_print_array,inner_adiabat,color='k',linestyle='solid',label='$a_1$',zorder=10)
            plot2.plot(t_print_array,outer_adiabat,color='k',linestyle='solid',label='$a_2$',zorder=10)
            
            handles,labels = plot1.get_legend_handles_labels()
            plot1.legend(handles,labels,loc="lower left",fontsize=18)

            plot1.set_ylabel(r"$a_1(m_1+m_2)/(\mathrm{au\,M_\odot})$")
            plot2.set_ylabel(r"$a_2(m_1+m_2+m_3)/(\mathrm{au\,M_\odot})$")
            plot2.set_xlabel("$t/\mathrm{Myr}$")

            pyplot.show()

    def test10(self,args):
        print("Test flybys module: instantaneous change -- SNe in binary")

        code = MSE()
        CONST_G = code.CONST_G
        CONST_km_per_s_to_AU_per_yr = code.CONST_KM_PER_S
                
        a1 = 10.0
        e1 = 0.1
        m1 = 1.0
        m2 = 0.8

        INCL1 = 0.1
        AP1 = 0.2
        LAN1 = 0.3
        f1 = 60.5*np.pi/180.0

        particles = Tools.create_fully_nested_multiple(2, [m1,m2], [a1], [e1], [INCL1], [AP1], [LAN1], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[1,1])

        binary = particles[2]
        binary.sample_orbital_phase_randomly = False
        binary.TA = f1
        
        delta_m1 = -0.5
        V_k_vec = np.array([0.0,1.0,2.0])*CONST_km_per_s_to_AU_per_yr
        #V_k_vec = np.array([0.0,0.0,0.0])

        particles[0].instantaneous_perturbation_delta_mass = delta_m1
        particles[0].instantaneous_perturbation_delta_VX = V_k_vec[0]
        particles[0].instantaneous_perturbation_delta_VY = V_k_vec[1]
        particles[0].instantaneous_perturbation_delta_VZ = V_k_vec[2]

        code.add_particles(particles)

        if args.verbose==True:
            print( '='*50)
            print( 'pre')
            print( 'a',binary.a,'e',binary.e,'INCL',binary.INCL*180.0/np.pi,'AP',binary.AP*180.0/np.pi,'LAN',binary.LAN*180.0/np.pi,'TA',binary.TA)

        code.apply_user_specified_instantaneous_perturbation()
        
        if args.verbose==True:
            print( '='*50)
            print( 'post')
            print( 'a',binary.a,'e',binary.e,'INCL',binary.INCL*180.0/np.pi,'AP',binary.AP*180.0/np.pi,'LAN',binary.LAN*180.0/np.pi,'TA',binary.TA)

        ### Compute analytic result (e.g., https://ui.adsabs.harvard.edu/abs/2016ComAC...3....6T/abstract) ###
        r1 = a1*(1.0-e1**2)/(1.0 + e1*np.cos(f1))
        v1_tilde = np.sqrt( CONST_G*(m1+m2)/(a1*(1.0-e1**2) ) )
        e1_vec_hat,j1_vec_hat = compute_e_and_j_hat_vectors(INCL1,AP1,LAN1)
        q1_vec_hat = np.cross(j1_vec_hat,e1_vec_hat)
        r1_vec = r1*( e1_vec_hat * np.cos(f1) + q1_vec_hat * np.sin(f1) )
        v1_vec = v1_tilde*( -e1_vec_hat * np.sin(f1) + q1_vec_hat * (e1 + np.cos(f1) ) )

        r1_dot_v1 = np.sum([x*y for x,y in zip(r1_vec,v1_vec)])
        r1_dot_V_k = np.sum([x*y for x,y in zip(r1_vec,V_k_vec)])
        v1_dot_V_k = np.sum([x*y for x,y in zip(v1_vec,V_k_vec)])
        V_k_dot_V_k = np.sum([x*y for x,y in zip(V_k_vec,V_k_vec)])
        v1c_sq = CONST_G*(m1+m2)/a1
        
        a1_p = a1*(1.0 + delta_m1/(m1+m2))*pow( 1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq, -1.0)
        j1_p = pow( (m1+m2)/(m1+m2+delta_m1), 2.0)*(1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq)*( 1.0 - e1**2 \
            + (1.0/(CONST_G*(m1+m2)*a1))*( r1**2*( 2.0*v1_dot_V_k + V_k_dot_V_k) - 2.0*r1_dot_v1*r1_dot_V_k - r1_dot_V_k**2) )
        e1_p = np.sqrt(1.0 - j1_p)
        
        if args.verbose==True:
            print( 'analytic results Toonen+16: ','new a1 = ',a1_p,'; new e1 = ',e1_p)
        
        N_r = 10
        assert round(binary.a,N_r) == round(a1_p,N_r)
        assert round(binary.e,N_r) == round(e1_p,N_r)
        
        print("Test passed")

        code.reset()
        
    def test11(self,args):
        print("Test flybys module: instantaneous change -- SNe in triple")

        code = MSE()
        CONST_G = code.CONST_G
        CONST_km_per_s_to_AU_per_yr = code.CONST_KM_PER_S        
        
        m1 = 1.0
        m2 = 0.8
        m3 = 1.2
        a1 = 10.0
        a2 = 100.0
        e1 = 0.1
        e2 = 0.3
        INCL1 = 0.1
        INCL2 = 0.5
        AP1 = 0.1
        AP2 = 1.0
        LAN1 = 0.1
        LAN2 = 2.0
        f1 = 60.5*np.pi/180.0
        f2 = 30.5*np.pi/180.0

        INCLs = [INCL1,INCL2]
        APs = [AP1,AP2]
        LANs = [LAN1,LAN2]
        masses = [m1,m2,m3]
        particles = Tools.create_fully_nested_multiple(3,masses, [a1,a2], [e1,e2], INCLs, APs, LANs, metallicities=[0.02,0.02,0.02],stellar_types=[1,1,1],object_types=[1,1,1])
        
        inner_binary = particles[3]
        outer_binary = particles[4]
        inner_binary.sample_orbital_phase_randomly = 0
        outer_binary.sample_orbital_phase_randomly = 0
        inner_binary.TA = f1
        outer_binary.TA = f2
        
        delta_m1 = -0.3
        
        km_p_s_to_AU_p_yr = 0.21094502112788768
        V_k_vec = np.array([1.0,2.0,2.0])*CONST_km_per_s_to_AU_per_yr
        V_k_vec = np.array([0.0,0.0,0.0])
        
        particles[0].instantaneous_perturbation_delta_mass = delta_m1
        particles[0].instantaneous_perturbation_delta_VX = V_k_vec[0]
        particles[0].instantaneous_perturbation_delta_VY = V_k_vec[1]
        particles[0].instantaneous_perturbation_delta_VZ = V_k_vec[2]

        code.add_particles(particles)
        
        if args.verbose==True:
            print( '='*50)
            print( 'pre')
            print( 'inner','a',inner_binary.a,'e',inner_binary.e,'INCL',inner_binary.INCL*180.0/np.pi,'AP',inner_binary.AP*180.0/np.pi,'LAN',inner_binary.LAN*180.0/np.pi)
            print( 'outer','a',outer_binary.a,'e',outer_binary.e,'INCL',outer_binary.INCL*180.0/np.pi,'AP',outer_binary.AP*180.0/np.pi,'LAN',outer_binary.LAN*180.0/np.pi)
        
        code.apply_user_specified_instantaneous_perturbation()
        
        if args.verbose==True:
            print( '='*50)
            print( 'post')
            print( 'inner','a',inner_binary.a,'e',inner_binary.e,'INCL',inner_binary.INCL*180.0/np.pi,'AP',inner_binary.AP*180.0/np.pi,'LAN',inner_binary.LAN*180.0/np.pi)
            print( 'outer','a',outer_binary.a,'e',outer_binary.e,'INCL',outer_binary.INCL*180.0/np.pi,'AP',outer_binary.AP*180.0/np.pi,'LAN',outer_binary.LAN*180.0/np.pi)

        
        ### Compute analytic result (e.g., https://ui.adsabs.harvard.edu/abs/2016ComAC...3....6T/abstract) ###
        r1 = a1*(1.0-e1**2)/(1.0 + e1*np.cos(f1))
        v1_tilde = np.sqrt( CONST_G*(m1+m2)/(a1*(1.0-e1**2) ) )
        e1_vec_hat,j1_vec_hat = compute_e_and_j_hat_vectors(INCL1,AP1,LAN1)
        q1_vec_hat = np.cross(j1_vec_hat,e1_vec_hat)
        r1_vec = r1*( e1_vec_hat * np.cos(f1) + q1_vec_hat * np.sin(f1) )
        v1_vec = v1_tilde*( -e1_vec_hat * np.sin(f1) + q1_vec_hat * (e1 + np.cos(f1) ) )

        r1_dot_v1 = np.sum([x*y for x,y in zip(r1_vec,v1_vec)])
        r1_dot_V_k = np.sum([x*y for x,y in zip(r1_vec,V_k_vec)])
        v1_dot_V_k = np.sum([x*y for x,y in zip(v1_vec,V_k_vec)])
        V_k_dot_V_k = np.sum([x*y for x,y in zip(V_k_vec,V_k_vec)])
        v1c_sq = CONST_G*(m1+m2)/a1

        r2 = a2*(1.0-e2**2)/(1.0 + e2*np.cos(f2))
        v2_tilde = np.sqrt( CONST_G*(m1+m2+m3)/(a2*(1.0-e2**2) ) )
        e2_vec_hat,j2_vec_hat = compute_e_and_j_hat_vectors(INCL2,AP2,LAN2)
        q2_vec_hat = np.cross(j2_vec_hat,e2_vec_hat)
        r2_vec = r2*( e2_vec_hat * np.cos(f2) + q2_vec_hat * np.sin(f2) )
        v2_vec = v2_tilde*( -e2_vec_hat * np.sin(f2) + q2_vec_hat * (e2 + np.cos(f2) ) )
    
        Delta_r2_vec = (delta_m1/( (m1+m2) + delta_m1) ) * (m2/(m1+m2)) * r1_vec
        Delta_v2_vec = (delta_m1/( (m1+m2) + delta_m1) ) * ( (m2/(m1+m2)) * v1_vec + V_k_vec*(1.0 + m1/delta_m1) )
        r2_vec_p = r2_vec + Delta_r2_vec
        r2p = np.sqrt(np.dot(r2_vec_p,r2_vec_p))


        r2_dot_v2 = np.sum([x*y for x,y in zip(r2_vec,v2_vec)])
        r2_dot_V_k = np.sum([x*y for x,y in zip(r2_vec,V_k_vec)])
        v2c_sq = CONST_G*(m1+m2+m3)/a2
        v2_dot_Delta_v2_vec = np.dot(v2_vec,Delta_v2_vec)
        Delta_v2_vec_dot_Delta_v2_vec = np.dot(Delta_v2_vec,Delta_v2_vec)
        
        a1_p = a1*(1.0 + delta_m1/(m1+m2))*pow( 1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq, -1.0)
        j1_p = pow( (m1+m2)/(m1+m2+delta_m1), 2.0)*(1.0 + 2.0*(a1/r1)*(delta_m1/(m1+m2)) - 2.0*v1_dot_V_k/v1c_sq - V_k_dot_V_k/v1c_sq)*( 1.0 - e1**2 \
            + (1.0/(CONST_G*(m1+m2)*a1))*( r1**2*( 2.0*v1_dot_V_k + V_k_dot_V_k) - 2.0*r1_dot_v1*r1_dot_V_k - r1_dot_V_k**2) )
        e1_p = np.sqrt(1.0 - j1_p)

        a2_p = a2*(1.0 + delta_m1/(m1+m2+m3))*pow( 1.0 + 2.0*(a2/r2p)*(delta_m1/(m1+m2+m3)) - 2.0*v2_dot_Delta_v2_vec/v2c_sq - Delta_v2_vec_dot_Delta_v2_vec/v2c_sq + 2.0*a2*(r2-r2p)/(r2*r2p), -1.0)

        alpha = (-delta_m1/(m1+m2+delta_m1)) * m2/(m1+m2)
        j2_p = ((m1+m2+m3)/(m1+m2+m3+delta_m1))**2 * (1.0 + 2.0*(a2/r2p)*(delta_m1/(m1+m2+m3)) + 2.0*a2*(r2-r2p)/(r2*r2p) - 2.0*np.dot(v2_vec,Delta_v2_vec)/v2c_sq - Delta_v2_vec_dot_Delta_v2_vec/v2c_sq)*( (1.0-e2**2) + (1.0/(CONST_G*(m1+m2+m3)*a2))*( r2**2*(2.0*np.dot(v2_vec,Delta_v2_vec) + np.dot(Delta_v2_vec,Delta_v2_vec)) + (-2.0*alpha*np.dot(r1_vec,r2_vec) + alpha**2*r1**2)*np.dot(v2_vec + Delta_v2_vec,v2_vec+Delta_v2_vec) + 2.0*np.dot(r2_vec,v2_vec)*( alpha*np.dot(r1_vec,v2_vec) - np.dot(r2_vec,Delta_v2_vec) + alpha*np.dot(r1_vec,Delta_v2_vec)) - (-alpha*np.dot(r1_vec,v2_vec) + np.dot(r2_vec,Delta_v2_vec) - alpha*np.dot(r1_vec,Delta_v2_vec))**2 ) ) 
        e2_p = np.sqrt(1.0 - j2_p)
        
        if args.verbose==True:
            print( 'analytic results Toonen+16: ','new a1 = ',a1_p,'; new e1 = ',e1_p)
            print( 'analytic results Toonen+16: ','new a2 = ',a2_p,'; new e2 = ',e2_p)
        
        N_r = 10
        assert round(inner_binary.a,N_r) == round(a1_p,N_r)
        assert round(inner_binary.e,N_r) == round(e1_p,N_r)

        assert round(outer_binary.a,N_r) == round(a2_p,N_r)
        assert round(outer_binary.e,N_r) == round(e2_p,N_r)

        print("Test passed")

        code.reset()
        
    def test12(self,args):
        print("Test flybys module: using analytic formulae")

        code = MSE()
        CONST_G = code.CONST_G
        
        ### binary orbit ###
        a = 1.0
        e = 0.1
        m1 = 1.0
        m2 = 0.8
        M_per = 1.0
        E = 2.0
        Q = 100.0
        INCL = 0.4*np.pi
        AP = 0.25*np.pi
        LAN = 0.25*np.pi
        
        masses = [m1,m2]
        m = m1 + m2
        particles = Tools.create_fully_nested_multiple(2,masses, [a], [e], [INCL], [AP], [LAN], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[1,1])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]
        binary = particles[2]

        for b in bodies:
            b.evolve_as_star = False
            b.include_mass_transfer_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False

        t_ref = 1.0e5 ### not actually used in this case, but needs to be specified for external particles

        external_particle = Particle(mass = M_per, is_binary=True, is_external=True, external_t_ref=t_ref, e=E, external_r_p = Q, INCL = 1.0e-10, AP = 1.0e-10, LAN = 1.0e-10)
        
        particles.append(external_particle)

        code = MSE()
        code.add_particles(particles)

        code.include_quadrupole_order_terms = True
        code.include_octupole_order_binary_pair_terms = True
        code.include_hexadecupole_order_binary_pair_terms = False
        code.include_dotriacontupole_order_binary_pair_terms = False

        code.include_flybys = False
        code.enable_tides = False
        code.enable_root_finding = False
        
        code.apply_external_perturbation_assuming_integrated_orbits()
        Delta_e = binary.e-e

        ### compute analytic result (https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.5630H/abstract) ###
        e_vec_hat,j_vec_hat = compute_e_and_j_hat_vectors(INCL,AP,LAN)
        e_vec = e_vec_hat*e
        j_vec = j_vec_hat*np.sqrt(1.0-e**2)
        ex = e_vec[0]
        ey = e_vec[1]
        ez = e_vec[2]
        jx = j_vec[0]
        jy = j_vec[1]
        jz = j_vec[2]

        eps_SA = (M_per/np.sqrt(m*(m+M_per)))*pow(a/Q,3.0/2.0)*pow(1.0 + E,-3.0/2.0)
        eps_oct = (a/Q)*(m1-m2)/((1.0+E)*m)
        Delta_e_an = (5*eps_SA*(np.sqrt(1 - E**(-2))*((1 + 2*E**2)*ey*ez*jx + (1 - 4*E**2)*ex*ez*jy + 2*(-1 + E**2)*ex*ey*jz) + 3*E*ez*(ey*jx - ex*jy)*np.arccos(-1.0/E)))/(2.*E*np.sqrt(ex**2 + ey**2 + ez**2))
        
        Delta_e_an += -(5*eps_oct*eps_SA*(np.sqrt(1 - E**(-2))*(ez*jy*(14*ey**2 + 6*jx**2 - 2*jy**2 + 8*E**4*(-1 + ey**2 + 8*ez**2 + 2*jx**2 + jy**2) + E**2*(-4 - 31*ey**2 + 32*ez**2 - 7*jx**2 + 9*jy**2)) - ey*(2*(7*ey**2 + jx**2 - jy**2) + 8*E**4*(-1 + ey**2 + 8*ez**2 + 4*jx**2 + jy**2) + E**2*(-4 - 31*ey**2 + 32*ez**2 + 11*jx**2 + 9*jy**2))*jz + ex**2*(-((14 + 45*E**2 + 160*E**4)*ez*jy) + 3*(14 - 27*E**2 + 16*E**4)*ey*jz) + 2*(-2 + 9*E**2 + 8*E**4)*ex*jx*(7*ey*ez + jy*jz)) + 3*E**3*(ez*jy*(-4 - 3*ey**2 + 32*ez**2 + 5*jx**2 + 5*jy**2) + ey*(4 + 3*ey**2 - 32*ez**2 - 15*jx**2 - 5*jy**2)*jz + ex**2*(-73*ez*jy + 3*ey*jz) + 10*ex*jx*(7*ey*ez + jy*jz))*np.arccos(-1.0/E)))/(32.*E**2*np.sqrt(ex**2 + ey**2 + ez**2))

        if args.verbose==True:
            print( 'SecularMultiple Delta e = ',Delta_e,'; analytic expression: Delta e = ',Delta_e_an)

        N_r = 8
        assert round(Delta_e,N_r) == round(Delta_e_an,N_r)

        print("Test passed")

        code.reset()

    def test13(self,args):
        print("Test flybys module: using analytic formulae and with single perturber")

        code = MSE()
        CONST_G = code.CONST_G
        
        ### binary orbit ###
        a = 1.0
        e = 0.1
        m1 = 1.0
        m2 = 0.8
        M_per = 1.0
        E = 2.0
        Q = 100.0
        INCL = 0.4*np.pi
        AP = 0.25*np.pi
        LAN = 0.25*np.pi
        
        masses = [m1,m2]
        m = m1 + m2
        particles = Tools.create_fully_nested_multiple(2,masses, [a], [e], [INCL], [AP], [LAN], metallicities=[0.02,0.02],stellar_types=[1,1],object_types=[1,1])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]
        binary = particles[2]

        for b in bodies:
            b.evolve_as_star = False
            b.include_mass_transfer_terms = False

        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = False

        t_ref = 1.0e5 ### not actually used in this case, but needs to be specified for external particles

        #external_particle = Particle(mass = M_per, is_binary=True, is_external=True, external_t_ref=t_ref, e=E, external_r_p = Q, INCL = 1.0e-10, AP = 1.0e-10, LAN = 1.0e-10)
        INCL_per = 1.0e-10
        AP_per = 1.0e-10
        LAN_per = 1.0e-10
        e_vec_hat,j_vec_hat = compute_e_and_j_hat_vectors(INCL_per,AP_per,LAN_per)
        
        code = MSE()
        code.add_particles(particles)

        code.include_quadrupole_order_terms = True
        code.include_octupole_order_binary_pair_terms = True
        code.include_hexadecupole_order_binary_pair_terms = False
        code.include_dotriacontupole_order_binary_pair_terms = False

        code.include_flybys = False
        code.enable_tides = False
        code.enable_root_finding = False
        
        code.apply_external_perturbation_assuming_integrated_orbits_single_perturber(M_per, E, Q, e_vec_hat[0], e_vec_hat[1], e_vec_hat[2], j_vec_hat[0], j_vec_hat[1], j_vec_hat[2])
        Delta_e = binary.e-e

        ### compute analytic result (https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.5630H/abstract) ###
        e_vec_hat,j_vec_hat = compute_e_and_j_hat_vectors(INCL,AP,LAN)
        e_vec = e_vec_hat*e
        j_vec = j_vec_hat*np.sqrt(1.0-e**2)
        ex = e_vec[0]
        ey = e_vec[1]
        ez = e_vec[2]
        jx = j_vec[0]
        jy = j_vec[1]
        jz = j_vec[2]

        eps_SA = (M_per/np.sqrt(m*(m+M_per)))*pow(a/Q,3.0/2.0)*pow(1.0 + E,-3.0/2.0)
        eps_oct = (a/Q)*(m1-m2)/((1.0+E)*m)
        Delta_e_an = (5*eps_SA*(np.sqrt(1 - E**(-2))*((1 + 2*E**2)*ey*ez*jx + (1 - 4*E**2)*ex*ez*jy + 2*(-1 + E**2)*ex*ey*jz) + 3*E*ez*(ey*jx - ex*jy)*np.arccos(-1.0/E)))/(2.*E*np.sqrt(ex**2 + ey**2 + ez**2))
        
        Delta_e_an += -(5*eps_oct*eps_SA*(np.sqrt(1 - E**(-2))*(ez*jy*(14*ey**2 + 6*jx**2 - 2*jy**2 + 8*E**4*(-1 + ey**2 + 8*ez**2 + 2*jx**2 + jy**2) + E**2*(-4 - 31*ey**2 + 32*ez**2 - 7*jx**2 + 9*jy**2)) - ey*(2*(7*ey**2 + jx**2 - jy**2) + 8*E**4*(-1 + ey**2 + 8*ez**2 + 4*jx**2 + jy**2) + E**2*(-4 - 31*ey**2 + 32*ez**2 + 11*jx**2 + 9*jy**2))*jz + ex**2*(-((14 + 45*E**2 + 160*E**4)*ez*jy) + 3*(14 - 27*E**2 + 16*E**4)*ey*jz) + 2*(-2 + 9*E**2 + 8*E**4)*ex*jx*(7*ey*ez + jy*jz)) + 3*E**3*(ez*jy*(-4 - 3*ey**2 + 32*ez**2 + 5*jx**2 + 5*jy**2) + ey*(4 + 3*ey**2 - 32*ez**2 - 15*jx**2 - 5*jy**2)*jz + ex**2*(-73*ez*jy + 3*ey*jz) + 10*ex*jx*(7*ey*ez + jy*jz))*np.arccos(-1.0/E)))/(32.*E**2*np.sqrt(ex**2 + ey**2 + ez**2))

        if args.verbose==True:
            print( 'SecularMultiple Delta e = ',Delta_e,'; analytic expression: Delta e = ',Delta_e_an)

        N_r = 8
        assert round(Delta_e,N_r) == round(Delta_e_an,N_r)

        print("Test passed")

        code.reset()


    def test14(self,args):
        print('Test compact object merger remnant properties')
        
        code = MSE()


        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        
        N = 20000
        
        m2 = 1.0
        h_vec_unit = np.array( [0.0,0.0,1.0] )
        e_vec_unit = np.array( [1.0,0.0,0.0] )
        
        seed = 0
        np.random.seed(seed)
        
        vs = []
        alphas = []
        delta_Ms = []
        for index in range(N):
            q = np.random.random()
            m1 = q * m2
            M = m1 + m2
            chi1 = np.random.random()
            chi2 = np.random.random()
                    
            spin_vec_1_unit = sample_random_vector_on_unit_sphere()
            spin_vec_2_unit = sample_random_vector_on_unit_sphere()
            v_recoil_vec,alpha_vec,M_final = code.determine_compact_object_merger_properties(m2,m1,chi2,chi1,spin_vec_2_unit,spin_vec_1_unit,h_vec_unit,e_vec_unit)
            v_recoil = np.linalg.norm(v_recoil_vec)/code.CONST_KM_PER_S ### convert AU/yr to km/s
            vs.append(v_recoil)
            
            alpha = np.linalg.norm(alpha_vec)

            alphas.append(alpha)
            delta_Ms.append( (M-M_final)/M )


        if args.verbose==True:
            print("mean vs/(km/s)",np.mean(np.array(vs)),round(np.mean(np.array(vs)),0))
            print("mean delta_Ms",np.mean(np.array(delta_Ms)))
            print("mean alphas",np.mean(np.array(alphas)))

        assert(round(np.mean(np.array(vs))/10.0,0) == 31)
        assert(round(np.mean(np.array(delta_Ms)),3) == 0.036)
        assert(round(np.mean(np.array(alphas)),3) == 0.631)
        
        code.reset()
        
        if args.plot == True:
            Nb=50
            fontsize=20
            from matplotlib import pyplot
            fig=pyplot.figure()
            plot=fig.add_subplot(1,1,1)
            plot.hist(vs,bins=np.linspace(0.0,1000.0,Nb),histtype='step')
            plot.set_xlabel("$v_\mathrm{recoil}/\mathrm{km/s}$",fontsize=fontsize)
            fig.savefig("v_recoil.pdf")
            
            fig=pyplot.figure()
            plot=fig.add_subplot(1,1,1)
            plot.hist(alphas,bins=np.linspace(0.0,1.0,Nb),histtype='step')
            plot.set_xlabel(r"$\chi_\mathrm{final}$",fontsize=fontsize)
            fig.savefig("chi_final.pdf")

            fig=pyplot.figure()
            plot=fig.add_subplot(1,1,1)
            plot.hist(delta_Ms,bins=np.linspace(0.0,0.2,Nb),histtype='step')
            plot.set_xlabel(r"$\Delta M/M$",fontsize=fontsize)
            fig.savefig("Delta_M.pdf")

            pyplot.show()
        
    def test15(self,args):
        print('Test sample elementary distributions')
        
        code = MSE()

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        
        N = 80000
        N2 = 20000

        seed = 0
        np.random.seed(seed)
       
        vs = []
        vs2 = []
        sigma_km_s = 265.0
        sigma = sigma_km_s*code.CONST_KM_PER_S
        mu_km_s = 20.0
        mu=mu_km_s*code.CONST_KM_PER_S

        r_hat_vecs = []
        theta_hat_vecs = []
        phi_hat_vecs = []

        for index in range(N):
            vx,vy,vz = code.test_sample_from_3d_maxwellian_distribution(sigma)
            v = np.array([vx,vy,vz])
            vs.append(np.linalg.norm(v))
            v2 = code.test_sample_from_normal_distribution(mu,sigma)
            vs2.append(v2)
            
        vs = np.array(vs)
        vs2 = np.array(vs2)

        for index in range(N2):
            r_hat_vec,theta_hat_vec,phi_hat_vec = code.test_sample_spherical_coordinates_unit_vectors_from_isotropic_distribution()
            r_hat_vecs.append(r_hat_vec)
            theta_hat_vecs.append(theta_hat_vec)
            phi_hat_vecs.append(phi_hat_vec)

        ### Check that r, theta, and phi hat vectors have zero components on average ###
        tol = 1e-2
        assert( np.mean( np.array( [x[0] for x in r_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[1] for x in r_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[2] for x in r_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[0] for x in theta_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[1] for x in theta_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[2] for x in theta_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[0] for x in phi_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[1] for x in phi_hat_vecs] )) <= tol )
        assert( np.mean( np.array( [x[2] for x in phi_hat_vecs] )) <= tol )

        ### Check orthogonality of r, theta, and phi hat vectors ###
        tol = 1e-12
        assert( np.sum( np.array([np.dot(x,y) for x,y in zip(r_hat_vecs,theta_hat_vecs)])) <= tol)
        assert( np.sum( np.array([np.dot(x,y) for x,y in zip(r_hat_vecs,phi_hat_vecs)])) <= tol)
        assert( np.sum( np.array([np.dot(x,y) for x,y in zip(theta_hat_vecs,phi_hat_vecs)])) <= tol)

        ### Check properties of Maxwellian and normal distributions ###
        N_r=0

        if args.verbose==True:
            print("Maxwellian mean vs/(km/s)",np.mean(np.array(vs))," an ", 2.0*sigma*np.sqrt(2.0/np.pi))
            print("Normal mean vs/(km/s)",np.mean(np.array(vs2))," an ", mu)
            print("Normal std vs/(km/s)",np.std(np.array(vs2))," an ", sigma)
        
        assert(round(np.mean(np.array(vs)),N_r) == round(2.0*sigma*np.sqrt(2.0/np.pi),N_r))
        assert(round(np.mean(np.array(vs2)),N_r) == round(mu,N_r))
        assert(round(np.std(np.array(vs2)),N_r) == round(sigma,N_r))

        code.reset()

        if args.plot == True:
            Nb=100
            fontsize=20
            from matplotlib import pyplot
            fig=pyplot.figure()
            plot=fig.add_subplot(1,1,1)
            plot.hist(vs/code.CONST_KM_PER_S,bins=np.linspace(0.0,1000.0,Nb),histtype='step',density=True)
            plot.set_xlabel("$v/(\mathrm{km/s})$",fontsize=fontsize)
            plot.set_title("Maxwellian")
            
            points=np.linspace(0.0,1000.0,1000)
            PDF_an = np.sqrt(2.0/np.pi) * (points**2/(sigma_km_s**3)) * np.exp( -points**2/(2.0*sigma_km_s**2) )
            plot.plot(points,PDF_an, color='tab:green')

            fig=pyplot.figure()
            plot=fig.add_subplot(1,1,1)
            plot.hist(vs2/code.CONST_KM_PER_S,bins=np.linspace(-1000,1000.0,Nb),histtype='step',density=True)
            plot.set_xlabel("$v/(\mathrm{km/s})$",fontsize=fontsize)
            plot.set_title("Normal")

            points=np.linspace(-1000,1000.0,1000)
            PDF_an = (1.0/(sigma_km_s*np.sqrt(2.0*np.pi))) * np.exp( - (points-mu_km_s)**2/(2.0*sigma_km_s**2))
            plot.plot(points,PDF_an, color='tab:green')

            pyplot.show()


    def test16(self,args):
        print('Test sample Kroupa 93 IMF')
        
        code = MSE()

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        
        N = 100000

        seed = 0
        np.random.seed(seed)

        ms = []
        for index in range(N):
            m = code.test_sample_from_kroupa_93_imf()
            ms.append(m)
        ms = np.array(ms)

        assert(round(np.mean(ms),1) == 0.50)

        code.reset()

        if args.verbose==True:
            print("mean ms/MSun",np.mean(ms))
        
        if args.plot == True:
            Nb=100
            fontsize=20
            from matplotlib import pyplot
            fig=pyplot.figure()
            plot=fig.add_subplot(1,1,1,yscale="log")
            plot.hist(np.log10(ms),bins=np.linspace(-1.0,2.0,Nb),histtype='step',density=True,color='tab:red')
            plot.set_xlabel("$\mathrm{log}_{10}(m/\mathrm{M}_\odot)$",fontsize=fontsize)
            plot.set_ylabel("$\mathrm{PDF}$",fontsize=fontsize)        

            
            points=np.linspace(-1.0,2.0,Nb)
            PDF_an = [np.log(10.0)*pow(10.0,log10m)*kroupa_93_imf(pow(10.0,log10m)) for log10m in points]
            plot.plot(points,PDF_an, color='tab:green')

            pyplot.show()

    def test17(self,args):
        print('Test kick velocity recipes')
        
        code = MSE()
        code.verbose_flag = 0
        
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_KM_PER_S = code.CONST_KM_PER_S
        N = 100
        if args.mode == 1:
            N = 500
        
        seed = 0
        np.random.seed(seed)

        kick_distributions = [1,2,3,4,5]
        N_k = len(kick_distributions)

        vs_NS = [[] for x in range(N_k)]
        vs_BH = [[] for x in range(N_k)]
        vs = [[] for x in range(N_k)]
        alpha = 2.7
        m1 = 8.0
        m2 = 100.0
        for index_kick,kick_distribution in enumerate(kick_distributions):
            for index in range(N):
                x = np.random.random()
                m = pow( x*(pow(m2,1.0-alpha) - pow(m1,1.0-alpha)) + pow(m1,1.0-alpha), 1.0/(1.0-alpha) )
                kw,v = code.test_kick_velocity(kick_distribution,m)
                
                if args.verbose == True:
                    print("index_kick",index_kick,"m",m,"kw",kw,"v/(km/s)",v/CONST_KM_PER_S)
                
                vs[index_kick].append(v/CONST_KM_PER_S)
                if kw==13:
                    vs_NS[index_kick].append(v/CONST_KM_PER_S)
                if kw==14:
                    vs_BH[index_kick].append(v/CONST_KM_PER_S)

            vs[index_kick] = np.array(vs[index_kick])
            vs_NS[index_kick] = np.array(vs_NS[index_kick])
            vs_BH[index_kick] = np.array(vs_BH[index_kick])

        code.reset()
       
        if args.plot == True:
            Nb=50
            fontsize=16
            labelsize=12
            from matplotlib import pyplot

            if args.fancy_plots == True:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)

            fig=pyplot.figure(figsize=(16,10))
            
            colors = ['k','tab:red','tab:blue','tab:green','tab:cyan']

            bins = np.linspace(0.0,1500.0,Nb)
            for index_kick in range(N_k):
                plot=fig.add_subplot(2,3,index_kick+1,yscale="log")
                color = 'k'

                plot.hist(vs[index_kick],bins=bins,histtype='step',density=True,color='k',linestyle='solid',label='$\mathrm{All}$')
                plot.hist(vs_BH[index_kick],bins=bins,histtype='step',density=True,color='tab:red',linestyle='dotted',label='$\mathrm{BH}$')
                plot.hist(vs_NS[index_kick],bins=bins,histtype='step',density=True,color='tab:blue',linestyle='dashed',label='$\mathrm{NS}$')
            
                plot.annotate("$\mathrm{Kick\,distribution\,%s}$"%(index_kick+1),xy=(0.1,0.9),xycoords='axes fraction',fontsize=fontsize)
                
                plot.set_xlabel("$V_\mathrm{kick}/(\mathrm{km\,s^{-1}})$",fontsize=fontsize)
                plot.set_ylabel("$\mathrm{PDF}$",fontsize=fontsize)
                    
                if 1==1:
                    points=np.linspace(10.0,1500.0,1000)
                    sigma_km_s = 265.0
                    PDF_an = np.sqrt(2.0/np.pi) * (points**2/(sigma_km_s**3)) * np.exp( -points**2/(2.0*sigma_km_s**2) )
                    plot.plot(points,PDF_an, color='tab:green',label='$\mathrm{Hobbs+05}$')
                handles,labels = plot.get_legend_handles_labels()
                plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
                
                plot.set_ylim(1.0e-5,1.0e-1)
                
            plot=fig.add_subplot(2,3,6,yscale="log")

            plot.legend(handles,labels,loc="upper left",fontsize=0.85*fontsize)
            plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            plot.axis('off')

            plot.set_xticks([])
            plot.set_yticks([])
            plot.set_xticklabels([])
            plot.set_yticklabels([])
            
            fig.savefig("kick_distributions.pdf")
            pyplot.show()

    def test18(self,args):
        print('Test flybys sampling')
        
        code = MSE()

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_KM_PER_S = code.CONST_KM_PER_S
        N = 20000
        
        seed = 0
        np.random.seed(seed)

        R_enc = 1.0e4
        n_star = 0.1*code.CONST_PER_PC3
        sigma_rel = 10.0*code.CONST_KM_PER_S
        M_int = 20.0

        b_vecs = []
        V_vecs = []
        M_pers = []
        bs = []
        Vs = []
        for index in range(N):
            M_per,b_vec,V_vec = code.test_flybys_perturber_sampling(R_enc,n_star,sigma_rel,M_int)
            M_pers.append(M_per)
            b_vecs.append(b_vec)
            V_vecs.append(V_vec)
            bs.append(np.linalg.norm(b_vec))
            Vs.append(np.linalg.norm(V_vec))

        bs = np.array(bs)
        Vs = np.array(Vs)
        
        b1 = 0.0
        b2 = R_enc
        b_mean_an = (2.0/3.0)*(b2**3 - b1**3)/(b2**2 - b1**2) ### mean b value assuming dN/db = 2b/(b2^2 - b1^2)
        
        tol = 1.0e-2
        assert( (b_mean_an - np.mean(bs))/b_mean_an <= tol)
        
        ### On average, individual components of b vec should be zero (isotropic orientations) ###
        assert( np.mean(np.array( [x[0] for x in b_vecs]))/b_mean_an <= tol ) 
        assert( np.mean(np.array( [x[1] for x in b_vecs]))/b_mean_an <= tol ) 
        assert( np.mean(np.array( [x[2] for x in b_vecs]))/b_mean_an <= tol ) 
        
        code.reset()

        if args.verbose==True:
            print("mean bs/au",np.mean(bs),' an ',b_mean_an)
        
        if args.plot == True:
            Nb=50
            fontsize=16
            labelsize=12
            from matplotlib import pyplot

            fig=pyplot.figure(figsize=(8,10))
            
            plot1=fig.add_subplot(2,1,1,yscale="log")
            plot2=fig.add_subplot(2,1,2,yscale="log")

            bins = np.linspace(2,np.log10(R_enc)+1,Nb)
            plot1.hist(np.log10(bs),bins=bins,histtype='step',density=True,color='k',linestyle='solid',label='$\mathrm{Num}$')
            bins = np.linspace(0,100.0,Nb)
            plot2.hist(Vs/code.CONST_KM_PER_S,bins=bins,histtype='step',density=True,color='k',linestyle='solid',label='$\mathrm{Num}$')
            
            points = np.linspace(2,np.log10(R_enc)+1,1000)
            plot1.plot(points, np.log(10.0)*2.0 * pow(10.0,2*points)/(b2**2 - b1**2), color='tab:green',label="$\mathrm{d}N/\mathrm{d}b \propto b$")

            plot1.set_xlabel("$\log_{10}(b/\mathrm{au})$",fontsize=fontsize)
            plot1.set_ylabel("$\mathrm{d} N/\mathrm{d} \log_{10} (b)$",fontsize=fontsize)

            plot2.set_xlabel("$V_\mathrm{enc}/(\mathrm{km \,s^{-1}}))$",fontsize=fontsize)
            plot2.set_ylabel("$\mathrm{d} N/\mathrm{d} V_\mathrm{enc}$",fontsize=fontsize)

            handles,labels = plot1.get_legend_handles_labels()
            plot1.legend(handles,labels,loc="best",fontsize=0.85*fontsize)
            
            plot1.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            plot2.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

            pyplot.show()

    def test19(self,args):
        print('Test white dwarf SNe single degenerate model 1')
        
        code = MSE()

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C

        luminosities = [0.03 * code.CONST_L_SUN,0.95 * code.CONST_L_SUN]

        ### He donor case ###
        #M_WDs = np.linspace(0.5,1.5,5)
        M_WDs = np.linspace(0.73,1.12,6)
        #M_WDs = [0.8,0.8,0.8,0.8,0.8,0.8]
        accretion_rates = pow(10.0,np.linspace(-9,-6,400))
        M_Hes = np.linspace(0.01,0.3,100)
        
        etas = [[[] for x in range(len(M_WDs))] for x in range(len(luminosities))]
        accretion_modes = [[[] for x in range(len(M_WDs))] for x in range(len(luminosities))]
        
        explosions = [[[[] for x in range(len(accretion_rates))] for x in range(len(M_WDs))] for x in range(len(luminosities))]
        all_etas = []
        for k,luminosity in enumerate(luminosities):
            for i,M_WD in enumerate(M_WDs):
                for j,m_dot in enumerate(accretion_rates):
                    eta,WD_accretion_mode = code.test_binary_evolution_SNe_Ia_single_degenerate_model_1_accumulation_efficiency(M_WD, m_dot, luminosity)
                    etas[k][i].append(eta)
                    all_etas.append(eta)
                    accretion_modes[k][i].append(WD_accretion_mode)


                    for l,M_He in enumerate(M_Hes):
                        explosion = code.test_binary_evolution_SNe_Ia_single_degenerate_model_1_explosion(M_WD, m_dot, M_He, luminosity)

                        explosions[k][i][j].append(int(explosion))

        ### Hydrogen donor case ###
        H_donor_m_dot_lower_boundaries = []
        H_donor_m_dot_upper_boundaries = []
        
        M_WDs2 = np.linspace(0.5,1.4,100)
        for i,M_WD in enumerate(M_WDs2):
            m_dot_lower,m_dot_upper = code.test_binary_evolution_SNe_Ia_single_degenerate_model_1_white_dwarf_hydrogen_accretion_boundaries(M_WD)
            H_donor_m_dot_lower_boundaries.append(m_dot_lower)
            H_donor_m_dot_upper_boundaries.append(m_dot_upper)

        ### Sanity checks ###
        all_etas = np.array(all_etas)
        assert( np.amin(all_etas) >= 0.0 and np.amax(all_etas) <= 1.0)
        
        if args.verbose==True:
            print("np.amin(all_etas)",np.amin(all_etas),"np.amax(all_etas)",np.amax(all_etas),"np.mean(all_etas)",np.mean(all_etas))

        code.reset()

        if args.plot == True:
            Nb=50
            fontsize=16
            labelsize=12
            from matplotlib import pyplot
            if args.fancy_plots == True:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)
                
            fig=pyplot.figure(figsize=(8,10))
            fige=pyplot.figure(figsize=(12,12))
            figh=pyplot.figure(figsize=(8,6))
            
            colors = ['k','tab:blue','tab:red','tab:green','tab:cyan','tab:brown']
            linestyles = ['solid','dotted','dashed','-.','solid','dotted']
            
            plot1=fig.add_subplot(2,1,1,xscale="log")
            plot2=fig.add_subplot(2,1,2,xscale="log")
            
            N_r=2
            plot1e=fige.add_subplot(N_r,N_r,1,xscale="linear",yscale="log")
            plot2e=fige.add_subplot(N_r,N_r,2,xscale="linear",yscale="log")
            plot3e=fige.add_subplot(N_r,N_r,3,xscale="linear",yscale="log")
            plot4e=fige.add_subplot(N_r,N_r,4,xscale="linear",yscale="log")
            
            plot1h=figh.add_subplot(1,1,1,xscale="linear",yscale="log")

            linewidth=1.5
            for i,M_WD in enumerate(M_WDs):
                plot1.plot(accretion_rates,etas[0][i],label="$M_\mathrm{WD}=%s \, \mathrm{M}_\odot$"%round(M_WD,2),color=colors[i],linestyle=linestyles[i],linewidth=linewidth)
                plot2.plot(accretion_rates,etas[1][i],label="$M_\mathrm{WD}=%s \, \mathrm{M}_\odot$"%round(M_WD,2),color=colors[i],linestyle=linestyles[i],linewidth=linewidth)
            
            plot1h.plot(M_WDs2,H_donor_m_dot_lower_boundaries,color='tab:red')
            plot1h.plot(M_WDs2,H_donor_m_dot_upper_boundaries,color='tab:red')
            
            xs = M_Hes
            ys = accretion_rates
            
            i=1
            j=0
            zs = explosions[j][i]
            plot1e.pcolormesh(xs,ys,zs,cmap = 'Reds')
            plot1e.set_title("$M_\mathrm{WD} = %s\,\mathrm{M}_\odot; \, L=%s\,\mathrm{L}_\odot$"%(round(M_WDs[i],2),round(luminosities[j]/code.CONST_L_SUN,2)))

            i=4
            j=0
            zs = explosions[j][i]
            plot2e.pcolormesh(xs,ys,zs,cmap = 'Reds')
            plot2e.set_title("$M_\mathrm{WD} = %s\,\mathrm{M}_\odot; \, L=%s\,\mathrm{L}_\odot$"%(round(M_WDs[i],2),round(luminosities[j]/code.CONST_L_SUN,2)))

            i=1
            j=1
            zs = explosions[j][i]
            plot3e.pcolormesh(xs,ys,zs,cmap = 'Reds')
            plot3e.set_title("$M_\mathrm{WD} = %s\,\mathrm{M}_\odot; \, L=%s\,\mathrm{L}_\odot$"%(round(M_WDs[i],2),round(luminosities[j]/code.CONST_L_SUN,2)))

            i=4
            j=1
            zs = explosions[j][i]
            plot4e.pcolormesh(xs,ys,zs,cmap = 'Reds')
            plot4e.set_title("$M_\mathrm{WD} = %s\,\mathrm{M}_\odot; \, L=%s\,\mathrm{L}_\odot$"%(round(M_WDs[i],2),round(luminosities[j]/code.CONST_L_SUN,2)))

            plot2.set_xlabel("$\log_{10}[\dot{M}/(\mathrm{M}_\odot/\mathrm{yr^{-1}})]$",fontsize=fontsize)
            plots = [plot1,plot2]
            for i,plot in enumerate(plots):
                plot.set_ylabel("$\eta$",fontsize=fontsize)

                handles,labels = plot.get_legend_handles_labels()
                plot.legend(handles,labels,loc="best",fontsize=0.85*fontsize)

                plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
                plot.set_title("$L=%s \, \mathrm{L}_\odot$"%(round(luminosities[i]/code.CONST_L_SUN,2)))


            plots = [plot1e,plot2e,plot3e,plot4e]
            for i,plot in enumerate(plots):
                plot.set_ylabel("$\log_{10}[\dot{M}/(\mathrm{M}_\odot/\mathrm{yr^{-1}})]$",fontsize=fontsize)

                plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            plot3e.set_xlabel("$M_\mathrm{He\,layer}/\mathrm{M}_\odot$",fontsize=fontsize)
            plot4e.set_xlabel("$M_\mathrm{He\,layer}/\mathrm{M}_\odot$",fontsize=fontsize)

            plot1h.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            plot1h.set_xlabel("$M_\mathrm{WD}/\mathrm{M}_\odot$",fontsize=fontsize)
            plot1h.set_ylabel("$\dot{M}/(\mathrm{M}_\odot/\mathrm{yr})$",fontsize=fontsize)

            plot1h.set_ylim(1e-8,1e-6)

            fig.savefig("white_dwarf_SNe_single_degenerate_model_1_eta.pdf")
            fige.savefig("white_dwarf_SNe_single_degenerate_model_1_explosion.pdf")
            figh.savefig("white_dwarf_SNe_single_degenerate_model_1_hydrogen_m_dot.pdf")
    
            pyplot.show()

    def test20(self,args):
        print("Test Helium star system")

        #particles = Tools.create_fully_nested_multiple(2,[5.958893384591632, 3.574547072971755], [0.03], [0.0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[11,7],object_types=[1,1])
        particles = Tools.create_fully_nested_multiple(2,[5.958893384591632, 1.0], [0.0043], [0.0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[11,7],object_types=[1,1])
        
        particles = Tools.create_fully_nested_multiple(2,[5.958893384591632, 1.0], [0.0015], [0.0], [0.01], [0.01], [0.01], metallicities=[0.02,0.02],stellar_types=[11,7],object_types=[1,1])
        
        binaries = [x for x in particles if x.is_binary == True]
        bodies = [x for x in particles if x.is_binary == False]

        code = MSE()
        code.add_particles(particles)

        t = 0.0
        N=1000
        tend = 1.0e9
        dt=tend/float(N)

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        code.evolve_model(0.0) ### make sure that stellar_evolution.cpp -- initialize_stars() is run
        code.particles[1].stellar_type = 10
        code.particles[0].age = 1.60722e+08
        code.particles[1].age = 0.0

        code.binary_evolution_SNe_Ia_single_degenerate_model = 0
        code.binary_evolution_SNe_Ia_double_degenerate_model = 0

        code.effective_radius_multiplication_factor_for_collisions_compact_objects = 1.0e0

        N_bodies = 2
        N_orbits = 1
        t_print = [[]]
        internal_indices_print = [[[] for x in range(N_bodies)]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
        L_print = [[[] for x in range(N_bodies)]]
        T_eff_print = [[[] for x in range(N_bodies)]]
        R_print = [[[] for x in range(N_bodies)]]
        X_print = [[[] for x in range(N_bodies)]]
        Y_print = [[[] for x in range(N_bodies)]]
        Z_print = [[[] for x in range(N_bodies)]]
        Rc_print = [[[] for x in range(N_bodies)]]
        R_L_print = [[[] for x in range(N_bodies)]]
        t_V_print = [[[] for x in range(N_bodies)]]
        mc_print = [[[] for x in range(N_bodies)]]
        a_print = [[[] for x in range(N_orbits)]]
        e_print = [[[] for x in range(N_orbits)]]
        rel_INCL_print = [[[] for x in range(N_orbits)]]
        spin_frequency_print = [[[] for x in range(N_bodies)]]
        
        i_status = 0
        integration_flags = [[]]
        code.verbose_flag = 1
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            
            #particles = code.particles
            #binaries = [x for x in particles if x.is_binary == True]
            #bodies = [x for x in particles if x.is_binary == False]


            particles = code.particles
            orbits = [x for x in particles if x.is_binary==True]
            bodies = [x for x in particles if x.is_binary==False]
            N_orbits = len(orbits)
            N_bodies = len(bodies)

            if args.verbose==True:
                print( 't/Myr',t,'smas',[x.a for x in orbits],'stellar types',[x.stellar_type for x in bodies],'masses',[x.mass for x in bodies],'ages',[x.age for x in bodies],'WD_He_layer_mass',[x.WD_He_layer_mass for x in bodies],'m_dot_accretion_SD',[x.m_dot_accretion_SD for x in bodies])
              
            if code.structure_change == True:
                #print("Python restruct")#,children1,children1_old,children2,children2_old)
                t_print.append([])
                integration_flags.append([])
                internal_indices_print.append([[] for x in range(N_bodies)])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
                L_print.append([[] for x in range(N_bodies)])
                T_eff_print.append([[] for x in range(N_bodies)])
                R_print.append([[] for x in range(N_bodies)])
                X_print.append([[] for x in range(N_bodies)])
                Y_print.append([[] for x in range(N_bodies)])
                Z_print.append([[] for x in range(N_bodies)])
                Rc_print.append([[] for x in range(N_bodies)])
                R_L_print.append([[] for x in range(N_bodies)])
                t_V_print.append([[] for x in range(N_bodies)])
                mc_print.append([[] for x in range(N_bodies)])
                a_print.append([[] for x in range(N_orbits)])
                e_print.append([[] for x in range(N_orbits)])
                rel_INCL_print.append([[] for x in range(N_orbits)])
                spin_frequency_print.append([[] for x in range(N_bodies)])
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
                
            for index in range(N_orbits):
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                internal_indices_print[i_status][index].append(bodies[index].index)
                m_print[i_status][index].append(bodies[index].mass)
                L_print[i_status][index].append(bodies[index].luminosity)
                k_print[i_status][index].append(bodies[index].stellar_type)
                R_print[i_status][index].append(bodies[index].radius)
                X_print[i_status][index].append(bodies[index].X)
                Y_print[i_status][index].append(bodies[index].Y)
                Z_print[i_status][index].append(bodies[index].Z)
                t_V_print[i_status][index].append(bodies[index].tides_viscous_time_scale)
                Rc_print[i_status][index].append(bodies[index].convective_envelope_radius)
                R_L_print[i_status][index].append(bodies[index].roche_lobe_radius_pericenter)
                mc_print[i_status][index].append(bodies[index].convective_envelope_mass)
                T_eff = Tools.compute_effective_temperature(bodies[index].luminosity, bodies[index].radius, code.CONST_L_SUN, code.CONST_R_SUN)
                T_eff_print[i_status][index].append(T_eff)
                spin_frequency_print[i_status][index].append( np.sqrt( bodies[index].spin_vec_x**2 + bodies[index].spin_vec_y**2 + bodies[index].spin_vec_z**2) )

            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)        

        code.write_final_log_entry() ### This has to be done within Python, since the C++ code does not know if the desired Python simulation end time has been reached!
        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])

        import copy
        error_code_copy = copy.deepcopy(code.error_code)
        log_copy = copy.deepcopy(code.log)

        code.reset()
                
        if HAS_MATPLOTLIB == True and args.plot==True:
            from matplotlib import lines
            
            parsec_in_AU = 206201.0
            CONST_L_SUN = 0.0002710404109745588

            plot_log = []
            previous_event_flag = -1
            for index_log,log in enumerate(log_copy):
                event_flag = log["event_flag"]
                if previous_event_flag == event_flag and (event_flag == 4 or event_flag == 10):
                    
                    continue
                plot_log.append(log)
                previous_event_flag = event_flag
                            
            N_l = len(plot_log)
            N_r = int(np.ceil(np.sqrt(N_l)))#+1
            N_c = N_r
            panel_length = 3

            fontsize=10
            fig=pyplot.figure(figsize=(N_r*panel_length,N_r*panel_length))
            
            legend_elements = []
            for k in range(15):
                color,s,description = Tools.get_color_and_size_and_description_for_star(k,1.0)
                legend_elements.append( lines.Line2D([0],[0], marker = 'o', markerfacecolor = color, color = 'w', markersize = 10 ,label="$\mathrm{%s}$"%description))

            for index_log,log in enumerate(plot_log):
                plot=fig.add_subplot(N_r,N_c,index_log+1)
                particles = log["particles"]
                event_flag = log["event_flag"]
                index1 = log["index1"]
                index2 = log["index2"]

                Tools.generate_mobile_diagram(particles,plot,fontsize=fontsize,index1=index1,index2=index2,event_flag=event_flag)

                text = Tools.get_description_for_event_flag(event_flag,log["SNe_type"])
                plot.set_title(text,fontsize=fontsize)
                plot.annotate("$t\simeq %s\,\mathrm{Myr}$"%round(log["time"]*1e-6,2),xy=(0.1,0.9),xycoords='axes fraction',fontsize=fontsize)

                if index_log == 0:
                    plot.legend(handles = legend_elements, bbox_to_anchor = (-0.05, 1.50), loc = 'upper left', ncol = 5,fontsize=0.85*fontsize)
                
            #fig.savefig(plot_filename + "_mobile.pdf")
        
            fig=pyplot.figure(figsize=(8,10))
            Np=4
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            plot4=fig.add_subplot(Np,1,4,yscale="log")
            
            fig_pos=pyplot.figure(figsize=(8,8))
            plot_pos=fig_pos.add_subplot(1,1,1)

            fig_HRD=pyplot.figure(figsize=(8,8))
            plot_HRD=fig_HRD.add_subplot(1,1,1)
            
            colors = ['k','tab:red','tab:green','tab:blue','y','k','tab:red','tab:green','tab:blue','y']
            linewidth=1.0
            for i_status in range(N_status):
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]
                
                #if i_status==0:
                #    plot1.plot(1.0e-6*t_print[i_status],np.array(m_print[i_status][0])**2*np.array(m_print[i_status][1])**2*np.array(a_print[i_status][0]),color='y',linewidth=3)
                #    plot1.plot(1.0e-6*t_print[i_status],(np.array(m_print[i_status][0])+np.array(m_print[i_status][1]))*np.array(a_print[i_status][0]),color='y',linewidth=3)
                for index in range(N_bodies):
                    color=colors[index]
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    
                    plot_pos.plot(np.array(X_print[i_status][index])/parsec_in_AU,np.array(Y_print[i_status][index])/parsec_in_AU,color=color,linestyle='solid',linewidth=linewidth)
                    
                    plot_HRD.scatter(np.log10(np.array(T_eff_print[i_status][index])), np.log10(np.array(L_print[i_status][index])/CONST_L_SUN),color=color,linestyle='solid',linewidth=linewidth,s=10)
                    #plot4.plot(1.0e-6*t_print[i_status],spin_frequency_print[i_status][index],color=color,linewidth=linewidth) 
                    plot4.plot(1.0e-6*t_print[i_status],3.15576e7*2.0*np.pi/np.array(spin_frequency_print[i_status][index]),color=color,linewidth=linewidth) 
                   
                linewidth=1.0
                for index in range(N_orbits):
                    color = colors[index]
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
                
                linewidth+=0.4

            fontsize=18
            labelsize=12

            plots = [plot1,plot2,plot3,plot4]
            for plot in plots:
                plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

            log_CEs = [x for x in plot_log if x["event_flag"] == 6]
            t_CEs_Myr = np.array([x["time"]*1e-6 for x in log_CEs])
            
            for k,t in enumerate(t_CEs_Myr):
                plot2.axvline(x=t,linestyle='dotted',color='tab:red',linewidth=0.5)
                plot2.annotate("$\mathrm{CE}$",xy=(1.02*t,1.0e3),fontsize=0.8*fontsize)
            
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{au}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            #plot4.set_ylabel("$\Omega_\mathrm{spin}/\mathrm{yr^{-1}}$",fontsize=fontsize)
            plot4.set_ylabel("$P_\mathrm{spin}/\mathrm{s}$",fontsize=fontsize)
            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            
            plot_pos.set_xlabel("$X/\mathrm{pc}$",fontsize=fontsize)
            plot_pos.set_ylabel("$Y/\mathrm{pc}$",fontsize=fontsize)
            plot_pos.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            

            plot_HRD.set_xlim(5.0,3.0)
            plot_HRD.set_ylim(-4.0,6.0)
            plot_HRD.set_xlabel("$\mathrm{log}_{10}(T_\mathrm{eff}/\mathrm{K})$",fontsize=fontsize)
            plot_HRD.set_ylabel(r"$\mathrm{log}_{10}(L/L_\odot)$",fontsize=fontsize)
            plot_HRD.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            
            #fig.savefig(plot_filename + ".pdf")
            #fig_pos.savefig(plot_filename + "_pos.pdf")
            #fig_HRD.savefig(plot_filename + "_HRD.pdf")
            
            pyplot.show()


    def test100(self,args):
        print('Unit tests')
        
        code = MSE()

        flag = code.unit_tests(args.mode)
        
        assert(flag == 0)
        
        print("Unit tests passed")
        

def kroupa_93_imf(m):

    alpha1 = -1.3
    alpha2 = -2.2
    alpha3 = -2.7
    m1 = 0.1
    m2 = 0.5
    m3 = 1.0
    m4 = 100.0
    C1 = 0.2905673356704877
    C2 = 0.155711179725752
    C3 = 0.155711179725752
    
    if (m>=m1 and m<m2):
        pdf = C1*pow(m,alpha1)
    elif (m>=m2 and m<m3):
        pdf = C2*pow(m,alpha2)
    elif (m>=m3 and m<=m4):
        pdf = C3*pow(m,alpha3)
    else:
        pdf = 0.0
        
    return pdf
    
def sample_random_vector_on_unit_sphere():
    INCL = np.arccos( 2.0*np.random.random() - 1.0)
    LAN = 2.0 * np.pi * np.random.random()
    return compute_unit_AM_vector(INCL,LAN)

def compute_unit_AM_vector(INCL,LAN):
    return np.array( [np.sin(LAN)*np.sin(INCL), -np.cos(LAN)*np.sin(INCL), np.cos(INCL)] )

def compute_e_and_j_hat_vectors(INCL,AP,LAN):
    sin_INCL = np.sin(INCL)
    cos_INCL = np.cos(INCL)
    sin_AP = np.sin(AP)
    cos_AP = np.cos(AP)
    sin_LAN = np.sin(LAN)
    cos_LAN = np.cos(LAN)
    
    e_hat_vec_x = (cos_LAN*cos_AP - sin_LAN*sin_AP*cos_INCL);
    e_hat_vec_y = (sin_LAN*cos_AP + cos_LAN*sin_AP*cos_INCL);
    e_hat_vec_z = (sin_AP*sin_INCL);
    
    j_hat_vec_x = sin_LAN*sin_INCL;
    j_hat_vec_y = -cos_LAN*sin_INCL;
    j_hat_vec_z = cos_INCL;

    e_hat_vec = np.array([e_hat_vec_x,e_hat_vec_y,e_hat_vec_z])
    j_hat_vec = np.array([j_hat_vec_x,j_hat_vec_y,j_hat_vec_z])

    return e_hat_vec,j_hat_vec
    
if __name__ == '__main__':
    args = parse_arguments()
    
    N_tests = 19
    if args.test==0:
        tests = list(range(1,N_tests+1)) + [100]
    else:
        tests = [args.test]

    t=test_mse()
    for i in tests:
        print( 'Running test number',i,'; verbose =',args.verbose,'; plot =',args.plot)
        function = getattr(t, 'test%s'%i)
        function(args)
    
    print("="*50)
    print("All tests passed!")
