"""
Some examples illustrating the usage of MSE
Adrian Hamers, August 2020
"""

import numpy as np
import numpy.random as randomf
import argparse

from mse import MSE,Particle,Tools


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
    
    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                         default=False,         help="verbose terminal output")
    add_bool_arg(parser, 'plot',                            default=False,         help="make plots")
    
    args = parser.parse_args()

    return args
    
class test_mse():

    def test1(self,args):
        print('Triple with stellar evolution')
        N_bodies = 3
        configuration="fully nested"
        masses = [1.0,1.0,20.5]
        #masses = [2.0,1.0,1.5]
        metallicities = [0.01,0.03,0.005]
        #semimajor_axes = [15.5,400.0]
        semimajor_axes = [15.5,120.0]
        eccentricities = [0.1,0.1]
        inclinations = [0.0001,15.0*np.pi/180.0]
        arguments_of_pericentre = [85.0*np.pi/180.0,0.01*np.pi/180.0]
        longitudes_of_ascending_node = [0.01,0.01]
        end_time = 3.0e7
        #end_time = 1.0e9
        N_steps = 500
        stellar_types = [1,1,1]
        Tools.evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,end_time,N_steps,stellar_types=stellar_types)

    def test2(self,args):
        print('Triple with stellar evolution')
        N_bodies = 3
        configuration="fully nested"
        masses = [50.0,40.0,7.5]
        #masses = [2.0,1.0,1.5]
        metallicities = [0.01,0.03,0.005]
        #semimajor_axes = [15.5,400.0]
        semimajor_axes = [0.01,400.0]
        eccentricities = [0.1,0.6]
        inclinations = [0.0001,85.0*np.pi/180.0]
        arguments_of_pericentre = [85.0*np.pi/180.0,0.01*np.pi/180.0]
        longitudes_of_ascending_node = [0.01,0.01]
        end_time = 5.0e7
        #end_time = 1.0e9
        N_steps = 1000
        stellar_types = [13,13,13]
        Tools.evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,end_time,N_steps,stellar_types=stellar_types)
        
    def test2b(self,args):
        print('Quadruple with stellar evolution')
        
        N_bodies=4
        
        #particles = Tools.create_nested_multiple(N_bodies, [40.0,14.8,8.5],[30.0,600.0],[0.2,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        
        #particles = Tools.create_nested_multiple(N_bodies, [34.0,25.8,8.5],[30.0,500.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        
        #particles = Tools.create_nested_multiple(N_bodies, [24.0,6.0,7.5],[15.5,600.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[85.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01]) ### promising
        #particles = Tools.create_nested_multiple(N_bodies, [11.0,1.0,7.5],[15.5,400.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[85.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01]) ### promising
        
        
        #particles = Tools.create_nested_multiple(N_bodies, [22.0,6.0,7.5],[12.5,600.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[15.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        particles = Tools.create_nested_multiple(N_bodies, [7.0,4.6,3.5,2.5],[10.0,1000.0,12000.0],[0.1,0.3,0.3],[0.0001,51.0*np.pi/180.0,123.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])

#        fig=pyplot.figure(figsize=(8,6))
#        plot=fig.add_subplot(1,1,1)
#        Tools.generate_mobile_diagram(particles,plot)
#        pyplot.show()

        #particles = Tools.create_nested_multiple(N_bodies, [4.0,2.8,1.5],[3000.0,400000.0],[0.1,0.3],[0.0001,89.9*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [6.0,4.6,3.5,2.5],[10.0,1000.0,12000.0],[0.1,0.3,0.3],[0.0001,51.0*np.pi/180.0,123.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
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


                
        code = MSE()
        code.enable_root_finding = True
        code.add_particles(particles)
        code.orbital_phases_random_seed = 1
        
        code.include_flybys = True
        #code.flybys_reference_binary = orbits[0].index
        #print 'a',orbits[0].a
        code.flybys_stellar_density = 0.1*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e5


        for b in bodies:
            b.kick_distribution_sigma = 0*265.0*code.CONST_KM_PER_S
       #     b.metallicity = 0.02
       #     b.include_tidal_friction_terms = True



        t_print = [[]]
        internal_indices_print = [[[] for x in range(N_bodies)]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
        R_print = [[[] for x in range(N_bodies)]]
        Rc_print = [[[] for x in range(N_bodies)]]
        R_L_print = [[[] for x in range(N_bodies)]]
        t_V_print = [[[] for x in range(N_bodies)]]
        mc_print = [[[] for x in range(N_bodies)]]
        a_print = [[[] for x in range(N_orbits)]]
        e_print = [[[] for x in range(N_orbits)]]
        rel_INCL_print = [[[] for x in range(N_orbits)]]
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        
        t = 0.0
        integration_flags = [[]]
        
        N = 10
        tend = 1e10

        N =1000
        tend = 5e8
        #tend = 1.4e10

        i_status = 0
        #print("rel_INCL_print",rel_INCL_print[0][2])
        seed=1
        dt = tend/float(N)
        i = 0
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state

            particles = code.particles
            orbits = [x for x in particles if x.is_binary==True]
            bodies = [x for x in particles if x.is_binary==False]
            N_orbits = len(orbits)
            N_bodies = len(bodies)

            if code.structure_change == True:
                print("Python restruct")#,children1,children1_old,children2,children2_old)
                t_print.append([])
                integration_flags.append([])
                internal_indices_print.append([[] for x in range(N_bodies)])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
                R_print.append([[] for x in range(N_bodies)])
                Rc_print.append([[] for x in range(N_bodies)])
                R_L_print.append([[] for x in range(N_bodies)])
                t_V_print.append([[] for x in range(N_bodies)])
                mc_print.append([[] for x in range(N_bodies)])
                a_print.append([[] for x in range(N_orbits)])
                e_print.append([[] for x in range(N_orbits)])
                rel_INCL_print.append([[] for x in range(N_orbits)])
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
                
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag,'integration_flag',code.integration_flag,'i_status',i_status)
           # print("bound",[x.is_bound for x in bodies])
            
            for index in range(N_orbits):
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                internal_indices_print[i_status][index].append(particles[index].index)
                #print("T",index,internal_index)
                m_print[i_status][index].append(particles[index].mass)
                k_print[i_status][index].append(particles[index].stellar_type)
                R_print[i_status][index].append(particles[index].radius)
                t_V_print[i_status][index].append(particles[index].tides_viscous_time_scale)
                Rc_print[i_status][index].append(particles[index].convective_envelope_radius)
                R_L_print[i_status][index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[i_status][index].append(particles[index].convective_envelope_mass)

            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)
            if flag==2:
                print("Root",t,'RLOF',[o.RLOF_at_pericentre_has_occurred for o in orbits],'col',[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits],'dyn inst',[o.dynamical_instability_has_occurred for o in orbits],'sec break',[o.secular_breakdown_has_occurred for o in orbits],'min peri',[o.minimum_periapse_distance_has_occurred for o in orbits],'GW',[o.GW_condition_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1
            
            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])

        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            
            colors = ['k','tab:red','tab:green','tab:blue','y','k','tab:red','tab:green','tab:blue','y']
            linewidth=1.0
            for i_status in range(N_status):
                #color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]

                
                
                
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color='k',linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    internal_index = internal_indices_print[i_status][index][0]
                    color = colors[internal_index]
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    #if index in [0,1]:
                        #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    #linewidth+=0.8
                    
                linewidth=1.0
                for index in range(N_orbits):
                    color = colors[index]
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,(180.0/np.pi)*rel_INCL_print[index],color='k',linestyle='solid',linewidth=linewidth)
                
                linewidth+=0.8



            fontsize=18
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            plot3.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            fig.savefig("test1.pdf")
            pyplot.show()


    def test2b(self,args):
        print('Dynamically unstable quadruple')
        
        N_bodies=4
        
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5],[30.0,400.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [4.0,2.8,1.5],[3000.0,400000.0],[0.1,0.3],[0.0001,89.9*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        particles = Tools.create_nested_multiple(N_bodies, [7.0,4.6,3.5,2.5],[10.0,1000.0,12000.0],[0.1,0.3,0.3],[0.0001,51.0*np.pi/180.0,123.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
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


        t_print = [[]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
        R_print = [[[] for x in range(N_bodies)]]
        Rc_print = [[[] for x in range(N_bodies)]]
        R_L_print = [[[] for x in range(N_bodies)]]
        t_V_print = [[[] for x in range(N_bodies)]]
        mc_print = [[[] for x in range(N_bodies)]]
        a_print = [[[] for x in range(N_orbits)]]
        e_print = [[[] for x in range(N_orbits)]]
        rel_INCL_print = [[[] for x in range(N_orbits)]]
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        
        t = 0.0
        integration_flags = [[]]

#        if 1==0:
#            for index in range(N_bodies):
#                m_print[index].append(particles[index].mass)
#                k_print[index].append(particles[index].stellar_type)
#                R_print[index].append(particles[index].radius)
#                t_V_print[index].append(particles[index].tides_viscous_time_scale)
#                Rc_print[index].append(particles[index].convective_envelope_radius)
#                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
        #t_print.append(t)
        N = 10
        tend = 1e10

        N =1000
        tend = 1e9
        #tend = 1e8

        i_status = 0
        #print("rel_INCL_print",rel_INCL_print[0][2])
        seed=1
        dt = tend/float(N)
        i = 0
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state

            particles = code.particles
            orbits = [x for x in particles if x.is_binary==True]
            bodies = [x for x in particles if x.is_binary==False]
            N_orbits = len(orbits)
            N_bodies = len(bodies)

            if code.structure_change == True:
                print("Python restruct")#,children1,children1_old,children2,children2_old)
                t_print.append([])
                integration_flags.append([])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
                R_print.append([[] for x in range(N_bodies)])
                Rc_print.append([[] for x in range(N_bodies)])
                R_L_print.append([[] for x in range(N_bodies)])
                t_V_print.append([[] for x in range(N_bodies)])
                mc_print.append([[] for x in range(N_bodies)])
                a_print.append([[] for x in range(N_orbits)])
                e_print.append([[] for x in range(N_orbits)])
                rel_INCL_print.append([[] for x in range(N_orbits)])
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
                
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag,'integration_flag',code.integration_flag,'i_status',i_status)
            print("bound",[x.is_bound for x in bodies])
            
            for index in range(N_orbits):
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[i_status][index].append(particles[index].mass)
                k_print[i_status][index].append(particles[index].stellar_type)
                R_print[i_status][index].append(particles[index].radius)
                t_V_print[i_status][index].append(particles[index].tides_viscous_time_scale)
                Rc_print[i_status][index].append(particles[index].convective_envelope_radius)
                R_L_print[i_status][index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[i_status][index].append(particles[index].convective_envelope_mass)

            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)
            if flag==2:
                print("Root",t,'RLOF',[o.RLOF_at_pericentre_has_occurred for o in orbits],'col',[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits],'dyn inst',[o.dynamical_instability_has_occurred for o in orbits],'sec break',[o.secular_breakdown_has_occurred for o in orbits],'min peri',[o.minimum_periapse_distance_has_occurred for o in orbits],'GW',[o.GW_condition_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1
            
            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])

        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])
        
        #if HAS_MATPLOTLIB==True and args.plot==True:
        if 1==0:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            
            colors = ['k','tab:red','tab:green','tab:blue','y','k','tab:red','tab:green','tab:blue','y']
            for i_status in range(N_status):
                color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]

                
                
                linewidth=1.0
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color=color,linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    #if index in [0,1]:
                        #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    linewidth+=0.8
                    
                linewidth=1.0
                for index in range(N_orbits):
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,(180.0/np.pi)*rel_INCL_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    linewidth+=0.8

            fontsize=18
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            plot3.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            fig.savefig("test2.pdf")
            pyplot.show()
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            
            colors = ['k','tab:red','tab:green','tab:blue','y','k','tab:red','tab:green','tab:blue','y']
            linewidth=1.0
            for i_status in range(N_status):
                #color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]

                
                
                
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color='k',linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    internal_index = internal_indices_print[i_status][index][0]
                    color = colors[internal_index]
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    #if index in [0,1]:
                        #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    #linewidth+=0.8
                    
                linewidth=1.0
                for index in range(N_orbits):
                    color = colors[index]
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,(180.0/np.pi)*rel_INCL_print[index],color='k',linestyle='solid',linewidth=linewidth)
                
                linewidth+=0.8



            fontsize=18
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{AU}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            plot3.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            fig.savefig("case2.pdf")
            pyplot.show()
        
    def test16_old(self,args):
        print('Dynamically unstable quadruple')
        
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



#        INCL_print = []
        t_print = [[]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
        R_print = [[[] for x in range(N_bodies)]]
        Rc_print = [[[] for x in range(N_bodies)]]
        R_L_print = [[[] for x in range(N_bodies)]]
        t_V_print = [[[] for x in range(N_bodies)]]
        mc_print = [[[] for x in range(N_bodies)]]
        a_print = [[[] for x in range(N_orbits)]]
        e_print = [[[] for x in range(N_orbits)]]
        rel_INCL_print = [[[] for x in range(N_orbits)]]
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        
        t = 0.0
        integration_flags = [[]]

#        if 1==0:
#            for index in range(N_bodies):
#                m_print[index].append(particles[index].mass)
#                k_print[index].append(particles[index].stellar_type)
#                R_print[index].append(particles[index].radius)
#                t_V_print[index].append(particles[index].tides_viscous_time_scale)
#                Rc_print[index].append(particles[index].convective_envelope_radius)
#                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
        #t_print.append(t)
        N = 10
        tend = 1e10

        N =1000
        tend = 1e9
        #tend = 1e8

        i_status = 0
        #print("rel_INCL_print",rel_INCL_print[0][2])
        seed=1
        dt = tend/float(N)
        i = 0
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state
        
            if state != 0:
                #print("Restructuring of system!")
                particles = code.particles
                orbits = [x for x in particles if x.is_binary==True]
                bodies = [x for x in particles if x.is_binary==False]
                N_orbits = len(orbits)
                N_bodies = len(bodies)
                
                t_print.append([])
                integration_flags.append([])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
                R_print.append([[] for x in range(N_bodies)])
                Rc_print.append([[] for x in range(N_bodies)])
                R_L_print.append([[] for x in range(N_bodies)])
                t_V_print.append([[] for x in range(N_bodies)])
                mc_print.append([[] for x in range(N_bodies)])
                a_print.append([[] for x in range(N_orbits)])
                e_print.append([[] for x in range(N_orbits)])
                rel_INCL_print.append([[] for x in range(N_orbits)])
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
            #if t > 50: i_status += 1
                
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag,'integration_flag',code.integration_flag,'i_status',i_status)
            print("bound",[x.is_bound for x in bodies])
            
            if state != 0:
                print("Restructuring of system!")
                
                #exit(0)

            for index in range(N_orbits):
                print(i_status,index)
                #print("in ",[o.index for o in orbits])
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[i_status][index].append(particles[index].mass)
                k_print[i_status][index].append(particles[index].stellar_type)
                R_print[i_status][index].append(particles[index].radius)
                t_V_print[i_status][index].append(particles[index].tides_viscous_time_scale)
                Rc_print[i_status][index].append(particles[index].convective_envelope_radius)
                R_L_print[i_status][index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[i_status][index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)
            if flag==2:
                print("Root",t,'RLOF',[o.RLOF_at_pericentre_has_occurred for o in orbits],'col',[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits],'dyn inst',[o.dynamical_instability_has_occurred for o in orbits],'sec break',[o.secular_breakdown_has_occurred for o in orbits],'min peri',[o.minimum_periapse_distance_has_occurred for o in orbits],'GW',[o.GW_condition_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1
            
            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])
        #for index in range(N_bodies):
        #        m_print[index] = np.array(m_print[index])
        #        t_V_print[index] = np.array(t_V_print[index])
        #for index in range(N_orbits):

#            rel_INCL_print[index] = np.array(rel_INCL_print[index])
#            e_print[index] = np.array(e_print[index])
#            a_print[index] = np.array(a_print[index])

        #rel_INCL_print = np.array(rel_INCL_print)
        #e_print = np.array(e_print)
        
        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            
            
            colors = ['k','tab:red','tab:green','tab:blue','y']
            for i_status in range(N_status):
                color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]

                
                
                linewidth=1.0
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color=color,linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    #if index in [0,1]:
                        #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    linewidth+=0.8
                    
                linewidth=1.0
                for index in range(N_orbits):
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
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
            fig.savefig("test16.pdf")
            pyplot.show()

    def test3(self,args):
        print('Mass transfer')
        
        N_bodies=3
        
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5],[30.0,400.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [4.0,2.8,1.5],[3000.0,400000.0],[0.1,0.3],[0.0001,89.9*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        
        #particles = Tools.create_nested_multiple(N_bodies, [6.0,4.6,3.5,2.5],[10.0,1000.0,12000.0],[0.1,0.3,0.3],[0.0001,51.0*np.pi/180.0,123.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01,0.01])
        particles = Tools.create_nested_multiple(N_bodies, [10.0,5.0,10.0],[6.0,1.0e2],[0.1,0.3],[0.0001,59.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01]) 
        
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



#        INCL_print = []
        t_print = [[]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
        R_print = [[[] for x in range(N_bodies)]]
        Rc_print = [[[] for x in range(N_bodies)]]
        R_L_print = [[[] for x in range(N_bodies)]]
        t_V_print = [[[] for x in range(N_bodies)]]
        mc_print = [[[] for x in range(N_bodies)]]
        a_print = [[[] for x in range(N_orbits)]]
        e_print = [[[] for x in range(N_orbits)]]
        rel_INCL_print = [[[] for x in range(N_orbits)]]
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        
        t = 0.0
        integration_flags = [[]]

#        if 1==0:
#            for index in range(N_bodies):
#                m_print[index].append(particles[index].mass)
#                k_print[index].append(particles[index].stellar_type)
#                R_print[index].append(particles[index].radius)
#                t_V_print[index].append(particles[index].tides_viscous_time_scale)
#                Rc_print[index].append(particles[index].convective_envelope_radius)
#                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
        #t_print.append(t)
        N = 10
        tend = 1e10

        N =1000
        tend = 1e9
        #tend = 1e8

        i_status = 0
        #print("rel_INCL_print",rel_INCL_print[0][2])
        seed=1
        dt = tend/float(N)
        i = 0
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state
        
            if state != 0:
                print("Restructuring of system!")
                particles = code.particles
                orbits = [x for x in particles if x.is_binary==True]
                bodies = [x for x in particles if x.is_binary==False]
                N_orbits = len(orbits)
                N_bodies = len(bodies)
                
                t_print.append([])
                integration_flags.append([])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
                R_print.append([[] for x in range(N_bodies)])
                Rc_print.append([[] for x in range(N_bodies)])
                R_L_print.append([[] for x in range(N_bodies)])
                t_V_print.append([[] for x in range(N_bodies)])
                mc_print.append([[] for x in range(N_bodies)])
                a_print.append([[] for x in range(N_orbits)])
                e_print.append([[] for x in range(N_orbits)])
                rel_INCL_print.append([[] for x in range(N_orbits)])
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
            #if t > 50: i_status += 1
                
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag,'integration_flag',code.integration_flag,'i_status',i_status)
            print("bound",[x.is_bound for x in bodies])
            
            if state != 0:
                print("Restructuring of system!")
                
                #exit(0)

            for index in range(N_orbits):
               # print(i_status,index)
                #print("in ",[o.index for o in orbits])
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[i_status][index].append(particles[index].mass)
                k_print[i_status][index].append(particles[index].stellar_type)
                R_print[i_status][index].append(particles[index].radius)
                t_V_print[i_status][index].append(particles[index].tides_viscous_time_scale)
                Rc_print[i_status][index].append(particles[index].convective_envelope_radius)
                R_L_print[i_status][index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[i_status][index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)
            if flag==2:
                print("Root",t,'RLOF',[o.RLOF_at_pericentre_has_occurred for o in orbits],'col',[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits],'dyn inst',[o.dynamical_instability_has_occurred for o in orbits],'sec break',[o.secular_breakdown_has_occurred for o in orbits],'min peri',[o.minimum_periapse_distance_has_occurred for o in orbits],'GW',[o.GW_condition_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1
            
            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])
        #for index in range(N_bodies):
        #        m_print[index] = np.array(m_print[index])
        #        t_V_print[index] = np.array(t_V_print[index])
        #for index in range(N_orbits):

#            rel_INCL_print[index] = np.array(rel_INCL_print[index])
#            e_print[index] = np.array(e_print[index])
#            a_print[index] = np.array(a_print[index])

        #rel_INCL_print = np.array(rel_INCL_print)
        #e_print = np.array(e_print)
        
        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            
            
            colors = ['k','tab:red','tab:green','tab:blue','y']
            for i_status in range(N_status):
                color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]

                
                
                linewidth=1.0
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color=color,linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    #if index in [0,1]:
                        #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    linewidth+=0.8
                    
                linewidth=1.0
                for index in range(N_orbits):
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
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
            fig.savefig("test16.pdf")
            pyplot.show()


    def test4(self,args):
        print('CE in binary')
        
        N_bodies=2
        
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,8.5],[30.0,400.0],[0.1,0.6],[0.0001,85.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [4.0,2.8,1.5],[3000.0,400000.0],[0.1,0.3],[0.0001,89.9*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        #particles = Tools.create_nested_multiple(N_bodies, [30.0,25.8,3.5,2.0],[35.0,400.0,5000.0],[0.1,0.6],[0.0001,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.01*np.pi/180.0],[0.01,0.01])
        particles = Tools.create_nested_multiple(N_bodies, [6.0,0.5],[1.0],[0.1],[0.0001],[45.0*np.pi/180.0],[0.01])
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
        
        code.include_flybys = False
        #code.flybys_reference_binary = orbits[0].index
        #print 'a',orbits[0].a
        code.flybys_stellar_density = 0.1*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e5



#        INCL_print = []
        t_print = [[]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
        R_print = [[[] for x in range(N_bodies)]]
        Rc_print = [[[] for x in range(N_bodies)]]
        R_L_print = [[[] for x in range(N_bodies)]]
        t_V_print = [[[] for x in range(N_bodies)]]
        mc_print = [[[] for x in range(N_bodies)]]
        a_print = [[[] for x in range(N_orbits)]]
        e_print = [[[] for x in range(N_orbits)]]
        rel_INCL_print = [[[] for x in range(N_orbits)]]
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        
        t = 0.0
        integration_flags = [[]]

#        if 1==0:
#            for index in range(N_bodies):
#                m_print[index].append(particles[index].mass)
#                k_print[index].append(particles[index].stellar_type)
#                R_print[index].append(particles[index].radius)
#                t_V_print[index].append(particles[index].tides_viscous_time_scale)
#                Rc_print[index].append(particles[index].convective_envelope_radius)
#                mc_print[index].append(particles[index].convective_envelope_mass)
        #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
        #rel_INCL_print.append(inner_binary.INCL_parent)
        #e_print.append(inner_binary.e)
        #INCL_print.append(inner_binary.INCL)
        #t_print.append(t)
        N = 10
        tend = 1e10

        N =1000
        tend = 1e9
        #tend = 1e8

        i_status = 0
        #print("rel_INCL_print",rel_INCL_print[0][2])
        seed=1
        dt = tend/float(N)
        i = 0
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            flag = code.CVODE_flag
            state = code.state
        
            if state != 0:
                print("Restructuring of system!")
                particles = code.particles
                orbits = [x for x in particles if x.is_binary==True]
                bodies = [x for x in particles if x.is_binary==False]
                N_orbits = len(orbits)
                N_bodies = len(bodies)
                
                t_print.append([])
                integration_flags.append([])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
                R_print.append([[] for x in range(N_bodies)])
                Rc_print.append([[] for x in range(N_bodies)])
                R_L_print.append([[] for x in range(N_bodies)])
                t_V_print.append([[] for x in range(N_bodies)])
                mc_print.append([[] for x in range(N_bodies)])
                a_print.append([[] for x in range(N_orbits)])
                e_print.append([[] for x in range(N_orbits)])
                rel_INCL_print.append([[] for x in range(N_orbits)])
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
            #if t > 50: i_status += 1
                
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag,'integration_flag',code.integration_flag,'i_status',i_status)
            print("bound",[x.is_bound for x in bodies])
            
            if state != 0:
                print("Restructuring of system!")
                
                #exit(0)

            for index in range(N_orbits):
                print(i_status,index)
                #print("in ",[o.index for o in orbits])
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                m_print[i_status][index].append(particles[index].mass)
                k_print[i_status][index].append(particles[index].stellar_type)
                R_print[i_status][index].append(particles[index].radius)
                t_V_print[i_status][index].append(particles[index].tides_viscous_time_scale)
                Rc_print[i_status][index].append(particles[index].convective_envelope_radius)
                R_L_print[i_status][index].append(particles[index].roche_lobe_radius_pericenter)
                mc_print[i_status][index].append(particles[index].convective_envelope_mass)
            #rel_INCL_print.append(compute_mutual_inclination(inner_binary.INCL,outer_binary.INCL,inner_binary.LAN,outer_binary.LAN))
            #rel_INCL_print.append(inner_binary.INCL_parent)
            #e_print.append(inner_binary.e)
            #INCL_print.append(inner_binary.INCL)
            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)
            if flag==2:
                print("Root",t,'RLOF',[o.RLOF_at_pericentre_has_occurred for o in orbits],'col',[o.physical_collision_or_orbit_crossing_has_occurred for o in orbits],'dyn inst',[o.dynamical_instability_has_occurred for o in orbits],'sec break',[o.secular_breakdown_has_occurred for o in orbits],'min peri',[o.minimum_periapse_distance_has_occurred for o in orbits],'GW',[o.GW_condition_has_occurred for o in orbits])
                break
            if flag==3:
                print("SNe")
                break
            if state==3:
                print("unbound")
                break
            code.random_seed = seed
            seed += 1
            
            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])
        #for index in range(N_bodies):
        #        m_print[index] = np.array(m_print[index])
        #        t_V_print[index] = np.array(t_V_print[index])
        #for index in range(N_orbits):

#            rel_INCL_print[index] = np.array(rel_INCL_print[index])
#            e_print[index] = np.array(e_print[index])
#            a_print[index] = np.array(a_print[index])

        #rel_INCL_print = np.array(rel_INCL_print)
        #e_print = np.array(e_print)
        
        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            if 1==1:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            #plot1.plot(t_print,np.array(INCL_print)*180.0/np.pi)
            
            
            colors = ['k','tab:red','tab:green','tab:blue','y']
            for i_status in range(N_status):
                color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]

                
                
                linewidth=1.0
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color=color,linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    #if index in [0,1]:
                        #plot2.plot(1.0e-6*t_print,R_L_print[index],color='g',linestyle='dotted',linewidth=linewidth)
                    #plot3.plot(1.0e-6*t_print,t_V_print[index],color='k',linestyle='solid',linewidth=linewidth)
                    linewidth+=0.8
                    
                linewidth=1.0
                for index in range(N_orbits):
                    smas = np.array(a_print[i_status][index])
                    es = np.array(e_print[i_status][index])
                    plot2.plot(1.0e-6*t_print[i_status],smas,color=color,linestyle='dotted',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],smas*(1.0-es),color=color,linestyle='solid',linewidth=linewidth)
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
            fig.savefig("test16.pdf")
            pyplot.show()



def reorientation_function(VRR_model,VRR_timescale,next_reorientation_time,orbit):

    print('='*50)
    print('reorientation_function')

    theta = 2.0*np.pi*randomf.random() - np.pi
    phi = 2.0*np.pi*randomf.random()
    r_hat_vec = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])

    if VRR_model == 1:

        Omega = 1.0/VRR_timescale
        
        orbit.VRR_Omega_vec_x = r_hat_vec[0]*Omega
        orbit.VRR_Omega_vec_y = r_hat_vec[1]*Omega
        orbit.VRR_Omega_vec_z = r_hat_vec[2]*Omega
    
    if VRR_model == 2:
        raise RuntimeError('VRR model 2 not supported')

    if VRR_model == 3:

        mu = 0.0
        sigma = (1.0/VRR_timescale)
        orbit.VRR_eta_20_init = orbit.VRR_eta_20_final
        orbit.VRR_eta_a_22_init = orbit.VRR_eta_a_22_final
        orbit.VRR_eta_b_22_init = orbit.VRR_eta_b_22_final
        orbit.VRR_eta_a_21_init = orbit.VRR_eta_a_21_final
        orbit.VRR_eta_b_21_init = orbit.VRR_eta_b_21_final

        orbit.VRR_eta_20_final = randomf.normal(mu,sigma)
        orbit.VRR_eta_a_22_final = randomf.normal(mu,sigma) 
        orbit.VRR_eta_b_22_final = randomf.normal(mu,sigma)
        orbit.VRR_eta_a_21_final = randomf.normal(mu,sigma)
        orbit.VRR_eta_b_21_final = randomf.normal(mu,sigma)
         
        orbit.VRR_initial_time = orbit.VRR_final_time
        orbit.VRR_final_time = next_reorientation_time

def compute_M_star_r(r,gamma,n_0,r_0,m_star):
    return 4.0*np.pi*(1.0/(3.0-gamma))*m_star*r_0**3*n_0*pow(r/r_0,3.0-gamma)

def compute_N_star_r(r,gamma,n_0,r_0,m_star):
    return compute_M_star_r(r,gamma,n_0,r_0,m_star)/m_star
    
def compute_n_star_r(r,gamma,n_0,r_0,m_star):
    return n_0*pow(r/r_0,-gamma)

def compute_rho_star_r(r,gamma,n_0,r_0,m_star):
    return m_star*compute_n_star_r(r,gamma,n_0,r_0,m_star)

def compute_sigma_r(r,gamma,n_0,r_0,m_star,M_MBH,CONST_G):
    return np.sqrt( CONST_G*M_MBH*(1.0/(r*(1.0+gamma)))*(1.0 + ((gamma+1.0)/(2.0*gamma-2.0))*compute_M_star_r(r,gamma,n_0,r_0,m_star)/M_MBH ) )
    
    
if __name__ == '__main__':
    args = parse_arguments()
    
    N_tests = 4
    if args.test==0:
        tests = range(1,N_tests+1)
    else:
        tests = [args.test]

    t=test_mse()
    for i in tests:
        print( 'Running test number',i,'; verbose =',args.verbose,'; plot =',args.plot)
        function = getattr(t, 'test%s'%i)
        function(args)
    
    print("="*50)
    print("All tests passed!")
