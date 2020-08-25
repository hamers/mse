import os
import numpy as np
import ctypes


class Tools(object):
 
    @staticmethod       
    def create_nested_multiple(N,masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None,metallicities=None,stellar_types=None):
        from mse import Particle

        """
        N is number of bodies
        masses should be N-sized array
        the other arguments should be (N-1)-sized arrays
        """

        N_bodies = N
        N_binaries = N-1

        particles = []

        for index in range(N_bodies):
            particle = Particle(is_binary=False,mass=masses[index])
            if radii is not None:
                particle.radius = radii[index]
            if metallicities is not None:
                particle.metallicity = metallicities[index]
            if stellar_types is not None:
                particle.stellar_type = stellar_types[index]
                
            particles.append(particle)
        
        for index in range(N_binaries):
            if index==0:
                child1 = particles[0]
                child2 = particles[1]
            else:
                child1 = previous_binary
                child2 = particles[index+1]
            particle = Particle(is_binary=True,child1=child1,child2=child2,a=semimajor_axes[index],e=eccentricities[index],INCL=inclinations[index],AP=arguments_of_pericentre[index],LAN=longitudes_of_ascending_node[index])

            previous_binary = particle
            particles.append(particle)
            
        return particles

    @staticmethod
    def create_2p2_quadruple_system(masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None,metallicities=None,stellar_types=None):
        from mse import Particle
        
        """
        Create a 2+2 quadruple system.
        Masses should contain the four masses.
        The other arguments should be length 3 arrays; first two entries: the two inner binaries; third entry: outer binary.
        """

        N_bodies = 4
        N_binaries = N_bodies-1

        particles = []

        ### Add the bodies ###
        for index in range(N_bodies):
            particle = Particle(is_binary=False,mass=masses[index])
            if radii is not None:
                particle.radius = radii[index]
            if metallicities is not None:
                particle.metallicity = metallicities[index]
            if stellar_types is not None:
                particle.stellar_type = stellar_types[index]
                
            particles.append(particle)

        ### Add the binaries ###
        for index in range(N_binaries):
            if index==0:
                child1 = particles[0]
                child2 = particles[1]
            elif index==1:
                child1 = particles[2]
                child2 = particles[3]
            elif index==2:
                child1 = particles[4]
                child2 = particles[5]

            particle = Particle(is_binary=True,child1=child1,child2=child2,a=semimajor_axes[index],e=eccentricities[index],INCL=inclinations[index],AP=arguments_of_pericentre[index],LAN=longitudes_of_ascending_node[index])
            particles.append(particle)
        
        return particles

    @staticmethod
    def compute_mutual_inclination(INCL_k,INCL_l,LAN_k,LAN_l):
        cos_INCL_rel = np.cos(INCL_k)*np.cos(INCL_l) + np.sin(INCL_k)*np.sin(INCL_l)*np.cos(LAN_k-LAN_l)
        return np.arccos(cos_INCL_rel)

    @staticmethod
    def determine_binary_masses(particles):
        ### set binary masses -- to ensure this happens correctly, do this from highest level to lowest level ###

        Tools.determine_binary_levels_in_particles(particles)

        max_level = np.amax([x.level for x in particles])
        level = max_level
        while (level > -1):
            for index,p in enumerate(particles):
                if (p.is_binary == True and p.level == level):
                    p.mass = p.child1.mass + p.child2.mass
            level -= 1

    @staticmethod
    def determine_binary_levels_in_particles(particles):
        for index,p in enumerate(particles):
            p.index_temp = index
            p.parent = None

        ### determine top binary ###
        for index_particle_1,particle_1 in enumerate(particles):
            if particle_1.is_binary == True:
                child1 = particle_1.child1
                child2 = particle_1.child2
                
                for index_particle_2,particle_2 in enumerate(particles):
                    if (index_particle_2 == child1.index_temp or index_particle_2 == child2.index_temp):
                        particle_2.parent = particle_1
                        
        for index_particle_1,particle_1 in enumerate(particles):
            particle_1.level = 0
            
            child = particle_1;
            parent = particle_1.parent

            if (parent != None): ### if parent == -1, P_p is the `top' binary, for which level=0 
                while (parent != None): ### search parents until reaching the top binary 
                    for index_particle_2,particle_2 in enumerate(particles):
                        
                        if parent == None: break
                        if (particle_2.index_temp == parent.index_temp):
                            particle_1.level += 1
                            
                            parent = particle_2.parent
                            
    @staticmethod
    def generate_mobile_diagram(particles,plot,line_width_horizontal=1.5,line_width_vertical = 0.2,line_color = 'k',line_width = 1.5,fontsize=12,use_default_colors=True):
        """
        Generate a Mobile diagram of a given multiple system.
        """
        
        try:
            import matplotlib
        except ImportError:
            print("mse.py -- generate_mobile_diagram -- unable to import Matplotlib which is needed to generate a Mobile diagram!")
            exit(0)

        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

#        print("N",len(bodies),len(binaries))

        if len(binaries)==0:
            if len(bodies)==0:
                print("mse.py -- generate_mobile_diagram -- zero bodies and zero binaries!")
                exit(0)
            else:
                Tools.draw_bodies(plot,bodies,fontsize)
                return

        Tools.determine_binary_levels_in_particles(particles)                    
        unbound_bodies = [x for x in particles if x.is_binary==False and x.parent == None]
        if len(unbound_bodies)>0:
            Tools.draw_bodies(plot,unbound_bodies,fontsize,y_ref = 1.4*line_width_vertical,dx=0.4*line_width_horizontal,dy=0.4*line_width_vertical)
            
        top_level_binary = [x for x in binaries if x.level==0][0]
        
        if use_default_colors==True:
            ### Assign some colors from mcolors to the orbits ###
            import matplotlib.colors as mcolors
            colors = mcolors.TABLEAU_COLORS
            color_names = list(colors)
            
            for index in range(len(binaries)):
                color_name = color_names[index]
                color=colors[color_name]
            
                o = binaries[index]
                o.color = color

        ### Make mobile diagram ###
        top_level_binary.x = 0.0
        top_level_binary.y = 0.0
        x_min = x_max = y_min = 0.0
        y_max = line_width_vertical
    
        plot.plot( [top_level_binary.x,top_level_binary.x], [top_level_binary.y,top_level_binary.y + line_width_vertical ], color=line_color,linewidth=line_width)
        x_min,x_max,y_min,y_max = Tools.draw_binary_node(plot,top_level_binary,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max)

        
        plot.set_xticks([])
        plot.set_yticks([])
        #print("minmax",x_min,x_max,y_min,y_max)
        beta = 0.7
        plot.set_xlim([x_min - beta*np.fabs(x_min),x_max + beta*np.fabs(x_max)])
        plot.set_ylim([y_min - beta*np.fabs(y_min),1.5*y_max + beta*np.fabs(y_max)])
        
        #plot.autoscale(enable=True,axis='both')
        
    @staticmethod
    def draw_binary_node(plot,particle,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max):
        x = particle.x
        y = particle.y
        
        child1 = particle.child1
        child2 = particle.child2

        #from decimal import Decimal
        #plot.annotate("$a = % \, \mathrm{AU}$"%(particle.semimajor_axis.value_in(units.AU)),xytext=(x,y))
        #text = "$a = %s \, \mathrm{AU}$"%(round(particle.semimajor_axis.value_in(units.AU),1))
        #text = "$a = \mathrm{%.1E}$"%(Decimal(particle.a))
        #text = "$a=\mathrm{%.1E\,au}$"%(Decimal(particle.a))
        text = "$a=\mathrm{%s\,au}$"%(round(particle.a,1))
        plot.annotate(text,xy=(x - 0.8*line_width_horizontal,y - 0.3*line_width_vertical),fontsize=fontsize,color=particle.color)
     
        #text = "$e = %.2f$"%(particle.e)
        text = "$e = %.2f$"%(particle.e)
        
        plot.annotate(text,xy=(x - 0.8*line_width_horizontal,y - 0.6*line_width_vertical),fontsize=fontsize,color=particle.color)

        alpha = 1.0
        if child1.is_binary == True and child2.is_binary == True:
            alpha = 3.5

        if child1.is_binary == True and child2.is_binary == False:
            alpha = 2.5
        if child1.is_binary == False and child2.is_binary == True:
            alpha = 2.5

        ### lines to child1 ###
        plot.plot( [x,x - alpha*line_width_horizontal],[y,y], color=line_color,linewidth=line_width)
        plot.plot( [x - alpha*line_width_horizontal,x - alpha*line_width_horizontal], [y,y - line_width_vertical], color=line_color,linewidth=line_width)
        
        ### lines to child2 ###
        plot.plot( [x,x + alpha*line_width_horizontal],[y,y], color=line_color,linewidth=line_width)
        plot.plot( [x + alpha*line_width_horizontal,x + alpha*line_width_horizontal], [y,y - line_width_vertical], color=line_color,linewidth=line_width)

        ### positions of children ###
        child1 = particle.child1
        child2 = particle.child2
        
        child1.x = particle.x - alpha*line_width_horizontal
        child2.x = particle.x + alpha*line_width_horizontal

        child1.y = particle.y - line_width_vertical
        child2.y = particle.y - line_width_vertical


        if (child1.x<x_min): x_min = child1.x
        if (child1.x>x_max): x_max = child1.x
        if (child2.x<x_min): x_min = child2.x
        if (child2.x>x_max): x_max = child2.x

        if (child1.y<y_min): y_min = child1.y
        if (child1.y>y_max): y_max = child1.y
        if (child2.y<y_min): y_min = child2.y
        if (child2.y>y_max): y_max = child2.y

        
        ### handle children ###
        if child1.is_binary == True:
            x_min,x_max,y_min,y_max = Tools.draw_binary_node(plot,child1,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max)
        else:
            color,s = Tools.get_color_and_size_for_star(child1.stellar_type,child1.radius)
            plot.scatter([child1.x],[child1.y],color=color,s=s,zorder=10)
            #text = "$%s\, M_\mathrm{J}$"%(round(child1.mass.value_in(units.MJupiter)))
            #text = "$m_i=\mathrm{%.1E}\,\mathrm{M}_\odot$"%(Decimal(child1.mass))
            text = "$\mathrm{%s}\,\mathrm{M}_\odot$"%(str(round(child1.mass,1)))
            #text = "$\mathrm{%s}\,\mathrm{M}_\odot\,(%d)$"%(str(round(child1.mass,1)),child1.index)
            plot.annotate(text,xy=(child1.x - 0.6*line_width_horizontal,child1.y - 0.5*line_width_vertical),color='k',fontsize=fontsize)
            text = "$%d; \,k=%d$"%(child1.index,child1.stellar_type)
            plot.annotate(text,xy=(child1.x - 0.6*line_width_horizontal,child1.y - 0.25*line_width_vertical),color='k',fontsize=0.5*fontsize)

        if child2.is_binary == True:
            x_min,x_max,y_min,y_max = Tools.draw_binary_node(plot,child2,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max)
        else:
            color,s = Tools.get_color_and_size_for_star(child2.stellar_type,child2.radius)
            plot.scatter([child2.x],[child2.y],color=color,s=s,zorder=10)
            #text = "$%s\, M_\mathrm{J}$"%(round(child2.mass.value_in(units.MJupiter)))
            #text = "$\mathrm{%.1E}\,\mathrm{M}_\odot$"%(Decimal(child2.mass))
            text = "$\mathrm{%s}\,\mathrm{M}_\odot$"%(str(round(child2.mass,1)))
            #text = "$\mathrm{%s}\,\mathrm{M}_\odot\,(%d)$"%(str(round(child2.mass,1)),child2.index)
            plot.annotate(text,xy=(child2.x - 0.3*line_width_horizontal,child2.y - 0.5*line_width_vertical),color='k',fontsize=fontsize)
            text = "$%d; \,k=%d$"%(child2.index,child2.stellar_type)
            plot.annotate(text,xy=(child2.x - 0.3*line_width_horizontal,child2.y - 0.25*line_width_vertical),color='k',fontsize=0.5*fontsize)

        return x_min,x_max,y_min,y_max

    @staticmethod
    def get_color_and_size_for_star(stellar_type,radius):
        if (stellar_type <= 1): color='gold'
        elif (stellar_type == 2): color='darkorange'
        elif (stellar_type == 3): color='firebrick'
        elif (stellar_type == 4): color='darkorange'
        elif (stellar_type == 5): color='orangered'
        elif (stellar_type == 6): color='crimson'
        elif (stellar_type == 7): color='royalblue'
        elif (stellar_type == 8): color='orangered'
        elif (stellar_type == 9): color='crimson'
        elif (stellar_type == 10): color='silver'
        elif (stellar_type == 11): color='silver'
        elif (stellar_type == 12): color='silver'
        elif (stellar_type == 13): color='gainsboro'
        elif (stellar_type == 14): color='k'
        else: color = 'k'

        CONST_R_SUN = 0.004649130343817401
        CONST_KM = 1.0/(1.4966e9)

        s = 5 + 50*np.log10(radius/CONST_R_SUN)
        if (stellar_type >= 7 and stellar_type <= 9):
            s = 5 + 5*np.log10(radius/CONST_KM)

        if (stellar_type >= 10):
            s = 5 + 20*np.log10(radius/CONST_KM)

        return color,s
        
    @staticmethod
    def draw_bodies(plot,bodies,fontsize,y_ref=1.0,dx=0.5,dy=0.5):
        #dx = 0.5
        #dy = 0.5
        for index,body in enumerate(bodies):
            color,s = Tools.get_color_and_size_for_star(body.stellar_type,body.radius)
            plot.scatter([index],[y_ref],color=color,s=s)
            #text = "$%s\, M_\mathrm{J}$"%(round(child1.mass.value_in(units.MJupiter)))
            #text = "$m_i=\mathrm{%.1E}\,\mathrm{M}_\odot$"%(Decimal(child1.mass))
            text = "$\mathrm{%s}\,\mathrm{M}_\odot$"%(str(round(body.mass,1)))
            plot.annotate(text,xy=(index - dx,y_ref-dy),color='k',fontsize=fontsize)

            #text = "$\mathrm{%d}$"%body.index
            text = "$%d; \,k=%d$"%(body.index,body.stellar_type)
            plot.annotate(text,xy=(index - dx,y_ref+0.5*dy),color='k',fontsize=fontsize)
            
            try:
                VX = body.VX
                VY = body.VY
                VZ = body.VZ
                V = np.sqrt(VX**2 + VY**2 + VZ**2)
                x = index
                y = y_ref
                Adx = 0.5*dx*VX/V
                Ady = 0.5*dx*VY/V
                plot.arrow(x, y, Adx, Ady,color=color,head_width=0.05, head_length=0.05)
            except AttributeError:
                pass

        plot.set_xlim([-2*dx,len(bodies)])
        plot.set_ylim([y_ref-2*dy,y_ref+2*dy])

        plot.set_xticks([])
        plot.set_yticks([])

    @staticmethod
    def get_description_from_event_flag(event_flag):
        if event_flag == 0:
            text = "$\mathrm{Initial\,system}$"
        elif event_flag == 1:
            text = "$\mathrm{Stellar\,type\,change}$"
        elif event_flag == 2:
            text = "$\mathrm{SNe}$"
        elif event_flag == 3:
            text = "$\mathrm{RLOF\,start}$"
        elif event_flag == 4:
            text = "$\mathrm{RLOF\,end}$"
        elif event_flag == 5:
            text = "$\mathrm{CE}$"
        elif event_flag == 6:
            text = "$\mathrm{Collision}$"
        elif event_flag == 7:
            text = "$\mathrm{Dyn.\,inst.}$"
        elif event_flag == 8:
            text = "$\mathrm{Sec.\,break.}$"
        else:
            text = ""
        return text

    @staticmethod
    def evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,tend,N_steps,stellar_types=None,make_plots=True,fancy_plots=False,plot_filename="test1",show_plots=True):

        if configuration == "fully nested":
            particles = Tools.create_nested_multiple(N_bodies, masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,metallicities=metallicities,stellar_types=stellar_types)
        elif configuration == "2+2 quadruple":
            particles = Tools.create_2p2_quadruple_system(masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,metallicities=metallicities,stellar_types=stellar_types)
        else:
            print("evolve_system.py: configuration ",configuration," currently not supported!")
            exit(-1)
        
        from mse import MSE

        orbits = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]

        N_orbits = len(orbits)
                
        code = MSE()
        code.enable_root_finding = True
        code.add_particles(particles)
        code.orbital_phases_random_seed = 1
        
        code.include_flybys = True
        code.flybys_stellar_density = 10*0.1*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e5

        t_print = [[]]
        internal_indices_print = [[[] for x in range(N_bodies)]]
        k_print = [[[] for x in range(N_bodies)]]
        m_print = [[[] for x in range(N_bodies)]]
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
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        
        t = 0.0
        integration_flags = [[]]

        i_status = 0

        seed=1
        dt = tend/float(N_steps)
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
    
            #i_rel = Tools.compute_mutual_inclination(code.log[-1]["particles"][3].INCL,code.log[-1]["particles"][4].INCL,code.log[-1]["particles"][3].LAN,code.log[-1]["particles"][4].LAN)
            #print("IREL",(180.0/np.pi)*i_rel)
            #print("LOG a",code.log[0]["particles"][4].a,code.log[-1]["time"],code.log[-1]["event_flag"])
            #print("LOG",code.log)
            #if code.log[-1]["event_flag"] != 0:
            #    print("LOG","event",code.log[-1]["event_flag"],"index1",code.log[-1]["index1"],"index2",code.log[-1]["index2"],"binary index",code.log[-1]["binary_index"],code.log[-1]["particles"])
            
            if code.structure_change == True:
                print("Python restruct")#,children1,children1_old,children2,children2_old)
                t_print.append([])
                integration_flags.append([])
                internal_indices_print.append([[] for x in range(N_bodies)])
                k_print.append([[] for x in range(N_bodies)])
                m_print.append([[] for x in range(N_bodies)])
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
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
                
            print( 't/Myr',t*1e-6,'es',[o.e for o in orbits],'as',[o.a for o in orbits],flag,'integration_flag',code.integration_flag,'i_status',i_status)
            
            for index in range(N_orbits):
                rel_INCL_print[i_status][index].append(orbits[index].INCL_parent)
                e_print[i_status][index].append(orbits[index].e)
                a_print[i_status][index].append(orbits[index].a)
            for index in range(N_bodies):
                internal_indices_print[i_status][index].append(particles[index].index)
                m_print[i_status][index].append(particles[index].mass)
                k_print[i_status][index].append(particles[index].stellar_type)
                R_print[i_status][index].append(particles[index].radius)
                X_print[i_status][index].append(particles[index].X)
                Y_print[i_status][index].append(particles[index].Y)
                Z_print[i_status][index].append(particles[index].Z)
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
            #code.random_seed = seed
            seed += 1
            
            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])

        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])

       
        if make_plots==True:
            
            try:
                from matplotlib import pyplot
            except ImportError:
                print("evolve_system.py -- ERROR: cannot import Matplotlib")
                exit(-1)
            
            if fancy_plots == True:
                print("Using LaTeX for plot text")
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            print("log",len(code.log))
            
            plot_log = []
            previous_event_flag = -1
            for index_log,log in enumerate(code.log):
                event_flag = log["event_flag"]
                if previous_event_flag == event_flag and event_flag == 3:
                    continue
                plot_log.append(log)
                previous_event_flag = event_flag
                            
            N_l = len(plot_log)
            fontsize=N_l
            fig=pyplot.figure(figsize=(N_l,N_l))
            N_r = int(np.sqrt(N_l))+1
            N_c = N_r
            for index_log,log in enumerate(plot_log):
                plot=fig.add_subplot(N_r,N_c,index_log+1)
                particles = log["particles"]
                event_flag = log["event_flag"]
                
                Tools.generate_mobile_diagram(particles,plot,fontsize=0.5*N_l)

                text = Tools.get_description_from_event_flag(event_flag)
                plot.set_title(text,fontsize=fontsize)
                plot.annotate("$t=%s\,\mathrm{Myr}$"%round(log["time"]*1e-6,1),xy=(0.1,0.9),xycoords='axes fraction',fontsize=fontsize)
                #if index_log>0: break
    
            fig.savefig(plot_filename + "_mobile.pdf")
        
            fig=pyplot.figure(figsize=(8,8))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            
            fig_pos=pyplot.figure(figsize=(8,8))
            plot_pos=fig_pos.add_subplot(1,1,1)
            
            colors = ['k','tab:red','tab:green','tab:blue','y','k','tab:red','tab:green','tab:blue','y']
            linewidth=1.0
            for i_status in range(N_status):
                #color = colors[i_status]
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]
                
                plot3.plot(1.0e-6*t_print[i_status],integration_flags[i_status],color='k',linestyle='dotted',linewidth=linewidth)
                
                for index in range(N_bodies):
                    #print("internal_indices_print[i_status][index]",internal_indices_print[i_status][index],index,internal_indices_print)
                    #internal_index = internal_indices_print[i_status][index][0]
                    #color = colors[internal_index]
                    color=colors[index]
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    
                    parsec_in_AU = code.CONST_PARSEC
                    plot_pos.plot(np.array(X_print[i_status][index])/parsec_in_AU,np.array(Y_print[i_status][index])/parsec_in_AU,color=color,linestyle='solid',linewidth=linewidth)
                    #print("i",index,"Final R",X_print[i_status][index][-1],Y_print[i_status][index][-1])
                    
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

            log_CEs = [x for x in code.log if x["event_flag"] == 5]
            t_CEs_Myr = np.array([x["time"]*1e-6 for x in log_CEs])
            
            for k,t in enumerate(t_CEs_Myr):
                plot2.axvline(x=t,linestyle='dashed',color='tab:red')
                plot2.annotate("$\mathrm{CE}$",xy=(t,1.0e3),fontsize=fontsize)
            
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{au}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            plot3.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            
            plot_pos.set_xlabel("$X/\mathrm{pc}$",fontsize=fontsize)
            plot_pos.set_ylabel("$Y/\mathrm{pc}$",fontsize=fontsize)
            
            fig.savefig(plot_filename + ".pdf")
            fig_pos.savefig(plot_filename + "_pos.pdf")
            
            if show_plots == True:
                pyplot.show()
    
        code.reset()
