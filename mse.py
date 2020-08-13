import os
import numpy as np
import ctypes

"""
SecularMultiple
    
A code to compute the secular (orbit-averaged) gravitational dynamics of hierarchical multiple systems
composed of nested binary orbits (simplex-type systems) with any configuration and any number of bodies.
A particle can repesent a binary (`is_binary = True') or a body (`is_binary = False').
The structure of the system is determined by linking to other particles with the attributes child1 and child2.
Tidal interactions and relativistic corrections are included in an ad hoc fashion
(tides: treating the companion as a single body, even if it is not; relativistic terms:
only including binary-binary interactions).
    
Includes routines for external perturbations (flybys & supernovae).

If you use this code for work in scientific publications, please cite:
https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2827H (the original paper)
https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.4139H (updates with external perturbations)

Make sure to first compile the code using `make'. The script `test_secularmultiple.py' can be used to test the
installation. See examples.py for some examples.

Adrian Hamers, June 2019
"""


class MSE(object):

    def __init__(self):
        self.__CONST_G = 4.0*np.pi**2
        self.__CONST_C = 63239.72638679138
        self.__CONST_M_SUN = 1.0
        self.__CONST_R_SUN = 0.004649130343817401
        self.__CONST_L_SUN = 0.0002710404109745588
        self.__CONST_KM_PER_S = 0.210862
        self.__CONST_PER_PC3 = 1.14059e-16


        self.__relative_tolerance = 1.0e-12
        self.__absolute_tolerance_eccentricity_vectors = 1.0e-10
        self.__include_quadrupole_order_terms = True
        self.__include_octupole_order_binary_pair_terms = True
        self.__include_octupole_order_binary_triplet_terms = True
        self.__include_hexadecupole_order_binary_pair_terms = True
        self.__include_dotriacontupole_order_binary_pair_terms = True
        self.__include_double_averaging_corrections = False

        self.__binary_evolution_CE_energy_flag = 0
        self.__binary_evolution_CE_spin_flag = 0

        self.__mstar_gbs_tolerance_default = 1.0e-10
        self.__mstar_gbs_tolerance_kick = 1.0e-6
        self.__mstar_collision_tolerance = 1.0e-10

        self.__particles_committed = False
        self.model_time = 0.0
        self.time_step = 0.0
        self.relative_energy_error = 0.0
        self.state = 0
        self.CVODE_flag = 0
        self.CVODE_error_code = 0
        self.integration_flag = 0 # start with secular approach 
        
        self.__random_seed = 0
        
        self.enable_tides = True
        self.enable_root_finding = True
        self.enable_VRR = False
        
        self.__include_flybys = True
        self.__flybys_correct_for_gravitational_focussing = True
        self.__flybys_reference_binary = -1
        self.__flybys_velocity_distribution = 0
        self.__flybys_mass_distribution = 0
        self.__flybys_mass_distribution_lower_value = 0.1
        self.__flybys_mass_distribution_upper_value = 100.0
        self.__flybys_encounter_sphere_radius = 1.0e5
        self.__flybys_stellar_density = 0.1*self.__CONST_PER_PC3 ### density at infinity
        self.__flybys_stellar_relative_velocity_dispersion = 30.0*self.__CONST_KM_PER_S

        __current_dir__ = os.path.dirname(os.path.realpath(__file__))
        lib_path = os.path.join(__current_dir__, 'libmse.so')

        #print 'lib_path',lib_path
#        if not os.path.isfile(lib_path):
            # try to find the library from the parent directory
#            lib_path = os.path.join(os.path.abspath(os.path.join(__current_dir__, os.pardir)), 'libmse.so')
            #print 'not fil'

        if not os.path.isfile(lib_path):
            print('Library libmse.so not exist -- trying to compile')
            os.system('make')
        
        self.lib = ctypes.cdll.LoadLibrary(lib_path)
        self.init_lib()
        self.particles = []

    def init_lib(self):
        self.lib.add_particle.argtypes = (ctypes.POINTER(ctypes.c_int),ctypes.c_bool,ctypes.c_bool)
        self.lib.add_particle.restype = ctypes.c_int
        
        self.lib.delete_particle.argtypes = (ctypes.c_int,)
        self.lib.delete_particle.restype = ctypes.c_int

        self.lib.set_children.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_int)
        self.lib.set_children.restype = ctypes.c_int

        self.lib.get_children.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int))
        self.lib.get_children.restype = ctypes.c_int

        self.lib.get_number_of_particles.argtypes = ()
        self.lib.get_number_of_particles.restype = ctypes.c_int

        self.lib.get_internal_index_in_particlesMap.argtypes = (ctypes.c_int,)
        self.lib.get_internal_index_in_particlesMap.restype = ctypes.c_int

        self.lib.get_is_binary.argtypes = (ctypes.c_int,)
        self.lib.get_is_binary.restype = ctypes.c_bool

        self.lib.get_is_bound.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_bool))
        self.lib.get_is_bound.restype = ctypes.c_int

        self.lib.set_mass.argtypes = (ctypes.c_int,ctypes.c_double)
        self.lib.set_mass.restype = ctypes.c_int

        self.lib.get_mass.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double))
        self.lib.get_mass.restype = ctypes.c_int

#        self.lib.set_mass_dot.argtypes = (ctypes.c_int,ctypes.c_double)
#        self.lib.set_mass_dot.restype = ctypes.c_int

        self.lib.set_mass_transfer_terms.argtypes = (ctypes.c_int,ctypes.c_bool)
        self.lib.set_mass_transfer_terms.restype = ctypes.c_int

        self.lib.get_mass_dot.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double))
        self.lib.get_mass_dot.restype = ctypes.c_int

        self.lib.set_radius.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.c_double)
        self.lib.set_radius.restype = ctypes.c_int

        self.lib.get_radius.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_radius.restype = ctypes.c_int

        self.lib.set_spin_vector.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_spin_vector.restype = ctypes.c_int

        self.lib.get_spin_vector.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_spin_vector.restype = ctypes.c_int

        self.lib.set_stellar_evolution_properties.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_bool,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int)
        self.lib.set_stellar_evolution_properties.restype = ctypes.c_int

        self.lib.get_stellar_evolution_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_stellar_evolution_properties.restype = ctypes.c_int

        ### kicks ###
        self.lib.set_kick_properties.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double)
        self.lib.set_kick_properties.restype = ctypes.c_int

        self.lib.get_kick_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_kick_properties.restype = ctypes.c_int

        ### orbital elements ###
        self.lib.set_orbital_elements.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int)
        self.lib.set_orbital_elements.restype = ctypes.c_int

        self.lib.get_orbital_elements.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),\
            ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_orbital_elements.restype = ctypes.c_int

        self.lib.get_inclination_relative_to_parent_interface.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double))
        self.lib.get_inclination_relative_to_parent_interface.restype = ctypes.c_int

        self.lib.get_relative_position_and_velocity.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_relative_position_and_velocity.restype = ctypes.c_int

        self.lib.get_absolute_position_and_velocity.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_absolute_position_and_velocity.restype = ctypes.c_int
        
        self.lib.set_integration_method.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_bool)
        self.lib.set_integration_method.restype = ctypes.c_int

        self.lib.set_PN_terms.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_int)
        self.lib.set_PN_terms.restype = ctypes.c_int

        self.lib.set_tides_terms.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_double)
        self.lib.set_tides_terms.restype = ctypes.c_int

        self.lib.set_root_finding_terms.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_bool,ctypes.c_bool,ctypes.c_double,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool);
        self.lib.set_root_finding_terms.restype = ctypes.c_int

        self.lib.set_root_finding_state.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool)
        self.lib.set_root_finding_state.restype = ctypes.c_int

        self.lib.get_root_finding_state.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool))
        self.lib.get_root_finding_state.restype = ctypes.c_int

        self.lib.set_constants.argtypes = (ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_constants.restype = ctypes.c_int

        self.__set_constants_in_code()

        self.lib.set_parameters.argtypes = (ctypes.c_double,ctypes.c_double,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_int,ctypes.c_bool, \
            ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double)
        self.lib.set_parameters.restype = ctypes.c_int
         
        self.__set_parameters_in_code() 
         
        self.lib.evolve_interface.argtypes = (ctypes.c_double,ctypes.c_double, \
            ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int))
        self.lib.evolve_interface.restype = ctypes.c_int

        self.lib.set_external_particle_properties.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_external_particle_properties.restype = ctypes.c_int

        self.lib.apply_external_perturbation_assuming_integrated_orbits_interface.argtypes = ()
        self.lib.apply_external_perturbation_assuming_integrated_orbits_interface.restype = ctypes.c_int

        self.lib.apply_user_specified_instantaneous_perturbation_interface.argtypes = ()
        self.lib.apply_user_specified_instantaneous_perturbation_interface.restype = ctypes.c_int

        self.lib.set_instantaneous_perturbation_properties.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_instantaneous_perturbation_properties.restype = ctypes.c_int

        self.lib.set_VRR_properties.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_VRR_properties.restype = ctypes.c_int

        self.lib.reset_interface.argtypes = ()
        self.lib.reset_interface.restype = ctypes.c_int
        
        self.lib.set_random_seed.argtypes = (ctypes.c_int,)
        self.lib.set_random_seed.restype = ctypes.c_int

        self.lib.initialize_code_interface.argtypes = ()
        self.lib.initialize_code_interface.restype = ctypes.c_int
        
        
        ### logging ###
        self.lib.get_size_of_log_data.argtypes = ()
        self.lib.get_size_of_log_data.restype = ctypes.c_int
        
        self.lib.get_log_entry_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int))
        self.lib.get_log_entry_properties.restype = ctypes.c_int
        
        self.lib.get_body_properties_from_log_entry.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_body_properties_from_log_entry.restype = ctypes.c_int
        
        
        ### tests ###
        self.lib.unit_tests_interface.argtypes = ()
        self.lib.unit_tests_interface.restype = ctypes.c_int

        self.lib.determine_compact_object_merger_properties_interface.argtypes = ( ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), \
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), \
            ctypes.POINTER(ctypes.c_double) )
        self.lib.determine_compact_object_merger_properties_interface.restype = ctypes.c_int

        self.lib.sample_from_3d_maxwellian_distribution_interface.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double) )
        self.lib.sample_from_3d_maxwellian_distribution_interface.restype = ctypes.c_int

        self.lib.sample_from_normal_distribution_interface.argtypes = (ctypes.c_double, ctypes.c_double)
        self.lib.sample_from_normal_distribution_interface.restype = ctypes.c_double

        self.lib.sample_from_kroupa_93_imf_interface.argtypes = ()
        self.lib.sample_from_kroupa_93_imf_interface.restype = ctypes.c_double

        self.lib.test_kick_velocity.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double))
        self.lib.test_kick_velocity.restype = ctypes.c_int

    ###############
    
#    def add_particle(self,particle):
#        index = ctypes.c_int(0)
#        self.lib.add_particle(ctypes.byref(index), particle.is_binary, particle.is_external)
#        particle.index = index.value
#        flag = self.__update_particle_in_code(particle)

#        self.particles.append(particle)

#    def add_particles(self,particles):
#        for index,particle in enumerate(particles):
#            self.add_particle(particle)
            
    def add_particle(self,particle):
        index = ctypes.c_int(0)

        self.lib.add_particle(ctypes.byref(index), particle.is_binary, particle.is_external)
        particle.index = index.value
        if particle.is_binary==False:
            flag = self.lib.set_mass(particle.index,particle.mass)

        self.particles.append(particle)

    def add_particles(self,particles):
        for index,particle in enumerate(particles):
            self.add_particle(particle)
        
        ### All particles need to be added individually to the code first before calling __update_particles_in_code, since the latter includes reference to particles' children ###
        flag = self.__update_particles_in_code(self.particles)


    def delete_particle(self,particle):
        flag = self.lib.delete_particle(particle.index)
        if flag==-1:
            raise RuntimeError('Could not delete particle with index {0}'.format(particle.index))
        self.particles.remove(particle)

    def commit_particles(self):
        self.__set_random_seed()
        
        flag = self.__update_particles_in_code()
        
        self.lib.initialize_code_interface()
        
        self.__update_particles_from_code()
        
        end_time,initial_hamiltonian,state,CVODE_flag,CVODE_error_code = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0)
        evolve_flag = self.lib.evolve_interface(0.0,0.0,ctypes.byref(end_time),ctypes.byref(initial_hamiltonian), \
            ctypes.byref(state),ctypes.byref(CVODE_flag),ctypes.byref(CVODE_error_code))

        self.initial_hamiltonian = initial_hamiltonian.value
        
        self.__particles_committed = True

    def evolve_model(self,end_time):
        
        if end_time is None:
            raise RuntimeError('End time not specified in evolve_model!')
        if self.__particles_committed == False:
            self.commit_particles()
        
        flag = self.__update_particles_in_code()

        ### get initial system structure ###
        orbits = [p for p in self.particles if p.is_binary == True]
        children1_old = [o.child1.index for o in orbits]
        children2_old = [o.child2.index for o in orbits]

        ### integrate system of ODEs ###
        start_time = self.model_time
#        time_step = end_time - start_time   

        output_time,hamiltonian,state,CVODE_flag,CVODE_error_code,integration_flag = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(self.integration_flag)
        evolve_flag = self.lib.evolve_interface(start_time,end_time,ctypes.byref(output_time),ctypes.byref(hamiltonian), \
            ctypes.byref(state),ctypes.byref(CVODE_flag),ctypes.byref(CVODE_error_code),ctypes.byref(integration_flag))
        output_time,hamiltonian,state,CVODE_flag,CVODE_error_code,integration_flag = output_time.value,hamiltonian.value,state.value,CVODE_flag.value,CVODE_error_code.value,integration_flag.value

        ### compute energy error ###
        self.hamiltonian = hamiltonian
        if self.initial_hamiltonian == 0.0:
            self.relative_energy_error = 0.0
        else:
            self.relative_energy_error = np.fabs( (self.initial_hamiltonian - self.hamiltonian)/self.initial_hamiltonian )

        ### update model time ###
        self.model_time = output_time

        if (flag==99):
            print('Error occurred during ODE integration; error code is {0}'.format(error_code))

        self.CVODE_flag = CVODE_flag
        self.CVODE_error_code = CVODE_error_code
        self.state = state
        self.integration_flag = integration_flag

        #if self.integration_flag == 0: 
        self.__copy_particle_structure_from_code()

        self.__update_particles_from_code()

        ### check if the system structure changed ###
        orbits = [p for p in self.particles if p.is_binary == True]
        children1 = [o.child1.index for o in orbits]
        children2 = [o.child2.index for o in orbits]
            
        self.structure_change = False
        if children1 != children1_old or children2 != children2_old:
            self.structure_change = True


        return self.state,self.structure_change,self.CVODE_flag,self.CVODE_error_code

    def apply_external_perturbation_assuming_integrated_orbits(self):
        self.__update_particles_in_code()
        self.lib.apply_external_perturbation_assuming_integrated_orbits_interface()
        self.__update_particles_from_code()
        
    def apply_user_specified_instantaneous_perturbation(self):
        self.__update_particles_in_code(set_instantaneous_perturbation_properties=True)
        self.lib.apply_user_specified_instantaneous_perturbation_interface()
        self.__update_particles_from_code()

    def __update_particle_in_code(self,particle,set_instantaneous_perturbation_properties=False):
        flag = 0
        if particle.is_binary == False:
            flag = self.lib.set_mass(particle.index,particle.mass)

#        if self.enable_tides == True:
#            particle.include_tidal_friction_terms = True
#            particle.tides_method = 1
#            particle.include_tidal_bulges_precession_terms = True
#            particle.include_rotation_precession_terms = True
        if self.enable_tides == False:
            particle.include_tidal_friction_terms = False
            particle.tides_method = 1
            particle.include_tidal_bulges_precession_terms = False
            particle.include_rotation_precession_terms = False
            
        flag += self.lib.set_tides_terms(particle.index,particle.include_tidal_friction_terms,particle.tides_method,particle.include_tidal_bulges_precession_terms,particle.include_rotation_precession_terms, \
            particle.minimum_eccentricity_for_tidal_precession)
            
        if self.enable_root_finding == False:
            particle.check_for_secular_breakdown = False
            particle.check_for_dynamical_instability = False
            particle.check_for_physical_collision_or_orbit_crossing = False
            particle.check_for_RLOF_at_pericentre = False
            particle.check_for_GW_condition = False
            
        flag += self.lib.set_root_finding_terms(particle.index,particle.check_for_secular_breakdown,particle.check_for_dynamical_instability,particle.dynamical_instability_criterion,particle.dynamical_instability_central_particle,particle.dynamical_instability_K_parameter, \
                particle.check_for_physical_collision_or_orbit_crossing,particle.check_for_minimum_periapse_distance,particle.check_for_minimum_periapse_distance_value,particle.check_for_RLOF_at_pericentre,particle.check_for_RLOF_at_pericentre_use_sepinsky_fit,particle.check_for_GW_condition)
        flag += self.lib.set_root_finding_state(particle.index,particle.secular_breakdown_has_occurred,particle.dynamical_instability_has_occurred, \
                particle.physical_collision_or_orbit_crossing_has_occurred,particle.minimum_periapse_distance_has_occurred,particle.RLOF_at_pericentre_has_occurred,particle.GW_condition_has_occurred)

        if self.enable_VRR == True:
            flag += self.lib.set_VRR_properties(particle.index,particle.VRR_model,particle.VRR_include_mass_precession,particle.VRR_mass_precession_rate, \
                particle.VRR_Omega_vec_x,particle.VRR_Omega_vec_y,particle.VRR_Omega_vec_z, \
                particle.VRR_eta_20_init,particle.VRR_eta_a_22_init,particle.VRR_eta_b_22_init,particle.VRR_eta_a_21_init,particle.VRR_eta_b_21_init, \
                particle.VRR_eta_20_final,particle.VRR_eta_a_22_final,particle.VRR_eta_b_22_final,particle.VRR_eta_a_21_final,particle.VRR_eta_b_21_final,particle.VRR_initial_time,particle.VRR_final_time)

        if particle.is_external==False:
            
            if particle.is_binary==True:
                flag += self.lib.set_children(particle.index,particle.child1.index,particle.child2.index)
                flag += self.lib.set_orbital_elements(particle.index,particle.a, particle.e, particle.TA, particle.INCL, particle.AP, particle.LAN, particle.sample_orbital_phase_randomly)
                flag += self.lib.set_PN_terms(particle.index,particle.include_pairwise_1PN_terms,particle.include_pairwise_25PN_terms)
                flag += self.lib.set_integration_method(particle.index,particle.integration_method,particle.KS_use_perturbing_potential)
            else:
                flag += self.lib.set_radius(particle.index,particle.radius,particle.radius_dot)
                #flag += self.lib.set_mass_dot(particle.index,particle.mass_dot)
                flag += self.lib.set_mass_transfer_terms(particle.index,particle.include_mass_transfer_terms)
                flag += self.lib.set_spin_vector(particle.index,particle.spin_vec_x,particle.spin_vec_y,particle.spin_vec_z)
                flag += self.lib.set_stellar_evolution_properties(particle.index,particle.stellar_type,particle.evolve_as_star,particle.sse_initial_mass,particle.metallicity,particle.sse_time_step,particle.epoch,particle.age, \
                    particle.convective_envelope_mass,particle.convective_envelope_radius,particle.core_mass,particle.core_radius,particle.luminosity,particle.apsidal_motion_constant,particle.gyration_radius,particle.tides_viscous_time_scale,particle.tides_viscous_time_scale_prescription)
                flag += self.lib.set_kick_properties(particle.index,particle.kick_distribution,particle.kick_distribution_sigma_km_s_NS,particle.kick_distribution_sigma_km_s_BH)

                if set_instantaneous_perturbation_properties==True:
                    flag += self.lib.set_instantaneous_perturbation_properties(particle.index,particle.instantaneous_perturbation_delta_mass, \
                        particle.instantaneous_perturbation_delta_X,particle.instantaneous_perturbation_delta_Y,particle.instantaneous_perturbation_delta_Z, \
                        particle.instantaneous_perturbation_delta_VX,particle.instantaneous_perturbation_delta_VY,particle.instantaneous_perturbation_delta_VZ)
                        
        else:
            flag += self.lib.set_external_particle_properties(particle.index, particle.external_t_ref, particle.e, particle.external_r_p, particle.INCL, particle.AP, particle.LAN)
    
        return flag

#    def __update_particles_in_code(self,set_instantaneous_perturbation_properties=False):
#        flag = 0
#        for index,particle in enumerate(self.particles):
#            flag += self.__update_particle_in_code(particle,set_instantaneous_perturbation_properties=set_instantaneous_perturbation_properties)
#        return flag
        
    def __update_particles_in_code(self,set_instantaneous_perturbation_properties=False):
        flag = 0
        for index,particle in enumerate(self.particles):
            if particle.is_binary==True:
                flag += self.lib.set_children(particle.index,particle.child1.index,particle.child2.index)
        
        flag = 0
        for index,particle in enumerate(self.particles):
            flag += self.__update_particle_in_code(particle,set_instantaneous_perturbation_properties=set_instantaneous_perturbation_properties)
        return flag

    def __update_particle_from_code(self,particle):
        mass = ctypes.c_double(0.0)
        flag = self.lib.get_mass(particle.index,ctypes.byref(mass))
        particle.mass = mass.value

        if self.enable_root_finding == True:
            secular_breakdown_has_occurred,dynamical_instability_has_occurred,physical_collision_or_orbit_crossing_has_occurred,minimum_periapse_distance_has_occurred,RLOF_at_pericentre_has_occurred,GW_condition_has_occurred = ctypes.c_bool(False),ctypes.c_bool(False),ctypes.c_bool(False),ctypes.c_bool(False),ctypes.c_bool(False),ctypes.c_bool(False)
            flag += self.lib.get_root_finding_state(particle.index,ctypes.byref(secular_breakdown_has_occurred),ctypes.byref(dynamical_instability_has_occurred), \
                ctypes.byref(physical_collision_or_orbit_crossing_has_occurred),ctypes.byref(minimum_periapse_distance_has_occurred),ctypes.byref(RLOF_at_pericentre_has_occurred),ctypes.byref(GW_condition_has_occurred))
            particle.secular_breakdown_has_occurred = secular_breakdown_has_occurred.value
            particle.dynamical_instability_has_occurred = dynamical_instability_has_occurred.value
            particle.physical_collision_or_orbit_crossing_has_occurred = physical_collision_or_orbit_crossing_has_occurred.value
            particle.minimum_periapse_distance_has_occurred = minimum_periapse_distance_has_occurred.value
            particle.RLOF_at_pericentre_has_occurred = RLOF_at_pericentre_has_occurred.value
            particle.GW_condition_has_occurred = GW_condition_has_occurred.value

        if particle.is_binary==True:
            a,e,TA,INCL,AP,LAN = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
            flag += self.lib.get_orbital_elements(particle.index,ctypes.byref(a),ctypes.byref(e),ctypes.byref(TA),ctypes.byref(INCL),ctypes.byref(AP),ctypes.byref(LAN))
            particle.a = a.value
            particle.e = e.value
            particle.TA = TA.value
            particle.INCL = INCL.value
            particle.AP = AP.value
            particle.LAN = LAN.value
            
            INCL_parent = ctypes.c_double(0.0)
            flag += self.lib.get_inclination_relative_to_parent_interface(particle.index,ctypes.byref(INCL_parent))
            particle.INCL_parent = INCL_parent.value
            
            x,y,z,vx,vy,vz = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
            flag = self.lib.get_relative_position_and_velocity(particle.index,ctypes.byref(x),ctypes.byref(y),ctypes.byref(z),ctypes.byref(vx),ctypes.byref(vy),ctypes.byref(vz))
            particle.x = x.value
            particle.y = y.value
            particle.z = z.value
            particle.vx = vx.value
            particle.vy = vy.value
            particle.vz = vz.value
        else:
            X,Y,Z,VX,VY,VZ = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
            flag = self.lib.get_absolute_position_and_velocity(particle.index,ctypes.byref(X),ctypes.byref(Y),ctypes.byref(Z),ctypes.byref(VX),ctypes.byref(VY),ctypes.byref(VZ))
            particle.X = X.value
            particle.Y = Y.value
            particle.Z = Z.value
            particle.VX = VX.value
            particle.VY = VY.value
            particle.VZ = VZ.value

            if particle.is_external==False:
                is_bound = ctypes.c_bool(True)
                flag += self.lib.get_is_bound(particle.index,ctypes.byref(is_bound))
                particle.is_bound = is_bound.value
                
                radius,radius_dot = ctypes.c_double(0.0),ctypes.c_double(0.0)
                flag += self.lib.get_radius(particle.index,ctypes.byref(radius),ctypes.byref(radius_dot))
                particle.radius = radius.value
                particle.radius_dot = radius_dot.value

                stellar_type,evolve_as_star,sse_initial_mass,metallicity,sse_time_step,epoch,age,convective_envelope_mass,convective_envelope_radius,core_mass,core_radius,luminosity,apsidal_motion_constant,gyration_radius,tides_viscous_time_scale,roche_lobe_radius_pericenter = ctypes.c_int(0),ctypes.c_bool(False), \
                    ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                flag += self.lib.get_stellar_evolution_properties(particle.index,ctypes.byref(stellar_type),ctypes.byref(evolve_as_star),ctypes.byref(sse_initial_mass),ctypes.byref(metallicity),ctypes.byref(sse_time_step), \
                    ctypes.byref(epoch),ctypes.byref(age),ctypes.byref(convective_envelope_mass),ctypes.byref(convective_envelope_radius),ctypes.byref(core_mass),ctypes.byref(core_radius),ctypes.byref(luminosity),ctypes.byref(apsidal_motion_constant),ctypes.byref(gyration_radius),ctypes.byref(tides_viscous_time_scale),ctypes.byref(roche_lobe_radius_pericenter))
                particle.stellar_type = stellar_type.value
                particle.evolve_as_star = evolve_as_star.value
                particle.sse_initial_mass = sse_initial_mass.value
                particle.metallicity = metallicity.value
                particle.sse_time_step = sse_time_step.value
                particle.epoch = epoch.value
                particle.age = age.value
                particle.convective_envelope_mass = convective_envelope_mass.value
                particle.convective_envelope_radius = convective_envelope_radius.value
                particle.core_mass = core_mass.value
                particle.core_radius = core_radius.value
                particle.luminosity = luminosity.value
                particle.apsidal_motion_constant = apsidal_motion_constant
                particle.gyration_radius = gyration_radius
                particle.tides_viscous_time_scale = tides_viscous_time_scale.value
                particle.roche_lobe_radius_pericenter = roche_lobe_radius_pericenter.value

                kick_distribution,kick_distribution_sigma_km_s_NS,kick_distribution_sigma_km_s_BH = ctypes.c_int(0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                flag += self.lib.get_kick_properties(particle.index,ctypes.byref(kick_distribution),ctypes.byref(kick_distribution_sigma_km_s_NS),ctypes.byref(kick_distribution_sigma_km_s_BH))
                particle.kick_distribution = kick_distribution.value
                particle.kick_distribution_sigma_km_s_NS = kick_distribution_sigma_km_s_NS.value
                particle.kick_distribution_sigma_km_s_BH = kick_distribution_sigma_km_s_BH.value

                mass_dot = ctypes.c_double(0.0)
                flag = self.lib.get_mass_dot(particle.index,ctypes.byref(mass_dot))
                particle.mass_dot = mass_dot.value

                spin_vec_x,spin_vec_y,spin_vec_z = ctypes.c_double(0.0), ctypes.c_double(0.0), ctypes.c_double(0.0)
                flag += self.lib.get_spin_vector(particle.index,ctypes.byref(spin_vec_x),ctypes.byref(spin_vec_y),ctypes.byref(spin_vec_z))
                particle.spin_vec_x = spin_vec_x.value
                particle.spin_vec_y = spin_vec_y.value
                particle.spin_vec_z = spin_vec_z.value

        return flag
        
    def __update_particles_from_code(self):
        
        flag = 0
        for index,particle in enumerate(self.particles):
            flag += self.__update_particle_from_code(particle)
        return flag

    def __copy_particle_structure_from_code(self):
        self.particles = []
        #N_particles = ctypes.c_int(0)
        #flag = self.lib.get_number_of_particles(ctypes.byref(N_particles))
        #N_particles = N_particles.value
        N_particles = self.lib.get_number_of_particles()
        
        for i in range(N_particles):
            
            internal_index = self.lib.get_internal_index_in_particlesMap(i)
            is_binary = self.lib.get_is_binary(internal_index)

            mass = ctypes.c_double(0.0)
            flag = self.lib.get_mass(internal_index,ctypes.byref(mass))
            mass = mass.value

            child1_index,child2_index = -1,-1
            if is_binary==True:
                child1,child2 = ctypes.c_int(0),ctypes.c_int(0)
                self.lib.get_children(internal_index,ctypes.byref(child1),ctypes.byref(child2))
                #child1 = self.particles[child1.value]
                #child2 = self.particles[child2.value]
                child1_index = child1.value
                child2_index = child2.value

            #print("isb",N_particles,"i",i,"internal_index",internal_index,"mass",mass,"is_binary",is_binary,"child1_index",child1_index,"child2_index",child2_index)#,child1,child2)
            p = Particle(is_binary=is_binary,mass=mass,child1=None,child2=None,a=0.0,e=0.0,INCL=0.0,AP=0.0,LAN=0.0) ### orbital elements should be updated later
            p.index = internal_index
                #a==None or e==None or INCL==None or LAN==None:
            #if is_binary==True:
            p.child1_index = child1_index
            p.child2_index = child2_index
            self.particles.append(p)
        #print("D1",[p.child1_index for p in self.particles])
        #print("D2",[p.child2_index for p in self.particles])
        binaries = [x for x in self.particles if x.is_binary == True]
        #if len(binaries)==0:
        #    exit(0)
        for i,p in enumerate(self.particles):
            if p.is_binary == True:
                i1 = [j for j in range(N_particles) if self.particles[j].index == p.child1_index][0]
                i2 = [j for j in range(N_particles) if self.particles[j].index == p.child2_index][0]
                #print("i1",i1,"i2",i2,"i",i,"p.child1_index",p.child1_index,"p.child2_index",p.child2_index)
                p.child1 = self.particles[i1]
                p.child2 = self.particles[i2]

    def __set_constants_in_code(self):
        self.lib.set_constants(self.__CONST_G,self.__CONST_C,self.__CONST_M_SUN,self.__CONST_R_SUN,self.__CONST_L_SUN,self.__CONST_KM_PER_S,self.__CONST_PER_PC3)


    def __set_parameters_in_code(self):
         self.lib.set_parameters(self.__relative_tolerance,self.__absolute_tolerance_eccentricity_vectors,self.__include_quadrupole_order_terms, \
             self.__include_octupole_order_binary_pair_terms,self.__include_octupole_order_binary_triplet_terms, \
             self.__include_hexadecupole_order_binary_pair_terms,self.__include_dotriacontupole_order_binary_pair_terms, self.__include_double_averaging_corrections, \
             self.__include_flybys, self.__flybys_reference_binary, self.__flybys_correct_for_gravitational_focussing, self.__flybys_velocity_distribution, self.__flybys_mass_distribution, \
             self.__flybys_mass_distribution_lower_value, self.__flybys_mass_distribution_upper_value, self.__flybys_encounter_sphere_radius, \
             self.__flybys_stellar_density, self.__flybys_stellar_relative_velocity_dispersion, \
             self.__binary_evolution_CE_energy_flag, self.__binary_evolution_CE_spin_flag, \
             self.__mstar_gbs_tolerance_default, self.__mstar_gbs_tolerance_kick, self.__mstar_collision_tolerance)

    def reset(self):
        self.__init__()
        self.lib.reset_interface()
        
    def __set_random_seed(self):
        self.lib.set_random_seed(self.random_seed)

    ### Logging ###
    def __get_log(self):
        N_log = self.lib.get_size_of_log_data()
        #print("get_log",N_log)
        log = []
        for index_log in range(N_log):
            entry = {}
            particles = []
            N_particles = self.lib.get_number_of_particles()
            
            time = ctypes.c_double(0.0)
            event_flag = ctypes.c_int(0)
            flag = self.lib.get_log_entry_properties(index_log,ctypes.byref(time),ctypes.byref(event_flag))
            entry.update({'time':time.value,'event_flag':event_flag.value})
            
            for index_particle in range(N_particles):
                internal_index = self.lib.get_internal_index_in_particlesMap(index_particle)

                parent,stellar_type = ctypes.c_int(0),ctypes.c_int(0)
                mass,radius,core_mass,sse_initial_mass,convective_envelope_mass,epoch,age,core_radius,convective_envelope_radius,luminosity,ospin = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                flag = self.lib.get_body_properties_from_log_entry(index_log,internal_index,ctypes.byref(parent),ctypes.byref(mass),ctypes.byref(radius),ctypes.byref(stellar_type), \
                    ctypes.byref(core_mass),ctypes.byref(sse_initial_mass),ctypes.byref(convective_envelope_mass), \
                    ctypes.byref(epoch),ctypes.byref(age), \
                    ctypes.byref(core_radius),ctypes.byref(convective_envelope_radius),ctypes.byref(luminosity),ctypes.byref(ospin))

                p = Particle(is_binary=False,mass=mass.value,radius=radius.value,stellar_type=stellar_type.value,core_mass=core_mass.value,sse_initial_mass=sse_initial_mass.value, \
                    convective_envelope_mass=convective_envelope_mass.value, epoch=epoch.value, age=age.value, core_radius=core_radius.value, convective_envelope_radius=convective_envelope_radius.value, luminosity=luminosity.value)
                p.index = internal_index
                p.parent = parent.value
                p.ospin = ospin.value

                particles.append(p)

            entry.update({'particles':particles})
            log.append(entry)
                
        return log

    @property
    def log(self):
        return self.__get_log()


    ### Tests ###
    def unit_tests(self):
        return self.lib.unit_tests_interface()

    def determine_compact_object_merger_properties(self,m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit):
        v_recoil_vec_x,v_recoil_vec_y,v_recoil_vec_z = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
        alpha_vec_final_x,alpha_vec_final_y,alpha_vec_final_z = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
        M_final = ctypes.c_double(0.0)
        self.lib.determine_compact_object_merger_properties_interface( m1,m2,chi1,chi2,spin_vec_1_unit[0],spin_vec_1_unit[1],spin_vec_1_unit[2],spin_vec_2_unit[0],spin_vec_2_unit[1],spin_vec_2_unit[2], \
            h_vec_unit[0],h_vec_unit[1],h_vec_unit[2],e_vec_unit[0], e_vec_unit[1], e_vec_unit[2], \
            ctypes.byref(v_recoil_vec_x),ctypes.byref(v_recoil_vec_y),ctypes.byref(v_recoil_vec_z), \
            ctypes.byref(alpha_vec_final_x),ctypes.byref(alpha_vec_final_y),ctypes.byref(alpha_vec_final_z), \
            ctypes.byref(M_final) )
        v_recoil_vec = np.array( [v_recoil_vec_x.value,v_recoil_vec_y.value,v_recoil_vec_z.value] )
        alpha_vec_final = np.array( [alpha_vec_final_x.value,alpha_vec_final_y.value,alpha_vec_final_z.value] )
        M_final = M_final.value

        return v_recoil_vec,alpha_vec_final,M_final
        
    def test_sample_from_3d_maxwellian_distribution(self,sigma):
        vx,vy,vz = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
        self.lib.sample_from_3d_maxwellian_distribution_interface(sigma,ctypes.byref(vx),ctypes.byref(vy),ctypes.byref(vz))
        vx,vy,vz = vx.value, vy.value, vz.value
        return vx,vy,vz

    def test_sample_from_normal_distribution(self,mu,sigma):
        v = self.lib.sample_from_normal_distribution_interface(mu,sigma)
        return v

    def test_sample_from_kroupa_93_imf(self):
        m = self.lib.sample_from_kroupa_93_imf_interface()
        return m

    def test_kick_velocity(self,kick_distribution,m):
        kw,v = ctypes.c_int(0),ctypes.c_double(0.0)
        self.lib.test_kick_velocity(kick_distribution,m,ctypes.byref(kw),ctypes.byref(v))
        return kw.value,v.value


    ### Constants ###
    @property
    def CONST_G(self):
        return self.__CONST_G

    @CONST_G.setter
    def CONST_G(self, value):
        self.__CONST_G = value
        self.__set_constants_in_code()
        
    @property
    def CONST_C(self):
        return self.__CONST_C

    @CONST_C.setter
    def CONST_C(self, value):
        self.__CONST_C = value
        self.__set_constants_in_code()

    @property
    def CONST_M_SUN(self):
        return self.__CONST_M_SUN

    @CONST_M_SUN.setter
    def CONST_M_SUN(self, value):
        self.__CONST_M_SUN = value
        self.__set_constants_in_code()

    @property
    def CONST_L_SUN(self):
        return self.__CONST_L_SUN

    @CONST_L_SUN.setter
    def CONST_L_SUN(self, value):
        self.__CONST_L_SUN = value
        self.__set_constants_in_code()

    @property
    def CONST_R_SUN(self):
        return self.__CONST_R_SUN

    @CONST_R_SUN.setter
    def CONST_R_SUN(self, value):
        self.__CONST_R_SUN = value
        self.__set_constants_in_code()

    @property
    def CONST_KM_PER_S(self):
        return self.__CONST_KM_PER_S

    @CONST_KM_PER_S.setter
    def CONST_KM_PER_S(self, value):
        self.__CONST_KM_PER_S = value
        self.__set_constants_in_code()

    @property
    def CONST_PER_PC3(self):
        return self.__CONST_PER_PC3

    @CONST_PER_PC3.setter
    def CONST_PER_PC3(self, value):
        self.__CONST_PER_PC3 = value
        self.__set_constants_in_code()

    ##################
    ### Parameters ###
    ##################
    
    @property
    def relative_tolerance(self):
        return self.__relative_tolerance

    @relative_tolerance.setter
    def relative_tolerance(self, value):
        self.__relative_tolerance = value
        self.__set_parameters_in_code()

    @property
    def absolute_tolerance_eccentricity_vectors(self):
        return self.__absolute_tolerance_eccentricity_vectors

    @absolute_tolerance_eccentricity_vectors.setter
    def absolute_tolerance_eccentricity_vectors(self, value):
        self.__absolute_tolerance_eccentricity_vectors = value
        self.__set_parameters_in_code()

    @property
    def include_quadrupole_order_terms(self):
        return self.__include_quadrupole_order_terms

    @include_quadrupole_order_terms.setter
    def include_quadrupole_order_terms(self, value):
        self.__include_quadrupole_order_terms = value
        self.__set_parameters_in_code()

    @property
    def include_octupole_order_binary_pair_terms(self):
        return self.__include_octupole_order_binary_pair_terms

    @include_octupole_order_binary_pair_terms.setter
    def include_octupole_order_binary_pair_terms(self, value):
        self.__include_octupole_order_binary_pair_terms = value
        self.__set_parameters_in_code()

    @property
    def include_octupole_order_binary_triplet_terms(self):
        return self.__include_octupole_order_binary_triplet_terms

    @include_octupole_order_binary_triplet_terms.setter
    def include_octupole_order_binary_triplet_terms(self, value):
        self.__include_octupole_order_binary_triplet_terms = value
        self.__set_parameters_in_code()

    @property
    def include_hexadecupole_order_binary_pair_terms(self):
        return self.__include_hexadecupole_order_binary_pair_terms

    @include_hexadecupole_order_binary_pair_terms.setter
    def include_hexadecupole_order_binary_pair_terms(self, value):
        self.__include_hexadecupole_order_binary_pair_terms = value
        self.__set_parameters_in_code()

    @property
    def include_dotriacontupole_order_binary_pair_terms(self):
        return self.__include_dotriacontupole_order_binary_pair_terms

    @include_dotriacontupole_order_binary_pair_terms.setter
    def include_dotriacontupole_order_binary_pair_terms(self, value):
        self.__include_dotriacontupole_order_binary_pair_terms = value
        self.__set_parameters_in_code()

    @property
    def include_double_averaging_corrections(self):
        return self.__include_double_averaging_corrections

    @include_double_averaging_corrections.setter
    def include_double_averaging_corrections(self, value):
        self.__include_double_averaging_corrections = value
        self.__set_parameters_in_code()

#    @property
#    def include_VRR(self):
#        return self.__include_VRR
    
#    @include_VRR.setter
#    def include_VRR(self, value):
#        self.__include_VRR = value
#        self.__set_parameters_in_code()


    @property
    def random_seed(self):
        return self.__random_seed

    @random_seed.setter
    def random_seed(self, value):
        self.__random_seed = value
        self.__set_random_seed()
        

    ### Flybys ###
    @property
    def include_flybys(self):
        return self.__include_flybys
    
    @include_flybys.setter
    def include_flybys(self, value):
        self.__include_flybys = value
        self.__set_parameters_in_code()

    @property
    def flybys_reference_binary(self):
        return self.__flybys_reference_binary
    
    @flybys_reference_binary.setter
    def flybys_reference_binary(self, value):
        self.__flybys_reference_binary = value
        self.__set_parameters_in_code()
        
    @property
    def flybys_correct_for_gravitational_focussing(self):
        return self.__flybys_correct_for_gravitational_focussing
    
    @flybys_correct_for_gravitational_focussing.setter
    def flybys_correct_for_gravitational_focussing(self, value):
        self.__flybys_correct_for_gravitational_focussing = value
        self.__set_parameters_in_code()

    @property
    def flybys_velocity_distribution(self):
        return self.__flybys_velocity_distribution
    
    @flybys_velocity_distribution.setter
    def flybys_velocity_distribution(self, value):
        self.__flybys_velocity_distribution = value
        self.__set_parameters_in_code()

    @property
    def flybys_mass_distribution(self):
        return self.__flybys_mass_distribution
    
    @flybys_mass_distribution.setter
    def flybys_mass_distribution(self, value):
        self.__flybys_mass_distribution = value
        self.__set_parameters_in_code()

    @property
    def flybys_mass_distribution_lower_value(self):
        return self.__flybys_mass_distribution_lower_value
    
    @flybys_mass_distribution_lower_value.setter
    def flybys_mass_distribution_lower_value(self, value):
        self.__flybys_mass_distribution_lower_value = value
        self.__set_parameters_in_code()

    @property
    def flybys_mass_distribution_upper_value(self):
        return self.__flybys_mass_distribution_upper_value
    
    @flybys_mass_distribution_upper_value.setter
    def flybys_mass_distribution_upper_value(self, value):
        self.__flybys_mass_distribution_upper_value = value
        self.__set_parameters_in_code()

    @property
    def flybys_encounter_sphere_radius(self):
        return self.__flybys_encounter_sphere_radius
    
    @flybys_encounter_sphere_radius.setter
    def flybys_encounter_sphere_radius(self, value):
        self.__flybys_encounter_sphere_radius = value
        self.__set_parameters_in_code()

    @property
    def flybys_stellar_density(self):
        return self.__flybys_stellar_density
    
    @flybys_stellar_density.setter
    def flybys_stellar_density(self, value):
        self.__flybys_stellar_density = value
        self.__set_parameters_in_code()
        
    @property
    def flybys_stellar_relative_velocity_dispersion(self):
        return self.__flybys_stellar_relative_velocity_dispersion
    
    @flybys_stellar_relative_velocity_dispersion.setter
    def flybys_stellar_relative_velocity_dispersion(self, value):
        self.__flybys_stellar_relative_velocity_dispersion = value
        self.__set_parameters_in_code()


    ### Binary evolution ###
    @property
    def binary_evolution_CE_energy_flag(self):
        return self.__binary_evolution_CE_energy_flag
    
    @binary_evolution_CE_energy_flag.setter
    def binary_evolution_CE_energy_flag(self, value):
        self.__binary_evolution_CE_energy_flag = value
        self.__set_parameters_in_code()
        
    @property
    def binary_evolution_CE_spin_flag(self):
        return self.__binary_evolution_CE_spin_flag
    
    @binary_evolution_CE_spin_flag.setter
    def binary_evolution_CE_spin_flag(self, value):
        self.__binary_evolution_CE_spin_flag = value
        self.__set_parameters_in_code()
        
    ### N-body ###
    @property
    def mstar_gbs_tolerance(self):
        return self.__mstar_gbs_tolerance
    
    @binary_evolution_CE_spin_flag.setter
    def mstar_gbs_tolerance(self, value):
        self.__mstar_gbs_tolerance = value
        self.__set_parameters_in_code()

    @property
    def mstar_collision_tolerance(self):
        return self.__mstar_collision_tolerance
    
    @binary_evolution_CE_spin_flag.setter
    def mstar_collision_tolerance(self, value):
        self.__mstar_collision_tolerance = value
        self.__set_parameters_in_code()
    

class Particle(object):
    def __init__(self, is_binary, mass=None, mass_dot=0.0, radius=1.0, radius_dot=0.0, child1=None, child2=None, a=None, e=None, TA=0.0, INCL=None, AP=None, LAN=None, \
            integration_method = 0, KS_use_perturbing_potential = True, \
            stellar_type=1, evolve_as_star=True, sse_initial_mass=None, metallicity=0.02, sse_time_step=1.0, epoch=0.0, age=0.0, core_mass=0.0, core_radius=0.0, \
            include_mass_transfer_terms=True, \
            kick_distribution = 1, kick_distribution_sigma_km_s_NS = 265.0, kick_distribution_sigma_km_s_BH=50.0, \
            spin_vec_x=0.0, spin_vec_y=0.0, spin_vec_z=1.0e-10, \
            include_pairwise_1PN_terms=True, include_pairwise_25PN_terms=True, \
            include_tidal_friction_terms=True, tides_method=1, include_tidal_bulges_precession_terms=True, include_rotation_precession_terms=True, \
            minimum_eccentricity_for_tidal_precession = 1.0e-3, apsidal_motion_constant=0.19, gyration_radius=0.08, tides_viscous_time_scale=1.0, tides_viscous_time_scale_prescription=1, \
            convective_envelope_mass=1.0, convective_envelope_radius=1.0, luminosity=1.0, \
            check_for_secular_breakdown=False,check_for_dynamical_instability=True,dynamical_instability_criterion=0,dynamical_instability_central_particle=0,dynamical_instability_K_parameter=0, \
            check_for_physical_collision_or_orbit_crossing=True,check_for_minimum_periapse_distance=False,check_for_minimum_periapse_distance_value=0.0,check_for_RLOF_at_pericentre=True,check_for_RLOF_at_pericentre_use_sepinsky_fit=False, check_for_GW_condition=False, \
            secular_breakdown_has_occurred=False, dynamical_instability_has_occurred=False, physical_collision_or_orbit_crossing_has_occurred=False, minimum_periapse_distance_has_occurred=False, RLOF_at_pericentre_has_occurred = False, GW_condition_has_occurred = False, \
            is_external=False, external_t_ref=0.0, external_r_p=0.0, \
            sample_orbital_phase_randomly=True, instantaneous_perturbation_delta_mass=0.0, instantaneous_perturbation_delta_X=0.0, instantaneous_perturbation_delta_Y=0.0, instantaneous_perturbation_delta_Z=0.0, \
            instantaneous_perturbation_delta_VX=0.0, instantaneous_perturbation_delta_VY=0.0, instantaneous_perturbation_delta_VZ=0.0, \
            VRR_model=0, VRR_include_mass_precession=0, VRR_mass_precession_rate=0.0, VRR_Omega_vec_x=0.0, VRR_Omega_vec_y=0.0, VRR_Omega_vec_z=0.0, \
            VRR_eta_20_init=0.0, VRR_eta_a_22_init=0.0, VRR_eta_b_22_init=0.0, VRR_eta_a_21_init=0.0, VRR_eta_b_21_init=0.0, \
            VRR_eta_20_final=0.0, VRR_eta_a_22_final=0.0, VRR_eta_b_22_final=0.0, VRR_eta_a_21_final=0.0, VRR_eta_b_21_final=0.0, \
            VRR_initial_time = 0.0, VRR_final_time = 1.0,roche_lobe_radius_pericenter=0.0):
                
                ### TO DO: remove default values for check_for_... here

        ### spin_vec: nonzero spin_vec_z: need to specify a finite initial direction 

        if is_binary==None:
            raise RuntimeError('Error when adding particle: particle should have property is_binary')

        self.is_external = is_external
        if is_external==True:
            is_binary = False ### for is_binary to true for external particles
            self.external_t_ref = external_t_ref
            self.external_r_p = external_r_p
            self.mass = mass
            self.e = e
            self.INCL = INCL
            self.AP = AP
            self.LAN = LAN


        self.index = None
        self.is_binary = is_binary

        self.include_tidal_friction_terms=include_tidal_friction_terms
        self.tides_method=tides_method
        self.include_tidal_bulges_precession_terms=include_tidal_bulges_precession_terms
        self.include_rotation_precession_terms=include_rotation_precession_terms
        self.minimum_eccentricity_for_tidal_precession=minimum_eccentricity_for_tidal_precession
        self.tides_viscous_time_scale=tides_viscous_time_scale
        self.tides_viscous_time_scale_prescription=tides_viscous_time_scale_prescription

        self.include_mass_transfer_terms = include_mass_transfer_terms

        self.kick_distribution = kick_distribution
        self.kick_distribution_sigma_km_s_NS = kick_distribution_sigma_km_s_NS
        self.kick_distribution_sigma_km_s_BH = kick_distribution_sigma_km_s_BH
          
        self.check_for_secular_breakdown=check_for_secular_breakdown
        self.check_for_dynamical_instability=check_for_dynamical_instability
        self.dynamical_instability_criterion=dynamical_instability_criterion
        self.dynamical_instability_central_particle=dynamical_instability_central_particle
        self.dynamical_instability_K_parameter=dynamical_instability_K_parameter
        self.check_for_physical_collision_or_orbit_crossing=check_for_physical_collision_or_orbit_crossing
        self.check_for_minimum_periapse_distance=check_for_minimum_periapse_distance
        self.check_for_minimum_periapse_distance_value=check_for_minimum_periapse_distance_value
        self.check_for_RLOF_at_pericentre=check_for_RLOF_at_pericentre
        self.check_for_RLOF_at_pericentre_use_sepinsky_fit=check_for_RLOF_at_pericentre_use_sepinsky_fit
        self.check_for_GW_condition=check_for_GW_condition

        self.secular_breakdown_has_occurred=secular_breakdown_has_occurred
        self.dynamical_instability_has_occurred=dynamical_instability_has_occurred
        self.physical_collision_or_orbit_crossing_has_occurred=physical_collision_or_orbit_crossing_has_occurred
        self.minimum_periapse_distance_has_occurred=minimum_periapse_distance_has_occurred
        self.RLOF_at_pericentre_has_occurred=RLOF_at_pericentre_has_occurred
        self.GW_condition_has_occurred=GW_condition_has_occurred

        self.sample_orbital_phase_randomly=sample_orbital_phase_randomly
        self.instantaneous_perturbation_delta_mass=instantaneous_perturbation_delta_mass
        self.instantaneous_perturbation_delta_X=instantaneous_perturbation_delta_X
        self.instantaneous_perturbation_delta_Y=instantaneous_perturbation_delta_Y
        self.instantaneous_perturbation_delta_Z=instantaneous_perturbation_delta_Z
        self.instantaneous_perturbation_delta_VX=instantaneous_perturbation_delta_VX
        self.instantaneous_perturbation_delta_VY=instantaneous_perturbation_delta_VY
        self.instantaneous_perturbation_delta_VZ=instantaneous_perturbation_delta_VZ

        self.VRR_model = VRR_model
        self.VRR_include_mass_precession = VRR_include_mass_precession
        self.VRR_mass_precession_rate = VRR_mass_precession_rate
        self.VRR_Omega_vec_x = VRR_Omega_vec_x
        self.VRR_Omega_vec_y = VRR_Omega_vec_y
        self.VRR_Omega_vec_z = VRR_Omega_vec_z
        self.VRR_eta_20_init = VRR_eta_20_init
        self.VRR_eta_a_22_init = VRR_eta_a_22_init
        self.VRR_eta_b_22_init = VRR_eta_b_22_init
        self.VRR_eta_a_21_init = VRR_eta_a_21_init
        self.VRR_eta_b_21_init = VRR_eta_b_21_init
        self.VRR_eta_20_final = VRR_eta_20_final
        self.VRR_eta_a_22_final = VRR_eta_a_22_final
        self.VRR_eta_b_22_final = VRR_eta_b_22_final
        self.VRR_eta_a_21_final = VRR_eta_a_21_final
        self.VRR_eta_b_21_final = VRR_eta_b_21_final
        self.VRR_initial_time = VRR_initial_time
        self.VRR_final_time = VRR_final_time

        if is_binary==False:
            if mass==None:
                raise RuntimeError('Error when adding particle: body should have mass specified') 
            self.mass = mass
            self.mass_dot = mass_dot
            self.stellar_type = stellar_type
            self.evolve_as_star = evolve_as_star
            self.sse_initial_mass = mass
            self.metallicity = metallicity
            self.sse_time_step = sse_time_step
            self.epoch = epoch
            self.age = age
            self.core_mass = core_mass
            self.core_radius = core_radius
            self.child1 = None
            self.child2 = None
            self.radius = radius
            self.radius_dot = radius_dot
            self.spin_vec_x = spin_vec_x
            self.spin_vec_y = spin_vec_y
            self.spin_vec_z = spin_vec_z
            self.apsidal_motion_constant=apsidal_motion_constant
            self.gyration_radius=gyration_radius
            self.convective_envelope_mass=convective_envelope_mass
            self.convective_envelope_radius=convective_envelope_radius
            self.luminosity=luminosity
            self.roche_lobe_radius_pericenter = roche_lobe_radius_pericenter
        
        else:
            if is_external==False:
                #if child1==None or child2==None:
                #    raise RuntimeError('Error when adding particle: a binary should have two children!')
                if a==None or e==None or INCL==None or LAN==None:
                    raise RuntimeError('Error when adding particle: a binary should have its orbital elements specified!')
                else:
                    self.child1 = child1
                    self.child2 = child2
                    #self.mass = child1.mass + child2.mass

                    self.a = a
                    self.e = e
                    self.TA = TA
                    self.INCL = INCL
                    self.AP = AP
                    self.LAN = LAN
                    
                    self.include_pairwise_1PN_terms = include_pairwise_1PN_terms
                    self.include_pairwise_25PN_terms = include_pairwise_25PN_terms

                    self.integration_method = integration_method
                    self.KS_use_perturbing_potential = KS_use_perturbing_potential
                    
    def __repr__(self):
        if self.index is None:
            if self.is_binary == False:
                return "Particle(is_binary={0}, mass={1:g})".format(self.is_binary,self.mass)
            else:
                return "Particle(is_binary={0}, mass={1:g}, child1={2:d}, child2={3:d}, a={4:g}, e={5:g}, INCL={6:g}, AP={7:g}, LAN={8:g})".format(self.is_binary,self.mass,self.child1.index,self.child2.index,self.a,self.e,self.INCL,self.AP,self.LAN)
        else:
            if self.is_binary == False:
                return "Particle(is_binary={0}, index={1:d}, mass={2:g})".format(self.is_binary,self.index,self.mass)
            else:
                return "Particle(is_binary={0}, index={1:d}, mass={2:g}, child1={3:d}, child2={4:d}, a={5:g}, e={6:g}, INCL={7:g}, AP={8:g}, LAN={9:g})".format(self.is_binary,self.index,self.mass,self.child1.index,self.child2.index,self.a,self.e,self.INCL,self.AP,self.LAN)

    @property
    def pos(self):
        return self.__pos

    @property
    def vel(self):
        return self.__vel

    @pos.setter
    def pos(self, pos_vec):
        if type(pos_vec).__module__ == np.__name__:
            if pos_vec.size == 3:
                self.x = pos_vec[0]
                self.y = pos_vec[1]
                self.z = pos_vec[2]
                self.__pos = pos_vec
            else:
                raise ValueError('Position vector must be len=3 vector.')
        else:
            raise TypeError('Position vector must be a np vector with len=3.')

    @vel.setter
    def vel(self, vel_vec):
        if type(vel_vec).__module__ == np.__name__:
            if vel_vec.size == 3:
                self.vx = vel_vec[0]
                self.vy = vel_vec[1]
                self.vz = vel_vec[2]
                self.__vel = vel_vec
            else:
                raise ValueError('Velocity vector must be len=3 vector.')
        else:
            raise TypeError('Velocity vector must be a np vector with len=3.')


class Tools(object):
 
    @staticmethod       
    def create_nested_multiple(N,masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None):
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
            particles.append(particle)

        
        #previous_binary = particles[-1]
        for index in range(N_binaries):
            if index==0:
                child1 = particles[0]
                child2 = particles[1]
            else:
                child1 = previous_binary
                child2 = particles[index+1]
            #print 'c',child1,child2
            particle = Particle(is_binary=True,child1=child1,child2=child2,a=semimajor_axes[index],e=eccentricities[index],INCL=inclinations[index],AP=arguments_of_pericentre[index],LAN=longitudes_of_ascending_node[index])

            previous_binary = particle
            particles.append(particle)
            
#            print 'p',particles
        
        return particles


    @staticmethod
    def compute_mutual_inclination(INCL_k,INCL_l,LAN_k,LAN_l):
        cos_INCL_rel = np.cos(INCL_k)*np.cos(INCL_l) + np.sin(INCL_k)*np.sin(INCL_l)*np.cos(LAN_k-LAN_l)
        return np.arccos(cos_INCL_rel)

