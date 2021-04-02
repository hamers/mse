import os
import numpy as np
import ctypes
import ast
import copy

"""
# Multiple Stellar Evolution (MSE) -- A Population Synthesis Code for Multiple-Star Systems #

MSE is code that models the long-term evolution of hierarchical multiple-star systems (binaries, triples, quadruples, and higher-order systems) from the main sequence until remnant phase. It takes into account gravitational dynamical evolution, stellar evolution (using the `sse` tracks), and binary interactions (such as mass transfer and common-envelope evolution).  It includes routines for external perturbations from flybys in the field, or (to limited extent) encounters in dense stellar systems such as galactic nuclei. 

C++ and Fortran compilers are required, as well as Python (2/3) for the Python interface. Make sure to first compile the code using `make`. Please modify the Makefile according to your installation (`CXX` and `FC` should be correctly assigned).  

The script `test_mse.py` can be used to test the installation. The script `run_system.py` is useful for quickly running a system. 

See the user guide (doc/doc.pdf) for more detailed information.

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
        self.__CONST_PARSEC = 206201.0
        self.__CONST_MJUP = 0.000954248

        self.__relative_tolerance = 1.0e-10
        self.__absolute_tolerance_eccentricity_vectors = 1.0e-8
        self.__absolute_tolerance_spin_vectors = 1.0e-3
        self.__absolute_tolerance_angular_momentum_vectors = 1.0e-2

        self.__wall_time_max_s = 1.8e4

        self.__include_quadrupole_order_terms = True
        self.__include_octupole_order_binary_pair_terms = True
        self.__include_octupole_order_binary_triplet_terms = True
        self.__include_hexadecupole_order_binary_pair_terms = True
        self.__include_dotriacontupole_order_binary_pair_terms = True
        self.__include_double_averaging_corrections = False

        self.__MSTAR_gbs_tolerance_default = 1.0e-10
        self.__MSTAR_gbs_tolerance_kick = 1.0e-8
        self.__MSTAR_collision_tolerance = 1.0e-10
        self.__MSTAR_output_time_tolerance = 1.0e-6
        self.__MSTAR_include_PN_acc_10 = True
        self.__MSTAR_include_PN_acc_20 = True
        self.__MSTAR_include_PN_acc_25 = True
        self.__MSTAR_include_PN_acc_30 = True
        self.__MSTAR_include_PN_acc_35 = True
        self.__MSTAR_include_PN_acc_SO = True
        self.__MSTAR_include_PN_acc_SS = True
        self.__MSTAR_include_PN_acc_Q = True
        self.__MSTAR_include_PN_spin_SO = True
        self.__MSTAR_include_PN_spin_SS = True
        self.__MSTAR_include_PN_spin_Q = True
        
        self.__nbody_analysis_fractional_semimajor_axis_change_parameter = 0.01
        self.__nbody_analysis_fractional_integration_time = 0.05
        self.__nbody_analysis_minimum_integration_time = 1.0e1
        self.__nbody_analysis_maximum_integration_time = 1.0e5
        self.__nbody_dynamical_instability_direct_integration_time_multiplier = 1.5
        self.__nbody_semisecular_direct_integration_time_multiplier = 1.0e2
        self.__nbody_supernovae_direct_integration_time_multiplier = 1.5
        self.__nbody_other_direct_integration_time_multiplier = 1.5
        
        self.__effective_radius_multiplication_factor_for_collisions_stars = 3.0
        self.__effective_radius_multiplication_factor_for_collisions_compact_objects = 1.0e3
        
        self.__binary_evolution_CE_energy_flag = 0
        self.__binary_evolution_CE_spin_flag = 1
        self.__binary_evolution_mass_transfer_timestep_parameter = 0.05
        self.__binary_evolution_CE_recombination_fraction = 1.0
        self.__binary_evolution_use_eCAML_model = False
        self.__chandrasekhar_mass = 1.44
        self.__eddington_accretion_factor = 10.0
        self.__nova_accretion_factor = 1.0e-3
        self.__alpha_wind_accretion = 1.5
        self.__beta_wind_accretion = 0.125

        self.__triple_mass_transfer_primary_star_accretion_efficiency_no_disk = 0.1
        self.__triple_mass_transfer_secondary_star_accretion_efficiency_no_disk = 0.1
        self.__triple_mass_transfer_primary_star_accretion_efficiency_disk = 0.9
        self.__triple_mass_transfer_secondary_star_accretion_efficiency_disk = 0.9
        self.__triple_mass_transfer_inner_binary_alpha_times_lambda = 5.0

        self.__particles_committed = False
        self.model_time = 0.0
        self.time_step = 0.0
        self.relative_energy_error = 0.0
        self.state = 0
        self.CVODE_flag = 0
        self.CVODE_error_code = 0
        self.integration_flag = 0
        self.__stop_after_root_found = False
        self.__random_seed = 0
        
        self.__verbose_flag = 0 ### 0: no verbose output in C++; > 0: verbose output, with increasing verbosity (>1 will slow down the code considerably)
        
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

        self.lib.set_mass_transfer_terms.argtypes = (ctypes.c_int,ctypes.c_bool)
        self.lib.set_mass_transfer_terms.restype = ctypes.c_int

        self.lib.get_mass_transfer_terms.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_bool))
        self.lib.get_mass_transfer_terms.restype = ctypes.c_int

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

        self.lib.set_stellar_evolution_properties.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int)
        self.lib.set_stellar_evolution_properties.restype = ctypes.c_int

        self.lib.get_stellar_evolution_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_stellar_evolution_properties.restype = ctypes.c_int


        ### kicks ###
        self.lib.set_kick_properties.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_bool,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_kick_properties.restype = ctypes.c_int

        self.lib.get_kick_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_kick_properties.restype = ctypes.c_int


        ### binary evolution ###
        self.lib.set_binary_evolution_properties.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_binary_evolution_properties.restype = ctypes.c_int

        self.lib.get_binary_evolution_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_binary_evolution_properties.restype = ctypes.c_int


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

        self.lib.set_PN_terms.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool)
        self.lib.set_PN_terms.restype = ctypes.c_int

        self.lib.get_PN_terms.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool))
        self.lib.get_PN_terms.restype = ctypes.c_int

        self.lib.set_tides_terms.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_double,ctypes.c_bool)
        self.lib.set_tides_terms.restype = ctypes.c_int

        self.lib.get_tides_terms.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_bool))
        self.lib.get_tides_terms.restype = ctypes.c_int

        self.lib.set_root_finding_terms.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_bool,ctypes.c_bool,ctypes.c_double,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool);
        self.lib.set_root_finding_terms.restype = ctypes.c_int

        self.lib.set_root_finding_state.argtypes = (ctypes.c_int,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,ctypes.c_bool)
        self.lib.set_root_finding_state.restype = ctypes.c_int

        self.lib.get_root_finding_state.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_bool))
        self.lib.get_root_finding_state.restype = ctypes.c_int

        self.lib.set_constants.argtypes = (ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        self.lib.set_constants.restype = ctypes.c_int

        self.__set_constants_in_code()

        self.lib.set_parameters.argtypes = (ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double, \
            ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,\
            ctypes.c_bool,ctypes.c_bool,ctypes.c_bool,\
            ctypes.c_bool,ctypes.c_int,ctypes.c_bool, ctypes.c_int,ctypes.c_int, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, \
            ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_bool, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
            ctypes.c_double, ctypes.c_double, \
            ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, \
            ctypes.c_bool, \
            ctypes.c_double)
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

        self.lib.set_verbose_flag.argtypes = (ctypes.c_int,)
        self.lib.set_verbose_flag.restype = ctypes.c_int

        self.lib.initialize_code_interface.argtypes = ()
        self.lib.initialize_code_interface.restype = ctypes.c_int
        
        
        ### logging ###
        self.lib.get_size_of_log_data.argtypes = ()
        self.lib.get_size_of_log_data.restype = ctypes.c_int
        
        self.lib.get_log_entry_properties.argtypes = (ctypes.c_int,ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double))
        self.lib.get_log_entry_properties.restype = ctypes.c_int
        
        self.lib.get_internal_index_in_particlesMap_log.argtypes = (ctypes.c_int,ctypes.c_int)
        self.lib.get_internal_index_in_particlesMap_log.restype = ctypes.c_int

        self.lib.get_is_binary_log.argtypes = (ctypes.c_int,ctypes.c_int)
        self.lib.get_is_binary_log.restype = ctypes.c_bool


        self.lib.get_body_properties_from_log_entry.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double), \
                    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_body_properties_from_log_entry.restype = ctypes.c_int
        
        self.lib.get_binary_properties_from_log_entry.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int), \
                        ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), \
                        ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
                        ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.get_binary_properties_from_log_entry.restype = ctypes.c_int


       
        ### tests ###
        self.lib.unit_tests_interface.argtypes = (ctypes.c_int,)
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

        self.lib.sample_spherical_coordinates_unit_vectors_from_isotropic_distribution_interface.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
            ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.sample_spherical_coordinates_unit_vectors_from_isotropic_distribution_interface.restype = ctypes.c_int

        self.lib.test_kick_velocity.argtypes = (ctypes.c_int,ctypes.c_double,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double))
        self.lib.test_kick_velocity.restype = ctypes.c_int

        self.lib.test_flybys_perturber_sampling.argtypes = (ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double, \
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), \
            ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
        self.lib.test_flybys_perturber_sampling.restype = ctypes.c_int

    ###############
    
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
        
        end_time,initial_hamiltonian,state,CVODE_flag,CVODE_error_code,integration_flag = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0)
        error_code = self.lib.evolve_interface(0.0,0.0,ctypes.byref(end_time),ctypes.byref(initial_hamiltonian), \
            ctypes.byref(state),ctypes.byref(CVODE_flag),ctypes.byref(CVODE_error_code),ctypes.byref(integration_flag))

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
        bodies_old = [b.index for b in self.particles if b.is_binary == False]
        children1_old = [o.child1.index for o in orbits]
        children2_old = [o.child2.index for o in orbits]

        ### integrate system of ODEs ###
        start_time = self.model_time

        output_time,hamiltonian,state,CVODE_flag,CVODE_error_code,integration_flag = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(self.integration_flag)
        error_code = self.lib.evolve_interface(start_time,end_time,ctypes.byref(output_time),ctypes.byref(hamiltonian), \
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

        self.error_code = error_code
        self.CVODE_flag = CVODE_flag
        self.CVODE_error_code = CVODE_error_code
        self.state = state
        self.integration_flag = integration_flag

        #if self.integration_flag == 0: 
        self.__copy_particle_structure_from_code()

        self.__update_particles_from_code()

        ### check if the system structure changed (including changes in bodies) ###
        orbits = [p for p in self.particles if p.is_binary == True]
        
        children1 = [o.child1.index for o in orbits]
        children2 = [o.child2.index for o in orbits]
        bodies = [b.index for b in self.particles if b.is_binary == False]            
        
        self.structure_change = False
        if bodies != bodies_old or children1 != children1_old or children2 != children2_old:
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

        if self.enable_tides == False:
            particle.include_tidal_friction_terms = False
            particle.tides_method = 1
            particle.include_tidal_bulges_precession_terms = False
            particle.include_rotation_precession_terms = False

        flag += self.lib.set_tides_terms(particle.index,particle.include_tidal_friction_terms,particle.tides_method,particle.include_tidal_bulges_precession_terms,particle.include_rotation_precession_terms, \
            particle.minimum_eccentricity_for_tidal_precession,particle.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary)
            
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

        flag += self.lib.set_binary_evolution_properties(particle.index,particle.dynamical_mass_transfer_low_mass_donor_timescale,particle.dynamical_mass_transfer_WD_donor_timescale,particle.compact_object_disruption_mass_loss_timescale, \
            particle.common_envelope_alpha, particle.common_envelope_lambda, particle.common_envelope_timescale, particle.triple_common_envelope_alpha)

        flag += self.lib.set_PN_terms(particle.index,particle.include_pairwise_1PN_terms,particle.include_pairwise_25PN_terms,particle.include_spin_orbit_1PN_terms,particle.exclude_1PN_precession_in_case_of_isolated_binary)
        
        if particle.is_external==False:
            
            if particle.is_binary==True:
                flag += self.lib.set_children(particle.index,particle.child1.index,particle.child2.index)
                flag += self.lib.set_orbital_elements(particle.index,particle.a, particle.e, particle.TA, particle.INCL, particle.AP, particle.LAN, particle.sample_orbital_phase_randomly)
                flag += self.lib.set_integration_method(particle.index,particle.integration_method,particle.KS_use_perturbing_potential)
            else:
                flag += self.lib.set_radius(particle.index,particle.radius,particle.radius_dot)
                flag += self.lib.set_mass_transfer_terms(particle.index,particle.include_mass_transfer_terms)
                flag += self.lib.set_spin_vector(particle.index,particle.spin_vec_x,particle.spin_vec_y,particle.spin_vec_z)
                flag += self.lib.set_stellar_evolution_properties(particle.index,particle.stellar_type,particle.object_type,particle.sse_initial_mass,particle.metallicity,particle.sse_time_step,particle.epoch,particle.age, \
                    particle.convective_envelope_mass,particle.convective_envelope_radius,particle.core_mass,particle.core_radius,particle.luminosity,particle.apsidal_motion_constant,particle.gyration_radius,particle.tides_viscous_time_scale,particle.tides_viscous_time_scale_prescription)
                flag += self.lib.set_kick_properties(particle.index,particle.kick_distribution,particle.include_WD_kicks,particle.kick_distribution_sigma_km_s_NS,particle.kick_distribution_sigma_km_s_BH,particle.kick_distribution_sigma_km_s_WD, \
                    particle.kick_distribution_2_m_NS,particle.kick_distribution_4_m_NS,particle.kick_distribution_4_m_ej,particle.kick_distribution_5_v_km_s_NS,particle.kick_distribution_5_v_km_s_BH,particle.kick_distribution_5_sigma)

                if set_instantaneous_perturbation_properties==True:
                    flag += self.lib.set_instantaneous_perturbation_properties(particle.index,particle.instantaneous_perturbation_delta_mass, \
                        particle.instantaneous_perturbation_delta_X,particle.instantaneous_perturbation_delta_Y,particle.instantaneous_perturbation_delta_Z, \
                        particle.instantaneous_perturbation_delta_VX,particle.instantaneous_perturbation_delta_VY,particle.instantaneous_perturbation_delta_VZ)
                        
        else:
            flag += self.lib.set_external_particle_properties(particle.index, particle.external_t_ref, particle.e, particle.external_r_p, particle.INCL, particle.AP, particle.LAN)
    
        return flag

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

        include_tidal_friction_terms,tides_method,include_tidal_bulges_precession_terms,include_rotation_precession_terms,minimum_eccentricity_for_tidal_precession,exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = ctypes.c_bool(True),ctypes.c_int(0),ctypes.c_bool(True),ctypes.c_bool(True),ctypes.c_double(0.0),ctypes.c_bool(True)
        flag += self.lib.get_tides_terms(particle.index,ctypes.byref(include_tidal_friction_terms),ctypes.byref(tides_method),ctypes.byref(include_tidal_bulges_precession_terms),ctypes.byref(include_rotation_precession_terms),ctypes.byref(minimum_eccentricity_for_tidal_precession),ctypes.byref(exclude_rotation_and_bulges_precession_in_case_of_isolated_binary))
        particle.include_tidal_friction_terms = include_tidal_friction_terms.value
        particle.tides_method = tides_method.value
        particle.include_tidal_bulges_precession_terms = include_tidal_bulges_precession_terms.value
        particle.include_rotation_precession_terms = include_rotation_precession_terms.value
        particle.minimum_eccentricity_for_tidal_precession = minimum_eccentricity_for_tidal_precession.value
        particle.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = exclude_rotation_and_bulges_precession_in_case_of_isolated_binary.value

        include_pairwise_1PN_terms,include_pairwise_25PN_terms, include_spin_orbit_1PN_terms, exclude_1PN_precession_in_case_of_isolated_binary = ctypes.c_bool(True),ctypes.c_bool(True),ctypes.c_bool(True),ctypes.c_bool(True)
        flag += self.lib.get_PN_terms(particle.index,ctypes.byref(include_pairwise_1PN_terms),ctypes.byref(include_pairwise_25PN_terms),ctypes.byref(include_spin_orbit_1PN_terms),ctypes.byref(exclude_1PN_precession_in_case_of_isolated_binary))
        particle.include_pairwise_1PN_terms = include_pairwise_1PN_terms.value
        particle.include_pairwise_25PN_terms = include_pairwise_25PN_terms.value
        particle.include_spin_orbit_1PN_terms = include_spin_orbit_1PN_terms.value
        particle.exclude_1PN_precession_in_case_of_isolated_binary = exclude_1PN_precession_in_case_of_isolated_binary.value

        dynamical_mass_transfer_low_mass_donor_timescale,dynamical_mass_transfer_WD_donor_timescale,compact_object_disruption_mass_loss_timescale,common_envelope_alpha,common_envelope_lambda,common_envelope_timescale,triple_common_envelope_alpha = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
        flag += self.lib.get_binary_evolution_properties(particle.index,ctypes.byref(dynamical_mass_transfer_low_mass_donor_timescale),ctypes.byref(dynamical_mass_transfer_WD_donor_timescale),ctypes.byref(compact_object_disruption_mass_loss_timescale), \
            ctypes.byref(common_envelope_alpha), ctypes.byref(common_envelope_lambda), ctypes.byref(common_envelope_timescale), ctypes.byref(triple_common_envelope_alpha))
        particle.dynamical_mass_transfer_low_mass_donor_timescale = dynamical_mass_transfer_low_mass_donor_timescale.value
        particle.dynamical_mass_transfer_WD_donor_timescale = dynamical_mass_transfer_WD_donor_timescale.value
        particle.compact_object_disruption_mass_loss_timescale = compact_object_disruption_mass_loss_timescale.value
        particle.common_envelope_alpha = common_envelope_alpha.value
        particle.common_envelope_lambda = common_envelope_lambda.value
        particle.common_envelope_timescale = common_envelope_timescale.value
        particle.triple_common_envelope_alpha = triple_common_envelope_alpha.value

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
            flag += self.lib.get_relative_position_and_velocity(particle.index,ctypes.byref(x),ctypes.byref(y),ctypes.byref(z),ctypes.byref(vx),ctypes.byref(vy),ctypes.byref(vz))
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

                stellar_type,object_type,sse_initial_mass,metallicity,sse_time_step,epoch,age,convective_envelope_mass,convective_envelope_radius,core_mass,core_radius,luminosity,apsidal_motion_constant,gyration_radius,tides_viscous_time_scale,roche_lobe_radius_pericenter = ctypes.c_int(0),ctypes.c_int(0), \
                    ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                flag += self.lib.get_stellar_evolution_properties(particle.index,ctypes.byref(stellar_type),ctypes.byref(object_type),ctypes.byref(sse_initial_mass),ctypes.byref(metallicity),ctypes.byref(sse_time_step), \
                    ctypes.byref(epoch),ctypes.byref(age),ctypes.byref(convective_envelope_mass),ctypes.byref(convective_envelope_radius),ctypes.byref(core_mass),ctypes.byref(core_radius),ctypes.byref(luminosity),ctypes.byref(apsidal_motion_constant),ctypes.byref(gyration_radius),ctypes.byref(tides_viscous_time_scale),ctypes.byref(roche_lobe_radius_pericenter))
                particle.stellar_type = stellar_type.value
                particle.object_type = object_type.value
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
                particle.apsidal_motion_constant = apsidal_motion_constant.value
                particle.gyration_radius = gyration_radius.value
                particle.tides_viscous_time_scale = tides_viscous_time_scale.value
                particle.roche_lobe_radius_pericenter = roche_lobe_radius_pericenter.value

                kick_distribution,include_WD_kicks,kick_distribution_sigma_km_s_NS,kick_distribution_sigma_km_s_BH,kick_distribution_sigma_km_s_WD,kick_distribution_2_m_NS,kick_distribution_4_m_NS,kick_distribution_4_m_ej,kick_distribution_5_v_km_s_NS,kick_distribution_5_v_km_s_BH,kick_distribution_5_sigma \
                    = ctypes.c_int(0),ctypes.c_bool(False),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                flag += self.lib.get_kick_properties(particle.index,ctypes.byref(kick_distribution),ctypes.byref(include_WD_kicks),ctypes.byref(kick_distribution_sigma_km_s_NS),ctypes.byref(kick_distribution_sigma_km_s_BH),ctypes.byref(kick_distribution_sigma_km_s_WD), \
                ctypes.byref(kick_distribution_2_m_NS),ctypes.byref(kick_distribution_4_m_NS),ctypes.byref(kick_distribution_4_m_ej),ctypes.byref(kick_distribution_5_v_km_s_NS),ctypes.byref(kick_distribution_5_v_km_s_BH),ctypes.byref(kick_distribution_5_sigma))
                particle.kick_distribution = kick_distribution.value
                particle.include_WD_kicks = include_WD_kicks.value
                particle.kick_distribution_sigma_km_s_NS = kick_distribution_sigma_km_s_NS.value
                particle.kick_distribution_sigma_km_s_BH = kick_distribution_sigma_km_s_BH.value
                particle.kick_distribution_sigma_km_s_WD = kick_distribution_sigma_km_s_WD.value
                particle.kick_distribution_2_m_NS = kick_distribution_2_m_NS.value
                particle.kick_distribution_4_m_NS = kick_distribution_4_m_NS.value
                particle.kick_distribution_4_m_ej = kick_distribution_4_m_ej.value
                particle.kick_distribution_5_v_km_s_NS = kick_distribution_5_v_km_s_NS.value
                particle.kick_distribution_5_v_km_s_BH = kick_distribution_5_v_km_s_BH.value
                particle.kick_distribution_5_sigma = kick_distribution_5_sigma.value

                mass_dot = ctypes.c_double(0.0)
                flag = self.lib.get_mass_dot(particle.index,ctypes.byref(mass_dot))
                particle.mass_dot = mass_dot.value

                spin_vec_x,spin_vec_y,spin_vec_z = ctypes.c_double(0.0), ctypes.c_double(0.0), ctypes.c_double(0.0)
                flag += self.lib.get_spin_vector(particle.index,ctypes.byref(spin_vec_x),ctypes.byref(spin_vec_y),ctypes.byref(spin_vec_z))
                particle.spin_vec_x = spin_vec_x.value
                particle.spin_vec_y = spin_vec_y.value
                particle.spin_vec_z = spin_vec_z.value
                
                include_mass_transfer_terms = ctypes.c_bool(True)
                flag += self.lib.get_mass_transfer_terms(particle.index,ctypes.byref(include_mass_transfer_terms))
                particle.include_mass_transfer_terms = include_mass_transfer_terms.value

        return flag
        
    def __update_particles_from_code(self):
        
        flag = 0
        for index,particle in enumerate(self.particles):
            flag += self.__update_particle_from_code(particle)
        return flag

    def __copy_particle_structure_from_code(self):
        self.particles = []
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
                child1_index = child1.value
                child2_index = child2.value

            #print("isb",N_particles,"i",i,"internal_index",internal_index,"mass",mass,"is_binary",is_binary,"child1_index",child1_index,"child2_index",child2_index)#,child1,child2)
            p = Particle(is_binary=is_binary,mass=mass,child1=None,child2=None,a=0.0,e=0.0,INCL=0.0,AP=0.0,LAN=0.0) ### orbital elements should be updated later
            p.index = internal_index

            p.child1_index = child1_index
            p.child2_index = child2_index
            self.particles.append(p)

        binaries = [x for x in self.particles if x.is_binary == True]

        for i,p in enumerate(self.particles):
            if p.is_binary == True:
                i1 = [j for j in range(N_particles) if self.particles[j].index == p.child1_index][0]
                i2 = [j for j in range(N_particles) if self.particles[j].index == p.child2_index][0]
                #print("i1",i1,"i2",i2,"i",i,"p.child1_index",p.child1_index,"p.child2_index",p.child2_index)
                p.child1 = self.particles[i1]
                p.child2 = self.particles[i2]

    def __set_constants_in_code(self):
        self.lib.set_constants(self.__CONST_G,self.__CONST_C,self.__CONST_M_SUN,self.__CONST_R_SUN,self.__CONST_L_SUN,self.__CONST_KM_PER_S,self.__CONST_PER_PC3,self.__CONST_MJUP)


    def __set_parameters_in_code(self):
        self.lib.set_parameters(self.__relative_tolerance,self.__absolute_tolerance_eccentricity_vectors,self.__absolute_tolerance_spin_vectors,self.__absolute_tolerance_angular_momentum_vectors,self.__include_quadrupole_order_terms, \
            self.__include_octupole_order_binary_pair_terms,self.__include_octupole_order_binary_triplet_terms, \
            self.__include_hexadecupole_order_binary_pair_terms,self.__include_dotriacontupole_order_binary_pair_terms, self.__include_double_averaging_corrections, \
            self.__include_flybys, self.__flybys_reference_binary, self.__flybys_correct_for_gravitational_focussing, self.__flybys_velocity_distribution, self.__flybys_mass_distribution, \
            self.__flybys_mass_distribution_lower_value, self.__flybys_mass_distribution_upper_value, self.__flybys_encounter_sphere_radius, \
            self.__flybys_stellar_density, self.__flybys_stellar_relative_velocity_dispersion, \
            self.__binary_evolution_CE_energy_flag, self.__binary_evolution_CE_spin_flag, self.__binary_evolution_mass_transfer_timestep_parameter, self.__binary_evolution_CE_recombination_fraction, self.__binary_evolution_use_eCAML_model, \
            self.__MSTAR_gbs_tolerance_default, self.__MSTAR_gbs_tolerance_kick, self.__MSTAR_collision_tolerance, self.__MSTAR_output_time_tolerance, \
            self.__nbody_analysis_fractional_semimajor_axis_change_parameter,self.__nbody_analysis_fractional_integration_time,self.__nbody_analysis_minimum_integration_time,self.__nbody_analysis_maximum_integration_time, \
            self.__nbody_dynamical_instability_direct_integration_time_multiplier,self.__nbody_semisecular_direct_integration_time_multiplier,self.__nbody_supernovae_direct_integration_time_multiplier,self.__nbody_other_direct_integration_time_multiplier, \
            self.__chandrasekhar_mass,self.__eddington_accretion_factor,self.__nova_accretion_factor,self.__alpha_wind_accretion,self.__beta_wind_accretion, \
            self.__triple_mass_transfer_primary_star_accretion_efficiency_no_disk,self.__triple_mass_transfer_secondary_star_accretion_efficiency_no_disk,self.__triple_mass_transfer_primary_star_accretion_efficiency_disk,self.__triple_mass_transfer_secondary_star_accretion_efficiency_disk,self.__triple_mass_transfer_inner_binary_alpha_times_lambda, \
            self.__effective_radius_multiplication_factor_for_collisions_stars, self.__effective_radius_multiplication_factor_for_collisions_compact_objects, \
            self.__MSTAR_include_PN_acc_10,self.__MSTAR_include_PN_acc_20,self.__MSTAR_include_PN_acc_25,self.__MSTAR_include_PN_acc_30,self.__MSTAR_include_PN_acc_35,self.__MSTAR_include_PN_acc_SO,self.__MSTAR_include_PN_acc_SS,self.__MSTAR_include_PN_acc_Q,self.__MSTAR_include_PN_spin_SO,self.__MSTAR_include_PN_spin_SS,self.__MSTAR_include_PN_spin_Q, \
            self.__stop_after_root_found, \
            self.__wall_time_max_s)

    def reset(self):
        self.__init__()
        self.lib.reset_interface()
        
    def __set_random_seed(self):
        self.lib.set_random_seed(self.random_seed)

    def __set_verbose_flag(self):
        self.lib.set_verbose_flag(self.verbose_flag)

    ### Logging ###
    def __get_log(self):
        N_log = self.lib.get_size_of_log_data()
        #print("get_log",N_log)
        log = []
        for index_log in range(N_log):
            entry = {}
            particles = []
            #N_particles = self.lib.get_number_of_particles()
            
            time = ctypes.c_double(0.0)
            N_particles,event_flag,integration_flag,index1,index2,binary_index,kick_speed_km_s = ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0),ctypes.c_double(0.0)
            flag = self.lib.get_log_entry_properties(index_log,ctypes.byref(time),ctypes.byref(event_flag),ctypes.byref(integration_flag),ctypes.byref(N_particles),ctypes.byref(index1),ctypes.byref(index2),ctypes.byref(binary_index),ctypes.byref(kick_speed_km_s))
            entry.update({'time':time.value,'event_flag':event_flag.value,'integration_flag':integration_flag.value,'index1':index1.value,'index2':index2.value,'binary_index':binary_index.value,'N_particles':N_particles.value,'kick_speed_km_s':kick_speed_km_s.value})
            #print("log i ",index_log,"N",N_log,"integration_flag",integration_flag.value)
            for index_particle in range(N_particles.value):
                internal_index = self.lib.get_internal_index_in_particlesMap_log(index_log,index_particle)

                is_binary = self.lib.get_is_binary_log(index_log,internal_index)
                
                append = False
                if is_binary == False:
                    parent,stellar_type,object_type = ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0)
                    mass,radius,core_mass,sse_initial_mass,convective_envelope_mass,epoch,age,core_radius,convective_envelope_radius,luminosity,ospin,X,Y,Z,VX,VY,VZ,metallicity,spin_vec_x,spin_vec_y,spin_vec_z = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                    flag = self.lib.get_body_properties_from_log_entry(index_log,internal_index,ctypes.byref(parent),ctypes.byref(mass),ctypes.byref(radius),ctypes.byref(stellar_type), \
                        ctypes.byref(core_mass),ctypes.byref(sse_initial_mass),ctypes.byref(convective_envelope_mass), \
                        ctypes.byref(epoch),ctypes.byref(age), \
                        ctypes.byref(core_radius),ctypes.byref(convective_envelope_radius),ctypes.byref(luminosity),ctypes.byref(ospin), \
                        ctypes.byref(X), ctypes.byref(Y), ctypes.byref(Z), ctypes.byref(VX), ctypes.byref(VY), ctypes.byref(VZ), \
                        ctypes.byref(object_type),ctypes.byref(metallicity), \
                        ctypes.byref(spin_vec_x),ctypes.byref(spin_vec_y),ctypes.byref(spin_vec_z))

                    p = Particle(is_binary=is_binary,mass=mass.value,radius=radius.value,stellar_type=stellar_type.value,core_mass=core_mass.value,sse_initial_mass=sse_initial_mass.value, \
                        convective_envelope_mass=convective_envelope_mass.value, epoch=epoch.value, age=age.value, core_radius=core_radius.value, convective_envelope_radius=convective_envelope_radius.value, luminosity=luminosity.value, metallicity=metallicity.value)
                    p.index = internal_index
                    p.parent = parent.value
                    p.ospin = ospin.value
                    p.object_type = object_type.value

                    p.X = X.value
                    p.Y = Y.value
                    p.Z = Z.value
                    p.VX = VX.value
                    p.VY = VY.value
                    p.VZ = VZ.value
                    
                    p.spin_vec_x = spin_vec_x.value
                    p.spin_vec_y = spin_vec_y.value
                    p.spin_vec_z = spin_vec_z.value
                    
                    append = True
                elif integration_flag.value == 0:
                    parent,child1,child2 = ctypes.c_int(0),ctypes.c_int(0),ctypes.c_int(0)
                    mass,a,e,TA,INCL,AP,LAN,h_vec_x,h_vec_y,h_vec_z,e_vec_x,e_vec_y,e_vec_z = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
                    flag = self.lib.get_binary_properties_from_log_entry(index_log,internal_index,ctypes.byref(parent),ctypes.byref(child1),ctypes.byref(child2), \
                        ctypes.byref(mass), ctypes.byref(a),ctypes.byref(e),ctypes.byref(TA),ctypes.byref(INCL), ctypes.byref(AP), ctypes.byref(LAN), \
                        ctypes.byref(h_vec_x),ctypes.byref(h_vec_y),ctypes.byref(h_vec_z), \
                        ctypes.byref(e_vec_x),ctypes.byref(e_vec_y),ctypes.byref(e_vec_z))

                    p = Particle(is_binary=is_binary,mass=mass.value,child1=child1.value,child2=child2.value,a=a.value,e=e.value,TA=TA.value,INCL=INCL.value,AP=AP.value,LAN=LAN.value)
                    p.index = internal_index
                    p.parent = parent.value
                    p.child1_index=child1.value
                    p.child2_index=child2.value
                    p.mass = mass.value
                    
                    p.h_vec_x = h_vec_x.value
                    p.h_vec_y = h_vec_y.value
                    p.h_vec_z = h_vec_z.value

                    p.e_vec_x = e_vec_x.value
                    p.e_vec_y = e_vec_y.value
                    p.e_vec_z = e_vec_z.value
                    
                    append = True

                if append==True:
                    particles.append(p)

            for i,p in enumerate(particles):
                if p.is_binary == True:
                    i1 = [j for j in range(N_particles.value) if particles[j].index == p.child1_index][0]
                    i2 = [j for j in range(N_particles.value) if particles[j].index == p.child2_index][0]
                    #print("i1",i1,"i2",i2,"i",i,"p.child1_index",p.child1_index,"p.child2_index",p.child2_index)
                    p.child1 = particles[i1]
                    p.child2 = particles[i2]

            entry.update({'N_particles':len(particles)})
            entry.update({'particles':particles})
            log.append(entry)
        #print("log done")
        return log

    @property
    def log(self):
        return self.__get_log()


    ### Tests ###
    def unit_tests(self,mode):
        return self.lib.unit_tests_interface(mode)

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

    def test_sample_spherical_coordinates_unit_vectors_from_isotropic_distribution(self):
        r_hat_vec_x,r_hat_vec_y,r_hat_vec_z,theta_hat_vec_x,theta_hat_vec_y,theta_hat_vec_z,phi_hat_vec_x,phi_hat_vec_y,phi_hat_vec_z = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0), \
            ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0), ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
        self.lib.sample_spherical_coordinates_unit_vectors_from_isotropic_distribution_interface(ctypes.byref(r_hat_vec_x),ctypes.byref(r_hat_vec_y),ctypes.byref(r_hat_vec_z), \
            ctypes.byref(theta_hat_vec_x),ctypes.byref(theta_hat_vec_y),ctypes.byref(theta_hat_vec_z), \
            ctypes.byref(phi_hat_vec_x),ctypes.byref(phi_hat_vec_y),ctypes.byref(phi_hat_vec_z))

        r_hat_vec = np.array( [r_hat_vec_x.value, r_hat_vec_y.value, r_hat_vec_z.value] )
        theta_hat_vec = np.array( [theta_hat_vec_x.value, theta_hat_vec_y.value, theta_hat_vec_z.value] )
        phi_hat_vec = np.array( [phi_hat_vec_x.value, phi_hat_vec_y.value, phi_hat_vec_z.value] )
        return r_hat_vec,theta_hat_vec,phi_hat_vec

    def test_kick_velocity(self,kick_distribution,m):
        kw,v = ctypes.c_int(0),ctypes.c_double(0.0)
        self.lib.test_kick_velocity(kick_distribution,m,ctypes.byref(kw),ctypes.byref(v))
        return kw.value,v.value

    def test_flybys_perturber_sampling(self,R_enc,n_star,sigma_rel,M_int):
        M_per,b_vec_x,b_vec_y,b_vec_z,V_vec_x,V_vec_y,V_vec_z = ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0),ctypes.c_double(0.0)
        self.lib.test_flybys_perturber_sampling(R_enc,n_star,sigma_rel,M_int,ctypes.byref(M_per),ctypes.byref(b_vec_x),ctypes.byref(b_vec_y),ctypes.byref(b_vec_z),ctypes.byref(V_vec_x),ctypes.byref(V_vec_y),ctypes.byref(V_vec_z)) 
        b_vec = np.array( [ b_vec_x.value,b_vec_y.value,b_vec_z.value] )
        V_vec = np.array( [ V_vec_x.value,V_vec_y.value,V_vec_z.value] )
        return M_per.value,b_vec,V_vec

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

    @property
    def CONST_PARSEC(self):
        return self.__CONST_PARSEC

    @CONST_PARSEC.setter
    def CONST_PARSEC(self, value):
        self.__CONST_PARSEC = value
        self.__set_constants_in_code()

    @property
    def CONST_MJUP(self):
        return self.__CONST_MJUP

    @CONST_MJUP.setter
    def CONST_MJUP(self, value):
        self.__CONST_MJUP = value
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
    def wall_time_max_s(self):
        return self.__wall_time_max_s
    @wall_time_max_s.setter
    def wall_time_max_s(self, value):
        self.__wall_time_max_s = value
        self.__set_parameters_in_code()

    @property
    def absolute_tolerance_eccentricity_vectors(self):
        return self.__absolute_tolerance_eccentricity_vectors
    @absolute_tolerance_eccentricity_vectors.setter
    def absolute_tolerance_eccentricity_vectors(self, value):
        self.__absolute_tolerance_eccentricity_vectors = value
        self.__set_parameters_in_code()

    @property
    def absolute_tolerance_spin_vectors(self):
        return self.__absolute_tolerance_spin_vectors
    @absolute_tolerance_spin_vectors.setter
    def absolute_tolerance_spin_vectors(self, value):
        self.__absolute_tolerance_spin_vectors = value
        self.__set_parameters_in_code()

    @property
    def absolute_tolerance_angular_momentum_vectors(self):
        return self.__absolute_tolerance_angular_momentum_vectors
    @absolute_tolerance_angular_momentum_vectors.setter
    def absolute_tolerance_angular_momentum_vectors(self, value):
        self.__absolute_tolerance_angular_momentum_vectors = value
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

    @property
    def verbose_flag(self):
        return self.__verbose_flag
    @verbose_flag.setter
    def verbose_flag(self, value):
        self.__verbose_flag = value
        self.__set_verbose_flag()

    @property
    def stop_after_root_found(self):
        return self.__stop_after_root_found
    @stop_after_root_found.setter
    def stop_after_root_found(self, value):
        self.__stop_after_root_found = value
        self.__set_parameters_in_code()


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


    ### N-body ###
    @property
    def MSTAR_gbs_tolerance_default(self):
        return self.__MSTAR_gbs_tolerance_default
    @MSTAR_gbs_tolerance_default.setter
    def MSTAR_gbs_tolerance_default(self, value):
        self.__MSTAR_gbs_tolerance_default = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_gbs_tolerance_kick(self):
        return self.__MSTAR_gbs_tolerance_kick
    @MSTAR_gbs_tolerance_kick.setter
    def MSTAR_gbs_tolerance_kick(self, value):
        self.__MSTAR_gbs_tolerance_kick = value
        self.__set_parameters_in_code()
        
    @property
    def MSTAR_collision_tolerance(self):
        return self.__MSTAR_collision_tolerance
    @MSTAR_collision_tolerance.setter
    def MSTAR_collision_tolerance(self, value):
        self.__MSTAR_collision_tolerance = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_output_time_tolerance(self):
        return self.__MSTAR_output_time_tolerance
    @MSTAR_output_time_tolerance.setter
    def MSTAR_output_time_tolerance(self, value):
        self.__MSTAR_output_time_tolerance = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_10(self):
        return self.__MSTAR_include_PN_acc_10
    @MSTAR_include_PN_acc_10.setter
    def MSTAR_include_PN_acc_10(self, value):
        self.__MSTAR_include_PN_acc_10 = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_20(self):
        return self.__MSTAR_include_PN_acc_20
    @MSTAR_include_PN_acc_20.setter
    def MSTAR_include_PN_acc_20(self, value):
        self.__MSTAR_include_PN_acc_20 = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_25(self):
        return self.__MSTAR_include_PN_acc_25
    @MSTAR_include_PN_acc_25.setter
    def MSTAR_include_PN_acc_25(self, value):
        self.__MSTAR_include_PN_acc_25 = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_30(self):
        return self.__MSTAR_include_PN_acc_30
    @MSTAR_include_PN_acc_30.setter
    def MSTAR_include_PN_acc_30(self, value):
        self.__MSTAR_include_PN_acc_30 = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_35(self):
        return self.__MSTAR_include_PN_acc_35
    @MSTAR_include_PN_acc_35.setter
    def MSTAR_include_PN_acc_35(self, value):
        self.__MSTAR_include_PN_acc_35 = value
        self.__set_parameters_in_code()
        
    @property
    def MSTAR_include_PN_acc_SO(self):
        return self.__MSTAR_include_PN_acc_SO
    @MSTAR_include_PN_acc_SO.setter
    def MSTAR_include_PN_acc_SO(self, value):
        self.__MSTAR_include_PN_acc_SO = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_SS(self):
        return self.__MSTAR_include_PN_acc_SS
    @MSTAR_include_PN_acc_SS.setter
    def MSTAR_include_PN_acc_SS(self, value):
        self.__MSTAR_include_PN_acc_SS = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_acc_Q(self):
        return self.__MSTAR_include_PN_acc_Q
    @MSTAR_include_PN_acc_Q.setter
    def MSTAR_include_PN_acc_Q(self, value):
        self.__MSTAR_include_PN_acc_Q = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_spin_SO(self):
        return self.__MSTAR_include_PN_spin_SO
    @MSTAR_include_PN_spin_SO.setter
    def MSTAR_include_PN_spin_SO(self, value):
        self.__MSTAR_include_PN_spin_SO = value
        self.__set_parameters_in_code()

    @property
    def MSTAR_include_PN_spin_SS(self):
        return self.__MSTAR_include_PN_spin_SS
    @MSTAR_include_PN_spin_SS.setter
    def MSTAR_include_PN_spin_SS(self, value):
        self.__MSTAR_include_PN_spin_SS = value
        self.__set_parameters_in_code()
        
    @property
    def MSTAR_include_PN_spin_Q(self):
        return self.__MSTAR_include_PN_spin_Q
    @MSTAR_include_PN_spin_Q.setter
    def MSTAR_include_PN_spin_Q(self, value):
        self.__MSTAR_include_PN_spin_Q = value
        self.__set_parameters_in_code()
        

    @property
    def nbody_analysis_fractional_semimajor_axis_change_parameter(self):
        return self.__nbody_analysis_fractional_semimajor_axis_change_parameter
    @nbody_analysis_fractional_semimajor_axis_change_parameter.setter
    def nbody_analysis_fractional_semimajor_axis_change_parameter(self, value):
        self.__nbody_analysis_fractional_semimajor_axis_change_parameter = value
        self.__set_parameters_in_code()

    @property
    def nbody_analysis_fractional_integration_time(self):
        return self.__nbody_analysis_fractional_integration_time
    @nbody_analysis_fractional_integration_time.setter
    def nbody_analysis_fractional_integration_time(self, value):
        self.__nbody_analysis_fractional_integration_time = value
        self.__set_parameters_in_code()

    @property
    def nbody_analysis_minimum_integration_time(self):
        return self.__nbody_analysis_minimum_integration_time
    @nbody_analysis_minimum_integration_time.setter
    def nbody_analysis_minimum_integration_time(self, value):
        self.__nbody_analysis_minimum_integration_time = value
        self.__set_parameters_in_code()

    @property
    def nbody_analysis_maximum_integration_time(self):
        return self.__nbody_analysis_maximum_integration_time
    @nbody_analysis_maximum_integration_time.setter
    def nbody_analysis_maximum_integration_time(self, value):
        self.__nbody_analysis_maximum_integration_time = value
        self.__set_parameters_in_code()

    @property
    def nbody_dynamical_instability_direct_integration_time_multiplier(self):
        return self.__nbody_dynamical_instability_direct_integration_time_multiplier
    @nbody_dynamical_instability_direct_integration_time_multiplier.setter
    def nbody_dynamical_instability_direct_integration_time_multiplier(self, value):
        self.__nbody_dynamical_instability_direct_integration_time_multiplier = value
        self.__set_parameters_in_code()

    @property
    def nbody_semisecular_direct_integration_time_multiplier(self):
        return self.__nbody_semisecular_direct_integration_time_multiplier
    @nbody_semisecular_direct_integration_time_multiplier.setter
    def nbody_semisecular_direct_integration_time_multiplier(self, value):
        self.__nbody_semisecular_direct_integration_time_multiplier = value
        self.__set_parameters_in_code()

    @property
    def nbody_supernovae_direct_integration_time_multiplier(self):
        return self.__nbody_supernovae_direct_integration_time_multiplier
    @nbody_supernovae_direct_integration_time_multiplier.setter
    def nbody_supernovae_direct_integration_time_multiplier(self, value):
        self.__nbody_supernovae_direct_integration_time_multiplier = value
        self.__set_parameters_in_code()

    @property
    def nbody_other_direct_integration_time_multiplier(self):
        return self.__nbody_other_direct_integration_time_multiplier
    @nbody_other_direct_integration_time_multiplier.setter
    def nbody_other_direct_integration_time_multiplier(self, value):
        self.__nbody_other_direct_integration_time_multiplier = value
        self.__set_parameters_in_code()

    @property
    def effective_radius_multiplication_factor_for_collisions_stars(self):
        return self.__effective_radius_multiplication_factor_for_collisions_stars
    @effective_radius_multiplication_factor_for_collisions_stars.setter
    def effective_radius_multiplication_factor_for_collisions_stars(self, value):
        self.__effective_radius_multiplication_factor_for_collisions_stars = value
        self.__set_parameters_in_code()

    @property
    def effective_radius_multiplication_factor_for_collisions_compact_objects(self):
        return self.__effective_radius_multiplication_factor_for_collisions_compact_objects
    @effective_radius_multiplication_factor_for_collisions_compact_objects.setter
    def effective_radius_multiplication_factor_for_collisions_compact_objects(self, value):
        self.__effective_radius_multiplication_factor_for_collisions_compact_objects = value
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
        
    @property
    def binary_evolution_CE_recombination_fraction(self):
        return self.__binary_evolution_CE_recombination_fraction
    @binary_evolution_CE_recombination_fraction.setter
    def binary_evolution_CE_recombination_fraction(self, value):
        self.__binary_evolution_CE_recombination_fraction = value
        self.__set_parameters_in_code()

    @property
    def binary_evolution_mass_transfer_timestep_parameter(self):
        return self.__binary_evolution_mass_transfer_timestep_parameter
    @binary_evolution_mass_transfer_timestep_parameter.setter
    def binary_evolution_mass_transfer_timestep_parameter(self, value):
        self.__binary_evolution_mass_transfer_timestep_parameter = value
        self.__set_parameters_in_code()

    @property
    def binary_evolution_use_eCAML_model(self):
        return self.__binary_evolution_use_eCAML_model
    @binary_evolution_use_eCAML_model.setter
    def binary_evolution_use_eCAML_model(self, value):
        self.__binary_evolution_use_eCAML_model = value
        self.__set_parameters_in_code()


    @property
    def chandrasekhar_mass(self):
        return self.__chandrasekhar_mass
    @chandrasekhar_mass.setter
    def chandrasekhar_mass(self, value):
        self.__chandrasekhar_mass = value
        self.__set_parameters_in_code()

    @property
    def eddington_accretion_factor(self):
        return self.__eddington_accretion_factor
    @eddington_accretion_factor.setter
    def eddington_accretion_factor(self, value):
        self.__eddington_accretion_factor = value
        self.__set_parameters_in_code()

    @property
    def nova_accretion_factor(self):
        return self.__nova_accretion_factor
    @nova_accretion_factor.setter
    def nova_accretion_factor(self, value):
        self.__nova_accretion_factor = value
        self.__set_parameters_in_code()

    @property
    def alpha_wind_accretion(self):
        return self.__alpha_wind_accretion
    @alpha_wind_accretion.setter
    def alpha_wind_accretion(self, value):
        self.__alpha_wind_accretion = value
        self.__set_parameters_in_code()

    @property
    def beta_wind_accretion(self):
        return self.__beta_wind_accretion
    @beta_wind_accretion.setter
    def beta_wind_accretion(self, value):
        self.__beta_wind_accretion = value
        self.__set_parameters_in_code()

    
    ### Triple evolution ###
    @property
    def triple_mass_transfer_primary_star_accretion_efficiency_no_disk(self):
        return self.__triple_mass_transfer_primary_star_accretion_efficiency_no_disk
    @triple_mass_transfer_primary_star_accretion_efficiency_no_disk.setter
    def triple_mass_transfer_primary_star_accretion_efficiency_no_disk(self, value):
        self.__triple_mass_transfer_primary_star_accretion_efficiency_no_disk = value
        self.__set_parameters_in_code()

    @property
    def triple_mass_transfer_secondary_star_accretion_efficiency_no_disk(self):
        return self.__triple_mass_transfer_secondary_star_accretion_efficiency_no_disk
    @triple_mass_transfer_secondary_star_accretion_efficiency_no_disk.setter
    def triple_mass_transfer_secondary_star_accretion_efficiency_no_disk(self, value):
        self.__triple_mass_transfer_secondary_star_accretion_efficiency_no_disk = value
        self.__set_parameters_in_code()

    @property
    def triple_mass_transfer_primary_star_accretion_efficiency_disk(self):
        return self.__triple_mass_transfer_primary_star_accretion_efficiency_disk
    @triple_mass_transfer_primary_star_accretion_efficiency_disk.setter
    def triple_mass_transfer_primary_star_accretion_efficiency_disk(self, value):
        self.__triple_mass_transfer_primary_star_accretion_efficiency_disk = value
        self.__set_parameters_in_code()

    @property
    def triple_mass_transfer_secondary_star_accretion_efficiency_disk(self):
        return self.__triple_mass_transfer_secondary_star_accretion_efficiency_disk
    @triple_mass_transfer_secondary_star_accretion_efficiency_disk.setter
    def triple_mass_transfer_secondary_star_accretion_efficiency_disk(self, value):
        self.__triple_mass_transfer_secondary_star_accretion_efficiency_disk = value
        self.__set_parameters_in_code()

    @property
    def triple_mass_transfer_inner_binary_alpha_times_lambda(self):
        return self.__triple_mass_transfer_inner_binary_alpha_times_lambda
    @triple_mass_transfer_inner_binary_alpha_times_lambda.setter
    def triple_mass_transfer_inner_binary_alpha_times_lambda(self, value):
        self.__triple_mass_transfer_inner_binary_alpha_times_lambda = value
        self.__set_parameters_in_code()


class Particle(object):
    def __init__(self, is_binary, mass=None, mass_dot=0.0, radius=1.0e-10, radius_dot=0.0, child1=None, child2=None, a=None, e=None, TA=0.0, INCL=None, AP=None, LAN=None, \
            integration_method = 0, KS_use_perturbing_potential = True, \
            stellar_type=1, object_type=1, sse_initial_mass=None, metallicity=0.02, sse_time_step=1.0, epoch=0.0, age=0.0, core_mass=0.0, core_radius=0.0, \
            include_mass_transfer_terms=True, \
            kick_distribution = 1, include_WD_kicks = False, kick_distribution_sigma_km_s_NS = 265.0, kick_distribution_sigma_km_s_BH=50.0, kick_distribution_sigma_km_s_WD = 1.0, \
            kick_distribution_2_m_NS=1.4, kick_distribution_4_m_NS=1.2, kick_distribution_4_m_ej=9.0, kick_distribution_5_v_km_s_NS=400.0,kick_distribution_5_v_km_s_BH=200.0, kick_distribution_5_sigma=0.3, \
            spin_vec_x=0.0, spin_vec_y=0.0, spin_vec_z=1.0e-10, \
            include_pairwise_1PN_terms=True, include_pairwise_25PN_terms=True, include_spin_orbit_1PN_terms=True, exclude_1PN_precession_in_case_of_isolated_binary=True, \
            include_tidal_friction_terms=True, tides_method=1, include_tidal_bulges_precession_terms=True, include_rotation_precession_terms=True, exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = True, \
            minimum_eccentricity_for_tidal_precession = 1.0e-3, apsidal_motion_constant=0.19, gyration_radius=0.08, tides_viscous_time_scale=1.0e100, tides_viscous_time_scale_prescription=1, \
            convective_envelope_mass=1.0e-10, convective_envelope_radius=1.0e-10, luminosity=1.0e-10, \
            check_for_secular_breakdown=True,check_for_dynamical_instability=True,dynamical_instability_criterion=0,dynamical_instability_central_particle=0,dynamical_instability_K_parameter=0, \
            check_for_physical_collision_or_orbit_crossing=True,check_for_minimum_periapse_distance=False,check_for_minimum_periapse_distance_value=0.0,check_for_RLOF_at_pericentre=True,check_for_RLOF_at_pericentre_use_sepinsky_fit=False, check_for_GW_condition=False, \
            secular_breakdown_has_occurred=False, dynamical_instability_has_occurred=False, physical_collision_or_orbit_crossing_has_occurred=False, minimum_periapse_distance_has_occurred=False, RLOF_at_pericentre_has_occurred = False, GW_condition_has_occurred = False, \
            is_external=False, external_t_ref=0.0, external_r_p=0.0, \
            sample_orbital_phase_randomly=True, instantaneous_perturbation_delta_mass=0.0, instantaneous_perturbation_delta_X=0.0, instantaneous_perturbation_delta_Y=0.0, instantaneous_perturbation_delta_Z=0.0, \
            instantaneous_perturbation_delta_VX=0.0, instantaneous_perturbation_delta_VY=0.0, instantaneous_perturbation_delta_VZ=0.0, \
            VRR_model=0, VRR_include_mass_precession=0, VRR_mass_precession_rate=0.0, VRR_Omega_vec_x=0.0, VRR_Omega_vec_y=0.0, VRR_Omega_vec_z=0.0, \
            VRR_eta_20_init=0.0, VRR_eta_a_22_init=0.0, VRR_eta_b_22_init=0.0, VRR_eta_a_21_init=0.0, VRR_eta_b_21_init=0.0, \
            VRR_eta_20_final=0.0, VRR_eta_a_22_final=0.0, VRR_eta_b_22_final=0.0, VRR_eta_a_21_final=0.0, VRR_eta_b_21_final=0.0, \
            VRR_initial_time = 0.0, VRR_final_time = 1.0,roche_lobe_radius_pericenter=0.0, \
            dynamical_mass_transfer_low_mass_donor_timescale=1.0e3, dynamical_mass_transfer_WD_donor_timescale=1.0e3, compact_object_disruption_mass_loss_timescale=1.0e3, common_envelope_alpha=1.0, common_envelope_lambda=1.0, common_envelope_timescale=1.0e3, triple_common_envelope_alpha=1.0):
                
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
        self.exclude_rotation_and_bulges_precession_in_case_of_isolated_binary = exclude_rotation_and_bulges_precession_in_case_of_isolated_binary

        self.include_pairwise_1PN_terms = include_pairwise_1PN_terms
        self.include_pairwise_25PN_terms = include_pairwise_25PN_terms
        self.include_spin_orbit_1PN_terms = include_spin_orbit_1PN_terms
        self.exclude_1PN_precession_in_case_of_isolated_binary = exclude_1PN_precession_in_case_of_isolated_binary
        
        self.include_mass_transfer_terms = include_mass_transfer_terms

        self.kick_distribution = kick_distribution
        self.include_WD_kicks = include_WD_kicks
        self.kick_distribution_sigma_km_s_NS = kick_distribution_sigma_km_s_NS
        self.kick_distribution_sigma_km_s_BH = kick_distribution_sigma_km_s_BH
        self.kick_distribution_sigma_km_s_WD = kick_distribution_sigma_km_s_WD
        self.kick_distribution_2_m_NS = kick_distribution_2_m_NS
        self.kick_distribution_4_m_NS = kick_distribution_4_m_NS
        self.kick_distribution_4_m_ej = kick_distribution_4_m_ej
        self.kick_distribution_5_v_km_s_NS = kick_distribution_5_v_km_s_NS
        self.kick_distribution_5_v_km_s_BH = kick_distribution_5_v_km_s_BH
        self.kick_distribution_5_sigma = kick_distribution_5_sigma
          
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

        self.dynamical_mass_transfer_low_mass_donor_timescale = dynamical_mass_transfer_low_mass_donor_timescale
        self.dynamical_mass_transfer_WD_donor_timescale = dynamical_mass_transfer_WD_donor_timescale
        self.compact_object_disruption_mass_loss_timescale = compact_object_disruption_mass_loss_timescale
        self.common_envelope_alpha = common_envelope_alpha
        self.common_envelope_lambda = common_envelope_lambda
        self.common_envelope_timescale = common_envelope_timescale
        self.triple_common_envelope_alpha = triple_common_envelope_alpha

        if is_binary==False:
            if mass==None:
                raise RuntimeError('Error when adding particle: body should have mass specified') 
            self.mass = mass
            self.mass_dot = mass_dot
            self.stellar_type = stellar_type
            self.object_type = object_type
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
                    
                    self.integration_method = integration_method
                    self.KS_use_perturbing_potential = KS_use_perturbing_potential
                    
    def __repr__(self):

        if self.index is None:
            if self.is_binary == False:
                return "Particle(is_binary={0}, mass={1:g})".format(self.is_binary,self.mass)
            else:
                #return "Particle(is_binary={0}, child1={1:d}, child2={2:d}, a={3:g}, e={4:g}, INCL={5:g}, AP={6:g}, LAN={7:g})".format(self.is_binary,self.child1,self.child2,self.a,self.e,self.INCL,self.AP,self.LAN)
                return "Particle(is_binary={0})".format(self.is_binary)
        else:
            if self.is_binary == False:
                return "Particle(is_binary={0}, index={1:d}, mass={2:g})".format(self.is_binary,self.index,self.mass)
            else:
                return "Particle(is_binary={0}, index={1:d}, child1={2:d}, child2={3:d}, a={4:g}, e={5:g}, INCL={6:g}, AP={7:g}, LAN={8:g})".format(self.is_binary,self.index,self.child1.index,self.child2.index,self.a,self.e,self.INCL,self.AP,self.LAN)


class Tools(object):

    @staticmethod
    def check_for_default_values(N_bodies,metallicities,stellar_types,object_types,inclinations,longitudes_of_ascending_node,arguments_of_pericentre):
        if stellar_types == []:
            for i in range(N_bodies):
                stellar_types.append(1)
            print("mse.py -- stellar_types not explicitly given -- setting initial stellar types to",stellar_types)

        if object_types == []:
            for i in range(N_bodies):
                object_types.append(1)
            print("mse.py -- object_types not explicitly given -- setting initial object types to",object_types)

        if metallicities == []:
            for i in range(N_bodies):
                metallicities.append(0.02)
            print("mse.py -- metallicities not explicitly given -- setting initial metallicities to",metallicities)
     
        if inclinations == []:
            for i in range(N_bodies-1):
                inclinations.append(np.arccos(np.random.random()))
            print("mse.py -- inclinations not explicitly given -- setting initial inclinations to",inclinations)
            
        if longitudes_of_ascending_node == []:
            for i in range(N_bodies-1):
                longitudes_of_ascending_node.append(2.0*np.pi * np.random.random())
            print("mse.py -- longitudes_of_ascending_node not explicitly given -- setting initial longitudes_of_ascending_node to",longitudes_of_ascending_node)
            
        if arguments_of_pericentre == []:
            for i in range(N_bodies-1):
                arguments_of_pericentre.append(2.0*np.pi * np.random.random())
            print("mse.py -- arguments_of_pericentre not explicitly given -- setting initial arguments_of_pericentre to",arguments_of_pericentre)
     
    @staticmethod       
    def create_fully_nested_multiple(N,masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None,metallicities=[],stellar_types=[],object_types=[]):

        """
        N is number of bodies
        masses should be N-sized array
        the other arguments should be (N-1)-sized arrays
        """

        N_bodies = N
        N_binaries = N-1

        Tools.check_for_default_values(N_bodies,metallicities,stellar_types,object_types,inclinations,longitudes_of_ascending_node,arguments_of_pericentre)

        particles = []

        for index in range(N_bodies):
            particle = Particle(is_binary=False,mass=masses[index])
            if radii is not None:
                particle.radius = radii[index]

            particle.metallicity = metallicities[index]
            particle.stellar_type = stellar_types[index]
            particle.object_type = object_types[index]
                
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
    def create_2p2_quadruple_system(masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None,metallicities=[],stellar_types=[],object_types=[]):
        
        """
        Create a 2+2 quadruple system.
        Masses should contain the four masses.
        The other arguments should be length 3 arrays; first two entries: the two inner binaries; third entry: outer binary.
        """

        N_bodies = 4
        N_binaries = N_bodies-1

        Tools.check_for_default_values(N_bodies,metallicities,stellar_types,object_types,inclinations,longitudes_of_ascending_node,arguments_of_pericentre)

        particles = []

        ### Add the bodies ###
        for index in range(N_bodies):
            particle = Particle(is_binary=False,mass=masses[index])
            if radii is not None:
                particle.radius = radii[index]

            particle.metallicity = metallicities[index]
            particle.stellar_type = stellar_types[index]
            particle.object_type = object_types[index]
            
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
    def parse_config(N_bodies,configuration):
        # convert input string '{num+[others]}' to the basic '[1,...[1,[others]]...]' list string 
        # for first encounter of '+' (and corresponding '{' and '}') only => call multiple times
        def first_plus_to_nested(old_str):
            err = 0
            # first occurances of '+' and '{'
            first_plus = old_str.find('+')
            first_left = old_str.find('{')
            # check ordering
            if first_left < first_plus:                
                chars = old_str[first_left+1:first_plus]
                # check if characters between '{' and '+' are int 
                try:
                    num = int(chars)                
                except:
                    err = 1
                    return '', err
                # variable to identify the '}' corresponding to '{'
                proper_brkt = -1
                for ind, char in enumerate(configuration):
                    if ind == first_left:
                        proper_brkt = 0
                        continue
                    # case of further nested '{' and '}'
                    if char == '{':
                        proper_brkt += 1
                    elif char == '}':
                        # error if '}' is found before '{'
                        if proper_brkt < 0:
                            err = 2
                            return '', err
                        # first corresponding '}' found
                        elif proper_brkt == 0:
                            first_right = ind
                            break
                        # case of further nested '{' and '}'
                        elif proper_brkt > 0:
                            proper_brkt -= 1
                # replacing '{num+' and '}' num times and writing in new string
                replace_left = ''
                replace_right = ''            
                for i in range(num):            
                    replace_left += '[1,'
                    replace_right += ']'
                new_str = old_str[:first_left] + replace_left + old_str[first_plus+1:first_right] + replace_right + old_str[first_right+1:] 
                return new_str, err 
            # error if placement of '{' and '+' is wrong                
            else:
                err = 3
                return '', err        

        # check if parsed list contains either list or int elements
        def list_is_int(lst):
            if type(lst) == int:            
                is_int = True
            # children of list need to be tested recursively
            else:
                is_int = True
                for elem in lst:
                    if type(elem) == list:  
                        is_int = list_is_int(elem)
                    elif type(elem) == int:
                        is_int = True
                    else:
                        is_int = False
                        break   
            return is_int

        # check if parsed list is a binary tree - each level has either two children (node) or one (body)
        def list_is_binary(lst):
            if type(lst) == int:            
                is_binary = True
            # children of list need to be tested recursively
            elif len(lst) == 2:
                is_binary = True
                for elem in lst:
                    if type(elem) == list:
                        if len(elem) == 2:
                            is_binary = list_is_binary(elem)
                        else:
                            is_binary = False
                            break
            else:
                is_binary = False 
            return is_binary

        # convert parsed list [num,others] to basic [[1,...[1,1]],others] containg only ones
        def convert_to_ones(lst):
            # fullly nested hierachy if input is int
            if type(lst) == int:            
                for i in range(lst):
                    if i == 0:
                        continue
                    elif i == 1:
                        sub_lst = [1,1]
                    else:
                        sub_lst = [1]+[sub_lst]
                lst = sub_lst
            else:
                for ind, elem in enumerate(lst):
                    if type(elem) == list:  
                        elem = convert_to_ones(elem)
                    elif elem == 1:
                        continue
                    elif type(elem) == int:
                        for i in range(elem):
                            if i == 0:
                                continue
                            # smallest child element with 2 bodies
                            elif i == 1:
                                sub_lst = [1,1]
                            # adding extra body to previous hierarchy
                            else:
                                sub_lst = [1]+[sub_lst]
                        lst[ind] = sub_lst
            return lst

        # find total number of bodies in a arbitrary nested list containing only ones
        def nested_list_size(lst):
            if type(lst) == int:            
                count = lst
            # iterate through children lists if they exist    
            else:
                count = 0            
                for elem in lst:
                    if type(elem) == list:  
                        count += nested_list_size(elem)
                    else:
                        count += 1    
            return count  

        # remove all spaces in string
        configuration = configuration.replace(' ','')
        # number of '+' should equal number of '{' and '}'
        num_plus = configuration.count('+')
        num_curl_left = configuration.count('{')
        num_curl_right = configuration.count('}')
        if num_curl_left == num_plus and num_curl_right == num_plus:
            while num_plus != 0:
                # convert sting notation with '+','{','}' to list string, and get error if any                    
                configuration, error = first_plus_to_nested(configuration)
                if configuration == '' and error == 1:
                    print("Value between '{' and '+' should be int. Input a valid configuration.")
                    print("Exiting...")
                    exit()
                elif configuration == '' and error == 2:
                    print("'{' should occur before '}'. Input a valid configuration.")
                    print("Exiting...")
                    exit()
                elif configuration == '' and error == 3:
                    print("'{' should occur before '+'. Input a valid configuration.")
                    print("Exiting...")
                    exit()
                num_plus -= 1
        else:
            print("Wrong number of curly brackets. Input a valid configuration.")
            print("Exiting...")
            exit()
        
        # check if string can be parsed to a meaningful list
        try:
            # function which parses string to Python expression
            lst = ast.literal_eval(configuration)
        except:
            print("String could not be parsed to meaningful list. Input a valid configuration.")
            print("Exiting...")
            exit()
        # required conditions to convert parsed list to contain only ones
        if list_is_binary(lst):
            if list_is_int(lst):
                lst = convert_to_ones(lst)   
                # check if configuration agrees with given number of bodies   
                if nested_list_size(lst) == N_bodies:
                    return lst
                else:
                    print("Number of elements incorrect. Input a valid configuration.")
                    print("Exiting...")
                    exit()
            else:
                print("Data types of elements should be int. Input a valid configuration.")
                print("Exiting...")
                exit()
        else:
            print("Configuration should be int or binary list. Input a valid configuration.")
            print("Exiting...")
            exit()

    @staticmethod
    def create_hierarchy(N_bodies,configuration,masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None,metallicities=[],stellar_types=[],object_types=[]):

        Tools.check_for_default_values(N_bodies,metallicities,stellar_types,object_types,inclinations,longitudes_of_ascending_node,arguments_of_pericentre)

        print("="*50)
        print("mse.py -- parsing given configuration into particles")

        config_list = Tools.parse_config(N_bodies, configuration)
        print("Verbose configuration :",config_list)
        print()

        N_binaries = N_bodies-1

        particles = []

        for index in range(N_bodies):
            particle = Particle(is_binary=False,mass=masses[index])
            if radii is not None:
                particle.radius = radii[index]

            particle.metallicity = metallicities[index]
            particle.stellar_type = stellar_types[index]
            particle.object_type = object_types[index]

            print("Particle id :",index)
            particles.append(particle)
        print()

        def nested_iteration(lst): 
            N_subbinary = 0   
            lst_tpl = ()
            N_subbinary_tpl = ()
            # iterate through binary list backward (rightmost first)
            for ind,elem in reversed(list(enumerate(lst))):        
                if type(elem) == list:     
                    ret_N, ret_lst_tpl, ret_N_tpl = nested_iteration(elem)
                    N_subbinary += ret_N 
                    lst_tpl += ret_lst_tpl
                    N_subbinary_tpl += ret_N_tpl
            N_subbinary += 1
            lst_tpl += (lst,)
            N_subbinary_tpl += (N_subbinary,) 
            return N_subbinary, lst_tpl, N_subbinary_tpl

        # lst_tpl : tuple of each of the N_binaries, starting from the rightmost in tree diagram [not used in the program; visual confirmation only]
        # N_subbinaries_tpl : tuple of number of subbinaries in each of the N_binaries (corresponding to lst_tpl) 
        _, lst_tpl, N_subbinary_tpl = nested_iteration(config_list)
        print("Tuple of each of the",N_binaries,"binaries (starting from rightmost, least hierarchy) :",lst_tpl)
        print()
        print("Tuple of number of subbinaries in each of the",N_binaries,"binaries (starting from rightmost, least hierarchy) :",N_subbinary_tpl)
        print()

        for index in range(N_binaries):
            if index == 0:
                child1 = particles[0]
                child2 = particles[1]
                # particle id only for printing; redundant otherwise
                id1 = 0
                id2 = 1

                particle_id = 2
                binary_id = N_bodies

            elif N_subbinary_tpl[index] == 1 and N_subbinary_tpl[index-1] >= 1:
                child1 = particles[particle_id]
                child2 = particles[particle_id+1]
                # particle id only for printing; redundant otherwise
                id1 = particle_id
                id2 = particle_id+1

                old_binary_id = binary_id
                particle_id += 2
                binary_id += 1

            elif N_subbinary_tpl[index] == N_subbinary_tpl[index-1]+1:
                child1 = particles[binary_id]
                child2 = particles[particle_id]
                # particle id only for printing; redundant otherwise
                id1 = binary_id
                id2 = particle_id

                particle_id += 1
                binary_id += 1

            elif N_subbinary_tpl[index] > N_subbinary_tpl[index-1]+1:
                child1 = particles[old_binary_id]
                child2 = particles[binary_id]
                # particle id only for printing; redundant otherwise
                id1 = old_binary_id
                id2 = binary_id

                binary_id += 1

            print("Particle id :",binary_id,"\tChild 1 id :",id1,"\tChild 2 id :",id2)
            particle = Particle(is_binary=True,child1=child1,child2=child2,a=semimajor_axes[index],e=eccentricities[index],INCL=inclinations[index],AP=arguments_of_pericentre[index],LAN=longitudes_of_ascending_node[index])
            particles.append(particle)

        return particles         


    @staticmethod
    def compute_mutual_inclination(INCL_k,INCL_l,LAN_k,LAN_l):
        cos_INCL_rel = np.cos(INCL_k)*np.cos(INCL_l) + np.sin(INCL_k)*np.sin(INCL_l)*np.cos(LAN_k-LAN_l)
        return np.arccos(cos_INCL_rel)

    @staticmethod
    def compute_effective_temperature(luminosity, radius, CONST_L_SUN, CONST_R_SUN):
        """
        Assumes black body radiation.
        Luminosity and radius should be in standard code units.
        Returns T_eff in K
        """
        
        T_Sun = 5770.0 ### (the Sun knows about ATI's old product stack)
        T_eff = T_Sun * pow(luminosity/CONST_L_SUN,0.25) * pow(radius/CONST_R_SUN,-0.5)
        return T_eff
    
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
    def evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,tend,N_steps,stellar_types=[],make_plots=True,fancy_plots=False,plot_filename="test1",show_plots=True,object_types=[],random_seed=0,verbose_flag=0,include_WD_kicks=False,kick_distribution_sigma_km_s_WD=1.0):

        np.random.seed(random_seed)
        
        if configuration == "fully_nested":
            particles = Tools.create_fully_nested_multiple(N_bodies, masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,metallicities=metallicities,stellar_types=stellar_types,object_types=object_types)
        elif configuration == "2+2_quadruple":
            particles = Tools.create_2p2_quadruple_system(masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,metallicities=metallicities,stellar_types=stellar_types,object_types=object_types)
        else:
            particles = Tools.create_hierarchy(N_bodies,configuration,masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,radii=None,metallicities=metallicities,stellar_types=stellar_types,object_types=object_types)

        print("="*50)
        print("mse.py -- evolve_system() -- running system with parameters:")
        print("Configuration: ",configuration)
        print("N_bodies: ",N_bodies)
        print("Object types:",object_types)
        print("Stellar types:",stellar_types)
        print("Masses/MSun: ",masses)
        print("Metallicities: ",metallicities)
        print("Semimajor axes (au): ",semimajor_axes)
        print("Eccentricities: ",eccentricities)
        print("Inclinations (rad): ",inclinations)
        print("Longitudes of the ascending node (rad): ",longitudes_of_ascending_node)
        print("Arguments of periapsis (rad): ",inclinations)
        print("Integration time (yr): ",tend)
        print("Number of plot output steps: ",N_steps)

        print("="*50)
        print("Starting evolution")
        
        from mse import MSE

        orbits = [x for x in particles if x.is_binary==True]
        bodies = [x for x in particles if x.is_binary==False]

        for b in bodies:
            b.include_WD_kicks = include_WD_kicks
            b.kick_distribution_sigma_km_s_WD = kick_distribution_sigma_km_s_WD
            b.common_envelope_timescale = 1.0e3
            
        N_bodies = len(bodies)
        N_orbits = len(orbits)

        code = MSE()
        code.enable_root_finding = True
        code.add_particles(particles)
        code.orbital_phases_random_seed = 1

        code.random_seed = random_seed
        code.verbose_flag = verbose_flag

        code.include_flybys = True
        code.flybys_stellar_density = 0.1*code.CONST_PER_PC3
        code.flybys_encounter_sphere_radius = 1.0e5
        code.binary_evolution_use_eCAML_model = False

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
        
        N_orbits_status = [N_orbits]
        N_bodies_status = [N_bodies]
        t = 0.0
        integration_flags = [[]]

        i_status = 0

        dt = tend/float(N_steps)
        i = 0
        
        while t<tend:

            t+=dt
            code.evolve_model(t)
            
            error_code = code.error_code
            if error_code != 0:
                print("mse.py -- Internal error with code ",error_code,"occurred -- stopping simulation.")
                return error_code,code.log
                
            CVODE_flag = code.CVODE_flag
            state = code.state

            particles = code.particles
            orbits = [x for x in particles if x.is_binary==True]
            bodies = [x for x in particles if x.is_binary==False]
            N_orbits = len(orbits)
            N_bodies = len(bodies)
               
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
                
                N_orbits_status.append(N_orbits)
                N_bodies_status.append(N_bodies)
                
                i_status += 1
                
            print( 't/Myr',t*1e-6,'masses/MSun',[b.mass for b in bodies],'es',[o.e for o in orbits],'smas/au',[o.a for o in orbits],'integration_flag',code.integration_flag)
            
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

            t_print[i_status].append(t)        
            integration_flags[i_status].append(code.integration_flag)

            i += 1

        N_status = i_status+1
        
        for i_status in range(N_status):
            t_print[i_status] = np.array(t_print[i_status])

        print("Final properties -- ","masses/MSun",[m_print[-1][i][-1] for i in range(N_bodies)])

        if verbose_flag > 0:
            print("log",code.log)
            
        if make_plots==True:
            
            try:
                from matplotlib import pyplot
                from matplotlib import lines
            except ImportError:
                print("evolve_system.py -- ERROR: cannot import Matplotlib")
                exit(-1)
            
            if fancy_plots == True:
                print("Using LaTeX for plot text")
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)  
        
            print("Number of log entries:",len(code.log))
            
            plot_log = []
            previous_event_flag = -1
            for index_log,log in enumerate(code.log):
                event_flag = log["event_flag"]
                if previous_event_flag == event_flag and (event_flag == 4 or event_flag == 10):
                    continue
                plot_log.append(log)
                previous_event_flag = event_flag
                            
            N_l = len(plot_log)
            N_r = int(np.ceil(np.sqrt(N_l)))#+1
            N_c = N_r
            panel_length = 3
            #fontsize=N_l
            fontsize=10
            fig=pyplot.figure(figsize=(N_r*panel_length,N_r*panel_length))
            
            legend_elements = []
            for k in range(15):
                color,s,description = Tools.get_color_and_size_and_description_for_star(k,1.0)
                legend_elements.append( lines.Line2D([0],[0], marker = 'o', markerfacecolor = color, color = 'w', markersize = 10 ,label="$\mathrm{%s}$"%description))#label = "$\mathrm{%s}$"%description) )

            for index_log,log in enumerate(plot_log):
                plot=fig.add_subplot(N_r,N_c,index_log+1)
                particles = log["particles"]
                event_flag = log["event_flag"]
                index1 = log["index1"]
                index2 = log["index2"]

                Tools.generate_mobile_diagram(particles,plot,fontsize=fontsize,index1=index1,index2=index2,event_flag=event_flag)

                text = Tools.get_description_for_event_flag(event_flag)
                plot.set_title(text,fontsize=fontsize)
                plot.annotate("$t\simeq %s\,\mathrm{Myr}$"%round(log["time"]*1e-6,2),xy=(0.1,0.9),xycoords='axes fraction',fontsize=fontsize)
                #if index_log>0: break

                if index_log == 0:
                    plot.legend(handles = legend_elements, bbox_to_anchor = (-0.05, 1.50), loc = 'upper left', ncol = 5,fontsize=0.85*fontsize)
                
            fig.savefig(plot_filename + "_mobile.pdf")
        
            fig=pyplot.figure(figsize=(8,10))
            Np=3
            plot1=fig.add_subplot(Np,1,1)
            plot2=fig.add_subplot(Np,1,2,yscale="log")
            plot3=fig.add_subplot(Np,1,3,yscale="linear")
            
            fig_pos=pyplot.figure(figsize=(8,8))
            plot_pos=fig_pos.add_subplot(1,1,1)

            fig_HRD=pyplot.figure(figsize=(8,8))
            plot_HRD=fig_HRD.add_subplot(1,1,1)
            
            colors = ['k','tab:red','tab:green','tab:blue','y','k','tab:red','tab:green','tab:blue','y']
            linewidth=1.0
            for i_status in range(N_status):
                N_bodies = N_bodies_status[i_status]
                N_orbits = N_orbits_status[i_status]
                
                for index in range(N_bodies):
                    color=colors[index]
                    plot1.plot(1.0e-6*t_print[i_status],m_print[i_status][index],color=color,linewidth=linewidth)
                    plot1.plot(1.0e-6*t_print[i_status],mc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    plot3.plot(1.0e-6*t_print[i_status],k_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],R_print[i_status][index],color=color,linestyle='solid',linewidth=linewidth)
                    plot2.plot(1.0e-6*t_print[i_status],Rc_print[i_status][index],color=color,linestyle='dotted',linewidth=linewidth)
                    
                    parsec_in_AU = code.CONST_PARSEC
                    plot_pos.plot(np.array(X_print[i_status][index])/parsec_in_AU,np.array(Y_print[i_status][index])/parsec_in_AU,color=color,linestyle='solid',linewidth=linewidth)
                    
                    plot_HRD.plot(np.log10(np.array(T_eff_print[i_status][index])), np.log10(np.array(L_print[i_status][index])/code.CONST_L_SUN),color=color,linestyle='solid',linewidth=linewidth)
                   
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

            plots = [plot1,plot2,plot3]
            for plot in plots:
                plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

            log_CEs = [x for x in code.log if x["event_flag"] == 6]
            t_CEs_Myr = np.array([x["time"]*1e-6 for x in log_CEs])
            
            for k,t in enumerate(t_CEs_Myr):
                plot2.axvline(x=t,linestyle='dotted',color='tab:red',linewidth=0.5)
                plot2.annotate("$\mathrm{CE}$",xy=(1.02*t,1.0e3),fontsize=0.8*fontsize)
            
            plot1.set_ylabel("$m/\mathrm{M}_\odot$",fontsize=fontsize)
            plot2.set_ylabel("$r/\mathrm{au}$",fontsize=fontsize)
            plot3.set_ylabel("$\mathrm{Stellar\,Type}$",fontsize=fontsize)
            plot3.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            plot2.set_ylim(1.0e-5,1.0e5)
            
            plot_pos.set_xlabel("$X/\mathrm{pc}$",fontsize=fontsize)
            plot_pos.set_ylabel("$Y/\mathrm{pc}$",fontsize=fontsize)
            plot_pos.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            
            
            plot_HRD.set_xlim(5.0,3.0)
            plot_HRD.set_ylim(-4.0,6.0)
            plot_HRD.set_xlabel("$\mathrm{log}_{10}(T_\mathrm{eff}/\mathrm{K})$",fontsize=fontsize)
            plot_HRD.set_ylabel(r"$\mathrm{log}_{10}(L/L_\odot)$",fontsize=fontsize)
            plot_HRD.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)
            
            fig.savefig(plot_filename + ".pdf")
            fig_pos.savefig(plot_filename + "_pos.pdf")
            fig_HRD.savefig(plot_filename + "_HRD.pdf")
            
            if show_plots == True:
                pyplot.show()
    
        error_code_copy = copy.deepcopy(code.error_code)
        log_copy = copy.deepcopy(code.log)

        code.reset()
      
        return error_code_copy, log_copy

       
    @staticmethod
    def generate_mobile_diagram(particles,plot,line_width_horizontal=1.5,line_width_vertical = 0.2,line_color = 'k',line_width = 1.5,fontsize=12,use_default_colors=True,index1=-1,index2=-1,event_flag=-1):
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

        if len(binaries)==0:
            if len(bodies)==0:
                print("mse.py -- generate_mobile_diagram -- zero bodies and zero binaries!")
                return
            else:
                Tools.draw_bodies(plot,bodies,fontsize,index1=index1,index2=index2,event_flag=event_flag)
                return

        Tools.determine_binary_levels_in_particles(particles)                    
        unbound_bodies = [x for x in particles if x.is_binary==False and x.parent == None]
        if len(unbound_bodies)>0:
            Tools.draw_bodies(plot,unbound_bodies,fontsize,y_ref = 1.5*line_width_vertical,dx=0.4*line_width_horizontal,dy=0.4*line_width_vertical,index1=index1,index2=index2,event_flag=event_flag)
            
        top_level_binaries = [x for x in binaries if x.level==0]
        
        for p in particles:
            p.color = 'k'
        
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
        top_x = 0.0
        top_y = 0.0
        
        if len(top_level_binaries)>1: ### adjustments for two unbound subsystems
            top_x = -3*line_width_horizontal
            top_y = 0.0
            
        for index,top_level_binary in enumerate(top_level_binaries):
        
            top_level_binary.x = top_x - 5*index * line_width_horizontal
            top_level_binary.y = top_y 

            x_min = x_max = y_min = 0.0
            y_max = line_width_vertical
    
            plot.plot( [top_level_binary.x,top_level_binary.x], [top_level_binary.y,top_level_binary.y + line_width_vertical ], color=line_color,linewidth=line_width)
            x_min,x_max,y_min,y_max = Tools.draw_binary_node(plot,top_level_binary,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max,index1=index1,index2=index2,event_flag=event_flag)

        plot.set_xticks([])
        plot.set_yticks([])
        #print("minmax",x_min,x_max,y_min,y_max)
        beta = 0.7
        plot.set_xlim([x_min - beta*np.fabs(x_min),x_max + beta*np.fabs(x_max)])
        plot.set_ylim([y_min - beta*np.fabs(y_min),1.5*y_max + beta*np.fabs(y_max)])
        
        #plot.autoscale(enable=True,axis='both')
        
    @staticmethod
    def draw_binary_node(plot,particle,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max,index1=-1,index2=-1,event_flag=-1):
        x = particle.x
        y = particle.y
        
        child1 = particle.child1
        child2 = particle.child2

        N_ra = 1
        if (particle.a<=0.1): N_ra = 2
        text = "$a=\mathrm{%s\,au}$\n $e= %.2f$"%(round(particle.a,N_ra),particle.e)
        x_plot = x - 0*0.8*line_width_horizontal
        y_plot = y - 0.3*line_width_vertical

        bbox_props = dict(boxstyle="round", pad=0.1, fc=particle.color, ec="k", lw=0.5,alpha=0.9)
        plot.text(x_plot, y_plot, text, ha="center", va="center", rotation=0,
            size=0.7*fontsize,
            bbox=bbox_props)
        
   
        #text = "$e = %.2f$"%(particle.e)
        #plot.annotate(text,xy=(x - 0.8*line_width_horizontal,y - 0.6*line_width_vertical),fontsize=fontsize,color=particle.color)

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

        CONST_MJUP = 0.000954248
        
        
       
        ### handle children ###
        if child1.is_binary == True:
            x_min,x_max,y_min,y_max = Tools.draw_binary_node(plot,child1,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max,index1=index1,index2=index2,event_flag=event_flag)
        else:
            color,s,description = Tools.get_color_and_size_and_description_for_star(child1.stellar_type,child1.radius)
            plot.scatter([child1.x],[child1.y],color=color,s=s,zorder=10)
            if child1.object_type == 2:
                text = "$\mathrm{%s}\,\mathrm{M_J}$"%(str(round(child1.mass/CONST_MJUP,1)))    
            else:
                text = "$\mathrm{%s}$"%(str(round(child1.mass,1)))
                
            plot.annotate(text,xy=(child1.x - 0.6*line_width_horizontal,child1.y - 0.5*line_width_vertical),color='k',fontsize=fontsize,zorder=10)

            if event_flag in [4,6,8]:
                if child1.index in [index1,index2] or child2.index in [index1,index2]:
                    plot.plot([child1.x,child2.x],[child1.y,child2.y],color='r',zorder=8)
                    if child1.index == index1:
                        plot.arrow(child1.x, child1.y, 0.5*(child2.x-child1.x), 0, head_width=0.05, head_length=0.1*np.fabs(child2.x-child1.x),zorder=9, color='r')
                    else:
                        plot.arrow(child2.x, child2.y, -0.5*np.fabs(child2.x-child1.x), 0, head_width=0.05, head_length=0.1*np.fabs(child2.x-child1.x),zorder=9, color='r')
            if event_flag in [2,12]:
                if child1.index == index1:
                    plot.scatter([child1.x],[child1.y],color=color,s=3*s,zorder=9,marker='*')
                if child2.index == index1:
                    plot.scatter([child2.x],[child2.y],color=color,s=3*s,zorder=9,marker='*')
            
            
        if child2.is_binary == True:
            x_min,x_max,y_min,y_max = Tools.draw_binary_node(plot,child2,line_width_horizontal,line_width_vertical,line_color,line_width,fontsize,x_min,x_max,y_min,y_max,index1=index1,index2=index2,event_flag=event_flag)
        else:
            color,s,description = Tools.get_color_and_size_and_description_for_star(child2.stellar_type,child2.radius)
            plot.scatter([child2.x],[child2.y],color=color,s=s,zorder=10)

            if child2.object_type == 2:
                text = "$\mathrm{%s}\,\mathrm{M_J}$"%(str(round(child2.mass/CONST_MJUP,1)))    
            else:
                text = "$\mathrm{%s}$"%(str(round(child2.mass,1)))

            plot.annotate(text,xy=(child2.x - 0.3*line_width_horizontal,child2.y - 0.5*line_width_vertical),color='k',fontsize=fontsize,zorder=10)

            if event_flag in [4,6,8]:
                if child1.index in [index1,index2] or child2.index in [index1,index2]:
                    plot.plot([child1.x,child2.x],[child1.y,child2.y],color='r',zorder=8)
                    if child1.index == index1:
                        plot.arrow(child1.x, child1.y, 0.5*(child2.x-child1.x), 0, head_width=0.05, head_length=0.1*np.fabs(child2.x-child1.x),zorder=9, color='r')
                    else:
                        plot.arrow(child2.x, child2.y, -0.5*np.fabs(child2.x-child1.x), 0, head_width=0.05, head_length=0.1*np.fabs(child2.x-child1.x),zorder=9, color='r')

            if event_flag in [2,12]:
                if child1.index == index1:
                    plot.scatter([child1.x],[child1.y],color=color,s=3*s,zorder=9,marker='*')
                if child2.index == index1:
                    plot.scatter([child2.x],[child2.y],color=color,s=3*s,zorder=9,marker='*')


        return x_min,x_max,y_min,y_max

    @staticmethod
    def draw_bodies(plot,bodies,fontsize,y_ref=1.0,dx=0.5,dy=0.5,index1=-1,index2=-1,event_flag=-1):

        for index,body in enumerate(bodies):
            color,s,description = Tools.get_color_and_size_and_description_for_star(body.stellar_type,body.radius)
            plot.scatter([index],[y_ref],color=color,s=s)
            text = "$\mathrm{%s}$"%(str(round(body.mass,1)))
            plot.annotate(text,xy=(index - dx,y_ref-dy),color='k',fontsize=fontsize)
            
            if body.index == index1:
                if event_flag in [2,12]:
                    plot.scatter([index],[y_ref],color=color,s=3*s,zorder=9,marker='*')

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
    def get_color_and_size_and_description_for_star(stellar_type,radius):
        if (stellar_type == 0): 
            color='gold'
            description = 'low-mass\, ' + 'MS'
        elif (stellar_type == 1): 
            color='gold'
            description = 'MS'
        elif (stellar_type == 2):
            color='darkorange'
            description = 'HG'
        elif (stellar_type == 3):
            color='firebrick'
            description = 'RGB'
        elif (stellar_type == 4):
            color='darkorange'
            description = 'CHeB'
        elif (stellar_type == 5):
            color='orangered'
            description = 'EAGB'
        elif (stellar_type == 6):
            color='crimson'
            description = 'TPAGB'
        elif (stellar_type == 7):
            color='royalblue'
            description = 'HeMS'
        elif (stellar_type == 8):
            color='orangered'
            description = 'HeHG'
        elif (stellar_type == 9):
            color='crimson'
            description = 'HeGB'
        elif (stellar_type == 10):
            color='silver'
            description = 'HeWD'
        elif (stellar_type == 11):
            color='silver'
            description = 'COWD'
        elif (stellar_type == 12):
            color='silver'
            description = 'ONeWD'
        elif (stellar_type == 13):
            color='gainsboro'
            description = 'NS'
        elif (stellar_type == 14):
            color='k'
            description = 'BH'
        else: 
            color = 'k'
            description = ''
        

        CONST_R_SUN = 0.004649130343817401
        CONST_KM = 1.0/(1.4966e9)

        s = 5 + 50*np.log10(radius/CONST_R_SUN)
        if (stellar_type >= 7 and stellar_type <= 9):
            s = 5 + 5*np.log10(radius/CONST_KM)

        if (stellar_type >= 10):
            s = 5 + 20*np.log10(radius/CONST_KM)

        return color,s,description
        
    @staticmethod
    def get_description_for_event_flag(event_flag):
        if event_flag == 0:
            text = "$\mathrm{Initial\,system}$"
        elif event_flag == 1:
            text = "$\mathrm{Stellar\,type\,change}$"
        elif event_flag == 2:
            text = "$\mathrm{SNe\,start}$"
        elif event_flag == 3:
            text = "$\mathrm{SNe\,end}$"
        elif event_flag == 4:
            text = "$\mathrm{RLOF\,start}$"
        elif event_flag == 5:
            text = "$\mathrm{RLOF\,end}$"
        elif event_flag == 6:
            text = "$\mathrm{CE\,start}$"
        elif event_flag == 7:
            text = "$\mathrm{CE\,end}$"
        elif event_flag == 8:
            text = "$\mathrm{Collision\,start}$"
        elif event_flag == 9:
            text = "$\mathrm{Collision\,end}$"
        elif event_flag == 10:
            text = "$\mathrm{Dyn.\,inst.}$"
        elif event_flag == 11:
            text = "$\mathrm{Sec.\,break.}$"
        elif event_flag == 12:
            text = "$\mathrm{WD\,kick\,start}$"
        elif event_flag == 13:
            text = "$\mathrm{WD\,kick\,end}$"
        elif event_flag == 14:
            text = "$\mathrm{Triple\,CE\,start}$"
        elif event_flag == 15:
            text = "$\mathrm{Triple\,CE\,end}$"
        else:
            text = ""
        return text
