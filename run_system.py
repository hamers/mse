"""
Run a system with Python directly from the command line.
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

    parser.add_argument("--configuration","--config",                 type=str,                   dest="configuration",                     default="fully_nested",             help="Configuration: `fully_nested', `2+2_quadruple', or general (see documentation) ")
    parser.add_argument("--masses", "--ms",                           type=float,  nargs="+",     dest="masses",                            default=[40.0,10.0,2.0],            help="Masses/MSun")
    parser.add_argument("--object_types", "--ots",                    type=int,    nargs="+",     dest="object_types",                      default=[],                         help="Object type (integer): 1: star, 2: planet")
    parser.add_argument("--stellar_types", "--sts",                   type=int,    nargs="+",     dest="stellar_types",                     default=[],                         help="Initial SSE stellar types (applies to stars only, i.e., object type = 1)")
    parser.add_argument("--metallicities", "--zs",                    type=float,  nargs="+",     dest="metallicities",                     default=[],                         help="Metallicities")
    parser.add_argument("--smas", "--semimajor_axes",                 type=float,  nargs="+",     dest="semimajor_axes",                    default=[15.0,120.0],               help="Semimajor axes (au)")
    parser.add_argument("--es",  "--eccentricities",                  type=float,  nargs="+",     dest="eccentricities",                    default=[0.1,0.1],                  help="Eccentricities")
    parser.add_argument("--is", "--inclinations",                     type=float,  nargs="+",     dest="inclinations",                      default=[0.0,1.0],                         help="Inclinations (rad)")
    parser.add_argument("--LANs", "--longitudes_of_ascending_node",   type=float,  nargs="+",     dest="longitudes_of_ascending_node",      default=[],                         help="Longitudes of the ascending node (rad)")
    parser.add_argument("--APs", "--arguments_of_pericentre",         type=float,  nargs="+",     dest="arguments_of_pericentre",           default=[],                         help="Arguments of periapsis (rad)")
    
    parser.add_argument("--tend",                           type=float,     dest="end_time",                    default=2.0e7,             help="Integration time (yr)")
    parser.add_argument("--Nsteps",                         type=int,       dest="N_steps",                     default=2000,              help="Number of plot output steps")
    
    parser.add_argument("--random_seed",                    type=int,       dest="random_seed",                 default=0,                 help="Random seed")
    parser.add_argument("--verbose_flag",                   type=int,       dest="verbose_flag",                default=1,                 help="0: no verbose output in C++; > 0: verbose output, with increasing verbosity (>1 will slow down the code considerably) ")
    
    parser.add_argument("--plot_filename",                  type=str,       dest="plot_filename",               default="test1",           help="Plot filename")
    
    parser.add_argument("--kick_distribution_sigma_km_s_NS",type=float,     dest="kick_distribution_sigma_km_s_NS", default=265.0,         help="NS kick sigma in km/s (assuming Maxwellian distribution)")

    parser.add_argument("--kick_distribution_sigma_km_s_BH",type=float,     dest="kick_distribution_sigma_km_s_BH", default=50.0,          help="BH kick sigma in km/s (assuming Maxwellian distribution)")

    parser.add_argument("--kick_distribution_sigma_km_s_WD",type=float,     dest="kick_distribution_sigma_km_s_WD", default=1.0,           help="WD kick sigma in km/s (assuming Maxwellian distribution); NOTE: include_WD_kicks should be toggled on to enable")
    parser.add_argument("--NS_model",                       type=int,       dest="NS_model",                    default=0,                 help="Model assumed for NSs: default (0), or Ye+19 (1)")
    parser.add_argument("--ECSNe_model",                    type=int,       dest="ECSNe_model",                 default=0,                 help="Model assumed for electron-capture SNe: default (0), or Ye+19 (1)")
    
    parser.add_argument("--flybys_stellar_density_per_cubic_pc",              type=float,       dest="flybys_stellar_density_per_cubic_pc",              default=0.1,        help="Density assumed for fly-bys (units: per cubic pc)")
    parser.add_argument("--flybys_encounter_sphere_radius_au",                type=float,       dest="flybys_encounter_sphere_radius_au",                default=1.0e5,      help="Encounter sphere radius assumed for fly-bys (au)")
    parser.add_argument("--flybys_stellar_relative_velocity_dispersion_km_s", type=float,       dest="flybys_stellar_relative_velocity_dispersion_km_s", default=30.0,       help="Relative velocity dispersion assumed for fly-bys (km/s)")
    
    parser.add_argument("--wall_time_max_s",                type=float,      dest="wall_time_max_s", default=1.8e4,       help="Maximum wall time allowed for the simulation (in s). Increase if needed.")

    
    
    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                         default=False,         help="Verbose terminal output (Python)")
    add_bool_arg(parser, 'make_plots',                      default=True,         help="Make plots")
    add_bool_arg(parser, 'plot_only',                       default=False,         help="Make plots by loading previously-generated data from disk; will not run the simulation")
    add_bool_arg(parser, 'save_data',                       default=True,          help="Save output data to pickle file for later replotting. Turn off if you do not use this feature and the pickle files are large.")
    add_bool_arg(parser, 'show_plots',                      default=True,          help="Display plots")
    add_bool_arg(parser, 'fancy_plots',                     default=False,         help="Use LaTeX fonts for plots (slower)")

    add_bool_arg(parser, 'include_WD_kicks',                default=False,         help="Let WDs receive a natal kick at birth. Assumes a Maxwellian distribution with \sigma = 1 km/s by default (can be changed).")
    add_bool_arg(parser, 'include_flybys',                  default=True,          help="Take into account fly-bys")
    add_bool_arg(parser, 'flybys_include_secular_encounters',default=False,        help="Include secular fly-bys")

    args = parser.parse_args()

    return args
    
  
if __name__ == '__main__':
    args = parse_arguments()
    
    N_bodies = len(args.masses)
    print("="*50)
    print("Parsed arguments:")
    print("Configuration: ",args.configuration)
    print("N_bodies: ",N_bodies)
    print("Object types:",args.object_types)
    print("Stellar types:",args.stellar_types)
    print("Masses/MSun: ",args.masses)
    print("Metallicities: ",args.metallicities)
    print("Semimajor axes (au): ",args.semimajor_axes)
    print("Eccentricities: ",args.eccentricities)
    print("Inclinations (rad): ",args.inclinations)
    print("Longitudes of the ascending node (rad): ",args.longitudes_of_ascending_node)
    print("Arguments of periapsis (rad): ",args.arguments_of_pericentre)
    print("Integration time (yr): ",args.end_time)
    print("Number of plot output steps: ",args.N_steps)

    print("="*50)
    
    error_code, log = Tools.evolve_system(args.configuration,N_bodies,args.masses,args.metallicities,args.semimajor_axes,args.eccentricities,args.inclinations,args.arguments_of_pericentre,args.longitudes_of_ascending_node,args.end_time,args.N_steps,stellar_types=args.stellar_types,plot_filename=args.plot_filename,object_types=args.object_types,fancy_plots=args.fancy_plots,show_plots=args.show_plots,random_seed=args.random_seed,verbose_flag=args.verbose_flag,include_WD_kicks=args.include_WD_kicks,kick_distribution_sigma_km_s_WD=args.kick_distribution_sigma_km_s_WD,NS_model=args.NS_model,ECSNe_model=args.ECSNe_model,kick_distribution_sigma_km_s_NS=args.kick_distribution_sigma_km_s_NS,kick_distribution_sigma_km_s_BH=args.kick_distribution_sigma_km_s_BH,flybys_stellar_density_per_cubic_pc=args.flybys_stellar_density_per_cubic_pc,flybys_encounter_sphere_radius_au=args.flybys_encounter_sphere_radius_au,flybys_stellar_relative_velocity_dispersion_km_s=args.flybys_stellar_relative_velocity_dispersion_km_s,flybys_include_secular_encounters=args.flybys_include_secular_encounters,include_flybys=args.include_flybys,plot_only=args.plot_only,save_data=args.save_data,make_plots=args.make_plots,wall_time_max_s=args.wall_time_max_s)
