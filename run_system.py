"""
Run a system from the command line
Adrian Hamers, September 2020
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
    parser.add_argument("--is", "--inclinations",                     type=float,  nargs="+",     dest="inclinations",                      default=[],                         help="Inclinations (rad)")
    parser.add_argument("--LANs", "--longitudes_of_ascending_node",   type=float,  nargs="+",     dest="longitudes_of_ascending_node",      default=[],                         help="Longitudes of the ascending node (rad)")
    parser.add_argument("--APs", "--arguments_of_pericentre",         type=float,  nargs="+",     dest="arguments_of_pericentre",           default=[],                         help="Arguments of periapsis (rad)")
    
    parser.add_argument("--tend",                           type=float,     dest="end_time",                    default=2.0e7,             help="Integration time (yr)")
    parser.add_argument("--Nsteps",                         type=int,       dest="N_steps",                     default=2000,              help="Number of plot output steps")
    
    parser.add_argument("--random_seed",                    type=int,       dest="random_seed",                 default=0,                 help="Random seed")
    
    parser.add_argument("--plot_filename",                  type=str,       dest="plot_filename",               default="test1",           help="Plot filename")
    
    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                         default=False,         help="Verbose terminal output")
    add_bool_arg(parser, 'plot',                            default=False,         help="Make plots")
    add_bool_arg(parser, 'show_plots',                      default=True,          help="Display plots")
    add_bool_arg(parser, 'fancy_plots',                     default=False,         help="Use LaTeX fonts for plots (slow)")
    
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
    print("Arguments of periapsis (rad): ",args.inclinations)
    print("Integration time (yr): ",args.end_time)
    print("Number of plot output steps: ",args.N_steps)

    print("="*50)
    
    Tools.evolve_system(args.configuration,N_bodies,args.masses,args.metallicities,args.semimajor_axes,args.eccentricities,args.inclinations,args.arguments_of_pericentre,args.longitudes_of_ascending_node,args.end_time,args.N_steps,stellar_types=args.stellar_types,plot_filename=args.plot_filename,object_types=args.object_types,fancy_plots=args.fancy_plots,show_plots=args.show_plots,random_seed=args.random_seed)
