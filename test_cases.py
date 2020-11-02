"""
A number of test systems
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
    add_bool_arg(parser, 'verbose',                         default=False,         help="Verbose terminal output")
    add_bool_arg(parser, 'plot',                            default=False,         help="Make plots")
    add_bool_arg(parser, 'show_plots',                      default=False,          help="Display plots")
    add_bool_arg(parser, 'fancy_plots',                     default=False,         help="Use LaTeX fonts for plots (slow)")
    parser.add_argument("--verbose_flag",                   type=int,       dest="verbose_flag",                default=1,                 help="0: no verbose output in C++; > 0: verbose output, with increasing verbosity (>1 will slow down the code considerably) ")

    
    args = parser.parse_args()

    return args
    
class test_mse():

    def test1(self,args):
        print('Massive triple with stellar evolution')
        N_bodies = 3
        configuration="fully_nested"
        masses = [30.0,10.0,2.0]
        metallicities = [0.01,0.03,0.005]
        semimajor_axes = [15.0,120.0]
        eccentricities = [0.1,0.1]
        inclinations = [0.0001,85.0*np.pi/180.0]
        arguments_of_pericentre = [85.0*np.pi/180.0,0.01*np.pi/180.0]
        longitudes_of_ascending_node = [0.01,0.01]
        end_time = 5.0e7
        N_steps = 5000
        stellar_types = [1,1,1]
        Tools.evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,end_time,N_steps,stellar_types=stellar_types,fancy_plots=args.fancy_plots,show_plots=args.show_plots,plot_filename="test1",verbose_flag=args.verbose_flag)

    def test2(self,args):
        print('Compact object merger (initial compact objects)')
        N_bodies = 3
        configuration="fully_nested"
        masses = [50.0,40.0,7.5]
        metallicities = [0.01,0.03,0.005]
        semimajor_axes = [0.01,400.0]
        eccentricities = [0.1,0.6]
        inclinations = [0.0001,85.0*np.pi/180.0]
        arguments_of_pericentre = [85.0*np.pi/180.0,0.01*np.pi/180.0]
        longitudes_of_ascending_node = [0.01,0.01]
        end_time = 5.0e7
        N_steps = 1000
        stellar_types = [13,13,13]
        Tools.evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,end_time,N_steps,stellar_types=stellar_types,fancy_plots=args.fancy_plots,show_plots=args.show_plots,plot_filename="test2",verbose_flag=args.verbose_flag)

    def test3(self,args):
        print('3+1 quadruple dynamically unstable due to secular evolution')
        N_bodies = 4
        configuration="fully_nested"
        masses = [3.0,2.0,1.5,2.0]

        metallicities = [0.02,0.02,0.02,0.02]
        semimajor_axes = [15.0,400.0,3000]
        eccentricities = [0.1,0.1,0.1]
        inclinations = [0.0001,85.0*np.pi/180.0,0.001]
        arguments_of_pericentre = [85.0*np.pi/180.0,0.01*np.pi/180.0,0.01*np.pi/180.0]
        longitudes_of_ascending_node = [0.01,0.01,0.01]
        end_time = 5.0e7
        N_steps = 5000
        stellar_types = [1,1,1,1]
        Tools.evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,end_time,N_steps,stellar_types=stellar_types,fancy_plots=args.fancy_plots,show_plots=args.show_plots,plot_filename="test3",verbose_flag=args.verbose_flag)
        
    def test4(self,args):
        print('Initially dynamically unstable triple')
        N_bodies = 3
        configuration="fully_nested"
        masses = [1.0,1.0,20.5]
        metallicities = [0.01,0.03,0.005]
        semimajor_axes = [15.5,120.0]
        eccentricities = [0.1,0.1]
        inclinations = [0.0001,15.0*np.pi/180.0]
        arguments_of_pericentre = [85.0*np.pi/180.0,0.01*np.pi/180.0]
        longitudes_of_ascending_node = [0.01,0.01]
        end_time = 3.0e7
        N_steps = 500
        stellar_types = [1,1,1]
        Tools.evolve_system(configuration,N_bodies,masses,metallicities,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node,end_time,N_steps,stellar_types=stellar_types,fancy_plots=args.fancy_plots,show_plots=args.show_plots,plot_filename="test4",verbose_flag=args.verbose_flag)

  
if __name__ == '__main__':
    args = parse_arguments()
    
    N_tests = 4
    if args.test==0:
        tests = range(1,N_tests+1)
    else:
        tests = [args.test]

    t=test_mse()
    for i in tests:
        print( 'Running system number',i,'; verbose =',args.verbose,'; plot =',args.plot)
        function = getattr(t, 'test%s'%i)
        function(args)
    
