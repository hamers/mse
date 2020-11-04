#!/bin/sh
mkdir -p figs
python3 run_system.py --plot_filename figs/ex1 --fancy_plots --no-show_plots --masses 3.0 2.0 1.0 --smas 15 500 --es 0.1 0.8 --is 0.01 1.5 --tend 5.0e8 --Nsteps 1000 &
python3 run_system.py --plot_filename figs/ex2 --fancy_plots --no-show_plots --masses 20.0 15.0 10.0 20.0 --smas 20 800 8000 --es 0.1 0.2 0.5 --is 0.01 1.5 0.01 --tend 1.0e7 --Nsteps 100 &
python3 run_system.py --plot_filename figs/ex3 --fancy_plots --no-show_plots --masses 4 0.0001 5 2 --smas 18 10 1000 --ots 1 2 1 1 --es 0.01 0.01 0.41 --is 0.01 0.01 1.0 --configuration "2+2_quadruple" --tend 5.0e8 --Nsteps 1000 --verbose_flag 1 &
