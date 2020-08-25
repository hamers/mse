# mse
Multiple Stellar Evolution (`MSE`) -- A Population Synthesis Code for Multiple-Star Systems

!!! This is a very early-stage version with still many bugs and missing features. !!!

A code to model the long-term evolution of hierarchical multiple-star systems (binaries, triples, quadruples, and higher-order systems) from the main sequence until  remnant phase. Takes into account gravitational dynamical evolution, stellar evolution (using the `sse` tracks), and binary interactions (such as mass transfer and common-envelope evolution). 
    
Includes routines for external perturbations from flybys in the field, or (to limited extent) encounters in dense stellar systems such as galactic nuclei. 

C++ compiler and Fortran compilers are required, as well as Python (2/3) for the Python interface. Make sure to first compile the code using `make`. Please modify the Makefile according to your installation (`CXX` and `FC` should be correctly assigned).  

The script `test_mse.py` can be used to test the
installation. 

Adrian Hamers, August 2020
