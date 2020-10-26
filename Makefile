CXX = c++-10
CXX = g++-10
CXX = g++

FC = gfortran

CXXSRC = interface.cpp src/types.cpp src/tools.cpp src/evolve.cpp src/structure.cpp src/ODE_system.cpp src/ODE_root_finding.cpp src/ODE_newtonian.cpp src/ODE_postnewtonian.cpp src/ODE_tides.cpp src/ODE_mass_changes.cpp src/ODE_VRR.cpp src/external.cpp src/stellar_evolution.cpp src/SNe.cpp src/flybys.cpp src/binary_evolution.cpp src/collision.cpp src/nbody_evolution.cpp src/test.cpp src/apsidal_motion_constant.cpp src/common_envelope_evolution.cpp src/logging.cpp src/parameters.cpp
CSRC = src/cvode/cvode.c src/cvode/cvode_dense.c src/cvode/cvode_direct.c src/cvode/cvode_io.c src/cvode/nvector_serial.c src/cvode/sundials_dense.c src/cvode/sundials_direct.c src/cvode/sundials_math.c src/cvode/sundials_nvector.c src/mstar/mst.c src/mstar/pn.c
FSRC = src/sse/evolv1.f src/sse/zcnsts.f src/sse/deltat.f src/sse/hrdiag.f src/sse/kick.f src/sse/mlwind.f src/sse/mrenv.f src/sse/ran3.f src/sse/star.f src/sse/zfuncs.f src/bse/dgcore.f src/bse/corerd.f src/bse/gntage.f

COBJ = $(CXXSRC:.cpp=.o) $(CSRC:.c=.o)
FOBJ = $(FSRC:.f=.o)
CXXHEADERS = $(CXXSRC:.cpp=.h)
CHEADERS = $(CSRC:.c=.h)
FHEADERS = src/sse/const_bse.h src/sse/zdata.h 

#CPPFLAGS = -fPIC -shared -O2 -lgfortran -Wno-comment -Wno-c++11-compat-deprecated-writable-strings -Wno-write-strings
CPPFLAGS = -fPIC -shared -lgfortran -Wno-comment -Wno-c++11-compat-deprecated-writable-strings -Wno-write-strings -g
FFLAGS = -fPIC

ifeq ($(DEBUG),1)
        CPPFLAGS += -DDEBUG
endif

all: $(COBJ) libmse.so

%.o: %.c $(CHEADERS)
	@echo "Compiling C source file $< ..."
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%.o: %.cpp $(CXXHEADERS)
	@echo "Compiling C++ source file $< ..."
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%.o: %.f $(FHEADERS)
	@echo "Compiling Fortran source file $< ..."
	$(FC) -c -o $@ $< $(FFLAGS)

#sse: $(FOBJ)
#	@echo $(FOBJ)
#	@echo "building sse"
#	$(FC) -c $(FFLAGS) $(FOBJ)
#	$(FC) -c -o $@ $(FOBJ) $(FFLAGS)
	
libmse.so: $(COBJ) $(FOBJ)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CXX) -o $@ $(COBJ) $(FOBJ) $(CPPFLAGS)
	@echo ""        
	@echo "The shared library $@ has been created successfully."

cleansse:
	@echo "Removing compiled sse/bse libraries"
	$(RM) src/sse/*.o*
	$(RM) src/bse/*.o*
cleanc:
	@echo "Removing compiled C/C++ libraries"
	$(RM) interface.o libmse.so src/*.o* src/cvode/*.o* src/mstar/*.o*
cleanlib:
	@echo "Removing libmse library"
	$(RM) libmse.so
clean: cleansse cleanc cleanlib
