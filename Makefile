CXX = /usr/local/bin/g++-9
CXX = /usr/local/bin/c++-9
CXX = c++-9
CXX = g++-9

FC = gfortran

CXXSRC = interface.cpp src/types.cpp src/evolve.cpp src/structure.cpp src/ODE_system.cpp src/root_finding.cpp src/newtonian.cpp src/postnewtonian.cpp src/tides.cpp src/external.cpp src/VRR.cpp src/stellar_evolution.cpp src/SNe.cpp src/flybys.cpp src/tools.cpp src/mass_changes.cpp src/binary_evolution.cpp src/collision.cpp src/nbody_evolution.cpp src/test.cpp
CSRC = src/cvode/cvode.c src/cvode/cvode_dense.c src/cvode/cvode_direct.c src/cvode/cvode_io.c src/cvode/nvector_serial.c src/cvode/sundials_dense.c src/cvode/sundials_direct.c src/cvode/sundials_math.c src/cvode/sundials_nvector.c src/mstar/mst.c src/mstar/pn.c
FSRC = src/sse/evolv1.f src/sse/zcnsts.f src/sse/deltat.f src/sse/hrdiag.f src/sse/kick.f src/sse/mlwind.f src/sse/mrenv.f src/sse/ran3.f src/sse/star.f src/sse/zfuncs.f src/bse/dgcore.f src/bse/corerd.f src/bse/gntage.f

#CXXHEADERS = interface.h src/types.h src/evolve.h src/structure.h src/ODE_system.h src/root_finding.h src/newtonian.h src/postnewtonian.h src/tides.h src/external.h src/VRR.h
#CHEADERS = src/cvode/cvode.h src/cvode/cvode_dense.h src/cvode/cvode_direct.h src/cvode/cvode_io.h src/cvode/nvector_serial.h src/cvode/sundials_dense.h src/cvode/sundials_direct.h src/cvode/sundials_math.h src/cvode/sundials_nvector.h
#COBJ = interface.o src/types.o src/evolve.o src/structure.o src/ODE_system.o src/root_finding.o src/newtonian.o src/postnewtonian.o src/tides.o src/external.o src/VRR.o src/cvode/cvode.o src/cvode/cvode_dense.o src/cvode/cvode_direct.o src/cvode/cvode_io.o src/cvode/nvector_serial.o src/cvode/sundials_dense.o src/cvode/sundials_direct.o src/cvode/sundials_math.o src/cvode/sundials_nvector.o

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
	$(RM) libmse.so src/*.o* src/cvode/*.o* src/mstar/*.o*
cleanlib:
	@echo "Removing libmse library"
	$(RM) libmse.so
clean: cleansse cleanc cleanlib
