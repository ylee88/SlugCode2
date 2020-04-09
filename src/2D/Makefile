OS = $(shell uname -s)

include Makefile.header

EXE_FILE = slugCode2

FC = gfortran
MC = mpif90

# # TAU
# MC = tau_f90.sh -tau_makefile=$(TAU_DIR)/lib/Makefile.tau-mpi -tau_options=-optCompInst

ifeq ($(OS),Darwin)
	LDFLAGS = -framework accelerate
endif
LDFLAGS += -lhdf5_fortran

# double precision
FFLAGS_OPT2 = -O3 -fdefault-real-8 -fdefault-double-8\
	-ffree-line-length-none

# quadruple precision
FFLAGS_OPT4 = -O3 -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none

# double precision
FFLAGS_DEBUG2 = -ggdb  -g -fdefault-real-8 -fdefault-double-8\
	-ffree-line-length-none -Wall

# quadruple precision
FFLAGS_DEBUG4 = -ggdb  -g -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wall

FFLAGS_DOUB = $(FFLAGS_DEBUG2)
FFLAGS_QUAD = $(FFLAGS_DEBUG4)

OBJS  = slug_driver.o         \
		abort_slug.o          \
		cfl.o                 \
		io.o                  \
		read_pars.o           \
		read_initFile.o       \
		bc.o                  \
		primconsflux.o        \
		bc_init.o             \
		block_data.o          \
		block_init.o          \
		block_finalize.o      \
		grid_data.o           \
		grid_init.o           \
		grid_finalize.o       \
		sim_data.o            \
		sim_initBlock.o       \
		WENO.o                \
		char_limiting.o       \
		eigensystem.o         \
		get_maxSpeed.o        \
		num_data.o            \
		num_init.o            \
		num_interface.o       \
		sim_init.o            \
		soln_RK2.o            \
		soln_RK3.o            \
		soln_RK4.o            \
		soln_WENO.o           \
		soln_getFlux.o        \
		soln_numeric.o        \
		soln_spatial.o        \
		soln_update.o


################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(MC) $(FFLAGS_DOUB) $(OBJS) $(LDFLAGS) -L$(HDF5_DIR)/lib -o $(EXE_FILE)
	@echo "code is now linking..."

%.o: %.F90
	$(MC) $(FFLAGS_DOUB) -c $<

io.o: io.F90
	$(MC) $(FFLAGS_DOUB) -I$(HDF5_DIR)/include -c $<


################################################################################
#FIND MODULE DEPENDENCIES
Makefile.dep:
	@echo "Making dependencies"
	../utils/finddep.py -r ./ -e F90 -o Makefile.dep

# include dependencies
# '-' suppress warning messages
# conditional statement avoid to making dependencies file while clean command
ifneq ($(MAKECMDGOALS),clean)
-include Makefile.dep
endif


################################################################################
#SOME USEFUL COMMANDS
.PHONY: clean

clean:
	rm -f *.o *.mod *~ $(EXE_FILE) Makefile.dep
