OS = $(shell uname -s)

include Makefile.header

EXE_FILE = slugCode2

# # TAU
# MC = tau_f90.sh -tau_makefile=$(TAU_DIR)/lib/Makefile.tau-mpi -tau_options=-optCompInst

ifeq ($(OS),Darwin)
    LDFLAGS += -framework accelerate
else
    LDFLAGS += -lblas -llapack
endif
LDFLAGS += -lhdf5_fortran

FFLAGS_DOUB = $(FFLAGS_DEBUG2)
FFLAGS_QUAD = $(FFLAGS_DEBUG4)

# opt: a fake target for optimization
ifeq ($(MAKECMDGOALS),opt)
    FFLAGS_DOUB = $(FFLAGS_OPT2)
    FFLAGS_QUAD = $(FFLAGS_OPT4)
endif
opt: $(EXE_FILE)

OBJS  = slug_driver.o         \
        abort_slug.o          \
        read_version.o        \
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
        eigensystem.o         \
        get_maxSpeed.o        \
        num_data.o            \
        num_init.o            \
        num_interface.o       \
        sim_init.o            \
        soln_RK2.o            \
        soln_RK3.o            \
        soln_RK4.o            \
        soln_sfPIF3.o         \
        soln_sfPIF4.o         \
        num_diff.o            \
        sfPIF.o               \
        WENO.o                \
        soln_WENO5.o          \
        gp_data.o             \
        gp_WENOinit.o         \
        gp_interface.o        \
        GP.o                  \
        linalg.o              \
        soln_gpWENO5.o        \
        soln_gpWENO7.o        \
        soln_numeric.o        \
        soln_spatial.o        \
        soln_update.o


################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(MC) $(FFLAGS_DOUB) $(OBJS) $(LDFLAGS) -L$(HDF5_LIB) -o $(EXE_FILE)
	@echo "code is now linking..."
	@git rev-parse HEAD >| GIT_VERSION

%.o: %.F90
	$(MC) $(FFLAGS_DOUB) -c $<

io.o: io.F90
	$(MC) $(FFLAGS_DOUB) -I$(HDF5_INC) -c $<


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
.PHONY: clean opt

clean:
	rm -f *.o *.mod *~ $(EXE_FILE) Makefile.dep GIT_VERSION
