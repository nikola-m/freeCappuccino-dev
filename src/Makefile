#
# Makefile for caffa3d program
#

# Compiler:
F90 = gfortran

# Compiler flags:
LFLAGS = -llapack
F90FLAGS = -Wall -O2

MOD_FILES=\
    modules_allocatable.f90 \
    utils.f90 \
    matrix.f90 \
    geometry.f90 \
    tensor_fields.f90 \
    field_initialization.f90\
    sparse_matrix.f90 \
    gradients.f90 \
    output.f90 \
    interpolation.f90 \
    scalar_fluxes.f90


LINEAR_SOLVER_FILES=\
    iccg.f90 \
    bicgstab.f90 \
    dpcg.f90 \
    mgmres.f90 \
    gauss-seidel.f90

TURBULENCE_FILES=\
    temperature.f90 \
    concentration.f90 \
    k_epsilon_std.f90\
    k_omega_sst.f90 \
    k_epsilon_rng.f90 \
    spalart_allmaras.f90 \
    k_eqn_eddy.f90
#   k_epsilon_zeta_f.f90 \

CAFFA_FILES=\
    allocate.f90 \
    asm_stress_terms.f90 \
    asm_heatflux_terms.f90 \
    adjustMassFlow.f90 \
    bcin.f90 \
    bpres.f90 \
    fieldManipulation.f90 \
    faceflux_velocity.f90 \
    facefluxmass.f90 \
    calcheatflux.f90 \
    calcp-multiple_correction_SIMPLE.f90 \
    calcstress.f90 \
    calc_strain_and_vorticity.f90 \
    calcuvw.f90 \
    calc_statistics.f90 \
    updateVelocityAtBoundary.f90 \
    correct_turbulence.f90 \
    correct_turbulence_inlet.f90 \
    fvm_laplacian.f90 \
    init.f90 \
    openfiles.f90 \
    PISO_multiple_correction.f90 \
    readfiles.f90 \
    read_input.f90 \
    time_shift.f90 \
    writefiles.f90 \
    write_restart_files.f90 \
    writehistory.f90 \
    main.f90 


POISSON_FILES=\
    fvm_laplacian.f90 \
    poisson.f90

STEADYHEAT_FILES=\
    fvm_laplacian.f90 \
    steadyHeatEq.f90


#
# How to create object files:
# 
MODS = ${MOD_FILES:.f90=.o}
TURBULENCE = ${TURBULENCE_FILES:.f90=.o}
LINEAR_SOLVERS = ${LINEAR_SOLVER_FILES:.f90=.o}
CAFFAOBJS = ${CAFFA_FILES:.f90=.o}
POISSONOBJS = ${POISSON_FILES:.f90=.o}
STEADYHEATOBJS = ${STEADYHEAT_FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: caffa3d poisson steadyHeatEq

caffa3d: ${MODS} ${TURBULENCE} ${LINEAR_SOLVERS} ${CAFFAOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${CAFFAOBJS} ${MODS} ${TURBULENCE} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o ../bin/caffa3d 

poisson: ${MODS} ${LINEAR_SOLVERS} ${POISSONOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${POISSONOBJS} ${MODS} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o ../bin/poisson 

steadyHeatEq: ${MODS} ${LINEAR_SOLVERS} ${STEADYHEATOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${STEADYHEATOBJS} ${MODS} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o ../bin/steadyHeatEq

.PHONY: clean
clean:
	@rm  *.o *.mod ../bin/caffa3d ../bin/poisson ../bin/steadyHeatEq

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90

.f90.o:
	${F90} ${F90FLAGS} -c ${INCS}  ${@:.o=.f90}
