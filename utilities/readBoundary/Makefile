#
# Makefile for readBoundaryFile program
#         

F90FLAGS = -O2 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      boundary_module.f90 \
      readBoundaryFile.f90 



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

readBoundaryFile: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o ../../bin/readBoundaryFile
clean:
	@rm  *.o *.mod ../../bin/readBoundaryFile

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
