#
# Makefile for foamToCappuccino program
#         

F90FLAGS = -O3 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      utils.f90 \
      foamToCappuccino.f90 



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

foamToCappuccino: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o ../../bin/foamToCappuccino 
clean:
	@rm  *.o *.mod ../../bin/foamToCappuccino 

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
