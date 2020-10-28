#
# Makefile for su2ToCappuccino program
#         

F90FLAGS = -O3 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      utils.f90 \
      su2ToCappuccino.f90 



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

su2ToCappuccino: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o ../../bin/su2ToCappuccino 
clean:
	@rm  *.o *.mod ../../bin/su2ToCappuccino 

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
