#
# Makefile for gambitToCappuccino program
#         

F90FLAGS = -O3 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      utils.f90 \
      qsort_c_module.f90 \
      gambitToCappuccino.f90 



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

gambitToCappuccino: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o ../../bin/gambitToCappuccino 
clean:
	@rm  *.o *.mod ../../bin/gambitToCappuccino 

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
