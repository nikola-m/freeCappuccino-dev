# This is a makefile for GNU make.

# This makefile builds "Cappuccino".

# 1. Source files

# Directory containing the main makefile, used by "VPATH":
makefile_dir = src/cappuccino

# Write paths relative to makefile_dir, separated by blanks (do not
# include ${makefile_dir} in the list):
# VPATH := ${makefile_dir} $(addprefix ${makefile_dir}/, dir)

# Write paths relative to makefile_dir in 'directories' file, if there are many source directories:
VPATH := ${makefile_dir} $(shell cat ${makefile_dir}/directories)

# Write source file names, without path, separated by blanks, not
# commas. Only .f* files, not included files. Also the file
# containinig the main program unit.
# sources := 

# Write source file names, without path, in a separate file,
# one file name per line:
sources := $(shell cat ${makefile_dir}/files)
# NOTE: In that file the file names are sorted so the file containing a module
# is listed before the file where the module has been used.

# 2. Objects and executable file

#objects := $(addsuffix .o, $(basename ${sources}))
# Or, if all the source files have the same suffix, more simply:
objects := $(sources:.f=.o)

execut = bin/cappuccino

# 3. Compiler-dependent part

# 4. Rules

FC = gfortran

# add -DLIS to define LIS library flag;
# add -DLAPACK to define Lapack library flag;
# The -cpp flag is here because of C style preprocessor directives.
#LDFLAGS = -O3 -cpp -DLAPACK -DLIS -Wall
LDFLAGS = -O3 -cpp -Wall

# Lapack is used for LSTSQ_QR gradients, LIS may be used for linear systems.
# If you don't use either LAPACK or LIS don't link these libraries below.
# OpenMP is used within LIS, only in a case you've compiled LIS with OpenMP support.
#LDLIBS = -llapack -llis -fopenmp
#LDLIBS = -llapack
LDLIBS =

# Extend known suffixes:
# (if the suffixes are not ".f" or ".F")

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

%.o: %.F90
	$(COMPILE.F) $(OUTPUT_OPTION) $<

.DELETE_ON_ERROR:
.PHONY: all clean

all: ${execut}

${execut}: ${objects}
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

clean:
	rm -f ${execut} *mod

clobber:
	rm -f *mod
