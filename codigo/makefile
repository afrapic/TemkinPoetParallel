#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
#	Modified Makefile from vasp, to compile code
#     	with scalapack bindings.
#
#
# 	Makefile for Intel Fortran compiler
#	See example in 
#	Intel Math Kernel Library User's guide version 10.0.xxx
#	
#	WARNING! WARNING! WARNING! WARNING! WARNING!
#
#	Run-time errors ocurr simply changing the order of linking
#	libraries. Hence check ***ALL**** options in this makefile very
#	carefully if you make changes.
#
#
#	To use this makefile, please first run
#       source /opt/intel/impi/3.1/bin64/mpivars.sh
#       to get the right variables in the environment.
#
#
#
#-----------------------------------------------------------------------

.SUFFIXES: .inc .f .f90 .F
SUFFIX=.f90

#-----------------------------------------------------------------------
# fortran compiler and linker
#-----------------------------------------------------------------------
#FC=/opt/mpich/intel/bin/mpif90 does not handle 64 bits pointers
FC=mpiifort

# fortran linker

FCL=$(FC)

#-----------------------------------------------------------------------
# Intel integer model
#  lp64 handles int as 32 bits
#  ilp64 handles int as 64 bits
#-----------------------------------------------------------------------

LP=lp64

#-----------------------------------------------------------------------
# general fortran flags  (there must a trailing blank on this line)
#-----------------------------------------------------------------------

FFLAGS = 
# Extensive runtime checks
#FFLAGS = -C -CB -CA -CU

OFLAG=-fast -heap-arrays
# Debug
#OFLAG= -g -heap-arrays

OFLAG_HIGH = $(OFLAG)
OBJ_HIGH =

OBJ_NOOPT =
#DEBUG  = -FR -O0
INLINE = $(OFLAG)
#
#	Where to look libraries. em64t is 64 bits with 
#	Pentium 4 - Core - Core 2 architecture. Xeon, Itanium, etc
#	are simply 64
#

LIBDIR = /opt/intel/mkl/10.0.1.014/lib/em64t

#-----------------------------------------------------------------------
#
#	BLAS
#
#-----------------------------------------------------------------------
#
# BLAS must be installed on the machine
# there are several options:
# 1) very slow but works:
#   retrieve the lapackage from ftp.netlib.org
#   and compile the blas routines (BLAS/SRC directory)
#   please use g77 or f77 for the compilation. When I tried to
#   use pgf77 or pgf90 for BLAS, VASP hang up when calling
#   ZHEEV  (however this was with lapack 1.1 now I use lapack 2.0)
# 2) most desirable: get an optimized BLAS
#
# the two most reliable packages around are presently:
# 3a) Intels own optimised BLAS (PIII, P4, Itanium)
#     http://developer.intel.com/software/products/mkl/
#   this is really excellent when you use Intel CPU's
#
# 3b) or obtain the atlas based BLAS routines
#     http://math-atlas.sourceforge.net/
#   you certainly need atlas on the Athlon, since the  mkl
#   routines are not optimal on the Athlon.
#   If you want to use atlas based BLAS, check the lines around LIB=
#
# 3c) mindblowing fast SSE2 (4 GFlops on P4, 2.53 GHz)
#   Kazushige Goto's BLAS
#   http://www.cs.utexas.edu/users/kgoto/signup_first.html
#
#-----------------------------------------------------------------------
# the following lines specify the position of BLAS  and LAPACK
#-----------------------------------------------------------------------
# faster Kazushige Goto's BLAS
# http://www.cs.utexas.edu/users/kgoto/signup_first.html
#BLAS=  /opt/libs/libgoto/libgoto_p4_512-r0.6.so

# mkl BLAS

BLAS= -L$(LIBDIR) -lmkl_intel_$(LP) -lmkl_sequential -lmkl_core -i-static
#BLAS= -L$(LIBDIR) -lmkl_intel_$(LP) -lmkl_core

#-----------------------------------------------------------------------
# options for linking 
#-----------------------------------------------------------------------
#LINK    = -lguide -lpthread -limf -lm 
#LINK    =  -liomp5 -limf -lm 


#=======================================================================
#	 	MPI section
#=======================================================================
#-----------------------------------------------------------------------
# fortran linker for mpi: if you use LAM and compiled it with the options
# suggested above,  you can use the following line
#-----------------------------------------------------------------------

FCL=$(FC)

#-----------------------------------------------------------------------
# location of SCALAPACK
#-----------------------------------------------------------------------
#
#	Always use mkl scalapack if available
#
SCA= -L$(LIBDIR) -lmkl_scalapack_$(LP) -lmkl_blacs_intelmpi_$(LP)

#-----------------------------------------------------------------------
# libraries for mpi
#-----------------------------------------------------------------------

MPILIB     = $(SCA)  $(LAPACK) $(BLAS) 

#-----------------------------------------------------------------------
# Project specific libraries
#-----------------------------------------------------------------------
FHOME   = /home/afrapic
ARCH = $(shell uname -m)
LDIR = $(FHOME)/lib
#
#   Sturmians base library
#
LIBNAME = base
BASELIB = $(LIBNAME)_$(ARCH).so
BASEINC = $(LDIR)/$(LIBNAME)/modules/

LOCALIB = $(LDIR)/$(BASELIB)

#-----------------------------------------------------------------------
# ALL Includes
#-----------------------------------------------------------------------

INC = -module $(MDIR) 
# -module $(BASEINC) 
INCDIR = -I$(MDIR) -I$(BASEINC)  

#-----------------------------------------------------------------------
# ALL libraries
#-----------------------------------------------------------------------

LIB =  $(MPILIB) $(LOCALIB)

#-----------------------------------------------------------------------
# general rules and compile lines
#-----------------------------------------------------------------------
BASIC=  
#
#  Simple mpi example
#
#SOURCE = hellompi.o
#
#  double precision example
#
#SOURCE = pdscaex.o pdscaexinfo.o pdlaread.o pdlawrite.o
#
#  original complex double precision example
#
#SOURCE = pzscaex.o pdscaexinfo.o pzlaread.o pzlawrite.o pzlaprnt.o
#
#  current complex double precision example v 0.2
#
#
#	 DIRS
#
ODIR = objs
SDIR = src
MDIR = modules
MAIN = main
#MODS := $(wildcard src/*module*) 
#MODULES = $(basename $(MODS))
VPATH = .:$(ODIR):$(SDIR):$(MDIR):$(BASEINC):$(MAIN)
#
#	Files
#
SRCS := $(wildcard src/*.f90 src/*.f src/*.for) 
BASE = $(notdir $(basename $(SRCS)))
OBJECTS =  $(addsuffix .o,$(BASE))
OBJS =  $(addprefix objs/,$(OBJECTS))

INCS=
EXE = ptp

all:  ptp 

ptp: $(OBJECTS)  $(INCS) main_ptp.o main/main_ptp.f90
	$(FCL) $(FFLAGS)  $(OBJS) $(ODIR)/main_ptp.o $(LIB) $(LINK) -o ptp

cs: $(OBJECTS)  $(INCS) main_cs.o main/main_cs.f90
	$(FCL) $(FFLAGS)  $(OBJS) $(ODIR)/main_cs.o $(LIB) $(LINK) -o cstp
clean: 
	rm -f *.g $(ODIR)/* *.L $(MDIR)/* src/*.*~ main/*.*~ *~

tpclean:
	rm -f *.o* *.po* tp*.dat

dep:    
	touch depend.mk depend.mk.old
	cd $(SDIR) ; ../sfmakedepend -f ../depend.mk  -I $(MDIR) -I $(BASEINC) *.f *.f90 ;	cd ..

echo:;
	@echo $(OBJECTS)
	@echo $(FFLAGS)
	@echo $(SRCS)
#
# MIND: I do not have a full dependency list for the include
# and MODULES: here are only the minimal basic dependencies
# if one strucuture is changed then touch_dep must be called
# with the corresponding name of the structure
#
#sc.o: sc.f
#hellompi.o: hellompi.f90
#problemdata_module.o: problemdata_module.f90
#problemdata_module.o: file_module.mod sturmians_data.mod
include depend.mk

#$(OBJ_HIGH):
#	$(FC) $(FFLAGS) $(OFLAG_HIGH) $(INCS) $*$(SUFFIX)
#$(OBJ_NOOPT):
#	$(FC) $(FFLAGS) $(INCS)  $*$(SUFFIX)


.F.o:
	$(FC) $(FFLAGS) $(OFLAG) $(INCS) $*$(SUFFIX)

.f.o:
	$(FC) $(FFLAGS) $(OFLAG) $(INCS) -c $< -o $(ODIR)/$*.o

%.o: %.f90
	$(FC)  $(FFLAGS) $(INC)  $(INCDIR)  $(OFLAG)  -c $< -o $(ODIR)/$*.o
