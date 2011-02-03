#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
# Project name, version, and POC
AC_INIT(towhee,6.2.15,towhee-bugs@lists.sourceforge.net)
AC_CONFIG_SRCDIR(Source)
AC_CONFIG_AUX_DIR(install)

# Do not require the Makefiles to be rebuilt using the autotools
# so the end user will not have to have the same version of autoconf
# and automake.
AM_MAINTAINER_MODE

# Initialize stuff for automake
AM_INIT_AUTOMAKE

# Define the argument to configure script to set the linkline for spirit
AC_ARG_ENABLE(link_spirit,
     [  --enable-link-spirit   Compile version using the library link line for spirit],
     [case "${enableval}" in
       yes) UseLINKSPIRIT=true ;;
       no)  UseLINKSPIRIT=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-link-spirit) ;;
     esac;echo "UseLINKSPIRIT is $UseLINKSPIRIT"],[UseLINKSPIRIT=false])

# Define the argument to configure script to set the linkline for srnsquall
AC_ARG_ENABLE(link_srnsquall,
     [  --enable-link-srnsquall   Compile version using the library link line for srnsquall],
     [case "${enableval}" in
       yes) UseLINKSQUALL=true ;;
       no)  UseLINKSQUALL=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-link-srnsquall) ;;
     esac;echo "UseLINKSQUALL is $UseLINKSQUALL"],[UseLINKSQUALL=false])

# Define the argument to configure script to enable lcao 
AC_ARG_ENABLE(lcao,
     [  --enable-lcao   Compile version using lcao libraries],
     [case "${enableval}" in
       yes) UseLCAO=true ;;
       no)  UseLCAO=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-lcao) ;;
     esac;echo "UseLCAO is $UseLCAO"],[UseLCAO=false])

# Define the argument to configure script to use internal erfs
AC_ARG_ENABLE(internal,
     [  --enable-interal   Compile version using internal erfs],
     [case "${enableval}" in
       yes) UseInternal=true ;;
       no)  UseInternal=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-internal) ;;
     esac;echo "UseInternal is $UseInternal"],[UseInternal=false])

# Define the argument to configure script to enable tramonto
AC_ARG_ENABLE(tramonto,
     [  --enable-tramonto   Compile version using tramonto libraries],
     [case "${enableval}" in
       yes) UseTramonto=true ;;
       no)  UseTramonto=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-tramonto) ;;
     esac;echo "UseTramonto is $UseTramonto"],[UseTramonto=false])

# Define the argument to enable Safe Compare
 AC_ARG_ENABLE(safe_compare,
      [  --enable-safe-compare    Emphasize cross-platform consistency 
           over performance.],
      [case "${enableval}" in
        yes) UseSafeCompare=true ;;
        no)  UseSafeCompare=false ;;
        *) AC_MSG_ERROR(bad value ${enableval} for --enable-safe-compare) ;;
      esac;echo "UseSafeCompare is $UseSafeCompare"],[UseSafeCompare=false])

# Define the argument to enable fix GNU
 AC_ARG_ENABLE(fix_GNU,
      [  --enable-fix-GNU   Fixes the problem with double underscore.],
      [case "${enableval}" in
        yes) UseFixGNU=true ;;
        no)  UseFixGNU=false ;;
        *) AC_MSG_ERROR(bad value ${enableval} for --enable-fix-GNU) ;;
      esac;echo "UseFixGNU is $UseFixGNU"],[UseFixGNU=false])

# Define the argument to configure script to enable MPI
AC_ARG_ENABLE(mpi,
     [  --enable-mpi    Compile parallel version using MPI],
     [case "${enableval}" in
       yes) UseMPI=true ;;
       no)  UseMPI=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
     esac;echo "UseMPI is $UseMPI"],[UseMPI=false])

# Turn off Command Line Arguments, if necessary
AC_ARG_ENABLE(command-line-args,
    [  --disable-command-line-args   Disable parsing of command line arguments by towhee executable],
    [case "${enableval}" in
      yes) DisableCLA=false ;;
      no)  DisableCLA=true ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --disable-command-line-args) ;;
    esac;echo "DisableCLA is $DisableCLA"],[DisableCLA=false])

# Tell automake about the conditionals we just defined
AM_CONDITIONAL(UseLINKSPIRIT, test x$UseLINKSPIRIT = xtrue)
AM_CONDITIONAL(UseLINKSQUALL, test x$UseLINKSQUALL = xtrue)
AM_CONDITIONAL(UseMPI, test x$UseMPI = xtrue)
AM_CONDITIONAL(UseTramonto, test x$UseTramonto = xtrue)
AM_CONDITIONAL(UseInternal, test x$UseInternal = xtrue)
AM_CONDITIONAL(UseLCAO, test x$UseLCAO = xtrue)
AM_CONDITIONAL(UseSafeCompare, test x$UseSafeCompare = xtrue)
AM_CONDITIONAL(UseFixGNU, test x$UseFixGNU = xtrue)
AM_CONDITIONAL(DisableCLA, test x$DisableCLA = xtrue)

# Checks for compilers
# It seems really silly to check for cc and when we know
#   where it is, and will set it, but otherwise configure
#   breaks.
# This has to be done oustide a conditional.  Silly, huh?
AC_PROG_CC

# Compiler Tests
# Already did CC test
AC_PROG_F77

# Link normally
AC_SUBST(F77LD,${CC})

#incorporate ADDLIB into our library
LIBS="${LIBS} ${ADDLIB}"

# extra items for UseLINKSPIRIT
if test x$UseLINKSPIRIT = xtrue ; then
  # Include the srnspirit link-line
  LIBS="${LIBS} -L/apps/intel/mkl701/lib/em64t -L/apps/scalapack -L/apps/blacs"
  LIBS="${LIBS} -llcao -lmkl -lguide -lmkl_lapack64 -lscalapack" 
  LIBS="${LIBS} -lblacsCinit_MPI-LINUX-0 -lblacsF77init_MPI-LINUX-0 -lblacs_MPI-LINUX-0" 
fi
# extra items for UseLINKSQUALL
if test x$UseLINKSQUALL = xtrue ; then
  # Include the srnsquall link-line
  LIBS="${LIBS} -L/opt/acml/pgi64/lib"
  LIBS="${LIBS} -llcao -lscalapack -lblacsCinit -lblacsF77init -lblacs -lacml" 
fi
# extra items for UseTramonto
if test x$UseTramonto = xtrue ; then
  # Include the Tramonto libraries
  LIBS="${LIBS} -L../lib -ldft -laztecoo -ltriutils -lifpack \
  -lepetraext -lepetra -llapack -lblas"
fi
# extra items for UseLCAO
if test x$UseLCAO = xtrue ; then
  # Include the LCAO libraries
  LIBS="${LIBS} -L../lib"
fi
# extra items for UseInternal
if test x$UseInternal = xtrue ; then
  # Nothing happens here
  LIBS="${LIBS}"
fi
# extra items for UseSafeCompare
if test x$UseSafeCompare = xtrue ; then
  # Nothing happens here
  LIBS="${LIBS}"
fi
# extra items for UseFixGNU
if test x$UseFixGNU = xtrue ; then
  # Nothing happens here
  LIBS="$(LIBS)"
fi
# If enable-mpi, check for MPI compile scripts
if test x${UseMPI} = xtrue ;then
  # Using default of C language
  ACX_MPI
  # Set default language to Fortran 77, to make ACX_MPI look for mpif77
  AC_LANG_PUSH(Fortran 77)
  ACX_MPI
  AC_LANG_POP(Fortran 77)

  # Use the MPI compilers and libraries
  CC="${MPICC}"
  F77="${MPIF77}"
  LIBS="${MPILIBS} ${LIBS}"
  FFLAGS="${FFLAGS}"
  CFLAGS="${CFLAGS}"
fi

# Checks for libraries.
# Find the F77 libraries needed for linking
ACX_F77_LIBRARY_LDFLAGS
# Then, link using cc, not f77
AC_SUBST(F77LD,${CC})


# Checks for header files.
# Only needed if using MPI (for control.c)
if test x${UseMPI} = xtrue ;then
  AC_HEADER_STDC
  AC_CHECK_HEADERS([stdlib.h string.h])
fi

# Checks for typedefs, structures, and compiler characteristics.
# Checks for library functions.
AC_CHECK_FUNCS([strstr])

# Generate the config
# Have to list all the Makefiles
AC_CONFIG_FILES([Makefile \
Source/Makefile \
ForceFields/Makefile \
Manual/Makefile \
Manual/projects/Makefile \
Manual/projects/towhee/Makefile \
Manual/projects/towhee/algorithm/Makefile \
Manual/projects/towhee/code/Makefile \
Manual/projects/towhee/database/Makefile \
Manual/projects/towhee/examples/Makefile \
Manual/projects/towhee/forcefields/Makefile \
Manual/projects/towhee/inpstyle/Makefile \
Manual/projects/towhee/input/Makefile \
Manual/projects/towhee/potentype/Makefile \
Manual/projects/towhee/utils/Makefile \
Utils/Makefile \
Examples/Makefile \
Examples/Amber_IsoPropanol/Makefile \
Examples/Amber_Villin/Makefile \
Examples/Au_Cu_Switch/Makefile \
Examples/AVB1_Methane/Makefile \
Examples/Benzene_AA_autofit/Makefile \
Examples/Catlow_Zeolite_4a/Makefile \
Examples/Charmm19_ubiquitin/Makefile \
Examples/Charmm27_benzene/Makefile \
Examples/Charmm27_Heme/Makefile \
Examples/Charmm27_NA/Makefile \
Examples/Charmm27_Polyalanine/Makefile \
Examples/Charmm22_Ethanethiol/Makefile \
Examples/Compass_Methanol/Makefile \
Examples/Convert/Makefile \
Examples/Convert/LAMMPS_class2/Makefile \
Examples/Convert/LAMMPS_decane/Makefile \
Examples/Convert/LAMMPS_lc/Makefile \
Examples/Cu_Pb_EAM/Makefile \
Examples/Cu_VLE/Makefile \
Examples/DFT_Field/Makefile \
Examples/Dick1994_PETN/Makefile \
Examples/DREIDING/Makefile \
Examples/Dubb_Zeolite/Makefile \
Examples/Energy_Biasing/Makefile \
Examples/EPM_VLCC/Makefile \
Examples/FENE_Hexamer/Makefile \
Examples/Formamide_Scaled/Makefile \
Examples/Fris_Walls/Makefile \
Examples/Gordon/Makefile \
Examples/Gromos_Isobutane/Makefile \
Examples/Gromos_Methylpropylsulfide/Makefile \
Examples/Hard_Sphere/Makefile \
Examples/Henry_Law/Makefile \
Examples/Histogram/Makefile \
Examples/Histogram/Phase/Makefile \
Examples/Histogram/PVT/Makefile \
Examples/Histogram/Weights/Makefile \
Examples/Ideal_Chain/Makefile \
Examples/MM2_Ethane/Makefile \
Examples/NaCl_1x1x1/Makefile \
Examples/OPLS_Propanamide/Makefile \
Examples/Parallel_Test/Makefile \
Examples/Potter_CF2H2/Makefile \
Examples/Shukla_Gasses/Makefile \
Examples/SKS_Pentane/Makefile \
Examples/Small_Peptide/Makefile \
Examples/SMMKmain_2244688nonane/Makefile \
Examples/SMMKnaip_Ethylpentane/Makefile \
Examples/Solid_LJium/Makefile \
Examples/Square_Well_Chain/Makefile \
Examples/Steele_Wall/Makefile \
Examples/TraPPE_Isomers/Makefile \
Examples/TraPPE_Molecules/Makefile \
Examples/TraPPE_Pentane/Makefile \
Examples/Triglycerol/Makefile \
Examples/TMMC/Makefile \
Examples/TMMC/LJ/Makefile \
Examples/TMMC/SPC_E_WATER/Makefile \
Examples/UFF/Makefile \
Examples/uVT_Ethane/Makefile \
Examples/Vink2001_Silicon/Makefile \
Examples/VLCC_Fit/Makefile \
Examples/Wall_Water/Makefile \
Examples/Walt2001_Nanotube/Makefile \
Examples/Weiner1984/Makefile \
])
AC_OUTPUT
