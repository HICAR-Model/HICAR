set(CMAKE_Fortran_COMPILER "/opt/cray/pe/craype/2.7.10/bin/ftn")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "9.3.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/cray/pe/petsc/3.13.3.0/real/GNU/8.2/haswell/include;/opt/cray/pe/libsci/20.09.1/GNU/8.1/x86_64/include;/opt/cray/pe/fftw/3.3.8.10/broadwell/include;/opt/cray/pe/mpt/7.7.18/gni/mpich-gnu/8.2/include;/opt/cray/pe/tpsl/20.03.2/GNU/8.2/haswell/include;/opt/cray/pe/hdf5-parallel/1.12.0.4/gnu/8.2/include;/opt/cray/pe/netcdf-hdf5parallel/4.7.4.4/gnu/8.2/include;/opt/cray/rca/2.2.20-7.0.3.1_3.15__g8e3fb5b.ari/include;/opt/cray/alps/6.6.67-7.0.3.1_3.18__gb91cd181.ari/include;/opt/cray/xpmem/default/include;/opt/cray/gni-headers/default/include;/opt/cray/pe/pmi/5.0.17/include;/opt/cray/ugni/default/include;/opt/cray/udreg/default/include;/opt/cray/pe/atp/3.14.5/include;/opt/cray/wlm_detect/1.3.3-7.0.3.1_3.6__g7109084.ari/include;/opt/cray/krca/2.2.8-7.0.3.1_3.14__g59af36e.ari/include;/opt/cray-hss-devel/9.0.0/include;/opt/gcc/9.3.0/snos/lib/gcc/x86_64-suse-linux/9.3.0/finclude;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/include;/opt/gcc/9.3.0/snos/lib/gcc/x86_64-suse-linux/9.3.0/include;/usr/local/include;/opt/gcc/9.3.0/snos/include;/opt/gcc/9.3.0/snos/lib/gcc/x86_64-suse-linux/9.3.0/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "AtpSigHandler;rca;hdf5_hl_parallel;hdf5_parallel;hdf5hl_fortran_parallel;hdf5_fortran_parallel;craypetsc_gnu_82_real;netcdf;netcdff;sci_gnu_82_mpi;sci_gnu_82;fftw3_mpi;fftw3_threads;fftw3;fftw3f_mpi;fftw3f_threads;fftw3f;mpich_gnu_82;mpichf90_gnu_82;gfortran;quadmath;pthread;gfortran;m;gomp;gcc_s;gcc;quadmath;m;gcc_s;gcc;pthread;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/petsc/3.13.3.0/real/GNU/8.2/haswell/lib;/opt/cray/pe/libsci/20.09.1/GNU/8.1/x86_64/lib;/opt/cray/pe/fftw/3.3.8.10/broadwell/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.18/gni/mpich-gnu/8.2/lib;/opt/cray/pe/hdf5-parallel/1.12.0.4/gnu/8.2/lib;/opt/cray/pe/netcdf-hdf5parallel/4.7.4.4/gnu/8.2/lib;/opt/cray/rca/2.2.20-7.0.3.1_3.15__g8e3fb5b.ari/lib64;/opt/cray/pe/atp/3.14.5/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/NCO/5.0.4-CrayGNU-21.09/lib64;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/lib64;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/lib64;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/lib64;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/lib64;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/lib64;/opt/gcc/9.3.0/snos/lib/gcc/x86_64-suse-linux/9.3.0;/opt/gcc/9.3.0/snos/lib64;/lib64;/usr/lib64;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/NCO/5.0.4-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/lib;/opt/gcc/9.3.0/snos/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
