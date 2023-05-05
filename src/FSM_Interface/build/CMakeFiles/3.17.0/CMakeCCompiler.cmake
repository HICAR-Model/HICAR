set(CMAKE_C_COMPILER "/opt/cray/pe/craype/2.7.10/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Clang")
set(CMAKE_C_COMPILER_VERSION "12.0.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/opt/cray/pe/cce/12.0.3/binutils/x86_64/x86_64-pc-linux-gnu/bin/ar")
set(CMAKE_C_COMPILER_AR "CMAKE_C_COMPILER_AR-NOTFOUND")
set(CMAKE_RANLIB "/opt/cray/pe/cce/12.0.3/binutils/x86_64/x86_64-pc-linux-gnu/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "CMAKE_C_COMPILER_RANLIB-NOTFOUND")
set(CMAKE_LINKER "/opt/cray/pe/cce/12.0.3/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/opt/cray/pe/petsc/3.13.3.0/real/CRAYCLANG/9.0/haswell/include;/opt/cray/pe/libsci/20.09.1/CRAYCLANG/9.0/x86_64/include;/opt/cray/pe/fftw/3.3.8.10/broadwell/include;/opt/cray/pe/mpt/7.7.18/gni/mpich-crayclang/10.0/include;/opt/cray/pe/tpsl/20.03.2/CRAYCLANG/9.0/haswell/include;/opt/cray/pe/hdf5-parallel/1.12.0.4/crayclang/10.0/include;/opt/cray/pe/netcdf-hdf5parallel/4.7.4.4/crayclang/10.0/include;/opt/cray/rca/2.2.20-7.0.3.1_3.15__g8e3fb5b.ari/include;/opt/cray/alps/6.6.67-7.0.3.1_3.18__gb91cd181.ari/include;/opt/cray/xpmem/default/include;/opt/cray/gni-headers/default/include;/opt/cray/dmapp/default/include;/opt/cray/pe/pmi/5.0.17/include;/opt/cray/ugni/default/include;/opt/cray/udreg/default/include;/opt/cray/pe/atp/3.14.5/include;/opt/cray/wlm_detect/1.3.3-7.0.3.1_3.6__g7109084.ari/include;/opt/cray/krca/2.2.8-7.0.3.1_3.14__g59af36e.ari/include;/opt/cray-hss-devel/9.0.0/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/include;/opt/cray/pe/cce/12.0.3/cce-clang/x86_64/lib/clang/12.0.0/include;/opt/cray/pe/cce/12.0.3/cce/x86_64/include/craylibs;/usr/local/include;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "AtpSigHandler;rca;hdf5_hl_parallel;hdf5_parallel;craypetsc_crayclang_real;netcdf;sci_cray_mpi;sci_cray;fftw3_mpi;fftw3_threads;fftw3;fftw3f_mpi;fftw3f_threads;fftw3f;mpich_crayclang;pgas-dmapp;quadmath;modules;fi;craymath;f;u;csup;pthread;atomic;m;craymp;-l:libunwind.so;pthread;c;-l:libunwind.so")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/petsc/3.13.3.0/real/CRAYCLANG/9.0/haswell/lib;/opt/cray/pe/libsci/20.09.1/CRAYCLANG/9.0/x86_64/lib;/opt/cray/pe/fftw/3.3.8.10/broadwell/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.18/gni/mpich-crayclang/10.0/lib;/opt/cray/pe/hdf5-parallel/1.12.0.4/crayclang/10.0/lib;/opt/cray/pe/netcdf-hdf5parallel/4.7.4.4/crayclang/10.0/lib;/opt/cray/rca/2.2.20-7.0.3.1_3.15__g8e3fb5b.ari/lib64;/opt/cray/pe/atp/3.14.5/lib;/opt/cray/pe/cce/12.0.3/cce/x86_64/lib;/opt/gcc/8.1.0/snos/lib/gcc/x86_64-suse-linux/8.1.0;/opt/gcc/8.1.0/snos/lib64;/lib64;/usr/lib64;/opt/gcc/8.1.0/snos/lib;/opt/cray/pe/cce/12.0.3/cce-clang/x86_64/lib;/lib;/usr/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/NCO/5.0.4-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
