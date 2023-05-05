set(CMAKE_CXX_COMPILER "/opt/cray/pe/craype/2.7.10/bin/CC")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Clang")
set(CMAKE_CXX_COMPILER_VERSION "12.0.0")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "14")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters;cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates;cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates;cxx_std_17;cxx_std_20")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates")
set(CMAKE_CXX17_COMPILE_FEATURES "cxx_std_17")
set(CMAKE_CXX20_COMPILE_FEATURES "cxx_std_20")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_CXX_SIMULATE_VERSION "")



set(CMAKE_AR "/opt/cray/pe/cce/12.0.3/binutils/x86_64/x86_64-pc-linux-gnu/bin/ar")
set(CMAKE_CXX_COMPILER_AR "CMAKE_CXX_COMPILER_AR-NOTFOUND")
set(CMAKE_RANLIB "/opt/cray/pe/cce/12.0.3/binutils/x86_64/x86_64-pc-linux-gnu/bin/ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "CMAKE_CXX_COMPILER_RANLIB-NOTFOUND")
set(CMAKE_LINKER "/opt/cray/pe/cce/12.0.3/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;CPP)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)

foreach (lang C OBJC OBJCXX)
  if (CMAKE_${lang}_COMPILER_ID_RUN)
    foreach(extension IN LISTS CMAKE_${lang}_SOURCE_FILE_EXTENSIONS)
      list(REMOVE_ITEM CMAKE_CXX_SOURCE_FILE_EXTENSIONS ${extension})
    endforeach()
  endif()
endforeach()

set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/opt/cray/pe/petsc/3.13.3.0/real/CRAYCLANG/9.0/haswell/include;/opt/cray/pe/libsci/20.09.1/CRAYCLANG/9.0/x86_64/include;/opt/cray/pe/fftw/3.3.8.10/broadwell/include;/opt/cray/pe/mpt/7.7.18/gni/mpich-crayclang/10.0/include;/opt/cray/pe/tpsl/20.03.2/CRAYCLANG/9.0/haswell/include;/opt/cray/pe/hdf5-parallel/1.12.0.4/crayclang/10.0/include;/opt/cray/pe/netcdf-hdf5parallel/4.7.4.4/crayclang/10.0/include;/opt/cray/rca/2.2.20-7.0.3.1_3.15__g8e3fb5b.ari/include;/opt/cray/alps/6.6.67-7.0.3.1_3.18__gb91cd181.ari/include;/opt/cray/xpmem/default/include;/opt/cray/gni-headers/default/include;/opt/cray/dmapp/default/include;/opt/cray/pe/pmi/5.0.17/include;/opt/cray/ugni/default/include;/opt/cray/udreg/default/include;/opt/cray/pe/atp/3.14.5/include;/opt/cray/wlm_detect/1.3.3-7.0.3.1_3.6__g7109084.ari/include;/opt/cray/krca/2.2.8-7.0.3.1_3.14__g59af36e.ari/include;/opt/cray-hss-devel/9.0.0/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/include;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/include;/opt/cray/pe/cce/12.0.3/cce-clang/x86_64/lib/clang/12.0.0/include;/opt/cray/pe/cce/12.0.3/cce/x86_64/include/craylibs;/opt/gcc/8.1.0/snos/include/g++;/opt/gcc/8.1.0/snos/include/g++/x86_64-suse-linux;/opt/gcc/8.1.0/snos/include/g++/backward;/usr/local/include;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "AtpSigHandler;rca;hdf5_hl_parallel;hdf5_parallel;craypetsc_crayclang_real;netcdf_c++4;netcdf;sci_cray_mpi;sci_cray;fftw3_mpi;fftw3_threads;fftw3;fftw3f_mpi;fftw3f_threads;fftw3f;mpich_crayclang;mpichcxx_crayclang;pgas-dmapp;quadmath;modules;fi;craymath;f;u;csup;pthread;atomic;m;stdc++;m;craymp;-l:libunwind.so;pthread;c;-l:libunwind.so")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/petsc/3.13.3.0/real/CRAYCLANG/9.0/haswell/lib;/opt/cray/pe/libsci/20.09.1/CRAYCLANG/9.0/x86_64/lib;/opt/cray/pe/fftw/3.3.8.10/broadwell/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.18/gni/mpich-crayclang/10.0/lib;/opt/cray/pe/hdf5-parallel/1.12.0.4/crayclang/10.0/lib;/opt/cray/pe/netcdf-hdf5parallel/4.7.4.4/crayclang/10.0/lib;/opt/cray/rca/2.2.20-7.0.3.1_3.15__g8e3fb5b.ari/lib64;/opt/cray/pe/atp/3.14.5/lib;/opt/cray/pe/cce/12.0.3/cce/x86_64/lib;/opt/gcc/8.1.0/snos/lib/gcc/x86_64-suse-linux/8.1.0;/opt/gcc/8.1.0/snos/lib64;/lib64;/usr/lib64;/opt/gcc/8.1.0/snos/lib;/opt/cray/pe/cce/12.0.3/cce-clang/x86_64/lib;/lib;/usr/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/NCO/5.0.4-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/JasPer/2.0.33-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/ANTLR/2.7.7-CrayGNU-21.09-python3/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/libpng/1.6.37-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/zlib/1.2.11-CrayGNU-21.09/lib;/apps/daint/UES/jenkins/7.0.UP03/21.09/daint-mc/software/UDUNITS/2.2.28-CrayGNU-21.09/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
