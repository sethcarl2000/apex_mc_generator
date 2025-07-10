
This is an example of a project which links the SIMC Fotran subroutines 'mc_hrsl' and 'mc_hrsr' to c++ code.

This will be patched into my build of G4MC, so that the SIMC routine may be directly used by the c++-code.

you can see how this works:

in the executable file, hrsl_exec.cc, we have the 'extern' statement which lets the cpp compiler 'see' the fortran subroutine; The actual subroutine is in 'simc_hrs/mc_hrsl.f'.

Also, there is an example of the opposite process; in 'simc_hrs/fortran_functs.hh/.cc', you can see how to give fortran subroutines access to cpp-written functions, albeit with some limitiations.

Check CMakeLists.txt for specific compiler flags used for compiling the Fortran77 code.

-Seth Hall; 2 July 2025

