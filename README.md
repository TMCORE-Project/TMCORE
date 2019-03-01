# tmcore
TRiSK-based Multiple-Conservation dynamical cORE

This model is a TRiSK-based shallow water dynamical core, which has 4 integral invariants including total mass, total absolute vorticity,

total energy and total potential enstrophy.

Here we provide 4 kinds of test cases, 

config_test_case = 2 ! Global Steady State Nolinear Zonal Geostrophic Flow, Williamson 1992, test case 2

config_test_case = 5 ! Zonal Flow Over an Isolated Moutain, Williamson 1992, test case 5

config_test_case = 6 ! Rossby Haurwitz Wave with 4 wavenumber, Williamson 1992, test case 6

config_test_case = 7 ! Linearized Rossby wave, Shamir et. al 2016

To obtain the energy conservation by setting namelist

config_energy_conservation = .true.

To obtain the potential enstrophy conservation by setting namelist

config_PV_scheme = 'APVM_Conservation'

To compile this model, you need netCDF Fortran as the I/O libaray, specifying the path to netCDF Fortran root directory by

export NETCDF_ROOT=$WHERE_YOUR_NETCDF_FORTRAN_IS

Next, configure the makefile info by cmake

cd tmcore

export tmcore=`pwd`

cd $WHERE_TO_RUN_TMCORE

cmake $tmcore

make

after successful compiling the model, tmcore_test.exe will appear, putting /run/namelist.sw in the same directory

as tmcore_test.exe, and modifying config_mesh_file in namelist to point at the mesh file, then

./tmcore_test.exe
