# Introduction

TMCORE stands for TRiSK-based Multiple-Conservation dynamical cORE. This dynamical core (shortly as dycore) has 4 integral invariants:

- total mass
- total absolute vorticity
- total energy
- total potential enstrophy

We are currently developing it with passion, and have set the target to delivery a simpler to understand and robust ensured version.

# Usage

Check environment for NetCDF library:

```
$ which ncdump
```

If it returns the path to `ncdump` command, set `NETCDF` or `NETCDF_ROOT` to the upper directory (e.g. `/usr`). Then change directory to TMCORE and use CMake to build:

```
$ cd <tmcore_root>/build
$ FC=<Fortran_compiler> cmake ..
$ make
```

Create or modify a `namelist.tmcore_sw`, then run `test_tmcore_sw.exe`