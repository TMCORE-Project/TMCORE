include ./configure.SBPEM

OBJS = module_para.o      \
       module_array.o     \
       B.o                \
       vary.o             \
       mesh.o             \
       haurwitz.o         \
       L.o                \
       output_netCDF.o    \
       main.o

all: EXE

EXE:  $(OBJS)
	$(F90) -o $(EXENAME) $(OBJS) \
	-L$(NETCDF)/lib -lnetcdf -lnetcdff $(OPT) $(FCFLAGS) $(CFLAGS)

module_para.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) module_para.f90

module_array.o : module_para.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) module_array.f90

haurwitz.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) haurwitz.f90

B.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) B.f90

mesh.o : module_para.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) mesh.f90

L.o : module_para.o mesh.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) L.f90

output_netCDF.o : module_para.o mesh.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) output_netCDF.f90
	
vary.o : module_para.o module_array.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) vary.f90

main.o : module_para.o module_array.o B.o vary.o mesh.o haurwitz.o L.o output_netCDF.o 
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) main.f90


clean:
	rm *.o *.mod *.exe *.a
