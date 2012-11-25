# Compiler
CC=nvcc

# Executable Name
EXECNAME=SEBS

# NetCDF Library
NETCDFDIR=/home/mabouali/mylib/netcdf-4.1.3
NETCDFINC=$(NETCDFDIR)/include
NETCDFLIB=$(NETCDFDIR)/lib

# Some Compile Options
CompArch=sm_20

CFLAGS=-arch=$(CompArch) -use_fast_math $(CompPath) -I$(NETCDFINC)
LINKFLAGS=-arch=$(CompArch)  $(CompPath) -L$(NETCDFLIB) -lnetcdf


all: SEBS

SEBS: AUXFunc.o SEBS_kb_1.o SEBS_EnergyBalance.o SEBS_Daily_Evapotranspiration.o Main.o
	$(CC) -o $(EXECNAME) $(LINKFLAGS) AUXFunc.o SEBS_kb_1.o SEBS_EnergyBalance.o SEBS_Daily_Evapotranspiration.o Main.o 

Main.o: Main.cu Constants.h
	$(CC) -c -O2 $(CFLAGS) Main.cu 

AUXFunc.o: AUXFunc.cu
	$(CC) -c -O2 $(CFLAGS) AUXFunc.cu 
	
SEBS_kb_1.o: SEBS_kb_1.cu Constants.h
	$(CC) -c -O2 $(CFLAGS) SEBS_kb_1.cu

SEBS_EnergyBalance.o: SEBS_EnergyBalance.cu Constants.h
	$(CC) -c -O2 $(CFLAGS) SEBS_EnergyBalance.cu

SEBS_Daily_Evapotranspiration.o: SEBS_Daily_Evapotranspiration.cu Constants.h
	$(CC) -c -O2 $(CFLAGS) SEBS_Daily_Evapotranspiration.cu

clean: cleanObj cleanBin cleanBackup

cleanObj:
	rm *.o

cleanBackup:
	rm *~

cleanBin:
	rm SEBS

cleanData:
	rm ./Output/*.dat



