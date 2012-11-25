# Compiler
TARGET=SEBS
CC=gcc

NETCDFDIR=/home/mabouali/mylib/netcdf-4.1.3
NETCDFINC=$(NETCDFDIR)/include
NETCDFLIB=$(NETCDFDIR)/lib

CFLAGS=-O3 -I$(NETCDFINC)

LINKFLAGS=-lm -L$(NETCDFLIB) -lnetcdf

all: SEBS

SEBS: AUXFunc.o SEBS_Daily_Evapotranspiration.o SEBS_EnergyBalance.o SEBS_kb_1.o SEBS_kb_1_yang.o SEBS_kb_1_yang_MOlength.o SEBS_kb_1_yang_bulktransfercoefficients.o Main.o
	$(CC) -o $(TARGET) AUXFunc.o SEBS_Daily_Evapotranspiration.o SEBS_EnergyBalance.o SEBS_kb_1.o SEBS_kb_1_yang.o SEBS_kb_1_yang_MOlength.o SEBS_kb_1_yang_bulktransfercoefficients.o Main.o $(LINKFLAGS)

Main.o: Main.c Constants.h
	$(CC) -c $(CFLAGS) Main.c

AUXFunc.o: AUXFunc.c
	$(CC) -c $(CFLAGS) AUXFunc.c

SEBS_Daily_Evapotranspiration.o: SEBS_Daily_Evapotranspiration.c  Constants.h
	$(CC) -c $(CFLAGS) SEBS_Daily_Evapotranspiration.c

SEBS_EnergyBalance.o: SEBS_EnergyBalance.c  Constants.h
	$(CC) -c $(CFLAGS) SEBS_EnergyBalance.c

SEBS_kb_1.o: SEBS_kb_1.c  Constants.h
	$(CC) -c $(CFLAGS) SEBS_kb_1.c

SEBS_kb_1_yang.o: SEBS_kb_1_yang.c  Constants.h
	$(CC) -c $(CFLAGS) SEBS_kb_1_yang.c

SEBS_kb_1_yang_MOlength.o: SEBS_kb_1_yang_MOlength.c  Constants.h
	$(CC) -c $(CFLAGS) SEBS_kb_1_yang_MOlength.c

SEBS_kb_1_yang_bulktransfercoefficients.o: SEBS_kb_1_yang_bulktransfercoefficients.c  Constants.h
	$(CC) -c $(CFLAGS) SEBS_kb_1_yang_bulktransfercoefficients.c

clean: cleanObj cleanBin cleanBackup

cleanObj:
	rm *.o

cleanBackup:
	rm *~

cleanBin:
	rm SEBS

cleanData:
	rm ./Output/*.dat

run:
	./$(TARGET)
	



