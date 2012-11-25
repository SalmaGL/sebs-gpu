# Compiler
CC=nvcc

#CompPath=-ccbin /home/abouali/MyCompiler
CompArch=sm_13

CFLAGS=-arch=$(CompArch) $(CompPath)

LINKFLAGS=-arch=$(CompArch)  $(CompPath)

all: SEBS

SEBS: SEBS_kb_1.ptx SEBS_EnergyBalance.ptx SEBS_Daily_Evapotranspiration.ptx
	cp -v *ptx ../

SEBS_kb_1.ptx: SEBS_kb_1.cu Constants.h
	$(CC) -ptx -O2 $(CFLAGS) SEBS_kb_1.cu

SEBS_EnergyBalance.ptx: SEBS_EnergyBalance.cu Constants.h
	$(CC) -ptx -O2 $(CFLAGS) SEBS_EnergyBalance.cu

SEBS_Daily_Evapotranspiration.ptx: SEBS_Daily_Evapotranspiration.cu Constants.h
	$(CC) -ptx -O2 $(CFLAGS) SEBS_Daily_Evapotranspiration.cu

clean: cleanObj cleanBackup

cleanObj:
	rm *.ptx

cleanBackup:
	rm *~



