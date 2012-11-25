#include <stdio.h>
#include <netcdf.h>

void checkCUDAError(const char *msg);

void CheckMemAllocationError(void *InVar, const char *InVarName);

void ReadNetCDF(float * InVar, char *InVarName, char *filePath);

void ReadNetCDF_int(int * InVar, char *InVarName, char *filePath);

void WriteOutput(float *InVar, const char *FileName, int NElem);
