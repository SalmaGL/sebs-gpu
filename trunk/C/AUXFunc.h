#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

float max(float a, float b);

float min(float a, float b);

void CheckMemAllocationError(void *InVar, const char *InVarName);

void ReadNetCDF(float * InVar, char *InVarName, char *filePath);

void ReadNetCDF_int(int * InVar, char *InVarName, char *filePath);

void WriteOutput(float *InVar, const char *FileName, int NElem);
