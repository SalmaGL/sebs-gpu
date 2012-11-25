#include "AUXFunc.h"

float max(float a, float b)
{
	return( a>=b? a : b );
}

float min(float a, float b)
{
	return ( a>=b? b : a );
}

void CheckMemAllocationError(void *InVar, const char *InVarName)
{
	if (InVar==NULL)
	{
		printf("Error Allocating Memory: %s\n",InVarName);
		exit(-1);
	}
}

void ReadNetCDF(float * InVar, char *InVarName, char *filePath)
{
	int ncid,varid;
	if (nc_open(filePath,NC_NOWRITE,&ncid)!=NC_NOERR)
	{
		printf("\n Error: Can not open %s.nc.\n",InVarName);
		exit(-1);
	}
	else
	{
		if ( nc_inq_varid(ncid,InVarName,&varid)!=NC_NOERR )
		{
			printf("\n Error: Can not find %s variable in the file.\n",InVarName);
			exit(-1);
		}
		else
		{
			printf("Reading %s:\n",InVarName);
			if(nc_get_var_float(ncid,varid,InVar)!=NC_NOERR)
			{
				printf("\n Fatal Error: Can not read %s.\n",InVarName);
				exit(-1);
			}			
		}
	}		
}

void ReadNetCDF_int(int * InVar, char *InVarName, char *filePath)
{
	int ncid,varid;
	if (nc_open(filePath,NC_NOWRITE,&ncid)!=NC_NOERR)
	{
		printf("\n Error: Can not open %s.nc.\n",InVarName);
		exit(-1);
	}
	else
	{
		if ( nc_inq_varid(ncid,InVarName,&varid)!=NC_NOERR )
		{
			printf("\n Error: Can not find %s variable in the file.\n",InVarName);
			exit(-1);
		}
		else
		{
			printf("Reading %s:\n",InVarName);
			if(nc_get_var_int(ncid,varid,InVar)!=NC_NOERR)
			{
				printf("\n Fatal Error: Can not read %s.\n",InVarName);
				exit(-1);
			}			
		}
	}		
}

void WriteOutput(float *InVar, const char *FileName, int NElem)
{
	FILE * fid;
	int i;
	printf("Writing %s\n",FileName);
	if ((fid=fopen(FileName,"w"))==NULL)
	{
		printf("FATAL ERROR: Can't open %s\n",FileName);
	}
	else
	{
		for(i=0;i<NElem;i++)
		{
			fprintf(fid,"%f\n",InVar[i]);
		}
		fclose(fid);
	}
}


