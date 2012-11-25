#include "SEBS_kb_1.h"
#include "SEBS_EnergyBalance.h"
#include "SEBS_Daily_Evapotranspiration.h"
#include "AUXFunc.h"
#include <time.h>
#include <stdio.h>

int main(void)
{
	/* Defining Variables *****************************************************/
	/* Input Variables */
	int NData=dimx*dimy;
	unsigned long int memsize=NData*sizeof(float);
	unsigned long int memsize_int=NData*sizeof(int);
	float *fc, *LAI, *NDVI, *LST_K, *hc, *Zref, *Uref, *Pref;
	float *P0, *Ps, *Tref_K, *qa_ref, *Albedo, *Emissivity;
	float *SWd, *LWd, *hpbl, *SWd24;
	float *lat_rad, *Ta_av_K, *Ns;
	int *ComptMask;
	float day_angle=194.0;
	
	/* Output Variables */
	float *d0, *z0h, *z0m;
	float *Rn, *G0, *H, *LE, *EF, *re_i, *ustar, *H_DL, *H_WL;
	float *Rndaily, *Edaily;

	/* General Variables */
	int i,j;
	clock_t stime,etime;	
	/* End of Defining Variables *********************************************/

	/* Allocating Memories ***************************************************/
	printf("Allocating Memory on host ... \n");
	CheckMemAllocationError(fc=(float *) malloc(memsize),"fc");
	CheckMemAllocationError(LAI=(float *) malloc(memsize),"LAI");
	CheckMemAllocationError(NDVI=(float *) malloc(memsize),"NDVI");
	CheckMemAllocationError(LST_K=(float *) malloc(memsize),"LST_K");
	CheckMemAllocationError(hc=(float *) malloc(memsize),"hc");
	CheckMemAllocationError(Zref=(float *) malloc(memsize),"Zref");
	CheckMemAllocationError(Uref=(float *) malloc(memsize),"Uref");
	CheckMemAllocationError(Pref=(float *) malloc(memsize),"Pref");
	CheckMemAllocationError(P0=(float *) malloc(memsize),"P0");
	CheckMemAllocationError(Ps=(float *) malloc(memsize),"Ps");
	CheckMemAllocationError(Tref_K=(float *) malloc(memsize),"Tref_K");
	CheckMemAllocationError(qa_ref=(float *) malloc(memsize),"qa_ref");
	CheckMemAllocationError(Albedo=(float *) malloc(memsize),"Albedo");
	CheckMemAllocationError(Emissivity=(float *) malloc(memsize),"Emissivity");
	CheckMemAllocationError(SWd=(float *) malloc(memsize),"SWd");
	CheckMemAllocationError(LWd=(float *) malloc(memsize),"LWd");
	CheckMemAllocationError(hpbl=(float *) malloc(memsize),"hpbl");
	CheckMemAllocationError(SWd24=(float *) malloc(memsize),"SWd24");
	CheckMemAllocationError(lat_rad=(float *) malloc(memsize),"lat_rad");
	CheckMemAllocationError(Ta_av_K=(float *) malloc(memsize),"Ta_av_K");
	CheckMemAllocationError(Ns=(float *) malloc(memsize),"Ns");
	CheckMemAllocationError(ComptMask=(int *) malloc(memsize_int),"ComptMask");

	CheckMemAllocationError(d0=(float *) malloc(memsize),"d0");
	CheckMemAllocationError(z0h=(float *) malloc(memsize),"z0h");
	CheckMemAllocationError(z0m=(float *) malloc(memsize),"z0m");
	CheckMemAllocationError(Rn=(float *) malloc(memsize),"Rn");
	CheckMemAllocationError(G0=(float *) malloc(memsize),"G0");
	CheckMemAllocationError(H=(float *) malloc(memsize),"H");
	CheckMemAllocationError(LE=(float *) malloc(memsize),"LE");
	CheckMemAllocationError(EF=(float *) malloc(memsize),"EF");
	CheckMemAllocationError(re_i=(float *) malloc(memsize),"re_i");
	CheckMemAllocationError(ustar=(float *) malloc(memsize),"ustar");
	CheckMemAllocationError(H_DL=(float *) malloc(memsize),"H_DL");
	CheckMemAllocationError(H_WL=(float *) malloc(memsize),"H_WL");
	CheckMemAllocationError(Rndaily=(float *) malloc(memsize),"Rndaily");
	CheckMemAllocationError(Edaily=(float *) malloc(memsize),"Edaily");
	printf("Done.\n");
	/* END of Allocating Memories ********************************************/
	
	/* Loading the input Variables *******************************************/
	ReadNetCDF(fc, "fc", "../InputData/sfc.nc");
	ReadNetCDF(LAI, "LAI", "../InputData/sLAI.nc");
	ReadNetCDF(NDVI, "NDVI", "../InputData/sNDVI.nc");
	ReadNetCDF(LST_K, "Ts", "../InputData/sTs.nc");
	ReadNetCDF(hc, "hc", "../InputData/shc.nc");
	ReadNetCDF(Zref, "Zref", "../InputData/Zref.nc");
	ReadNetCDF(Uref, "Uref", "../InputData/sUref.nc");
	ReadNetCDF(Pref, "Pressure", "../InputData/sPressure.nc");
	ReadNetCDF(P0, "prmsl", "../InputData/sprmsl.nc");
	ReadNetCDF(Ps, "Pressure", "../InputData/sPressure.nc");
	ReadNetCDF(Tref_K, "Ta_ref", "../InputData/sTa_ref.nc");
	ReadNetCDF(qa_ref, "qa_ref", "../InputData/sqa_ref.nc");
	ReadNetCDF(Albedo, "Albedo", "../InputData/sAlbedo.nc");
	ReadNetCDF(Emissivity, "Emissivity", "../InputData/sEmissivity.nc");
	ReadNetCDF(SWd, "SWd", "../InputData/sSWd.nc");
	ReadNetCDF(LWd, "LWd", "../InputData/sLWd.nc");
	ReadNetCDF(hpbl, "hpbl", "../InputData/shpbl.nc");
	ReadNetCDF(SWd24, "SWd24", "../InputData/sSWd24.nc");
	ReadNetCDF(lat_rad, "lat", "../InputData/slat.nc");
	ReadNetCDF(Ta_av_K, "Ta_avg", "../InputData/sTa_avg.nc");
	ReadNetCDF(Ns, "Ns", "../InputData/sNs.nc");
	ReadNetCDF_int(ComptMask, "ComptMask", "../InputData/sComptMask.nc");
	/* END OF Loading the input Variables ************************************/	
	
	stime=clock();
	for(i=0;i<NData;i++)
	{
		if (ComptMask[i]==1.0)
		{
			SEBS_kb_1(fc[i],LAI[i],NDVI[i],LST_K[i],
					  hc[i],Zref[i],Uref[i],Pref[i],
					  P0[i],Ps[i],Tref_K[i],qa_ref[i],
					  (z0m+i),(d0+i),(z0h+i) );
	
			SEBS_EnergyBalance(d0[i],z0m[i],z0h[i],
			           fc[i],LAI[i],hc[i],Albedo[i],Emissivity[i],LST_K[i],
			           NDVI[i],SWd[i],LWd[i],hpbl[i],Zref[i],Tref_K[i],
			           Uref[i],qa_ref[i],Pref[i],Ps[i],P0[i],				
			           (Rn+i),(G0+i),(H+i),(LE+i),(EF+i),(re_i+i),(ustar+i),(H_DL+i),(H_WL+i) );
			SEBS_Daily_Evapotranspiration(day_angle,
					    lat_rad[i],Albedo[i],SWd24[i],Ta_av_K[i],EF[i], Ns[i],
					    (Rndaily+i),(Edaily+i) );
		}
		else
		{
			z0m[i]=-9999.0;
			d0[i]=-9999.0;
			z0h[i]=-9999.0;
			Rn[i]=-9999.0;
			G0[i]=-9999.0;
			H[i]=-9999.0;
			LE[i]=-9999.0;
			EF[i]=-9999.0;
			re_i[i]=-9999.0;
			ustar[i]=-9999.0;
			H_DL[i]=-9999.0;
			H_WL[i]=-9999.0;
			Rndaily[i]=-9999.0;
			Edaily[i]=-9999.0;
		}
	}
	etime=clock();
	printf("Calculation time: %f\n",(double) (etime-stime)/ (double) CLOCKS_PER_SEC);
	
	/* Writting the Output ***************************************************/
	WriteOutput(z0m,"./Output/z0m.dat",NData);
	WriteOutput(z0h,"./Output/z0h.dat",NData);
	WriteOutput(d0,"./Output/d0.dat",NData);
	WriteOutput(Rn,"./Output/Rn.dat",NData);
	WriteOutput(G0,"./Output/G0.dat",NData);
	WriteOutput(H,"./Output/H.dat",NData);
	WriteOutput(LE,"./Output/LE.dat",NData);
	WriteOutput(EF,"./Output/EF.dat",NData);
	WriteOutput(re_i,"./Output/re_i.dat",NData);
	WriteOutput(ustar,"./Output/ustar.dat",NData);
	WriteOutput(H_DL,"./Output/H_DL.dat",NData);
	WriteOutput(H_WL,"./Output/H_WL.dat",NData);
	WriteOutput(Rndaily,"./Output/Rndaily.dat",NData);
	WriteOutput(Edaily,"./Output/Edaily.dat",NData);
	/* END OF Writting the Output ********************************************/

	printf("All done!\n");
	return(0);
}

