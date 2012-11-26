#include "Constants.h"
#include "SEBS_kb_1.cuh"
#include "SEBS_EnergyBalance.cuh"
#include "SEBS_Daily_Evapotranspiration.cuh"
#include "AUXFunc.cuh"
#include <stdio.h>

/* NOTE: Skip to "Main Calculation Loop" to see how to use the SEBS functions. 
         Most of the code here is just to show how to read data and feed it to
         the calculation part. Depending on your program this part can change.
*/

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

	/* Input Variables on Device */
	float *fc_d, *LAI_d, *NDVI_d, *LST_K_d, *hc_d, *Zref_d, *Uref_d, *Pref_d;
	float *P0_d, *Ps_d, *Tref_K_d, *qa_ref_d, *Albedo_d, *Emissivity_d;
	float *SWd_d, *LWd_d, *hpbl_d, *SWd24_d;
	float *lat_rad_d, *Ta_av_K_d, *Ns_d;
	int *ComptMask_d;
	float *day_angle_d;
	
	/* Output Variables on Device*/
	float *d0_d, *z0h_d, *z0m_d;
	float *Rn_d, *G0_d, *H_d, *LE_d, *EF_d, *re_i_d, *ustar_d, *H_DL_d, *H_WL_d;
	float *Rndaily_d, *Edaily_d;

	/* General Variables */
	dim3 dimGrid((int)NData/nThreadsPerBlock+1); 
    dim3 dimBlock(nThreadsPerBlock);
    
    int DeviceCount;
    cudaDeviceProp DeviceProp;

	/* Timing Variables */
	cudaEvent_t start_event, stop_event;
	float ElapTime;
	/* End of Defining Variables *********************************************/

	/* Getting Some info regarding the GPU Device ****************************/
	cudaGetDeviceCount(&DeviceCount);
	if (DeviceCount==0)
	{
		printf("FATAL ERROR: There is no CUDA-Enabled Device available!!!\n");
		return(-1);
	}
	else
	{
			printf("Number of GPU device available: %d\n",DeviceCount);
	}
	if (SelectedDevice>DeviceCount-1)
	{
		printf("FATAL ERROR: The selected device does not exist.\n");
		return(-1);
	}
	cudaSetDevice(SelectedDevice);
	printf("Selected Device:\n");
	cudaGetDeviceProperties(&DeviceProp,SelectedDevice);
	printf("Device Name: %s\n",DeviceProp.name);
	printf("Compute Capability: %d.%d\n",DeviceProp.major,DeviceProp.minor);
	printf("Number of MultiProcessors: %d\n",DeviceProp.multiProcessorCount);
	/* END OF Getting Some info regarding the GPU Device *********************/

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

	printf("Allocating Memory on device ... \n");
	cudaMalloc( (void **) &fc_d, memsize);
	checkCUDAError("allocating memory for fc_d");
	cudaMalloc( (void **) &LAI_d, memsize);
	checkCUDAError("allocating memory for LAI_d");
	cudaMalloc( (void **) &NDVI_d, memsize);
	checkCUDAError("allocating memory for NDVI_d");
	cudaMalloc( (void **) &LST_K_d, memsize);
	checkCUDAError("allocating memory for LST_K_d");
	cudaMalloc( (void **) &hc_d, memsize);
	checkCUDAError("allocating memory for hc_d");
	cudaMalloc( (void **) &Zref_d, memsize);
	checkCUDAError("allocating memory for Zref_d");
	cudaMalloc( (void **) &Uref_d, memsize);
	checkCUDAError("allocating memory for Uref_d");
	cudaMalloc( (void **) &Pref_d, memsize);
	checkCUDAError("allocating memory for Pref_d");
	cudaMalloc( (void **) &P0_d, memsize);
	checkCUDAError("allocating memory for P0_d");
	cudaMalloc( (void **) &Ps_d, memsize);
	checkCUDAError("allocating memory for Ps_d");
	cudaMalloc( (void **) &Tref_K_d, memsize);
	checkCUDAError("allocating memory for Tref_K_d");
	cudaMalloc( (void **) &qa_ref_d, memsize);
	checkCUDAError("allocating memory for qa_ref_d");
	cudaMalloc( (void **) &Albedo_d, memsize);
	checkCUDAError("allocating memory for Albedo_d");
	cudaMalloc( (void **) &Emissivity_d, memsize);
	checkCUDAError("allocating memory for Emissivity_d");
	cudaMalloc( (void **) &SWd_d, memsize);
	checkCUDAError("allocating memory for SWd_d");
	cudaMalloc( (void **) &LWd_d, memsize);
	checkCUDAError("allocating memory for LWd_d");
	cudaMalloc( (void **) &hpbl_d, memsize);
	checkCUDAError("allocating memory for hpbl_d");
	cudaMalloc( (void **) &SWd24_d, memsize);
	checkCUDAError("allocating memory for SWd24_d");
	cudaMalloc( (void **) &lat_rad_d, memsize);
	checkCUDAError("allocating memory for lat_rad_d");
	cudaMalloc( (void **) &Ta_av_K_d, memsize);
	checkCUDAError("allocating memory for Ta_av_K_d");
	cudaMalloc( (void **) &Ns_d, memsize);
	checkCUDAError("allocating memory for Ns_d");
	cudaMalloc( (void **) &day_angle_d, sizeof(float));
	checkCUDAError("allocating memory for day_angle_d");
	cudaMalloc( (void **) &ComptMask_d, memsize_int);
	checkCUDAError("allocating memory for ComptMask_d");

	cudaMalloc( (void **) &d0_d, memsize);
	checkCUDAError("allocating memory for d0_d");
	cudaMalloc( (void **) &z0h_d, memsize);
	checkCUDAError("allocating memory for z0h_d");
	cudaMalloc( (void **) &z0m_d, memsize);
	checkCUDAError("allocating memory for z0m_d");
	cudaMalloc( (void **) &Rn_d, memsize);
	checkCUDAError("allocating memory for Rn_d");
	cudaMalloc( (void **) &G0_d, memsize);
	checkCUDAError("allocating memory for G0_d");
	cudaMalloc( (void **) &H_d, memsize);
	checkCUDAError("allocating memory for H_d");
	cudaMalloc( (void **) &LE_d, memsize);
	checkCUDAError("allocating memory for LE_d");
	cudaMalloc( (void **) &EF_d, memsize);
	checkCUDAError("allocating memory for EF_d");
	cudaMalloc( (void **) &re_i_d, memsize);
	checkCUDAError("allocating memory for re_i_d");
	cudaMalloc( (void **) &ustar_d, memsize);
	checkCUDAError("allocating memory for ustar_d");
	cudaMalloc( (void **) &H_DL_d, memsize);
	checkCUDAError("allocating memory for H_DL_d");
	cudaMalloc( (void **) &H_WL_d, memsize);
	checkCUDAError("allocating memory for H_WL_d");
	cudaMalloc( (void **) &Rndaily_d, memsize);
	checkCUDAError("allocating memory for Rndaily_d");
	cudaMalloc( (void **) &Edaily_d, memsize);
	checkCUDAError("allocating memory for Edaily_d");	
	printf("Done.\n");
	/* END of Allocating Memories ********************************************/
	
	/* Creating CUDA event to time the code **********************************/
	cudaEventCreate(&start_event);
	cudaEventCreate(&stop_event);
	/* END OF Creating CUDA event to time the code ***************************/

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

	/* Starting the Timer ****************************************************/
	printf("Starting the Timer ...");
	cudaEventRecord( start_event, 0 );
	checkCUDAError("Starting the Timer.");
	printf(" done.\n");
	
	/* Uploading Data to Device **********************************************/
	printf("Uploading Data to Device Memory ...\n");
	cudaMemcpy( fc_d, fc, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading fc");
	cudaMemcpy( LAI_d, LAI, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading LAI");
	cudaMemcpy( NDVI_d, NDVI, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading NDVI");
	cudaMemcpy( LST_K_d, LST_K, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading LST_K");
	cudaMemcpy( hc_d, hc, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading hc");
	cudaMemcpy( Zref_d, Zref, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Zref");
	cudaMemcpy( Uref_d, Uref, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Uref");
	cudaMemcpy( Pref_d, Pref, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Pref");
	cudaMemcpy( P0_d, P0, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading P0");
	cudaMemcpy( Ps_d, Ps, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Ps");
	cudaMemcpy( Tref_K_d, Tref_K, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Tref_K");
	cudaMemcpy( qa_ref_d, qa_ref, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading qa_ref");
	cudaMemcpy( Albedo_d, Albedo, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Albedo");
	cudaMemcpy( Emissivity_d, Emissivity, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Emissivity");
	cudaMemcpy( SWd_d, SWd, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading SWd");
	cudaMemcpy( LWd_d, LWd, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading LWd");
	cudaMemcpy( hpbl_d, hpbl, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading hpbl");
	cudaMemcpy( SWd24_d, SWd24, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading SWd24");
	cudaMemcpy( lat_rad_d, lat_rad, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading lat_rad");
	cudaMemcpy( Ta_av_K_d, Ta_av_K, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Ta_av_K");
	cudaMemcpy( Ns_d, Ns, memsize, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading Ns");
	cudaMemcpy( ComptMask_d, ComptMask, memsize_int, cudaMemcpyHostToDevice );
	checkCUDAError("Uploading ComptMask");
	cudaMemcpy( day_angle_d, &day_angle, sizeof(float), cudaMemcpyHostToDevice );
	checkCUDAError("Uploading day_angle");
	printf(" done.\n");
	/* END OF Uploading Data to Device ***************************************/
	
/*********************************************************************************************/
/*********************************** Main Calculation Loop ***********************************/
/*********************************************************************************************/

	printf("-- running SEBS_kb_1 \n");
	SEBS_kb_1_Kernel<<< dimGrid, dimBlock >>>(ComptMask_d,fc_d,LAI_d,NDVI_d,
											  LST_K_d,hc_d,Zref_d,Uref_d,
											  Pref_d,P0_d,Ps_d,Tref_K_d,
											  qa_ref_d,z0m_d,d0_d,z0h_d,NData);
	checkCUDAError("Running SEBS_kb_1_Kernel");
	cudaThreadSynchronize();
		
	printf("-- SEBS_EnergyBalance ...\n");
	SEBS_EnergyBalance_Kernel<<< dimGrid, dimBlock >>>(ComptMask_d,d0_d,z0m_d,z0h_d,
									fc_d,LAI_d,hc_d,Albedo_d,Emissivity_d,LST_K_d,NDVI_d,
									SWd_d,LWd_d,hpbl_d,Zref_d,Tref_K_d,Uref_d,qa_ref_d,
									Pref_d,Ps_d,P0_d,
									Rn_d,G0_d,H_d,LE_d,EF_d,re_i_d,ustar_d,H_DL_d,H_WL_d,NData);		
	checkCUDAError("Running SEBS_EnergyBalance_Kernel");
	cudaThreadSynchronize();
	
	printf("-- SEBS_Daily_Evapotranspiration\n");
	SEBS_Daily_Evapotranspiration_Kernel<<< dimGrid, dimBlock >>>(
							ComptMask_d,day_angle_d,lat_rad_d,
							Albedo_d,SWd24_d,Ta_av_K_d,EF_d, Ns_d,
							Rndaily_d,Edaily_d,NData); 
	checkCUDAError("Running SEBS_Daily_Evapotranspiration_Kernel");
	cudaThreadSynchronize();
	printf("-- Calculationg is done.\n");
/*********************************************************************************************/
/******************************** END OF Main Calculation Loop *******************************/
/*********************************************************************************************/

	/* Downloading Data to Device ********************************************/
	printf("-- Transfering data from Device to Host ...\n");
	cudaMemcpy( d0, d0_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading d0");
	cudaMemcpy( z0m, z0m_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading z0m");
	cudaMemcpy( z0h, z0h_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading z0h");
	cudaMemcpy( Rn, Rn_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading Rn");
	cudaMemcpy( G0, G0_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading G0");
	cudaMemcpy( H, H_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading H");
	cudaMemcpy( LE, LE_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading LE");
	cudaMemcpy( EF, EF_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading EF");
	cudaMemcpy( re_i, re_i_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading re_i");
	cudaMemcpy( ustar, ustar_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading ustar");
	cudaMemcpy( H_DL, H_DL_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading H_DL");
	cudaMemcpy( H_WL, H_WL_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading H_WL");
	cudaMemcpy( Rndaily, Rndaily_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading Rndaily");
	cudaMemcpy( Edaily, Edaily_d, memsize, cudaMemcpyDeviceToHost );
	checkCUDAError("downloading Edaily");
	/* END OF Downloading Data to Device *************************************/

	/* Stoping the Timer *****************************************************/
	cudaEventRecord( stop_event, 0 );
	cudaEventSynchronize( stop_event );
	cudaEventElapsedTime( &ElapTime, start_event, stop_event );
	printf("Total Calculation Time including memory transfer [ms]: %f\n",ElapTime);

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

	/* Freeing Memory ********************************************************/	
	cudaFree(ComptMask);
	cudaFree(fc_d);
	cudaFree(LAI_d);
	cudaFree(NDVI_d);
	cudaFree(LST_K_d);
	cudaFree(hc_d);
	cudaFree(Zref_d);
	cudaFree(Uref_d);
	cudaFree(Pref_d);
	cudaFree(P0_d);
	cudaFree(Ps_d);
	cudaFree(Tref_K_d);
	cudaFree(qa_ref_d);
	cudaFree(Albedo_d);
	cudaFree(Emissivity_d);
	cudaFree(SWd_d);
	cudaFree(LWd_d);
	cudaFree(hpbl_d);
	cudaFree(day_angle_d);
	cudaFree(lat_rad_d);
	cudaFree(SWd24_d);
	cudaFree(Ta_av_K_d);
	cudaFree(Ns_d); 

	cudaFree(d0_d);
	cudaFree(z0h_d);
	cudaFree(z0m_d);
	cudaFree(Rn_d);
	cudaFree(G0_d);
	cudaFree(H_d);
	cudaFree(LE_d);
	cudaFree(EF_d);
	cudaFree(re_i_d);
	cudaFree(ustar_d);
	cudaFree(H_DL_d);
	cudaFree(H_WL_d);
	cudaFree(Rndaily_d);
	cudaFree(Edaily_d);
	/* END OF Freeing Memory ********************************************************/	

	cudaEventDestroy(start_event);
	cudaEventDestroy(stop_event);

	printf("All done!\n");
	return(0);
}


