%% Inivtializing the SEBS_GPU
disp('Initializing SEBS GPU Modules and Kernels ...')
SEBSGPU_WORLD=InitSEBSGPU_WORLD();

%% Loading Input Data Set
% Depending on how you have stored your data, you must modify this part of
% the code.

disp('Loading sample data set ...')
albedo=netcdf.getVar(netcdf.open('../InputData/sAlbedo.nc','NC_NOWRITE'),0);
ComptMask=netcdf.getVar(netcdf.open('../InputData/sComptMask.nc','NC_NOWRITE'),0);
emissivity=netcdf.getVar(netcdf.open('../InputData/sEmissivity.nc','NC_NOWRITE'),0);
fc=netcdf.getVar(netcdf.open('../InputData/sfc.nc','NC_NOWRITE'),0);
hc=netcdf.getVar(netcdf.open('../InputData/shc.nc','NC_NOWRITE'),0);
hpbl=netcdf.getVar(netcdf.open('../InputData/shpbl.nc','NC_NOWRITE'),0);
LAI=netcdf.getVar(netcdf.open('../InputData/sLAI.nc','NC_NOWRITE'),0);
lat=netcdf.getVar(netcdf.open('../InputData/slat.nc','NC_NOWRITE'),0);
lon=netcdf.getVar(netcdf.open('../InputData/slon.nc','NC_NOWRITE'),0);
LWd=netcdf.getVar(netcdf.open('../InputData/sLWd.nc','NC_NOWRITE'),0);
NDVI=netcdf.getVar(netcdf.open('../InputData/sNDVI.nc','NC_NOWRITE'),0);
Pref=netcdf.getVar(netcdf.open('../InputData/sPressure.nc','NC_NOWRITE'),0);
P0=netcdf.getVar(netcdf.open('../InputData/sprmsl.nc','NC_NOWRITE'),0);
qa_ref=netcdf.getVar(netcdf.open('../InputData/sqa_ref.nc','NC_NOWRITE'),0);
SWd24=netcdf.getVar(netcdf.open('../InputData/sSWd24.nc','NC_NOWRITE'),0);
SWd=netcdf.getVar(netcdf.open('../InputData/sSWd.nc','NC_NOWRITE'),0);
Ta_av_K=netcdf.getVar(netcdf.open('../InputData/sTa_avg.nc','NC_NOWRITE'),0);
Tref_K=netcdf.getVar(netcdf.open('../InputData/sTa_ref.nc','NC_NOWRITE'),0);
Ns=netcdf.getVar(netcdf.open('../InputData/sNs.nc','NC_NOWRITE'),0);
LST_K=netcdf.getVar(netcdf.open('../InputData/sTs.nc','NC_NOWRITE'),0);
Uref=netcdf.getVar(netcdf.open('../InputData/sUref.nc','NC_NOWRITE'),0);
Zref=netcdf.getVar(netcdf.open('../InputData/Zref.nc','NC_NOWRITE'),0);

Ps=Pref;
day_angle=ones(120,120)*194;
%% Calling SEBS

disp('Applying SEBS Algorithm ...')
% tic
% [z0m, d0, z0h]=SEBS_kb_1(...
%                 fc, LAI, NDVI, LST_K, ...
%                 hc, Zref, Uref, Pref, ...
%                 P0, Ps, Tref_K, qa_ref,...
%                 ComptMask, SEBSGPU_WORLD);
% 
% [Rn, G0, H_i, LE, ...
%  evap_fr, re_i, ustar, ...
%  H_DL, H_WL]=SEBS_EnergyBalance(ComptMask, ...
%                 d0, z0m, z0h, fc, ...
%                 LAI, hc, albedo, emissivity, ...
%                 LST_K, NDVI, SWd, LWd, ...
%                 hpbl, Zref, Tref_K, ...
%                 Uref, qa_ref, Pref, Ps, P0, ...
%                 SEBSGPU_WORLD);
% 
% [Rndaily, Edaily]=  ...
%     SEBS_Daily_Evapotranspiration(ComptMask, day_angle, lat, ...
%                 albedo, SWd24, Ta_av_K, evap_fr, Ns,SEBSGPU_WORLD);
% 
% toc
% 
% disp('Plotting ...')
% close
% surface(lon*180.0/pi,lat*180/pi,double(Edaily),'EdgeColor','none');
% xlim([min(min(lon*180.0/pi)) max(max(lon*180.0/pi))]);
% ylim([min(min(lat*180.0/pi)) max(max(lat*180.0/pi))]);
% title('Edaily');
   

tic
[z0m, d0, z0h, ...
    Rn, G0, H_i, LE, evap_fr, re_i, ustar, H_DL, H_WL, ...
    Rndaily, Edaily2]=SEBS (ComptMask, ...
    fc, LAI, NDVI, LST_K, hc, Zref, Uref, Pref, P0, Ps, Tref_K, qa_ref, ...
    albedo, emissivity, SWd, LWd, hpbl, ...
    day_angle, lat, SWd24, Ta_av_K, Ns,SEBSGPU_WORLD);
toc
%%
% figure
% disp('Plotting ...')
% close
% surface(lon*180.0/pi,lat*180/pi,double(Edaily2),'EdgeColor','none');
% xlim([min(min(lon*180.0/pi)) max(max(lon*180.0/pi))]);
% ylim([min(min(lat*180.0/pi)) max(max(lat*180.0/pi))]);
% title('Edaily2');
