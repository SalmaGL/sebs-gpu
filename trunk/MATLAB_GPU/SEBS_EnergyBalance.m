function [Rn, G0, H_i, LE, evap_fr, re_i, ustar, H_DL, H_WL]=SEBS_EnergyBalance(ComptMask, d0, z0m, z0h, fc, LAI, hc, albedo, emissivity, LST_K, NDVI, SWd, LWd, hpbl, Zref, Tref_K, Uref, qaref, Pref, Ps, P0, SEBSGPU_WORLD)
InputVarSize=size(ComptMask);
if (ndims(ComptMask)>3) 
    error('Input variables can be at most 3 dimensional')
end
DataSize    = numel(ComptMask);
DataSize_D  = gpuArray(int32(DataSize));
SEBSGPU_WORLD.SEBS_EnergyBalance.GridSize=ceil(DataSize/SEBSGPU_WORLD.SEBS_EnergyBalance.ThreadBlockSize(1));
% Linearizing the Input Variables
ComptMask_D 	= gpuArray(single(reshape(ComptMask,DataSize,1))); 
d0_D 			= gpuArray(single(reshape(d0,DataSize,1))); 
z0m_D 			= gpuArray(single(reshape(z0m,DataSize,1))); 
z0h_D 			= gpuArray(single(reshape(z0h,DataSize,1))); 
fc_D 			= gpuArray(single(reshape(fc,DataSize,1))); 
LAI_D 			= gpuArray(single(reshape(LAI,DataSize,1))); 
hc_D 			= gpuArray(single(reshape(hc,DataSize,1))); 
albedo_D 		= gpuArray(single(reshape(albedo,DataSize,1))); 
emissivity_D 	= gpuArray(single(reshape(emissivity,DataSize,1))); 
LST_K_D 		= gpuArray(single(reshape(LST_K,DataSize,1))); 
NDVI_D 			= gpuArray(single(reshape(NDVI,DataSize,1))); 
SWd_D 			= gpuArray(single(reshape(SWd,DataSize,1))); 
LWd_D 			= gpuArray(single(reshape(LWd,DataSize,1))); 
hpbl_D 			= gpuArray(single(reshape(hpbl,DataSize,1))); 
Zref_D 			= gpuArray(single(reshape(Zref,DataSize,1))); 
Tref_K_D 		= gpuArray(single(reshape(Tref_K,DataSize,1))); 
Uref_D 			= gpuArray(single(reshape(Uref,DataSize,1))); 
qaref_D 		= gpuArray(single(reshape(qaref,DataSize,1))); 
Pref_D 			= gpuArray(single(reshape(Pref,DataSize,1))); 
Ps_D 			= gpuArray(single(reshape(Ps,DataSize,1))); 
P0_D 			= gpuArray(single(reshape(P0,DataSize,1))); 

Rn_D        = gpuArray(single(zeros(DataSize,1)));
G0_D		= gpuArray(single(zeros(DataSize,1)));
H_i_D		= gpuArray(single(zeros(DataSize,1)));
LE_D		= gpuArray(single(zeros(DataSize,1)));
evap_fr_D	= gpuArray(single(zeros(DataSize,1)));
re_i_D      = gpuArray(single(zeros(DataSize,1)));
ustar_D     = gpuArray(single(zeros(DataSize,1)));
H_DL_D      = gpuArray(single(zeros(DataSize,1)));
H_WL_D      = gpuArray(single(zeros(DataSize,1)));

[Rn_D, G0_D, H_i_D, LE_D, ...
 evap_fr_D, re_i_D, ustar_D, ...
 H_DL_D, H_WL_D]= feval(SEBSGPU_WORLD.SEBS_EnergyBalance, ...
                    Rn_D, G0_D, H_i_D, LE_D, ...
                    evap_fr_D, re_i_D, ustar_D, ...
                    H_DL_D, H_WL_D, ...
                    ComptMask_D, d0_D, z0m_D, z0h_D, ...
                    fc_D, LAI_D, hc_D, albedo_D, ...
                    emissivity_D, LST_K_D, NDVI_D, ...
                    SWd_D, LWd_D, hpbl_D, Zref_D, ...
                    Tref_K_D, Uref_D, qaref_D, Pref_D, ...
                    Ps_D, P0_D, DataSize_D);

if (ndims(ComptMask)==2)
    Rn=reshape(gather(Rn_D),InputVarSize(1),InputVarSize(2)); 
    G0=reshape(gather(G0_D),InputVarSize(1),InputVarSize(2)); 
    H_i=reshape(gather(H_i_D),InputVarSize(1),InputVarSize(2)); 
    LE=reshape(gather(LE_D),InputVarSize(1),InputVarSize(2)); 
    evap_fr=reshape(gather(evap_fr_D),InputVarSize(1),InputVarSize(2)); 
    re_i=reshape(gather(re_i_D),InputVarSize(1),InputVarSize(2)); 
    ustar=reshape(gather(ustar_D),InputVarSize(1),InputVarSize(2)); 
    H_DL=reshape(gather(H_DL_D),InputVarSize(1),InputVarSize(2)); 
    H_WL=reshape(gather(H_WL_D),InputVarSize(1),InputVarSize(2)); 
else
    Rn=reshape(gather(Rn_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    G0=reshape(gather(G0_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    H_i=reshape(gather(H_i_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    LE=reshape(gather(LE_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    evap_fr=reshape(gather(evap_fr_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    re_i=reshape(gather(re_i_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    ustar=reshape(gather(ustar_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    H_DL=reshape(gather(H_DL_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
    H_WL=reshape(gather(H_WL_D),InputVarSize(1),InputVarSize(2),InputVarSize(3)); 
end

Rn(Rn==-9999)=NaN;
G0(G0==-9999)=NaN;
H_i(H_i==-9999)=NaN;
LE(LE==-9999)=NaN;
evap_fr(evap_fr==-9999)=NaN;
re_i(re_i==-9999)=NaN;
ustar(ustar==-9999)=NaN;
H_DL(H_DL==-9999)=NaN;
H_WL(H_WL==-9999)=NaN;

end



