function [Rndaily, Edaily]=SEBS_Daily_Evapotranspiration(ComptMask, day_angle, lat_rad, albedo, SWd24_Wm, Ta_av_K, evapfr, Ns,SEBSGPU_WORLD)
InputVarSize=size(ComptMask);
if (ndims(ComptMask)>3) 
    error('Input variables can be at most 3 dimensional')
end
DataSize    = numel(ComptMask);
DataSize_D  = gpuArray(int32(DataSize));
SEBSGPU_WORLD.SEBS_Daily_Evapotranspiration.GridSize=ceil(DataSize/SEBSGPU_WORLD.SEBS_Daily_Evapotranspiration.ThreadBlockSize(1));
% Linearizing the Input Variables
ComptMask_D    = gpuArray(single(reshape(ComptMask,DataSize,1)));
day_angle_D     = gpuArray(single(reshape(day_angle,DataSize,1)));
lat_rad_D       = gpuArray(single(reshape(lat_rad,DataSize,1)));
albedo_D        = gpuArray(single(reshape(albedo,DataSize,1)));
SWd24_Wm_D      = gpuArray(single(reshape(SWd24_Wm,DataSize,1)));
Ta_av_K_D       = gpuArray(single(reshape(Ta_av_K,DataSize,1)));
evapfr_D        = gpuArray(single(reshape(evapfr,DataSize,1)));
Ns_D          = gpuArray(single(reshape(Ns,DataSize,1)));

Rndaily_D       = gpuArray(single(zeros(DataSize,1)));
Edaily_D        = gpuArray(single(zeros(DataSize,1)));

[Rndaily_D, Edaily_D]=feval(SEBSGPU_WORLD.SEBS_Daily_Evapotranspiration, ...
                            Rndaily_D, Edaily_D, ...
                            ComptMask_D, day_angle_D, ...
                            lat_rad_D, albedo_D, ...
                            SWd24_Wm_D, Ta_av_K_D, ...
                            evapfr_D, Ns_D, DataSize_D);

if (ndims(ComptMask)==2)
    Rndaily = reshape(gather(Rndaily_D),InputVarSize(1),InputVarSize(2));
    Edaily  = reshape(gather(Edaily_D),InputVarSize(1),InputVarSize(2));
else
    Rndaily = reshape(gather(Rndaily_D),InputVarSize(1),InputVarSize(2),InputVarSize(3));
    Edaily  = reshape(gather(Edaily_D),InputVarSize(1),InputVarSize(2),InputVarSize(3));
end

Rndaily(Rndaily==-9999)=NaN;
Edaily(Edaily==-9999)=NaN;

end