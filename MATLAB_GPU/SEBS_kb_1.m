function [z0m, d0, z0h]=SEBS_kb_1(ComptMask, fc, LAI, NDVI, LST_K, hc, Zref, Uref, Pref, P0, Ps, Tref_K, qa_ref,SEBSGPU_WORLD)

InputVarSize=size(ComptMask);
if (ndims(ComptMask)>3) 
    error('Input variables can be at most 3 dimensional')
end
DataSize    = numel(ComptMask);
DataSize_D  = gpuArray(int32(DataSize));
SEBSGPU_WORLD.SEBS_kb_1.GridSize=ceil(DataSize/SEBSGPU_WORLD.SEBS_kb_1.ThreadBlockSize(1));
% Linearizing the Input Variables
ComptMask_D = gpuArray(single(reshape(ComptMask,DataSize,1)));
fc_D        = gpuArray(single(reshape(fc,DataSize,1)));
LAI_D       = gpuArray(single(reshape(LAI,DataSize,1)));
NDVI_D      = gpuArray(single(reshape(NDVI,DataSize,1)));
LST_K_D     = gpuArray(single(reshape(LST_K,DataSize,1)));
hc_D        = gpuArray(single(reshape(hc,DataSize,1)));
Zref_D      = gpuArray(single(reshape(Zref,DataSize,1)));
Uref_D      = gpuArray(single(reshape(Uref,DataSize,1)));
Pref_D      = gpuArray(single(reshape(Pref,DataSize,1)));
P0_D        = gpuArray(single(reshape(P0,DataSize,1)));
Ps_D        = gpuArray(single(reshape(Ps,DataSize,1)));
Tref_K_D    = gpuArray(single(reshape(Tref_K,DataSize,1)));
qa_ref_D    = gpuArray(single(reshape(qa_ref,DataSize,1)));

z0m_D   = gpuArray(single(zeros(DataSize,1)));
d0_D	= gpuArray(single(zeros(DataSize,1)));
z0h_D   = gpuArray(single(zeros(DataSize,1)));

[z0m_D, d0_D, z0h_D]=feval(SEBSGPU_WORLD.SEBS_kb_1, ...
                            z0m_D, d0_D, z0h_D, ...
                            fc_D, LAI_D, NDVI_D, ...
                            LST_K_D, hc_D, Zref_D, ...
                            Uref_D, Pref_D, P0_D, Ps_D, ...
                            Tref_K_D, qa_ref_D, ...
                            ComptMask_D, DataSize_D);


if (ndims(ComptMask)==2)
    z0m = reshape(gather(z0m_D),InputVarSize(1),InputVarSize(2));
    d0  = reshape(gather(d0_D),InputVarSize(1),InputVarSize(2));
    z0h = reshape(gather(z0h_D),InputVarSize(1),InputVarSize(2));
else
    z0m = reshape(gather(z0m_D),InputVarSize(1),InputVarSize(2),InputVarSize(3));
    d0  = reshape(gather(d0_D),InputVarSize(1),InputVarSize(2),InputVarSize(3));
    z0h = reshape(gather(z0h_D),InputVarSize(1),InputVarSize(2),InputVarSize(3));
end

z0m(z0m==-9999)=NaN;
d0(d0==-9999)=NaN;
z0h(z0h==-9999)=NaN;

end
