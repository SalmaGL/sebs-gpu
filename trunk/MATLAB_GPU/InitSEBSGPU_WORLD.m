function [SEBSGPU_WORLD]=InitSEBSGPU_WORLD()
% Preparing SEBS GPU Kernels
SEBS_kb_1_Interface='float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int';
SEBSGPU_WORLD.SEBS_kb_1=parallel.gpu.CUDAKernel('SEBS_kb_1.ptx',SEBS_kb_1_Interface);
SEBSGPU_WORLD.SEBS_kb_1.ThreadBlockSize=64;

SEBS_EnergyBalance_Interface='float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int DataSize';
SEBSGPU_WORLD.SEBS_EnergyBalance=parallel.gpu.CUDAKernel('SEBS_EnergyBalance.ptx',SEBS_EnergyBalance_Interface);
SEBSGPU_WORLD.SEBS_EnergyBalance.ThreadBlockSize=64;

SEBS_Daily_Evapotranspiration_Interface='float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int';
SEBSGPU_WORLD.SEBS_Daily_Evapotranspiration=parallel.gpu.CUDAKernel('SEBS_Daily_Evapotranspiration.ptx',SEBS_Daily_Evapotranspiration_Interface);
SEBSGPU_WORLD.SEBS_Daily_Evapotranspiration.ThreadBlockSize=64;
end