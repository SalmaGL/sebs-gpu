# Introduction #

The Surface Energy Balance System (SEBS) for estimation of
turbulent heat fluxes;

A High Performance GPU Implementation of
Surface Energy Balance System (SEBS) based on
CUDA-C


# Details #

**The Surface Energy Balance System (SEBS) for estimation of
turbulent heat fluxes**

Abstract:
A Surface Energy Balance System (SEBS) is proposed for the estimation of atmospheric turbulent fluxes and evaporative fraction using
satellite earth observation data, in combination with meteorological information at proper scales. SEBS consists of: a set of tools for the
determination of the land surface physical parameters, such as albedo, emissivity, temperature, vegetation coverage etc., from spectral reflectance and radiance measurements; a model for the determination of the roughness length for heat transfer; and a new formulation for the determination of the evaporative fraction on the basis of energy balance at limiting cases. Four experimental data sets are used to assess the reliabilities of SEBS. Based on these case studies, SEBS has proven to be capable to estimate turbulent heat fluxes and evaporative fraction at various scales
with acceptable accuracy. The uncertainties in the estimated heat fluxes are comparable to in-situ measurement uncertainties.


_Access the full publicaiton at:_ http://www.hydrol-earth-syst-sci.net/6/85/2002/hess-6-85-2002.pdf;

**A High Performance GPU Implementation of
Surface Energy Balance System (SEBS) based on
CUDA-C**

Abstract:
This paper introduces a new implementation of the Surface Energy Balance System (SEBS) algorithm harnessing the many cores available on Graphics Processing Units (GPUs). This new implementation uses Compute Unified Device Architecture C (CUDA-C) programming model and is designed to be executed on a system equipped with NVIDIA®'s graphic cards. The output of the new implementation is compared to a MATLAB code that has already been fully tested in the Water Cycle Multimission Observation Strategy (WACMOS) project. The code is timed against both MATLAB and a purely high-performance C implementation of the same algorithm. The code has been tested on several different NVIDIA® cards, with different compute capabilities. The authors have decided to provide the entire source code to the scientific community free of charge; hence, at the end, the instruction on how to obtain the code is also presented.

_Access the full publicaiton at:_
http://www.sciencedirect.com/science/article/pii/S1364815212003106;

_Access the code at:_
https://code.google.com/p/sebs-gpu/