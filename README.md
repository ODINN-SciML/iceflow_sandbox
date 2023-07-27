# ODINN.jl iceflow sandbox

Ice flow model experiments for OGGM-ODINN. Currently under development by Lucille Gimenes.

Current developments:

Accelerating ice flow simulations in OGGM with Julia
- Coupling a 1D Shallow Ice Approximation model in Julia to OGGM.
- Coupling a 2D Shallow Ice Approximation from ODINN.jl to OGGM.

At the moment, the options available are :
- Run a glacier simulation with any type of climate with the Julia SIA model (by chosing `cfg.PARAMS["evolution_model"]="SIA_1D"` in OGGM) 
- Run a glacier spin up (fixed geometry or dynamic) 

## Overview 

### Ice flow model 1D 

The new OGGM class (available in the forked repository of OGGM on this user profile) `IceflowJuliaModel` and its child `SIA_1D` calls Julia code (availble in oggm.core > `SIA_1D.jl`) and a Julia solver from `DifferentialEquations.jl` to solve the SIA equation on glaciers. 

`IceFlowJuliaModel` is a copy of the `FlowlineModel` class, with the methods `run_until` and `run_until_and_store` adapted to be used with the functions `glacier_evolution` and `glacier_evolution_store` from `SIA_1D.jl`. `SIA_1D` has no method `step` compared to `FluxBasedModel` or `SemiImplicitModel`. The new class isn ow functionnal only works elevation bands (only one flowline) and doesn't have claving implemented. 


<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/glacier_plots/Siachen_profile.png" width="850"></center>

> Example with Siachen glacier : profile at the end of a GCM climate data run (GCM MRI-ESM2-0 with ssp126). In green is the initial surface of the glacier, and in blue is the solution given by the SemiImplicitModel for the year 2300. In dashed orange is the surface given by the SIA_1D model.

<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/glacier_plots/Siachen_evolution.png" width="900"></center>

> Evolution of volume, area and length of Siachen glacier, under the SSP126 of GCM MRI-ESM2-0

### Benchmarking 

To find which solver was the most efficient, a benchmark was conducted on 12 different glaciers (see notebook Benchmark_from_OGGM)

<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/benchmark_OGGM/random_500y.png" width="900"></center>

> Benchmark conducted with saving annual variables, as well as model geometry and flowline diagnostic variables.

It appears that the RDPK3Sp35() solver ("5-stage, third order low-storage scheme with embedded error estimator, optimized for compressible fluid mechanics")
is the most efficient solver. 
