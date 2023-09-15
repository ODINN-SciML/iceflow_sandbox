# ODINN.jl iceflow sandbox

Ice flow model experiments for OGGM-ODINN. Currently under development by Lucille Gimenes.

Current developments:

Accelerating ice flow simulations in OGGM with Julia
- Coupling a 1D Shallow Ice Approximation model in Julia to OGGM.
- Coupling a 2D Shallow Ice Approximation from ODINN.jl to OGGM.

At the moment, the options available are :
- Run a glacier simulation with any type of climate with the Julia SIA model (simply by chosing `cfg.PARAMS["evolution_model"]="SIA1D"` in OGGM) 
- Run a glacier spin up (fixed geometry or dynamic) 

## Overview 

### Ice flow model 1D 

The new OGGM class (available in the forked repository of OGGM on this user profile : https://github.com/lucillegimenes/oggm ) `IceflowJuliaModel` and its child `SIA1D` calls Julia code (availble in oggm.core > `SIA1D.jl`& `SIA1D_utils.jl`) and a Julia solver from `DifferentialEquations.jl` to solve the SIA equation on glaciers. 

`IceFlowJuliaModel` is a copy of the `FlowlineModel` class, with the methods `run_until` and `run_until_and_store` adapted to be used with the functions `glacier_evolution` and `glacier_evolution_store` from `SIA_1D.jl`. `SIA_1D` has no method `step` compared to `FluxBasedModel` or `SemiImplicitModel`. The new class is now functionnal only with elevation bands (only one flowline) and doesn't have calving implemented. 


<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/OGGM_plots/Siachen_profile.png" width="850"></center>

> Example with Siachen glacier : profile at the end of a GCM climate data run (GCM MRI-ESM2-0 with ssp126). In green is the initial surface of the glacier, and in blue is the solution given by the SemiImplicitModel for the year 2300. In dashed orange is the surface given by the SIA_1D model.

<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/OGGM_plots/Siachen_evolution.png" width="900"></center>

> Evolution of volume, area and length of Siachen glacier, under the SSP126 of GCM MRI-ESM2-0

### Benchmarking 

To find which solver was the most efficient, a benchmark was conducted on 12 different glaciers (see notebook Benchmark_from_OGGM)

<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/OGGM_benchmarks/random_300y_20r.png" width="900"></center>

> Benchmark conducted with saving annual variables, as well as model geometry and flowline diagnostic variables.

It appears that the RDPK3Sp35() solver ("5-stage, third order low-storage scheme with embedded error estimator, optimized for compressible fluid mechanics")
is the most efficient solver. 

### Assessment of the new model SIA1D with climate experiment 

<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots/OGGM_plots/Experiments/sinus/khumbu_exp_sinus_spin_up_20r_length.png" width="800"></center>

> Evolution of volume, area and length of glacier RGI60-15.03733 (also called Khumbu) under experimental sinusoidal climate represented in bottom subplot.

## In this repository, you will find : 

- Notebooks illustrating the results from the iceflow model `SIA1D` : `OGGM_SIA1D_test_(Python).ipynb` is an exemple with simulated glacier Mer de glace, and `Benchmark_from_OGGM_(Python).ipynb` aims to benchmark overall the new iceflow model. Climate experiments are done with `Report_figures.ipynb`. These notebooks will work with `SIA1D` only if the used OGGM version is the one available at https://github.com/lucillegimenes/oggm.
- Scripts : most of the Julia files available in this directory was used during the early developpements of the iceflow model, but is not useful anymore. In scripts > AD is available a version of `SIA1D` that is compatible with Julia automatic differentiation. 
- Plots obtained from the notebooks

