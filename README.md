# iceflow_sandbox
Ice flow model experiments for OGGM-ODINN 

At the moment, the options available are :
- Retrieve flowline parameters and mass balance (for a random climate) from OGGM, to then solve the Shallow Ice Approximation (SIA) equation on the glacier flowline in Julia. 

## Overview 

### Ice flow model

The function `glacier_evolution` implemented in `1D_SIA.jl` solves the SIA equation along the flowline of a given glacier by using a solver from the `DifferentialEquations.jl` package. Useful packages and functions from OGGM are called with `oggm_access.jl`.

<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots_rgi_ids/100y/Merdeglace.svg" width="700"></center>

> Example with Mer de Glace glacier : In green is the initial surface of the glacier, and in dashed black is the solution given by OGGM after 100 years of simulation under a random climate. In blue is the surface computed by a Julia solver and in red, the 'raw' solution that is computed using a basic numerical scheme (function `glacier_evolution_optim` in `1D_SIA_raw.jl`). 

### Benchmarking 

To find which solver was the most efficient, a benchmark was conducted on 12 different glaciers (in `benchmark.jl`)


<center><img src="https://github.com/lucillegimenes/iceflow_sandbox/blob/main/plots_rgi_ids/benchmark100y_glena_reltol-8.png" width="700"></center>

It appears that the RDPK3Sp35() solver ("5-stage, third order low-storage scheme with embedded error estimator, optimized for compressible fluid mechanics")
is the most efficient solver, being overall 12 times quicker than the OGGM ice flow model. 

