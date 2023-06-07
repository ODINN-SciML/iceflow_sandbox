###################################################################
#### This script calls the functions 'galcier_evolution' and   ####
#### 'glacier_evolution_optim' and plots the resulting         ####
#### thicknesses                                              ####
###################################################################


using Revise
using Distributed
using BenchmarkTools
using Plots
using ProgressMeter
using OrdinaryDiffEq
using Statistics: median
# using AbbreviatedStackTraces
using Logging: global_logger
using LinearAlgebra
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#Importing necessary functions
include("oggm_access.jl")
include("1D_SIA.jl")
include("1D_SIA_raw.jl")

#Choose one glacier
rgi_ids=["RGI60-11.03638"]
gdirs=init_gdirs(rgi_ids)
gdir=gdirs[1]
gla_name=gdir.name

#Getting necessary parameters from OGGM
PARAMS["evolution_model"] = "FluxBased" #or SemiImplicit
tasks.init_present_time_glacier(gdir)

fls=gdir.read_pickle("model_flowlines")
bed_o = fls[end].bed_h #Glacier bed elevation
surface_o = (fls[end].surface_h) #Initial surface elevation
widths_o = fls[end].widths_m  #Widths along the flowline
nx_o = length(bed_o) #grid size
dx_o = fls[end].dx_meter #grid resolution in meters

diag = gdir.get_diagnostics()
glen_a_o = diag["inversion_glen_a"] #creep parameter A calibrated by OGGM 
println(glen_a_o)


#Numerical parameters for the simulation
n_years=100.0
tspan = (0.0, n_years*sec_in_year)
y0=2003.0

#using the solver
iceflow_sol =glacier_evolution(gdir=gdir, 
                                dx=dx_o,
                                nx=nx_o, 
                                width=widths_o,
                                glen_a= glen_a_o, 
                                solver = RDPK3Sp35(), #CKLLSRK54_3C()
                                reltol=1e-8,
                                bed_hs=bed_o,
                                surface_ini=surface_o,
                                y0=y0,
                                tspan=tspan)


plot(bed_o, c="brown",label="bed",title="$gla_name (for $n_years years)",ylabel="Elevation (m.a.s.l.)")
plot!(iceflow_sol[end] .+ bed_o, c="blue",label="surface solver")

#using a numerical scheme ("raw")
xc, bed_h, surface_h, years, volume, long =glacier_evolution_optim(gdir=gdir,
                                                                    dx=dx_o,nx=length(bed_o),
                                                                    width=widths_o,
                                                                    glen_a= glen_a_o,
                                                                    bed_h=bed_o,
                                                                    surface_ini=surface_o,
                                                                    n_years=n_years,
                                                                    y0=y0)



#Comparing now with oggm
workflow.execute_entity_task(tasks.run_random_climate, gdir, y0 = y0, nyears=n_years,
                                seed=1,store_fl_diagnostics=true)



plot!(surface_h, color="red", label="surface raw")
plot!(surface_o,color="green",label="surface initiale")

f = gdir.get_filepath("fl_diagnostics")
ds = xr.open_dataset(f,group="fl_0")
plot!(ds["bed_h"].data+ds["thickness_m"].data[Int(n_years+1),:],linestyle=:dash,color="black",label="surface oggm")


#println(mean(surface_h.-iceflow_sol[end].-bed_o))
