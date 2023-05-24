
using Revise
using Distributed
using NCDatasets
using BenchmarkTools
using Plots

using ProgressMeter
using OrdinaryDiffEq
using Statistics: median
# using AbbreviatedStackTraces
using Logging: global_logger
using LinearAlgebra
using DataFrames
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("oggm_access.jl")
include("1D_SIA.jl")
include("1D_SIA_raw.jl")

#Choose one glacier
rgi_ids=["RGI60-01.05355"]
gdirs=init_gdirs(rgi_ids)
gdir=gdirs[1]
gla_name=gdir.name

#Getting the flowlines 
PARAMS["evolution_model"] = "FluxBased"
tasks.init_present_time_glacier(gdir)


fls=gdir.read_pickle("model_flowlines")
bed_o = fls[end].bed_h
surface_o = fls[end].surface_h
widths_o = fls[end].widths_m
dx_o = fls[end].dx_meter


diag = gdir.get_diagnostics()
glen_a_o = diag["inversion_glen_a"]
println(glen_a_o)


plot(bed_o, c="brown",label="bed",title="$gla_name",ylabel="Elevation (m.a.s.l.)")
#using the solver

iceflow_sol = glacier_evolution(gdir=gdir, 
                                dx=dx_o, # grid resolution in m
                                nx=length(bed_o),  # grid size
                                width=widths_o,  # glacier width in m 
                                glen_a= glen_a_o,  # ice stiffness 2.4e-24
                                n_years=100.0,  # simulation time in years
                                solver = CKLLSRK54_3C(),
                                reltol=1e-6,
                                bed_hs=bed_o,
                                surface_ini=surface_o)


plot(bed_o, c="brown",label="bed",title="$gla_name",ylabel="Elevation (m.a.s.l.)")
plot!(iceflow_sol[end] .+ bed_o, c="blue",label="surface solver")

#using a numerical scheme ("raw")


xc, bed_h, surface_h, years, volume, long = glacier_evolution_optim(gdir=gdir,
                                                                    dx=dx_o,nx=length(bed_o),
                                                                    width=widths_o,
                                                                    glen_a= 2.4e-24,
                                                                    bed_h=bed_o,
                                                                    surface_ini=surface_o,
                                                                    n_years=100.0)



#Comparing now with oggm



#workflow.execute_entity_task(tasks.run_from_climate_data, gdir,climate_filename="climate_historical",ys=2004, ye=2019,store_fl_diagnostics=true)
workflow.execute_entity_task(tasks.run_random_climate, gdir, y0 = 2003, nyears=100,
                                seed=1,store_fl_diagnostics=true)



plot!(surface_h, color="red", label="surface raw")
plot!(surface_o,color="green",label="surface initiale")

f = gdir.get_filepath("fl_diagnostics")
ds = NCDataset(f)
fl_id=0
ds2=ds.group["fl_$fl_id"]
plot!(ds2["bed_h"][:,end]+ds2["thickness_m"][:,end],linestyle=:dash,color="black",label="surface oggm")