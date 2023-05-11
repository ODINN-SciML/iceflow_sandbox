include("oggm_access.jl")
include("1D_SIA.jl")
include("1D_SIA_raw.jl")
using NCDatasets
using BenchmarkTools
using Distributed
using Plots

#Choose one glacier
rgi_ids=["RGI60-11.03671"]
gdirs=init_gdirs(rgi_ids)
gdir=gdirs[1]
gla_name=gdir.name

#Getting the flowlines 
PARAMS["evolution_model"] = "FluxBased"
tasks.init_present_time_glacier(gdir)

#=
f = gdir.get_filepath("fl_diagnostics")
ds = NCDataset(f)
fl_id=ds.dim["flowlines"]-1
ds2=ds.group["fl_$fl_id"]
bed_oggm=ds2["bed_h"][:]
surface_t0=bed_oggm .+ ds2["thickness_m"][:,1]
=#

fls=gdir.read_pickle("model_flowlines")
bed_o = fls[end].bed_h
surface_o = fls[end].surface_h
widths_o = fls[end].widths_m
dx_o = fls[end].dx_meter


diag = gdir.get_diagnostics()
glen_a_o = diag["inversion_glen_a"]
println(glen_a_o)



#using the solver

iceflow_sol =@btime glacier_evolution(gdir=gdir, 
                                dx=dx_o, # grid resolution in m
                                nx=length(bed_o),  # grid size
                                width=widths_o,  # glacier width in m 
                                glen_a= glen_a_o,  # ice stiffness 2.4e-24
                                n_years=15.0,  # simulation time in years
                                solver = Ralston(),
                                reltol=1e-6,
                                bed_hs=bed_o,
                                surface_ini=surface_o)


plot(bed_o, c="brown",label="bed",title="$gla_name",ylabel="Elevation (m.a.s.l.)")
display(plot!(iceflow_sol[end] .+ bed_o, c="blue",label="surface solver"))

#using a numerical scheme ("raw")

xc, bed_h, surface_h, years, volume, long =@btime glacier_evolution_optim(gdir=gdir,
                                                                    dx=dx_o,nx=length(bed_o),
                                                                    width=widths_o,
                                                                    glen_a= 2.4e-14,
                                                                    bed_h=bed_o,
                                                                    surface_ini=surface_o,
                                                                    n_years=15.0)


#Comparing now with oggm

@btime begin 
workflow.execute_entity_task(tasks.run_from_climate_data, gdir,
                    climate_filename="climate_historical",
                    ys=2004, ye=2019,store_fl_diagnostics=true)

end 


#plot!(bed_h, color="black",label="Bedrock raw")
p_flowline = plot!(surface_h, color="red", label="surface raw")
display(p_flowline)
plot!(surface_o,color="green",label="surface initiale")

f = gdir.get_filepath("fl_diagnostics")
ds = NCDataset(f)
fl_id=0
ds2=ds.group["fl_$fl_id"]
plot!(ds2["bed_h"][:,end]+ds2["thickness_m"][:,end],linestyle=:dash,color="black",label="surface oggm")