include("oggm_access.jl")
include("1D_SIA_1.jl")
include("1D_SIA_raw_1.jl")
using NCDatasets
using BenchmarkTools

#Choose one glacier
rgi_ids=["RGI60-15.03591"]
gdirs=init_gdirs(rgi_ids)
gdir=gdirs[1]


#Getting the flowlines 
#PARAMS["evolution_model"] = "FluxBased"
tasks.init_present_time_glacier(gdir)
#tasks.run_constant_climate(gdir, nyears=10, y0=2000, store_fl_diagnostics=1)

#Getting oggm parameters for the linear MB gradient
tasks.apparent_mb_from_linear_mb(gdir)
lin_mb = gdir.read_pickle("linear_mb_params")

grad_mb_o=lin_mb["grad"]
ela_o=lin_mb["ela_h"]

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
@time begin

iceflow_sol = glacier_evolution(dx=dx_o,  # grid resolution in m
                                nx=length(bed_o),  # grid size
                                width=widths_o,  # glacier width in m 
                                glen_a= glen_a_o,  # ice stiffness 2.4e-24
                                ela_h=ela_o, # mass balance model Equilibrium Line Altitude
                                mb_grad=grad_mb_o,  # linear mass balance gradient (unit: [mm w.e. yr-1 m-1])
                                n_years=200.0,  # simulation time in years
                                solver = nothing,bed_hs=bed_o,surface_ini=surface_o)

end 

plot(bed_o, c="brown",label="bed",title="Glacier geometry at the end of the simulation",ylabel="Elevation (m.a.s.l.)")
display(plot!(iceflow_sol[end] .+ bed_o, c="blue",label="surface solver"))

#using a numerical scheme ("raw")
@time begin 
xc, bed_h, surface_h, years, volume, long = glacier_evolution_optim(dx=dx_o,nx=length(bed_o),
                                                                    width=widths_o,
                                                                    glen_a= glen_a_o,
                                                                    mb_grad=grad_mb_o,
                                                                    ela_h=ela_o,
                                                                    bed_h=bed_o,
                                                                    surface_ini=surface_o)

end 

#plot!(bed_h, color="black",label="Bedrock raw")
p_flowline = plot!(surface_h, color="red", label="surface raw")
display(p_flowline)
plot!(surface_o,color="green",label="surface initiale")
