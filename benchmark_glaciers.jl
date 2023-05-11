include("oggm_access.jl")
include("1D_SIA.jl")
include("1D_SIA_raw.jl")
using NCDatasets
using BenchmarkTools
using Distributed

@everywhere begin
using DifferentialEquations

using JLD2
using Random
using SugarBLAS
using Statistics: median
# using AbbreviatedStackTraces
using Logging: global_logger
using Revise, BenchmarkTools
using LinearAlgebra
using TerminalLoggers: TerminalLogger
using DataFrames
using Plots
global_logger(TerminalLogger())


PARAMS["evolution_model"] = "FluxBased"


function benchmark_setting_gla(rgi)
    rgi_id = [rgi]
    gdirs=init_gdirs(rgi_id)
    gdir=gdirs[1]
    gla_name=gdir.name
    gla_nb = findall(x->x==rgi,rgi_ids)[1]

    #Getting the flowlines 
    tasks.init_present_time_glacier(gdir)

    fls=gdir.read_pickle("model_flowlines")
    bed_o = fls[end].bed_h
    surface_o = fls[end].surface_h
    widths_o = fls[end].widths_m
    dx_o = fls[end].dx_meter

    diag = gdir.get_diagnostics()
    glen_a_o = diag["inversion_glen_a"]



    ude_benchmark = Dict("glacier_id"=>[], "time_stats"=>[],"time_stats_oggm"=>[])

    println("Benchmarking glacier settings: ", rgi)
    push!(ude_benchmark["glacier_id"], rgi)

    try
        t_stats = @timed glacier_evolution(gdir=gdir, 
        dx=dx_o, # grid resolution in m
        nx=length(bed_o),  # grid size
        width=widths_o,  # glacier width  in years
        glen_a= 2.4e-24,
        n_years=15.0,
        solver = Ralston(),
        reltol=1e-6,
        bed_hs=bed_o,
        surface_ini=surface_o)

        push!(ude_benchmark["time_stats"], t_stats)

        t_stats_o =@timed workflow.execute_entity_task(tasks.run_from_climate_data, gdir,
                                climate_filename="climate_historical",
                                ys=2004, ye=2019,store_fl_diagnostics=true)
            
        push!(ude_benchmark["time_stats_oggm"], t_stats_o)
        


    catch error
        println("ERROR: ", error)
        @warn "Solver not working. Skipping..."
    end

    iceflow_sol = glacier_evolution(gdir=gdir, dx=dx_o, nx=length(bed_o), width=widths_o,glen_a= 2.4e-24,n_years=15.0,solver = Ralston(),
    reltol=1e-6,bed_hs=bed_o,surface_ini=surface_o)

    plot(bed_o, c="brown",label="bed",title="$gla_name",ylabel="Elevation (m.a.s.l.)")
    plot!(iceflow_sol[end] .+ bed_o, c="blue",label="surface solver")
    plot!(surface_o,color="green",label="surface initiale")

    f = gdir.get_filepath("fl_diagnostics")
    ds = NCDataset(f)
    fl_id=0
    ds2=ds.group["fl_$fl_id"]
    plot!(ds2["bed_h"][:,end]+ds2["thickness_m"][:,end],linestyle=:dash,color="black",label="surface oggm")

    savefig("myplot_$gla_nb.png")   

    GC.gc()
    return ude_benchmark
end


end #everywhere 


### Main ###

rgi_ids=["RGI60-11.03638","RGI60-11.03671","RGI60-11.03643","RGI60-11.03674","RGI60-11.03756", #Argentière, Gébroulaz, Mer de Glace,St-Sorlin, Sarennes
        "RGI60-16.00543","RGI60-16.01339", #Zongo, Antizana
        "RGI60-11.03232", #Ossoue
        "RGI60-15.03591"] #Mera 

ude_benchmarks=zeros(length(rgi_ids))

# Benchmark every solver in parallel
ude_benchmarks = pmap(glacier -> benchmark_setting_gla(glacier), rgi_ids) 
                
save_object("benchmark_gla.jld2",ude_benchmarks)


### A look at the results ###
ude = load("benchmark_gla.jld2")

#=
for ude_benchmark in ude["single_stored_object"]
    if length(ude_benchmark["ude_settings"]) > 0
            println(ude_benchmark["ude_settings"][1]["solver"], " - ", ude_benchmark["time_stats"][1].time)
    end
end
=#

n=length(rgi_ids)
stats = DataFrame(:id => 1:n, :name => rgi_ids, 
                  :t => [ude["single_stored_object"][i]["time_stats"][1].time for i =1:n],
                  :bytes => [ude["single_stored_object"][i]["time_stats"][1].bytes for i =1:n],
                  :t_oggm=> [ude["single_stored_object"][i]["time_stats_oggm"][1].time for i =1:n],
                  :bytes_oggm => [ude["single_stored_object"][i]["time_stats_oggm"][1].bytes for i =1:n])
