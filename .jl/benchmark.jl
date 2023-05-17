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
global_logger(TerminalLogger())


#Choose one glacier
rgi_ids=["RGI60-11.01450"]
gdirs=init_gdirs(rgi_ids)
gdir=gdirs[1]


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


reltol = 1e-6

function benchmark_setting(setting)

    ude_benchmark = Dict("ude_settings"=>[], "time_stats"=>[])

    UDE_settings = Dict("reltol"=>reltol,"solver"=>[])
    UDE_settings["solver"] = setting

    println("Benchmarking UDE settings: ", UDE_settings)
    push!(ude_benchmark["ude_settings"], UDE_settings)

    try
        t_stats = @timed glacier_evolution(gdir=gdir, 
        dx=dx_o, # grid resolution in m
        nx=length(bed_o),  # grid size
        width=widths_o,  # glacier width in m 
        glen_a= 2.4e-24,  # ice stiffness 2.4e-24
        n_years=15.0,  # simulation time in years
        solver = setting,
        reltol=UDE_settings["reltol"],
        bed_hs=bed_o,
        surface_ini=surface_o)

        # Save stats for each solver

        push!(ude_benchmark["time_stats"], t_stats)

  
    catch error
        println("ERROR: ", error)
        @warn "Solver not working. Skipping..."
    end

    GC.gc()
    return ude_benchmark
end

ude_solvers = [Euler(),BS3(),OwrenZen3(), Ralston(), RDPK3Sp35(), CKLLSRK54_3C()]
#Tried but glacier exceeds boundaries : 
#VCABM(), Vern6(), AN5(),AB3(), KenCarp3(autodiff=false),TRBDF2(autodiff=false),ROCK4(),QNDF(autodiff=false), Tsit5(),Rodas4P()
#ImplicitEuler, ImplicitMidpoint,RK4, DP5

#Vern8,9 , Feagin10,MSRK5, KuttaPRK2p5,SSPRK22, TanYam7 instables

end #everywhere 


### Main ###

ude_benchmarks=zeros(length(ude_solvers))

# Benchmark every solver in parallel
ude_benchmarks = pmap(ude_solver -> benchmark_setting(ude_solver), ude_solvers) 
                
save_object("benchmark.jld2",ude_benchmarks)


### A look at the results ###
ude = load("benchmark.jld2")

#=
for ude_benchmark in ude["single_stored_object"]
    if length(ude_benchmark["ude_settings"]) > 0
            println(ude_benchmark["ude_settings"][1]["solver"], " - ", ude_benchmark["time_stats"][1].time)
    end
end
=#

n=length(ude_solvers)
stats = DataFrame(:id => 1:n, :name => ude_solvers, 
                  :reltol => ones(n)*reltol,:t => [ude["single_stored_object"][i]["time_stats"][1].time for i =1:n],
                  :bytes => [ude["single_stored_object"][i]["time_stats"][1].bytes for i =1:n],
                  :gctime => [ude["single_stored_object"][i]["time_stats"][1].gctime for i =1:n])
