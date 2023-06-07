###################################################################
#### This script defines and runs functions useful to benchmark ####
#### various Julia solver on a set of 12 different glaciers    ####
###################################################################


using Revise
using Distributed
using ProgressMeter
using OrdinaryDiffEq
using Random
using SugarBLAS
using Statistics: median
# using AbbreviatedStackTraces
using Logging: global_logger
using LinearAlgebra
using BenchmarkTools
using JLD2
using DataFrames
using TerminalLoggers: TerminalLogger
using Plots 
using StatsPlots
using CategoricalArrays 
global_logger(TerminalLogger())

include("oggm_access.jl")
include("1D_SIA.jl")
include("1D_SIA_raw.jl") 


### FUNCTIONS ### 

function bench_solv(setting, gdir, reltol) #to benchmark the different solvers

    fls::PyObject = gdir.read_pickle("model_flowlines")
    bed_o::Vector{Float64} = fls[end].bed_h
    surface_o::Vector{Float64} = fls[end].surface_h
    widths_o::Vector{Float64} = fls[end].widths_m
    dx_o::Float64 = fls[end].dx_meter

    nx_o = length(bed_o)

    diag = gdir.get_diagnostics()
    glen_a_o = diag["inversion_glen_a"]

    n_years=100.0
    tspan = (0.0, n_years*sec_in_year)
    y0=2003.0

    ude_benchmark = Dict("ude_settings"=>[], "time_stats"=>[],"time_stats_oggm"=>[])

    UDE_settings = Dict("reltol"=>reltol,"solver"=>[])
    UDE_settings["solver"] = setting

    println("Benchmarking UDE settings: ", UDE_settings)
    push!(ude_benchmark["ude_settings"], UDE_settings)
           
    t_stats = @benchmark glacier_evolution(gdir=$gdir, 
                                        dx=$dx_o, # grid resolution in m
                                        nx=$nx_o,  # grid size
                                        width=$widths_o,  # glacier width in m 
                                        glen_a=$glen_a_o,  # ice stiffness 2.4e-24
                                        solver = $setting,
                                        reltol=$reltol,
                                        bed_hs=$bed_o,
                                        surface_ini=$surface_o,
                                        y0=$y0,
                                        tspan=$tspan)


    # Save stats for each solver
    push!(ude_benchmark["time_stats"], t_stats)

    #also benchmarking oggm solution
    t_stats_o =@benchmark workflow.execute_entity_task(tasks.run_random_climate, $gdir, y0 = 2003, nyears=100,seed=1,store_fl_diagnostics=true)

            
    push!(ude_benchmark["time_stats_oggm"], t_stats_o)


    GC.gc()
    return ude_benchmark
end

function bench_gla(rgi_id,reltol) #to benchmark the solvers on different glaciers
    
    reltol = reltol

    #Initialize glacier
    rgi = [rgi_id]
    gdirs=init_gdirs(rgi)
    gdir=gdirs[1]

    PARAMS["evolution_model"] = "SemiImplicit"
    tasks.init_present_time_glacier(gdir)

   #Benchmark every solver for one given glacier
    ude_solvers =[BS3(),CKLLSRK54_3C(),OwrenZen3(),RDPK3Sp35(),Ralston()] #[BS3(),OwrenZen3(), Ralston(), RDPK3Sp35(), CKLLSRK54_3C()]

    ude_benchmarks = pmap(ude_solver -> bench_solv(ude_solver, gdir, reltol), ude_solvers) 


    gla_benchmark= Dict("id"=>[], "solvers"=>[])
    push!(gla_benchmark["id"], rgi_id)
    push!(gla_benchmark["solvers"], ude_benchmarks)

    return gla_benchmark

end


function BenchmarkingAll(filename) #to benchmark all the glaciers from the list on different solvers
    rgi_ids=["RGI60-11.03638","RGI60-11.03671","RGI60-11.03643","RGI60-11.03674","RGI60-11.03756", #Argentière, Gébroulaz, Mer de Glace,St-Sorlin, Sarennes
    "RGI60-16.00543","RGI60-16.01339", #Zongo, Antizana
    "RGI60-11.03232", #Ossoue
    "RGI60-15.03591", #Mera
    "RGI60-11.03646",
    "RGI60-14.07524", #Siachen
    "RGI60-01.05355"] #Alexander (Alaska)
    reltol = 1e-8
    results = pmap(r -> bench_gla(r,reltol), rgi_ids)

    save(filename, "data", results)

end 


### MAIN ###
filename= "data/benchmark100y_glena_reltol-8_semiimp.jld2"
#BenchmarkingAll(filename)

rgi_ids=["RGI60-11.03638","RGI60-11.03671","RGI60-11.03643","RGI60-11.03674","RGI60-11.03756", #Argentière, Gébroulaz, Mer de Glace,St-Sorlin, Sarennes
    "RGI60-16.00543","RGI60-16.01339", #Zongo, Antizana
    "RGI60-11.03232", #Ossoue
    "RGI60-15.03591", #Mera
    "RGI60-11.03646",
    "RGI60-14.07524", #Siachen
    "RGI60-01.05355"] #Alexander (Alaska)

ude = load(filename)
ude_solvers_str=["BS3","CKLLSRK54_3C","OwrenZen3","RDPK3Sp35","Ralston"]

ng=length(rgi_ids)
ns= length(ude_solvers_str)

#Creating an adequate Dataframe with the results 
t_ms=[]
t_o =[]
mem=[]
mem_o=[]
for r=1:ng
    append!(t_ms,[mean(ude["data"][r]["solvers"][1][i]["time_stats"][1]).time*10^(-6) for i =1:ns]) #getting the data in the right order 
    append!(mem,[mean(ude["data"][r]["solvers"][1][i]["time_stats"][1]).memory*10^(-6) for i =1:ns])
    append!(t_o,[mean(ude["data"][r]["solvers"][1][i]["time_stats_oggm"][1]).time*10^(-6) for i =1:ns]) 
    append!(mem_o,[mean(ude["data"][r]["solvers"][1][i]["time_stats_oggm"][1]).memory*10^(-6) for i =1:ns])
        
end 

df=DataFrame(:id => repeat(rgi_ids,inner=ns),
                    :solv_name => repeat(ude_solvers_str,outer=ng),
                    :t_ms => t_ms,
                    :t_oggm_ms => t_o,
                    :memory_MiB => mem,
                    :memory_oggm_MiB => mem_o )

### Plotting 
m_oggm = round(median(df[:,:t_oggm_ms]),digits=1)
df_group = groupby(df,:solv_name)
label_solvers=ude_solvers_str
for k=1:ns
    med = round(median(df_group[k][:,:t_ms]),digits=1)
    label_solvers[k] = ude_solvers_str[k]*" (median = $med ms)"
end 

df.solv_name = categorical(df.solv_name)

gr(bg = :ghostwhite)
@df df scatter(:id,:t_ms,group=:solv_name,size=(850, 750),title="Benchmarking solvers on 12 glaciers (glen A calibrated)",
                                            xlabel="RGI ID",xrotation=90,xtickfont=8, ylabel ="time [ms]",yscale=:log10,
                                            ylimits=(10,10000),minorgrid=true,legend=:top, legendcolumns=2,
                                            label=reshape(label_solvers, 1, :),legendfontsize=7,palette=:Paired_9,
                                                markershape=:circle,ma=0.9,markerstrokecolor = :white,ms=5) #palette=cgrad(:matter, 5, categorical = true)


df_group_o = groupby(df,:id)
df_t_o=[]
for l=1:ng
    med = round(median(df_group_o[l][:,:t_oggm_ms]),digits=1)
    append!(df_t_o,med)
end 

scatter!(rgi_ids, df_t_o,markershape=:cross,ma=0.9,markerstrokecolor = :white,label="OGGM (median = $m_oggm ms)",
                            ms=5,markerclor=:red)

