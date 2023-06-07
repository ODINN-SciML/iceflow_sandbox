####################################################################
#### This script defines and runs functions useful to benchmark ####
#### the differences in ice thickness depending on the simulation ##
#### time on a set of 12 different glaciers                     ####
####################################################################

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
function bench_time_diff(n_y, gdir, reltol) #to benchmark the different simulation time

    fls=gdir.read_pickle("model_flowlines")
    bed_o = fls[end].bed_h
    surface_o = fls[end].surface_h
    widths_o = fls[end].widths_m
    dx_o = fls[end].dx_meter

    nx_o = length(bed_o)

    diag = gdir.get_diagnostics()
    glen_a_o = diag["inversion_glen_a"]


    #time parameters
    tspan = (0.0, n_y*sec_in_year)
    y0=2003.0

    ude_benchmark = Dict("ude_settings"=>[], "difference_m"=>[])

    UDE_settings = Dict("reltol"=>reltol,"n_years"=>[])
    UDE_settings["n_years"] = n_y

    println("Benchmarking UDE settings: ", UDE_settings)
    push!(ude_benchmark["ude_settings"], UDE_settings)
           
    iceflow_sol = glacier_evolution(gdir=gdir, 
                                        dx=dx_o, # grid resolution in m
                                        nx=nx_o,  # grid size
                                        width=widths_o,  # glacier width in m 
                                        glen_a=glen_a_o,  # ice stiffness 2.4e-24
                                        solver = RDPK3Sp35(),
                                        reltol=reltol,
                                        bed_hs=bed_o,
                                        surface_ini=surface_o,
                                        y0=y0,
                                        tspan=tspan)


    #also benchmarking oggm solution
    
    workflow.execute_entity_task(tasks.run_random_climate, gdir, y0 = y0, nyears=n_y,seed=1,store_fl_diagnostics=true)

    f = gdir.get_filepath("fl_diagnostics")
    ds = xr.open_dataset(f,group="fl_0")
            
    push!(ude_benchmark["difference_m"], mean(ds["thickness_m"].data[Int(n_y+1),:] .- iceflow_sol[end]))

    GC.gc()
    return ude_benchmark
end

function bench_gla(rgi_id,reltol) #to benchmark on different glaciers
    
    reltol = reltol

    #Initialize glacier
    rgi = [rgi_id]
    gdirs=init_gdirs(rgi)
    gdir=gdirs[1]

    PARAMS["evolution_model"] = "FluxBased"
    tasks.init_present_time_glacier(gdir)

   #Benchmark every simulation time for one given glacier
    ude_time =[50,100,150,200,400,600,800,1000]

    ude_benchmarks = pmap(ny -> bench_time_diff(ny, gdir, reltol), ude_time) 


    gla_benchmark= Dict("id"=>[], "sim_time"=>[])
    push!(gla_benchmark["id"], rgi_id)
    push!(gla_benchmark["sim_time"], ude_benchmarks)

    return gla_benchmark

end


function BenchmarkingAll(filename) #to benchmark all the glaciers from the list

    rgi_ids = ["RGI60-01.05355","RGI60-11.03232","RGI60-11.03638","RGI60-11.03643",
            "RGI60-11.03646","RGI60-11.03671","RGI60-11.03674","RGI60-11.03756",
            "RGI60-14.07524","RGI60-15.03591","RGI60-16.00543","RGI60-16.01339"]

    reltol = 1e-8
    results = pmap(r -> bench_gla(r,reltol), rgi_ids)

    save(filename, "data", results)

end 


### MAIN ###
filename= "data/benchmark_time_glena_reltol-8_diff.jld2"
#BenchmarkingAll(filename)


rgi_ids = ["RGI60-01.05355","RGI60-11.03232","RGI60-11.03638","RGI60-11.03643",
"RGI60-11.03646","RGI60-11.03671","RGI60-11.03674","RGI60-11.03756",
"RGI60-14.07524","RGI60-15.03591","RGI60-16.00543","RGI60-16.01339"]

ude = load(filename)
#ude_solvers_str=["BS3","OwrenZen3","Ralston","RDPK3Sp35","CKLLSRK54_3C"]
ude_time_str=["50","100","150","200","400","600","800","1000"]
ude_time=[50,100,150,200,400,600,800,1000]

ng=length(rgi_ids)
ns= length(ude_time_str)

#Creating an adequate Dataframe with the results 
diff=[]
for r=1:ng
    append!(diff,[ude["data"][r]["sim_time"][1][i]["difference_m"][1] for i =1:ns]) #getting the data in the right order 
        
end 

df=DataFrame(:id => repeat(rgi_ids,inner=ns),
                    :sim_time => repeat(ude_time,outer=ng),
                    :difference => diff)


### Plotting 
m_oggm = round(median(df[:,:difference]),digits=1)
df_group = groupby(df,:sim_time)
df_gla = groupby(df,:id)



df.sim_time= categorical(df.sim_time)
df.id= categorical(df.id)

label_gla=rgi_ids
for k=1:ng
    med = round(median(df_gla[k][:,:difference]),digits=1)
    label_gla[k] = rgi_ids[k]*" (median = $med m)"
end 


gr(bg = :ghostwhite)
@df df plot(:sim_time,:difference,group=:id,size=(850, 750),title="Benchmarking ice thickness difference depending on simulation time",
                                            xlabel="Simulation time [years]",xrotation=90,xtickfont=8, ylabel ="Thickness difference [m] (oggm - julia)",linealpha=0.3,
                                            minorgrid=true,legend=:bottom,legendcolumns=2,
                                            label=reshape(rgi_ids, 1, :),legendfontsize=7,palette=:Paired_12,
                                                markershape=:circle,ma=0.9,markerstrokecolor = :white,ms=5) #palette=cgrad(:matter, 5, categorical = true)
