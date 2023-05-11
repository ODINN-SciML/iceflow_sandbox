#%%

# -*- coding: utf-8 -*-
using DifferentialEquations

using Random
using SugarBLAS
using Statistics: median
# using AbbreviatedStackTraces
using Logging: global_logger
using Revise, BenchmarkTools
using LinearAlgebra
using TerminalLoggers: TerminalLogger
using NetCDF
using Plots
include("oggm_access.jl")
global_logger(TerminalLogger())


const sec_in_day = 24.0 * 60.0 * 60.0
const sec_in_year = sec_in_day * 365.0
const glen_n = 3.0
const ρ = 900.0
const g = 9.81

#@views diff1(A) = (A[2:end] .- A[1:end - 1])
#@views diff2(A) = (A[3:end] .- A[1:end - 2])
diff1!(O, I, dx) = @views @. O[2:end] = (I[2:end] - I[1:end - 1]) / dx
diff2!(O, I, dx) = @views @. O[2:end-1] = (I[3:end] - I[1:end - 2]) / dx
avg!(O, I) = @views @. O[2:end] = (I[2:end] + I[1:end - 1])./2.0

function glacier_evolution(;
    gdir,
    dx,
    nx,  # grid size e.g 200
    width,  # glacier width in m
    glen_a,  # ice stiffness 2.4e-24
    n_years,  # simulation time in years
    solver,
    reltol,
    bed_hs,
    surface_ini
)
    
    mbmod=massbalance.MultipleFlowlineMassBalance(gdir,
                                                mb_model_class=massbalance.MonthlyTIModel,
                                                filename="climate_historical")



    function get_mb2(heights,y)
        mb = mbmod.get_annual_mb(heights, year=y, fl_id=0)
        return mb 
    end
    
    #bed_hs = collect(LinRange(top_h, bottom_h, nx))

    let
    #surface_h = copy(bed_hs)
    surface_h = copy(surface_ini)
    
    #H = bed_hs .* 0.0
    H = surface_h .- bed_hs
    
    #pre-allocation
    surface_gradient = zeros(nx)
    surface_gradient_s = zeros(nx)

    diffusivity = zeros(nx)
    diffusivity_s = zeros(nx)

    grad_x_diff = zeros(nx)
    flux_div = zeros(nx)

    mb=zeros(nx-1)
    dH=zeros(nx)


    p = (nx, dx, width, glen_a, bed_hs, surface_gradient,surface_gradient_s,
            diffusivity,diffusivity_s,grad_x_diff)

    Γ = (2 .* glen_a .* width .* (ρ * g).^glen_n)/(glen_n .+ 2)
    y=2004

    function iceflow!(dH, H, p, t)
        # Retrieve model parameters
        nx, dx, width, glen_a, bed_hs,surface_gradient, surface_gradient_s, diffusivity,
             diffusivity_s, grad_x_diff= p

        # Clip negative ice thickness values
        H .= ifelse.(H.<0.0, 0.0, H)
        @assert H[end-2] .== 0.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"

        # Surface gradient
        #surface_gradient[2:end-1] .= diff2(surface_h) ./ (2*dx)
        diff2!(surface_gradient, surface_h, 2.0*dx)

        # Diffusivity
        diffusivity .= Γ .* H.^(glen_n.+2) .* surface_gradient.^(glen_n.-1)

        # Ice flux in a staggered grid
        avg!(diffusivity_s,diffusivity)

        #surface_gradient_s[2:end] .= diff1(surface_h) ./ dx
        diff1!(surface_gradient_s, surface_h, dx)



        grad_x_diff .= surface_gradient_s .* diffusivity_s

        #flux_div = diff1(grad_x_diff) ./ dx
        diff1!(flux_div, grad_x_diff, dx)

        # Mass balance
        mb .= get_mb2(surface_h[begin:end-1],y)

        
        # # Ice thickness update: old + flux div + mb
        #dH .= zeros(size(dH))
        # @show median(mb)
        dH[begin:end-1] .= (flux_div[2:end] ./ width[2:end]) .+ mb


        

        surface_h .= bed_hs .+ H

        annee = t ÷ sec_in_year
        y=2004+annee


    end

    tspan = (0.0, n_years*sec_in_year)

    iceflow_prob = ODEProblem(iceflow!,H,tspan,p)

    iceflow_sol = solve(iceflow_prob,solver,dt=10*sec_in_day, reltol=reltol,save_everystep=false,dense=false)

   #= 
    if solver == nothing
        iceflow_sol = solve(iceflow_prob,alg=Euler(),#alg=ImplicitMidpoint(autodiff=false), #alg=ImplicitEuler(autodiff=false), 
                            dt=10*sec_in_day,reltol=1e-6,adaptive=false, save_everystep=false, progress=true, 
                            progress_steps = 10)
    else
        iceflow_sol = solve(iceflow_prob, solver,reltol=1e-6,dense=false,dt=10*sec_in_day,
            save_everystep=false, progress=true, progress_steps = 10)
    end

    #plot(bed_hs, c="brown",label="bed",title="Glacier geometry at the end of the simulation",ylabel="Elevation (m.a.s.l.)")
    #display(plot!(iceflow_sol[end] .+ bed_hs, c="blue",label="surface"))
    =#
    return iceflow_sol

    end #let

end

