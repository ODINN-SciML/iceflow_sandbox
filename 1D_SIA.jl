# -*- coding: utf-8 -*-
using DifferentialEquations
using Plots
using Statistics: median
# using AbbreviatedStackTraces
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

const sec_in_day = 24.0 * 60.0 * 60.0
const sec_in_year = sec_in_day * 365.0
const glen_n = 3.0
const ρ = 900.0
const g = 9.81

@views diff1(A) = (A[2:end] .- A[1:end - 1])
@views diff2(A) = (A[3:end] .- A[1:end - 2])
@views avg(A) = (A[2:end] .+ A[1:end - 1])./2.0

function glacier_evolution(;
    dx=100.0,  # grid resolution in m
    nx=200,  # grid size
    width=600.0,  # glacier width in m
    top_h=3000.0,  # bed top altitude
    bottom_h=1200.0,  # bed bottom altitude
    glen_a=2.4e-24,  # ice stiffness
    ela_h=2600.0,  # mass balance model Equilibrium Line Altitude
    mb_grad=3.0,  # linear mass balance gradient (unit: [mm w.e. yr-1 m-1])
    n_years=200.0,  # simulation time in years
    solver = Ralston()
)
    bed_h = collect(LinRange(top_h, bottom_h, nx))
    
    # surface_h = copy(bed_h)
    H = bed_h .* 0.0

    p = (nx, dx, width, glen_a, bed_h)

    function get_mb(heights)
        mb = (heights .- ela_h) .* mb_grad
        return mb ./ sec_in_year ./ ρ
    end

    function iceflow!(dH, H, p, t)
        # Retrieve model parameters
        nx, dx, width, glen_a, bed_h = p

        surface_h = bed_h .+ H

        # Clip negative ice thickness values
        H .= ifelse.(H.<0.0, 0.0, H)
        @assert H[end-2] == 0.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"

        # Surface gradient
        surface_gradient = zeros(nx)
        surface_gradient[2:end-1] = diff2(surface_h) ./ (2*dx)

        # Diffusivity
        Γ = (2 * glen_a * width * (ρ * g)^glen_n)/(glen_n + 2)
        diffusivity = Γ .* H.^(glen_n+2) .* surface_gradient.^(glen_n-1)

        # Ice flux in a staggered grid
        diffusivity_s = zeros(nx)
        diffusivity_s[2:end] = avg(diffusivity)

        surface_gradient_s = zeros(nx)
        surface_gradient_s[2:end] = diff1(surface_h) ./ dx

        grad_x_diff = surface_gradient_s .* diffusivity_s
        flux_div = diff1(grad_x_diff) ./ dx

        # Mass balance
        mb = get_mb(surface_h[begin:end-1])

        # # Ice thickness update: old + flux div + mb
        dH .= zeros(size(dH))
        # @show median(mb)
        dH[begin:end-1] .= (flux_div ./ width) .+ mb

    end

    tspan = (0.0, n_years*sec_in_year)

    iceflow_prob = ODEProblem(iceflow!,H,tspan,p)
    if solver == nothing
        iceflow_sol = solve(iceflow_prob, 
                            alg=Euler(),
                            #alg=ImplicitMidpoint(autodiff=false), 
                            #alg=ImplicitEuler(autodiff=false), 
                            dt=10*sec_in_day,
                            adaptive=false,
                            reltol=1e-7, 
                            save_everystep=false, 
                            progress=true, 
                            progress_steps = 10)
    else
        iceflow_sol = solve(iceflow_prob, solver,
                            reltol=1e-6, save_everystep=false, 
                            progress=true, progress_steps = 10)
    end


    bed_h = collect(LinRange(top_h, bottom_h, nx))
    plot(bed_h, c="brown")
    display(plot!(iceflow_sol[end] .+ bed_h, c="blue"))

    return iceflow_sol

end

@time iceflow_sol = glacier_evolution(dx=100.0,  # grid resolution in m
                                nx=200,  # grid size
                                width=600.0,  # glacier width in m
                                top_h=3000.0,  # bed top altitude

                                bottom_h=1200.0,  # bed bottom altitude
                                glen_a=2.4e-24,  # ice stiffness
                                ela_h=2600.0,  # mass balance model Equilibrium Line Altitude
                                mb_grad=3.0,  # linear mass balance gradient (unit: [mm w.e. yr-1 m-1])
                                n_years=700.0,  # simulation time in years
                                solver = nothing
                                )

# plot(iceflow_sol.u[end])

