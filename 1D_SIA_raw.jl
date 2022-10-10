# -*- coding: utf-8 -*-
# using AbbreviatedStackTraces
# using Distributed

# processes = 10
# addprocs(processes - nprocs(); exeflags="--project")

# @everywhere begin
using Revise
using DifferentialEquations
using ProgressMeter, Infiltrator
using Plots
using Statistics: median
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

function glacier_evolution_optim(;
    dx=100.0,  # grid resolution in m
    nx=200,  # grid size
    width=600.0,  # glacier width in m
    top_h=3000.0,  # bed top altitude
    bottom_h=1200.0,  # bed bottom altitude
    glen_a=2.4e-24,  # ice stiffness
    ela_h=2600.0,  # mass balance model Equilibrium Line Altitude
    mb_grad=3.0,  # linear mass balance gradient (unit: [mm w.e. yr-1 m-1])
    n_years=700  # simulation time in years
)
    function get_mb(heights)
        mb = (heights .- ela_h) .* mb_grad
        return mb ./ sec_in_year ./ ρ
    end

    let 
    bed_h = collect(LinRange(top_h, bottom_h, nx))
    surface_h = copy(bed_h)
    thick = bed_h .* 0.0

    t = 0.0
    dt = sec_in_day * 10.0

    years = collect(0:(n_years+1))
    volume = zeros(size(years))
    length = zeros(size(years))

    new_thick = zeros(nx)
    diffusivity = zeros(nx)
    diffusivity_s = zeros(nx)
    surface_gradient = zeros(nx)
    surface_gradient_s = zeros(nx)
    grad_x_diff = zeros(nx)
    flux_div = zeros(nx-1)
    mb = zeros(nx-1)

    for (i, y) in enumerate(years)
        let end_t = y * sec_in_year
        # Time integration
        while t < end_t
            # This is to guarantee a precise arrival on a specific date if asked
            remaining_t = end_t - t
            if remaining_t < dt
                dt = remaining_t
            end

            # Surface gradient
            surface_gradient[2:end-1] .= diff2(surface_h) ./ (2.0*dx)

            # Diffusivity
            diffusivity .= width * ((ρ*g)^3.0) .* (thick.^3.0) .* surface_gradient.^2.0
            diffusivity .*= 2.0/(glen_a+2.0) * glen_a .* thick.^2.0

            # Ice flux in a staggered grid
            diffusivity_s[2:end] .= avg(diffusivity)

            surface_gradient_s[2:end] .= diff1(surface_h) ./ dx

            grad_x_diff .= surface_gradient_s .* diffusivity_s
            flux_div .= diff1(grad_x_diff) ./ dx

            # Mass balance
            mb .= get_mb(surface_h[begin:end-1])

            # Ice thickness update: old + flux div + mb
            new_thick[begin:end-1] .= thick[begin:end-1] .+ (dt/width) .* flux_div .+ dt.*mb

            if maximum(new_thick) > 400.0
                # xcoordinates
                xc = collect(0:nx-1) .* dx
                thick_dummy = ifelse.(new_thick.<0.0, 0.0, new_thick)
                surface_h .= bed_h .+ thick_dummy
                plot(xc, bed_h, color="black", title="Glacier geometry at the end of the simulation", label="Bedrock", ylabel="Elevation (m.a.s.l.)")
                p_flowline = plot!(xc, surface_h, color="slateblue", label="Ice")
                display(p_flowline)
            end

             # We can have negative thickness because of MB - correct here
            thick .= ifelse.(new_thick.<0.0, 0.0, new_thick)
            
            @assert !isnan(maximum(thick)) "NaN values!"
            @assert thick[end] == 0.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"

            # Prepare for next step 
            surface_h .= bed_h .+ thick
            t += dt
        end
        end # let

        volume[i] = sum(thick .* width .* dx)
        length[i] = sum(thick .> 0.0) .* dx
    end

    # xcoordinates
    xc = collect(0:nx-1) .* dx

    plot(xc, bed_h, color="black", title="Glacier geometry at the end of the simulation", label="Bedrock", ylabel="Elevation (m.a.s.l.)")
    p_flowline = plot!(xc, surface_h, color="slateblue", label="Ice")
    display(p_flowline)


    return xc, bed_h, surface_h, years, volume, length
    end # let
    
end

function wrapper_grad(grad)
    return glacier_evolution_optim(mb_grad=grad)
end

# end # @everywhere


#######  MAIN ########

@time xc, bed_h, surface_h, years, volume, length = glacier_evolution_optim()

# plot(xc, bed_h, color="black", title="Glacier geometry at the end of the simulation", label="Bedrock", ylabel="Elevation (m.a.s.l.)")
# p_flowline = plot!(xc, surface_h, color="slateblue", label="Ice")
# display(p_flowline)


# multi_grad = false

# if multi_grad

#     grads = collect(LinRange(1, 5, 20))

#     out = @showprogress pmap(grad -> wrapper_grad(grad), grads)

#     volumes = []
#     for o in out
#         xc, bed_h, surface_h, years, volume, length = o
#         push!(volumes, volume)
#     end

#     Plots.plot(grads, volumes[:, end-1] .* 1e-9, 
#             title="Final volume as function of MB gradient",
#             xlabel="MB grad",
#             ylabel = "Volume (km³)")
# end


