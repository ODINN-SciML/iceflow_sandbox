# -*- coding: utf-8 -*-
# using AbbreviatedStackTraces
# using Distributed

# processes = 10
# addprocs(processes - nprocs(); exeflags="--project")

# @everywhere begin
using Revise, BenchmarkTools
using ProgressMeter, Infiltrator
using Plots

const sec_in_day = 24.0f0 * 60.0f0 * 60.0f0
const sec_in_year = sec_in_day * 365.0f0
const glen_n = 3.0f0
const ρ = 900.0f0
const g = 9.81f0

@views diff1(A) = (A[2:end] .- A[1:end - 1])
@views diff2(A) = (A[3:end] .- A[1:end - 2])
@views avg(A) = (A[2:end] .+ A[1:end - 1])./2.0f0

function glacier_evolution_optim(;
    dx=100.0f0,  # grid resolution in m
    nx=200,  # grid size
    width=600.0f0,  # glacier width in m
    top_h=3000.0f0,  # bed top altitude
    bottom_h=1200.0f0,  # bed bottom altitude
    glen_a=2.4f-24,  # ice stiffness
    ela_h=2600.0f0,  # mass balance model Equilibrium Line Altitude
    mb_grad=3.0f0,  # linear mass balance gradient (unit: [mm w.e. yr-1 m-1])
    n_years=700  # simulation time in years
)
    function get_mb(heights)
        mb = (heights .- ela_h) .* mb_grad
        return mb ./ sec_in_year ./ ρ
    end

    let 
    bed_h = collect(Float32,LinRange(top_h, bottom_h, nx))
    surface_h = copy(bed_h)
    thick = bed_h .* 0.0f0

    t = 0.0f0
    dt = sec_in_day * 10.0f0

    years = collect(Int32,0:(n_years+1))
    volume = zeros(Float32,size(years))
    length = zeros(Float32,size(years))

    new_thick = zeros(Float32,nx)
    diffusivity = zeros(Float32,nx)
    diffusivity_s = zeros(Float32,nx)
    surface_gradient = zeros(Float32,nx)
    surface_gradient_s = zeros(Float32,nx)
    grad_x_diff = zeros(Float32,nx)
    flux_div = zeros(Float32,nx-1)
    mb = zeros(Float32,nx-1)

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
            surface_gradient[2:end-1] .= diff2(surface_h) ./ (2.0f0*dx)

            # Diffusivity
            diffusivity .= width * ((ρ*g)^3.0f0) .* (thick.^3.0f0) .* surface_gradient.^2.0f0
            diffusivity .*= 2.0f0/(glen_a+2.0f0) * glen_a .* thick.^2.0f0

            # Ice flux in a staggered grid
            diffusivity_s[2:end] .= avg(diffusivity)

            surface_gradient_s[2:end] .= diff1(surface_h) ./ dx

            grad_x_diff .= surface_gradient_s .* diffusivity_s
            flux_div .= diff1(grad_x_diff) ./ dx

            # Mass balance
            mb .= get_mb(surface_h[begin:end-1])

            # Ice thickness update: old + flux div + mb
            new_thick[begin:end-1] .= thick[begin:end-1] .+ (dt/width) .* flux_div .+ dt.*mb

            # We can have negative thickness because of MB - correct here
            thick .= ifelse.(new_thick.<0.0f0, 0.0f0, new_thick)
            
            @assert thick[end] == 0.0f0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"

            # Prepare for next step 
            surface_h .= bed_h .+ thick
            t += dt
        end
        end # let

        volume[i] = sum(thick .* width .* dx)
        length[i] = sum(thick .> 0.0f0) .* dx
    end

    # xcoordinates
    xc = collect(0:nx-1) .* dx

    return xc, bed_h, surface_h, years, volume, length
    end # let
    
end

function wrapper_grad(grad)
    return glacier_evolution_optim(mb_grad=grad)
end

# end # @everywhere


#######  MAIN ########

@btime xc, bed_h, surface_h, years, volume, length = glacier_evolution_optim()

plot(xc, bed_h, color="black", title="Glacier geometry at the end of the simulation", label="Bedrock", ylabel="Elevation (m.a.s.l.)")
p_flowline = plot!(xc, surface_h, color="slateblue", label="Ice")
display(p_flowline)


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


