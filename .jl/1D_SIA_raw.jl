# -*- coding: utf-8 -*-
# using AbbreviatedStackTraces
# using Distributed

# processes = 10
# addprocs(processes - nprocs(); exeflags="--project")
include("oggm_access.jl")

# @everywhere begin
    using Revise, BenchmarkTools
    using ProgressMeter, Infiltrator
    using Plots
    using NetCDF
    
    const sec_in_day = 24.0 * 60.0 * 60.0
    const sec_in_year = sec_in_day * 365.0
    const glen_n = 3.0
    const ρ = 900.0
    const g = 9.81
    
    diff1!(O, I, dx) = @views @. O[2:end] = (I[2:end] - I[1:end - 1]) / dx
    diff2!(O, I, dx) = @views @. O[2:end-1] = (I[3:end] - I[1:end - 2]) / dx
    avg!(O, I) = @views @. O[2:end] = (I[2:end] + I[1:end - 1])./2.0
    
    function glacier_evolution_optim(;
        gdir,
        dx,  # grid resolution in m
        nx,  # grid size
        width,  # glacier width in m
        glen_a,  # ice stiffness 2.4e-24
        n_years,
        bed_h,
        surface_ini # simulation time in years
        
    )

        mbmod=massbalance.MultipleFlowlineMassBalance(gdir)    

        function get_mb2(heights,y)
            mb = mbmod.get_annual_mb(heights, year=y, fl_id=0)
            return mb 
        end
    
        let 
        #bed_h = collect(LinRange(top_h, bottom_h, nx))
        surface_h = copy(surface_ini)
        thick = surface_h .- bed_h
    
        t = 0.0
        dt = sec_in_day * 10.0
    
        years = collect(0:(n_years+1))
        volume = zeros(size(years))
        long = zeros(size(years))
    
        new_thick = zeros(nx)
        diffusivity = zeros(nx)
        diffusivity_s = zeros(nx)
        surface_gradient = zeros(nx)
        surface_gradient_s = zeros(nx)
        grad_x_diff = zeros(nx)
        flux_div = zeros(nx)
        mb = zeros(nx-1)
        yy=2004
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
                diff2!(surface_gradient, surface_h, 2.0*dx)
    
                # Diffusivity
                diffusivity .= 2 .*glen_a .*width .* ((ρ*g).^glen_n) ./(glen_n .+2) 

                diffusivity .*= (surface_gradient.^(glen_n .-1)) .*(thick.^(glen_n .+2))

                # Ice flux in a staggered grid
                avg!(diffusivity_s, diffusivity)
    
                diff1!(surface_gradient_s, surface_h, dx)
    
                grad_x_diff .= surface_gradient_s .* diffusivity_s
                diff1!(flux_div, grad_x_diff, dx)
    
                # Mass balance
                mb .= get_mb2(surface_h[begin:end-1],yy)
                # Ice thickness update: old + flux div + mb
                new_thick[begin:end-1] .= thick[begin:end-1] .+ (dt./width[2:end]) .* flux_div[2:end] .+ dt.*mb
    
                # We can have negative thickness because of MB - correct here
                thick .= ifelse.(new_thick.<0.0, 0.0, new_thick)
                
                @assert thick[end] == 0.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"
    
                # Prepare for next step 
                surface_h .= bed_h .+ thick
                t += dt

                annee = t ÷ sec_in_year
                yy=2004+annee

            end
            end # let
    
            volume[i] = sum(thick .* width .* dx)
            long[i] = sum(thick .> 0.0) .* dx
        end
    
        # xcoordinates
        xc = collect(0:nx-1) .* dx
    
        return xc, bed_h, surface_h, years, volume, long
        end # let

       
        
    end
    
    function wrapper_grad(grad)
        return glacier_evolution_optim(mb_grad=grad)
    end
    
    # end # @everywhere
    
    
    #######  MAIN ########
    #x=ncread("saved_11-00897.nc","bed_h")
    # Let's test the performance

    #=
    xc, bed_h, surface_h, years, volume, long = @time glacier_evolution_optim(nx=length(x),bed_h=x)
    =#
    # Plot the results of the simulation
    #=
    plot(bed_h, color="black", title="Glacier geometry at the end of the simulation", label="Bedrock", ylabel="Elevation (m.a.s.l.)")
    p_flowline = plot!(surface_h, color="slateblue", label="Ice")
    display(p_flowline)
    =#

    
    
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
    
    
    
