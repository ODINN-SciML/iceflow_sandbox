###################################################################
####  Solving the SIA equation along the glacier flowline with ####
####  a numerical scheme (function 'glacier_evolution_optim') #####
###################################################################



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
    surface_ini, # simulation time in years,
    y0)


    #Mass balance model from OGGM
    mbmod = massbalance.MultipleFlowlineMassBalance(gdir,mb_model_class=functools.partial(massbalance.RandomMassBalance,
    mb_model_class=massbalance.MonthlyTIModel),y0=y0, halfsize=15,bias=0, seed=1,
    filename="climate_historical",input_filesuffix="",unique_samples=false)
 

    function get_mb_r!(mb,heights,y,mbmod)
        mb .= (mbmod.get_annual_mb(heights, year=Int(y), fl_id=0)) .*sec_in_year
    end
    
    let 

    
    #thickness at t=0
    thick = surface_ini .- bed_h
    
    t = 0.0
    dt = sec_in_day * 10.0
    
    #Pre-allocation
    years = collect(0:(n_years+1))
    volume = zeros(size(years))
    long = zeros(size(years))

    surface_h = zeros(nx)
    new_thick = zeros(nx)
    diffusivity = zeros(nx)
    diffusivity_s = zeros(nx)
    surface_gradient = zeros(nx)
    surface_gradient_s = zeros(nx)
    grad_x_diff = zeros(nx)
    flux_div = zeros(nx)
    mb = zeros(nx)
    
    #To add the annual mass balance every year during the loop
    tstops= collect(range(0,n_years,step=1))
        
    for (i, y) in enumerate(years)
        let end_t = y * sec_in_year
        # Time integration
        while t < end_t
            # This is to guarantee a precise arrival on a specific date if asked
            remaining_t = end_t - t
            if remaining_t < dt
                dt = remaining_t
            end

            surface_h .= bed_h .+ thick

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

            #ice thickness update with the annual mass balance
            if t/sec_in_year in tstops 
                y = y0 + t/sec_in_year
                get_mb_r!(mb,surface_h,y,mbmod)
                thick .= thick .+mb   

            end 

            # Ice thickness update: old + flux div + mb     
            new_thick[begin:end-1] .= thick[begin:end-1] .+ (dt./width[2:end]) .* flux_div[2:end]

    
            # We can have negative thickness because of MB - correct here
            thick .= ifelse.(new_thick.<0.0, 0.0, new_thick)
                
            @assert thick[end] == 0.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"
                
            # Prepare for next step 
            t += dt

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

    
