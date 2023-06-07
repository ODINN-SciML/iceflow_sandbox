#%%
# -*- coding: utf-8 -*-

###################################################################
#### Solving the SIA equation along the glacier flowline with  ####
#### a Julia solver (function 'glacier_evolution'). The actual ####
#### coding of the SIA equation is in 'iceflow!' function      ####
###################################################################

#Constants definition
const sec_in_day = 24.0 * 60.0 * 60.0
const sec_in_year = sec_in_day * 365.0
const glen_n = 3.0
const ρ = 900.0
const g = 9.81

diff1!(O, I, dx) = @views @. O[2:end] = (I[2:end] - I[1:end - 1]) / dx
diff2!(O, I, dx) = @views @. O[2:end-1] = (I[3:end] - I[1:end - 2]) / dx
avg!(O, I) = @views @. O[2:end] = (I[2:end] + I[1:end - 1])./2.0

function glacier_evolution(;
                            gdir::PyObject,
                            dx::Float64,
                            nx::Int64,  # grid size e.g 200
                            width::Vector{Float64},  # glacier width in m
                            glen_a::Float64,  # ice stiffness 2.4e-24
                            solver,
                            reltol::Float64,
                            bed_hs::Vector{Float64},
                            surface_ini::Vector{Float64},
                            y0::Float64,
                            tspan::Tuple{Float64, Float64})
    
    #Getting the accurate mass balance model from OGGM  
    mbmod = massbalance.MultipleFlowlineMassBalance(gdir,mb_model_class=functools.partial(massbalance.RandomMassBalance,
                            mb_model_class=massbalance.MonthlyTIModel),
                            y0=y0, halfsize=15,bias=0, seed=1,
                            filename="climate_historical",input_filesuffix="",unique_samples=false)
    

    #H at t_0
    H0 = (surface_ini .- bed_hs)
    
    #pre-allocation
    surface = zeros(Float64,nx)
    surface_gradient = zeros(Float64,nx)
    surface_gradient_s = zeros(Float64,nx)

    diffusivity = zeros(Float64,nx)
    diffusivity_s = zeros(Float64,nx)

    grad_x_diff = zeros(Float64,nx)
    flux_div = zeros(Float64,nx)

    mb=zeros(Float64,nx)

    Γ = (2 .* glen_a .* width .* (ρ * g).^glen_n)/(glen_n .+ 2)

    p = (dx, width, bed_hs, surface_gradient,surface_gradient_s,
            diffusivity,diffusivity_s,grad_x_diff, Γ,surface, flux_div)

    #Callback functions for the mas balance 
    tstops, _ = define_callback_steps(tspan)
    stop_condition(u,t,integrator) = stop_condition_tstops(u,t,integrator, tstops) #closure
    
    #the function to be called buring the callbacks
    function action!(integrator)
        y = y0 .+ (integrator.t ./ sec_in_year)
        get_mb!(mb, bed_hs .+ integrator.u,y,mbmod)
        integrator.u .+= mb
        integrator.u[integrator.u.<0.0] .= 0.0

    end        
    
    #Defining the callback
    cb_MB = DiscreteCallback(stop_condition, action!)

    #Solving the problem         
    iceflow_prob = ODEProblem(iceflow!,H0,tspan,tstops=tstops,p)
    iceflow_sol = solve(iceflow_prob,solver,callback=cb_MB, tstops=tstops,  reltol=reltol,save_everystep=false,dense=false)

    return iceflow_sol

end

function get_mb!(mb, heights,y,mbmod)
    mb .= (mbmod.get_annual_mb(heights, year=Int(y), fl_id=0)) .*sec_in_year
end

function define_callback_steps(tspan; step=sec_in_year)
    tmin_int = Int(tspan[1])
    tmax_int = Int(tspan[2])+1
    tstops = range(tmin_int+step, tmax_int, step=step) |> collect
    tstops = filter(x->( (Int(tspan[1])<x) & (x<=Int(tspan[2])) ), tstops)
    return tstops, step
end

function stop_condition_tstops(u,t,integrator, tstops) 
    t in tstops
end

function iceflow!(dH, H, p, t)
    # Retrieve model parameters
   dx::Float64, width::Vector{Float64}, bed_hs::Vector{Float64},surface_gradient::Vector{Float64}, surface_gradient_s::Vector{Float64},
    diffusivity::Vector{Float64},diffusivity_s::Vector{Float64}, grad_x_diff::Vector{Float64}, Γ::Vector{Float64}, surface::Vector{Float64} ,
    flux_div::Vector{Float64} = p

    surface .= bed_hs .+ H 

    # Clip negative ice thickness values
    @views H[H.<0.0] .= 0.0
    @assert H[end-2] .== 0.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"

    # Surface gradient
    diff2!(surface_gradient, surface, 2.0*dx)

    # Diffusivity
    diffusivity .= Γ .* H.^(glen_n.+2) .* surface_gradient.^(glen_n.-1)

    # Ice flux in a staggered grid
    avg!(diffusivity_s,diffusivity)

    diff1!(surface_gradient_s, surface, dx)

    grad_x_diff .= surface_gradient_s .* diffusivity_s

    diff1!(flux_div, grad_x_diff, dx)

    
    # # Ice thickness update: old + flux div + mb
    #dH .= zeros(size(dH))
    # @show median(mb)
    dH[begin:end-1] .= (flux_div[2:end] ./ width[2:end]) 


end
