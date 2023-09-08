
#For differential computation in the numerical scheme
diff_x!(O, I, Δx) = @views @. O = (I[2:end] - I[1:end - 1])./ Δx
avg!(O, I) = @views @. O = (I[2:end] + I[1:end - 1])./2.0

include("SIA1D.jl")
using LinearSolve


"""
SIA1D!(dH, H, SIA1Dmodel,t)

Compute an in-place step of the Shallow Ice Approximation PDE in a forward model
"""

function SIA1D!(dH, H, iceflow_model, t) #where {F <: AbstractFloat, I <: Integer}
    # Retrieve model parameters
    Δx = iceflow_model.Δx
    Γ = iceflow_model.Γ
    Γs = iceflow_model.Γs
    #dSdx = iceflow_model.dSdx
    #∇Sx = iceflow_model.∇Sx
    w = iceflow_model.w 
    #w̅ = iceflow_model.w̅
    w0 = iceflow_model.w0
    w0_stag = iceflow_model.w0_stag
    λ = iceflow_model.λ
    #S = iceflow_model.S
    #Sexp = iceflow_model.Sexp
    B = iceflow_model.B
    #H̅ = iceflow_model.H̅
    D = iceflow_model.D
    #Fx = iceflow_model.Fx 
    #Fxx = iceflow_model.Fxx


    # First, enforce values to be positive
    map!(x -> ifelse(x>0.0,x,0.0), H, H)
    
    S_tmp = get_tmp(iceflow_model.S_tmp,H)
    @. S_tmp = B + H

    w_tmp = get_tmp(iceflow_model.w_tmp,H)
    @. w_tmp = w0 + (λ* H)

    H̅_tmp = get_tmp(iceflow_model.H̅_tmp,H)
    #Staggered thickness
    avg!(H̅_tmp,H)

    w̅_tmp = get_tmp(iceflow_model.w̅_tmp,H)
    #Staggered width
    avg!(w̅_tmp,w_tmp)

    dSdx_tmp = get_tmp(iceflow_model.dSdx_tmp, H)
    #(Staggered) slope
    diff_x!(dSdx_tmp,S_tmp,Δx)

    # Diffusivity
    D_tmp = get_tmp(iceflow_model.D_tmp,H)
    @views @. D_tmp[2:end-1] = ((dSdx_tmp^2.0)^((n-1)/2.0)) * ((w̅_tmp + w0_stag)/2.0) * (Γ * H̅_tmp^(n.+2) + Γs*(H̅_tmp^n))
    D_tmp[1] = D[1] 
    D_tmp[end] = D[end] 

    Sexp_tmp = get_tmp(iceflow_model.Sexp_tmp,H)
    @views @. Sexp_tmp[2:end-1] = S_tmp
    Sexp_tmp[1]=S_tmp[1]
    Sexp_tmp[end]=S_tmp[end]

    ∇Sx_tmp = get_tmp(iceflow_model.∇Sx_tmp, H) 
    diff_x!(∇Sx_tmp,Sexp_tmp,Δx)

    Fx_tmp = get_tmp(iceflow_model.Fx_tmp, H)
    
    @. Fx_tmp = -D_tmp * ∇Sx_tmp 

    Fxx_tmp = get_tmp(iceflow_model.Fxx_tmp, H)
    diff_x!(Fxx_tmp,Fx_tmp,Δx)

    dH .=  .-Fxx_tmp./w

    # Clip negative ice thickness values
    #@views H[H.<0.0] .= 0.0
    # @assert H[end] .< 10.0 "Glacier exceeding boundaries! at time $(t/sec_in_year)"

end

"""
glacier_evolution_store(SIA1D, solver, reltol, y0, y1, l_var, do_geom, do_fl_diag, mb_step, dyn_spin_thick)

Runs the flowline dynamics up to given year date y1
This function runs the model for the time difference y1-self.y0 and stores diagnostic variables to give back to OGGM
"""

function glacier_evolution_store(;SIA1D::Py, 
                                solver, 
                                reltol::F,
                                y0::I, 
                                y1::I, 
                                l_var::PyList{Any}, 
                                do_geom::I, 
                                do_fl_diag::I,
                                mb_step::String,
                                dyn_spin_thick::F) where {F <: AbstractFloat, I <: Integer}

    # Initialize SIA1D model with variables from current flowline 
    iceflow_model::SIA1Dmodel = SIA1Dmodel()
    initialize_iceflow_model!(iceflow_model, SIA1D)

    #Time parameters 
    sim_params::SimulationParameters = SimulationParameters()
    initialize_sim_params!(sim_params, SIA1D, y0, y1, mb_step)

    
    #### Initialize necessary storage structs     
    #Diagnostic
    diag_jl::diag = diag()
    initialize_diag!(diag_jl,sim_params.nₜm,sim_params.fact_m,
                    sim_params.diff, SIA1D,l_var,dyn_spin_thick)

    ##Geometry and flowline diagnostic
    do_geom = Bool(do_geom)
    do_fl_diag = Bool(do_fl_diag)
    is_tidewater = pyconvert(Bool,SIA1D.is_tidewater)
    geom_var::geom = geom()
    fl_var::fl_diag = fl_diag()
    initialize_geom_diag!(geom_var,fl_var,do_geom,do_fl_diag,
                        sim_params.nₜ,iceflow_model.nx, SIA1D,sim_params.diff, is_tidewater)


    #Callback functions for the mass balance and update the SIA1D object 
    tstops, _ = define_callback_steps(sim_params.tspan)
    stop_condition(u,t,integrator) = stop_condition_tstops(u,t,integrator, tstops) #closure

    #To be called during the callbacks
    function action!(integrator)
        #Adding the mass balance at each time stops
        sim_params.yr = sim_params.y₀slf + ((integrator.t - sec_in_month)/sec_in_year)
        
        get_mb!(iceflow_model,integrator.u,sim_params, SIA1D)
            
        integrator.u .+= (iceflow_model.mb .* sec_in_month)
        integrator.u[integrator.u.<0.0] .= 0.0

        #writing the thickness solution every year as a property of the SIA1D object
        setproperty!(SIA1D, :t, integrator.t)
        setproperty!(SIA1D.fls[0], :thick, np.array(integrator.u))

        #Yearly Index
        tᵢ = floor(Int,integrator.t ./ sec_in_year) + 1 
        #Monthly index
        tᵢm = floor(Int,(integrator.t ./sec_in_year)./(1/12)) + 1

        if sim_params.yr >= sim_params.y₀
            if sim_params.mb_step == "annual" 
                if tᵢ > 1 && tᵢm%12 == 1 #first year is already saved
                    tᵢ = tᵢ + sim_params.diff
                    #Writing the diagnostic variables
                    update_diag!(diag_jl,SIA1D,l_var,dyn_spin_thick,tᵢ)    

                    #Writing for the geometry (and flowline diagnostic)
                    if do_geom || do_fl_diag
                        update_geom!(geom_var,SIA1D,is_tidewater,tᵢ)
                        if do_fl_diag
                            update_fl_diag!(fl_var,SIA1D,tᵢ)
                        end 
                    end 
                end 
            else #store monthly steps
                #Writing the diagnostic variables
                update_diag!(diag_jl,SIA1D,l_var,dyn_spin_thick,tᵢm+(sim_params.diff*12))
                if tᵢm%12 == 1
                    tᵢ = tᵢ + sim_params.diff
                    #Writing for the geometry (and flowline diagnostic) only once a year 
                    if do_geom || do_fl_diag
                        update_geom!(geom_var,SIA1D,is_tidewater,tᵢ)
                        if do_fl_diag
                            update_fl_diag!(fl_var,SIA1D,tᵢ)
                        end 
                    end    
                end 

            end 
        end 


    end   

    #Defining the callback
    cb_MB = DiscreteCallback(stop_condition, action!)

     
    u0= iceflow_model.H
    du0 = copy(u0)
    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> SIA1D!(du, u, iceflow_model, 0.0), du0,u0)

    f = ODEFunction(SIA1D!; jac_prototype = float.(jac_sparsity))
    #Solving the problem         
    iceflow_prob = ODEProblem(f,iceflow_model.H,sim_params.tspan,tstops=tstops,iceflow_model)

    solve(iceflow_prob, solver,callback=cb_MB, tstops=tstops, reltol=reltol,save_everystep=false,dense=false)

    return SIA1D,diag_jl, geom_var, fl_var
    
end 



#Basically like glacier_evolution_store but without storing anything
function glacier_evolution(;
    SIA1D::Py, #IceflowJuliaModel object
    solver, #Solver to use
    reltol::F, #Relative tolerance
    y0::I, #starting year of the simulation
    y1::I) where {F <: AbstractFloat, I <: Integer} #ending year of the simulation

    mb_step::String ="annual"

    # Initialize SIA1D model with variables from current flowline 
    iceflow_model::SIA1Dmodel = SIA1Dmodel()
    initialize_iceflow_model!(iceflow_model, SIA1D)
    
    #Time parameters 
    sim_params::SimulationParameters = SimulationParameters()
    initialize_sim_params!(sim_params, SIA1D, y0, y1, mb_step)


    #Callback functions for the mass balance and update the SIA1D object 
    tstops, _ = define_callback_steps(sim_params.tspan)
    stop_condition(u,t,integrator) = stop_condition_tstops(u,t,integrator, tstops) #closure

    #To be called during the callbacks
    function action!(integrator)
        sim_params.yr = sim_params.y₀slf + ((integrator.t - sec_in_month)/sec_in_year)
        get_mb!(iceflow_model,integrator.u,sim_params, SIA1D)
        integrator.u .+= (iceflow_model.mb .* sec_in_month)
        integrator.u[integrator.u.<0.0] .= 0.0

        #writing the thickness solution every year as a property of the SIA1D object
        setproperty!(SIA1D, :t, integrator.t)
        setproperty!(SIA1D.fls[0], :thick, np.array(integrator.u))

    end   

    #Defining the callback
    cb_MB = DiscreteCallback(stop_condition, action!)

    #Solving the problem         
    iceflow_prob = ODEProblem(SIA1D!,iceflow_model.H,sim_params.tspan,tstops=tstops,iceflow_model)
    iceflow_sol = solve(iceflow_prob,solver,callback=cb_MB, tstops=tstops, reltol=reltol,save_everystep=false,dense=false)

    return SIA1D

end


##### CALLBACK FUNCTIONS #####
#=
function get_mb!(mb::Vector{Float64}, heights::Vector{Float64},y::Float64,SIA1D::Py)
    mb .= pyconvert(Vector{Float64},SIA1D.get_mb(heights,y,fl_id=0))     
end
=#

function get_mb!(iceflow_model::SIA1Dmodel{F,I},
                 H::Vector{F}, 
                 sim_params::SimulationParameters{F,I},
                 SIA1D::Py) where {F <: AbstractFloat, I <: Integer}

    iceflow_model.mb .= pyconvert(Vector{F}, SIA1D.get_mb(iceflow_model.B + H, sim_params.yr,fl_id=0))                
end 

function define_callback_steps(tspan::Tuple)
    step = sec_in_month
    tmin_int = Int(tspan[1])
    tmax_int = Int(tspan[2])+1
    tstops = range(tmin_int+step, tmax_int, step=step) |> collect
    tstops = filter(x->( (Int(tspan[1])<x) & (x<=Int(tspan[2])) ), tstops)
    return tstops, step
end

function stop_condition_tstops(u,t,integrator, tstops) 
    t in tstops
end


##############################