np=pyimport("numpy")
using PreallocationTools

#Constants definition
const sec_in_day = 24.0 * 60.0 * 60.0
const sec_in_year = sec_in_day * 365.0
const sec_in_month = 2628000.0
const n = 3.0 #glen n
const ρ = 900.0
const g = 9.81

###############################################
###### SHALLOW ICE APPROXIMATION MODELS #######
###############################################


mutable struct SIA1Dmodel{F <: AbstractFloat, I <: Integer}
    nx::Union{I,Nothing}
    A::Union{F, Nothing}
    Γ::Union{F, Nothing}
    Γs::Union{F, Nothing}
    Δx::Union{I, Nothing}
    H::Union{Vector{F}, Nothing}
    H̅::Union{Vector{F}, Nothing}
    H̅_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    w::Union{Vector{F}, Nothing}
    w_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    w̅ ::Union{Vector{F}, Nothing}
    w̅_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    w0::Union{Vector{F}, Nothing}
    w0_stag::Union{Vector{F}, Nothing}
    λ::Union{Vector{F}, Nothing}
    S::Union{Vector{F}, Nothing} #Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    S_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    Sexp::Union{Vector{F}, Nothing}
    Sexp_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    B::Union{Vector{F}, Nothing} #Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    dSdx::Union{Vector{F}, Nothing}
    dSdx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    ∇Sx::Union{Vector{F}, Nothing}
    ∇Sx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    D::Union{Vector{F}, Nothing}
    D_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    Fx::Union{Vector{F}, Nothing}
    Fx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    Fxx::Union{Vector{F}, Nothing}
    Fxx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing}
    mb::Union{Vector{F}, Nothing}
end 


function SIA1Dmodel(
    nx::Union{I,Nothing} = nothing,
    A::Union{F, Nothing} = nothing,
    Γ::Union{F, Nothing} = nothing,
    Γs::Union{F, Nothing} = nothing,
    Δx::Union{I, Nothing} = nothing,
    H::Union{Vector{F}, Nothing} = nothing, 
    H̅::Union{Vector{F}, Nothing} = nothing,
    H̅_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    w::Union{Vector{F}, Nothing} = nothing,
    w_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    w̅ ::Union{Vector{F}, Nothing} = nothing, 
    w̅_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    w0::Union{Vector{F}, Nothing} = nothing, 
    w0_stag::Union{Vector{F}, Nothing} = nothing, 
    λ::Union{Vector{F}, Nothing} = nothing,
    S::Union{Vector{F}, Nothing} = nothing,
    S_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    Sexp::Union{Vector{F}, Nothing}=nothing, 
    Sexp_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    B::Union{Vector{F}, Nothing} = nothing, 
    dSdx::Union{Vector{F}, Nothing} = nothing, 
    dSdx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    ∇Sx::Union{Vector{F}, Nothing} = nothing, 
    ∇Sx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    D::Union{Vector{F}, Nothing} = nothing, 
    D_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    Fx::Union{Vector{F}, Nothing}= nothing, 
    Fx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    Fxx::Union{Vector{F}, Nothing}= nothing,
    Fxx_tmp::Union{DiffCache{Vector{F},Vector{F}}, Nothing} = nothing,
    mb::Union{Vector{F}, Nothing}= nothing) where {F <: AbstractFloat, I <: Integer}

    SIA1D_model = SIA1Dmodel{Float64,Int64}(nx, A, Γ, Γs, Δx, H, H̅, H̅_tmp, w, w_tmp, w̅, w̅_tmp, w0, w0_stag,
     λ, S, S_tmp, Sexp, Sexp_tmp, B, dSdx, dSdx_tmp, ∇Sx, ∇Sx_tmp, D, D_tmp, Fx, Fx_tmp, Fxx, Fxx_tmp, mb)

    return SIA1D_model
end

function initialize_iceflow_model!(iceflow_model::SIA1Dmodel{F, I}, SIA1D::Py) where {F <: AbstractFloat, I <: Integer}
    fl = SIA1D.fls[0]
    nx = pyconvert(I,fl.nx)
    iceflow_model.nx = pyconvert(I,fl.nx)
    iceflow_model.A = pyconvert(F,SIA1D.glen_a)
    iceflow_model.Γs = pyconvert(F,SIA1D.fs)*((ρ .* g).^n)
    iceflow_model.Γ = ((2.0 * iceflow_model.A * ((ρ .* g).^n) )/(n + 2.0))
    iceflow_model.Δx = pyconvert(I,fl.dx_meter)
    iceflow_model.H = pyconvert(Vector{F},fl.thick)
    iceflow_model.H̅ = zeros(F,nx-1)
    iceflow_model.H̅_tmp = DiffCache(iceflow_model.H̅)
    iceflow_model.w = pyconvert(Vector{F},fl.widths_m)
    iceflow_model.w_tmp = DiffCache(iceflow_model.w)
    iceflow_model.w̅ = zeros(F, nx-1)
    iceflow_model.w̅_tmp = DiffCache(iceflow_model.w̅)
    iceflow_model.w0 = pyconvert(Vector{F},fl._w0_m)
    iceflow_model.w0_stag = zeros(F, nx-1)
    avg!(iceflow_model.w0_stag, iceflow_model.w0)
    iceflow_model.λ = pyconvert(Vector{F},fl._lambdas)
    iceflow_model.S = pyconvert(Vector{F},fl.surface_h)
    iceflow_model.S_tmp = DiffCache(iceflow_model.S)
    iceflow_model.Sexp = zeros(F,nx + 2)
    iceflow_model.Sexp_tmp = DiffCache(iceflow_model.Sexp)
    iceflow_model.B = pyconvert(Vector{F},fl.bed_h)
    iceflow_model.dSdx = zeros(F, nx-1)
    iceflow_model.dSdx_tmp = DiffCache(iceflow_model.dSdx)
    iceflow_model.∇Sx = zeros(F, nx + 1)
    iceflow_model.∇Sx_tmp = DiffCache(iceflow_model.∇Sx) 
    iceflow_model.D = zeros(F, nx + 1)
    iceflow_model.D_tmp = DiffCache(iceflow_model.D)
    iceflow_model.Fx = zeros(F,nx + 1)
    iceflow_model.Fx_tmp = DiffCache(iceflow_model.Fx)
    iceflow_model.Fxx = zeros(F,nx)
    iceflow_model.Fxx_tmp = DiffCache(iceflow_model.Fxx)
    iceflow_model.mb = zeros(F,nx)
end

###############################################
######    VARIABLE STORAGE FOR OGGM     #######
###############################################

mutable struct diag{F <: AbstractFloat}
    volume_m3::Union{Vector{Union{F,Py}}, Nothing}
    area_m2::Union{Vector{Union{F,Py}}, Nothing}
    length_m::Union{Vector{Union{F,Py}}, Nothing}
    calving_m3::Union{Vector{Union{F,Py}}, Nothing}
    calving_rate_myr::Union{Vector{Union{F,Py}}, Nothing}
    volume_bsl_m3::Union{Vector{Union{F,Py}}, Nothing}
    volume_bwl_m3::Union{Vector{Union{F,Py}}, Nothing}
    area_m2_min_h::Union{Vector{Union{F,Py}}, Nothing}
end

mutable struct geom{F <: AbstractFloat}
    sects::Union{Vector{Matrix{Union{F,Py}}}, Nothing}
    widths::Union{Vector{Matrix{Union{F,Py}}}, Nothing}
    buckets::Union{Vector{Vector{Union{F,Py}}}, Nothing} 
end 

mutable struct fl_diag{F <: AbstractFloat, I <: Integer}
    thick_fl::Union{Vector{Matrix{Union{F,Py}}}, Nothing}
    volume_bsl_fl::Union{Vector{Vector{Union{F,Py}}}, Nothing} 
    volume_bwl_fl::Union{Vector{Vector{Union{F,Py}}}, Nothing} 
end 

function diag(volume_m3::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    area_m2::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    length_m::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    calving_m3::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    calving_rate_myr::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    volume_bsl_m3::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    volume_bwl_m3::Union{Vector{Union{F,Py}}, Nothing} = nothing,
    area_m2_min_h::Union{Vector{Union{F,Py}}, Nothing} = nothing) where {F <: AbstractFloat}

    diag_jl = diag{Float64}(volume_m3,area_m2,length_m,calving_m3,calving_rate_myr,
                                volume_bsl_m3, volume_bwl_m3, area_m2_min_h)
    return diag_jl

end 

function initialize_diag!(diag_jl::diag{F},nₜm::I,fac::I,
                            diff::I, SIA1D::Py,l_var::PyList{Any},
                            dyn_spin_thick::F) where {F <: AbstractFloat, I <: Integer}
    
    diag_jl.volume_m3 = Vector{F}(zeros(nₜm))
    diag_jl.area_m2 = Vector{F}(zeros(nₜm))
    diag_jl.length_m = Vector{F}(zeros(nₜm))
    diag_jl.calving_m3 = Vector{F}(zeros(nₜm))
    diag_jl.calving_rate_myr = Vector{F}(zeros(nₜm))
    diag_jl.volume_bsl_m3 = Vector{F}(zeros(nₜm))
    diag_jl.volume_bwl_m3 = Vector{F}(zeros(nₜm))
    diag_jl.area_m2_min_h = Vector{F}(zeros(nₜm))

    if diff >= 0               
        #Diagnostic at year 0 (and more if there is a spinup) 
        for i in 1:(diff*fac)+1
            update_diag!(diag_jl,SIA1D,l_var,dyn_spin_thick,i)
        end 
    end

    return diag_jl
end  

function update_diag!(diag_jl::diag{F},SIA1D::Py,
    l_var::PyList{Any},dyn_spin_thick::F,j::I) where {F <: AbstractFloat, I <: Integer}
    diag_jl.volume_m3[j] = SIA1D.volume_m3
    diag_jl.area_m2[j] = SIA1D.area_m2
    diag_jl.length_m[j] = SIA1D.length_m
    diag_jl.calving_m3[j] = SIA1D.calving_m3_since_y0
    diag_jl.calving_rate_myr[j] = SIA1D.calving_rate_myr
    diag_jl.volume_bsl_m3[j] = SIA1D.volume_bsl_m3
    diag_jl.volume_bwl_m3[j] = SIA1D.volume_bwl_m3
    if "area_min_h" in l_var
        diag_jl.area_m2_min_h[j] = sum([sum(fl.bin_area_m2[fl.thick > dyn_spin_thick]) for fl in SIA1D.fls])
    end 
end 


function geom(sects::Union{Vector{Matrix{Union{F,Py}}}, Nothing} = nothing,
    widths::Union{Vector{Matrix{Union{F,Py}}}, Nothing} = nothing ,
    buckets::Union{Vector{Vector{Union{F,Py}}}, Nothing}= nothing) where {F <: AbstractFloat}

    geom_var= geom{Float64}(sects, widths, buckets)

    return geom_var

end

function fl_diag(thick_fl::Union{Vector{Matrix{Union{F,Py}}}, Nothing} = nothing,
    volume_bsl_fl::Union{Vector{Matrix{Union{F,Py}}}, Nothing} = nothing,
    volume_bwl_fl::Union{Vector{Matrix{Union{F,Py}}}, Nothing} = nothing) where {F <: AbstractFloat, I <: Integer}

    fl_var = fl_diag{Float64,Int64}(thick_fl,volume_bsl_fl,volume_bwl_fl)

    return fl_var

end 

function initialize_geom_diag!(geom_var::geom{F},
                            fl_var::fl_diag{F,I},
                            do_geom::Bool,
                            do_fl_diag::Bool,
                            nₜ::I, nx::I, 
                            SIA1D::Py, diff::I,
                            is_tidewater::Bool) where {F <: AbstractFloat, I <: Integer}

    geom_var.sects = [zeros(F,nₜ,nx) for fl in SIA1D.fls]
    geom_var.widths = [zeros(F,nₜ,nx) for fl in SIA1D.fls]
    geom_var.buckets = [zeros(F,nₜ) for fl in SIA1D.fls]
    
    fl_var.thick_fl = [zeros(F,nₜ,nx) for fl in SIA1D.fls]
    fl_var.volume_bsl_fl = [zeros(F,nₜ) for fl in SIA1D.fls]
    fl_var.volume_bwl_fl = [zeros(F,nₜ) for fl in SIA1D.fls]

    #Storage for geometry
    if do_geom || do_fl_diag
        if diff >= 0               
            for i in 1:diff+1
                update_geom!(geom_var,SIA1D,is_tidewater,i)
            end 
        end 

        if do_fl_diag
            if diff >= 0               
                for i in 1:diff+1
                    update_fl_diag!(fl_var,SIA1D,i)
                end 
            end     

        end 
    end 

    return geom_var, fl_var
end 

function update_geom!(geom_var::geom{F},SIA1D::Py,
    is_tidewater::Bool,j::I) where {F <: AbstractFloat, I <: Integer}
    for (s,w,b,fl) in zip(geom_var.sects,geom_var.widths,geom_var.buckets,SIA1D.fls)
        @views s[j,:] .= fl.section
        @views w[j,:] .= fl.widths_m
        if is_tidewater
            try 
                b[j] = fl.calving_bucket_m3
            catch 
                println("AttributeError")
            end 
        end
    end 
end 

function update_fl_diag!(fl_var::fl_diag{F,I},SIA1D::Py,j::I) where {F <: AbstractFloat, I <: Integer}
    for (t,vs,vw,fl) in zip(fl_var.thick_fl,fl_var.volume_bsl_fl,fl_var.volume_bwl_fl,SIA1D.fls)
        @views t[j,:] .= fl.thick
        vs[j] = fl.volume_bsl_m3
        vw[j] = fl.volume_bwl_m3
    end 
end 

###############################################
######    SIMULATION PARAMETERS         #######
###############################################

mutable struct SimulationParameters{F <: AbstractFloat, I <: Integer} 
    y₀::Union{I,Nothing}
    y₀slf::Union{I,Nothing}
    y₁::Union{I,Nothing}
    yr::Union{F,Nothing}
    nₜ::Union{I, Nothing}
    diff::Union{I, Nothing}
    tspan::Union{Tuple, Nothing}
    mb_step::Union{String, Nothing}
    nₜm::Union{I, Nothing}
    fact_m::Union{I, Nothing}
end

function SimulationParameters(
    y₀::Union{I,Nothing} = nothing, 
    y₀slf::Union{I,Nothing} = nothing, 
    y₁::Union{I,Nothing} = nothing, 
    yr::Union{F,Nothing} = nothing, 
    nₜ::Union{I, Nothing} = nothing,
    diff::Union{I, Nothing} = nothing,
    tspan::Union{Tuple, Nothing} = nothing,
    mb_step::Union{String, Nothing} = nothing,
    nₜm::Union{I, Nothing} = nothing,
    fact_m::Union{I, Nothing} = nothing) where {F <: AbstractFloat, I <: Integer}

    Simulation_Parameters = SimulationParameters{Float64,Int64}(y₀, y₀slf, y₁, yr, nₜ, diff, tspan, mb_step, nₜm, fact_m)
                                                    
    return Simulation_Parameters
end

function initialize_sim_params!(Simulation_Parameters::SimulationParameters{F,I},
                             SIA1D::Py, y0::I, y1::I,
                              mb_step::String) where {F <: AbstractFloat, I <: Integer}
    Simulation_Parameters.y₀ = pyconvert(I,y0)
    Simulation_Parameters.y₀slf = pyconvert(I,SIA1D.y0)
    Simulation_Parameters.y₁ = pyconvert(I,y1)
    Simulation_Parameters.yr = Simulation_Parameters.y₀slf
    Simulation_Parameters.nₜ = Int(y1-y0) + 1
    Simulation_Parameters.diff = Simulation_Parameters.y₀slf - Simulation_Parameters.y₀
    Simulation_Parameters.tspan = (0.0 , (Simulation_Parameters.y₁ - Simulation_Parameters.y₀slf)*sec_in_year)
    Simulation_Parameters.mb_step = mb_step
    if Simulation_Parameters.mb_step =="annual"
        Simulation_Parameters.nₜm = Simulation_Parameters.nₜ
        Simulation_Parameters.fact_m = 1
    else 
        Simulation_Parameters.nₜm = (Simulation_Parameters.nₜ-1)*12 +1
        Simulation_Parameters.fact_m = 12
    end 

end 


