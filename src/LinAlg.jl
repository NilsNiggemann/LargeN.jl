function M_inverse(T::Real,J̃::AbstractMatrix,λ::Real)
    inv(J̃ ./T +λ *I)
end

function getNCell(Jfunc::Function)
    try
        return maximum(Jfunc.alpha_vec)
    catch
        error("Function does not contain sufficient information, please provide dimension argument explicitly")
    end
end

getNCell(JMat::AbstractMatrix) = size(JMat)[1]


function X_q(T::Real,J̃::AbstractMatrix,λ::Real)
    real(sum(M_inverse(T,J̃,λ)))
end

function constraint_cont(J̃::Function,λ::Real,T::Real,BZextent::Real,d::Integer = getDimfromFunc(J̃);kwargs...)
    f(q) = real(tr(M_inverse(T,J̃(q),λ)))
    return 1/(BZextent)^d * BZIntegral(f,d,BZextent;kwargs...) -1/3
end

function optimizeConstraint(J̃::Function,T::Real,BZextent::Real,d::Integer= getDimfromFunc(J̃);min = 2/T+0.01, max=2/T+10,guess = (min,max),method = Roots.Bisection(),kwargs... )
    f(λ) = constraint_cont(J̃,λ,T,BZextent,d;kwargs...)
    find_zero(f,guess,method;atol = 1e-15, rtol = 1e-15,kwargs...)
end

function BZIntegral(f::Function,d::Integer,BZextent::Real; kwargs...)
    val,err = hcubature(f, zeros(d), BZextent .*ones(d); reltol=1e-10, abstol=1e-10,maxevals =5000,kwargs...)
    return val
end

function X_q_eig(T::Real,J̃::AbstractMatrix,λ::Real)
    evals, evecs = eigen(J̃)
    sum = 0. +0im
    for alpha in eachindex(evals),beta in eachindex(evals),gamma in eachindex(evals)
        sum += (evecs[alpha,beta] *conj(evecs[gamma,beta] ))/(evals[beta]/T + λ)
    end
    return real(sum)
end

function zeroTweight(J̃::AbstractMatrix)
    evals, evecs = eigen(J̃)
    Emin = minimum(evals)
    minEnergies = findall(x -> isapprox(x,Emin,atol = 1e-12),evals)

    weight = zero(real(eltype(evecs)))
    for beta in minEnergies
        @views weight += abs2(sum(evecs[:,beta]))
    end
    weight
end

function ComputeEig2D(Jfunc::Function,nk::Integer,NCell = getNCell(Jfunc);ext = 2pi,min = -ext,max = ext,)
    karray = range(min,max,length = nk)
    values = Array{Float64}(undef,NCell,nk,nk)
    vectors = Array{ComplexF64}(undef,NCell,NCell,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray)
        sys = LargeN.eigen(Jfunc(SA[kx,ky]))
        values[:,i,j] .= sys.values
        vectors[:,:,i,j] .= sys.vectors
    end
    return (;values,vectors)
end

function ComputeEig3D(Jfunc::Function,nk::Integer,NCell = getNCell(Jfunc);ext = 2pi,min = -ext,max = ext)
    karray = range(min,max,length = nk)
    eig = Array{Float64}(undef,NCell,nk,nk,nk)
    vec = Array{ComplexF64}(undef,NCell,NCell,nk,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray),(k,kz) in enumerate(karray)
        sys = eigen(Jfunc(SA[kx,ky,kz]))
        eig[:,i,j,k] .= sys.values
        vec[:,:,i,j,k] .= sys.vectors
    end
    return eig,vec
end


function ComputeEigvals2D(Jfunc::Function,nk::Integer,ext::T;min = -ext,max = ext) where T <: AbstractFloat
    karray = range(min,max,length = nk)
    k0 = first(karray)
    NCell = size(Jfunc(SVector(k0,k0)),1)

    eig = Array{T}(undef,NCell,nk,nk)

    Threads.@threads for i in 1:nk
        JArr = Array{Complex{T}}(undef,NCell,NCell) # eigvals! allocates much more anyway so this is fine
        kx = karray[i]
        for j in 1:nk
            ky = karray[j]
            JArr .= Jfunc(SA[kx,ky])
            EV = eigvals!(JArr)
            eig[:,i,j] .= EV
        end
    end
    return eig
end

function ComputeEigvals3D(Jfunc::Function,nk::Integer,ext::T;min = -ext,max = ext) where T <: AbstractFloat
    karray = range(min,max,length = nk)
    k0 = first(karray)
    NCell = size(Jfunc(SVector(k0,k0,k0)),1)
    eig = Array{T}(undef,NCell,nk,nk,nk)

    Threads.@threads for i in 1:nk
        JArr = Array{Complex{T}}(undef,NCell,NCell) # eigvals! allocates much more anyway so this is fine
        kx = karray[i]
        for j in 1:nk
            ky = karray[j]
            for k in 1:nk
                kz = karray[k]
                JArr .= Jfunc(SA[kx,ky,kz])
                EV = eigvals!(JArr)
                eig[:,i,j,k] .= EV
            end
        end
    end
    return eig
end


@inline function constraint(eval::AbstractArray,Lam::Real,T::FloatType) where FloatType
    s = zero(FloatType)
    LamT = Lam*T
    for e in eval
        s += 1/(e+LamT)
    end
    T*s/length(eval) - 1/3*one(FloatType)
end

function optimizeConstraint(eval::AbstractArray,T::Real;guess =3/T,kwargs...)
    @inline constr(Lam) = constraint(eval,Lam,T)
    Lambda = find_zero(constr,guess,atol = 1e-15, rtol = 1e-15;kwargs...)
end

function getEvals(JFunc, Dim ,BZextent = 4pi;nk = 20)
    evalFunc = (nothing,ComputeEigvals2D,ComputeEigvals3D)[Dim]
    EV= evalFunc(JFunc,nk,BZextent)
end

function analyzeSpectrum(T,JFunc::Function, NCell = getNCell(JFunc),Dim = getDimfromFunc(JFunc);BZextent = 4pi,nk = 20,verbose = true)
    EV = getEvals(JFunc,Dim,BZextent,nk = nk)

    Emin = minimum(EV)
    degeneracy = length(filter(x -> isapprox(x,Emin,atol = 1e-8),EV))#/length(EV)
    

    verbose && println("ground state degeneracy: $degeneracy. Degeneracy per unit cell: $(degeneracy*NCell/length(EV))")

    LamSing = -Emin/T
    verbose && println("-Emin/T: $LamSing")

    LeftBound = LamSing
    RightBound = LamSing+10+T
    constr = x->constraint(EV,x,T)
    return (JFunc = JFunc, EV = EV, LamSing = LamSing, LeftBound = LeftBound, RightBound = RightBound,constraint = constr)
end

function analyzeSpectrum(T,Sys::Geometry,Basis::Basis_Struct,pairToInequiv::Function;kwargs...)
    JFunc = constructJ(Sys,Basis,pairToInequiv)
    return analyzeSpectrum(T,JFunc,kwargs...)
end

analyzeSpectrum(T,Sys::Geometry,Mod::Module;kwargs...) = analyzeSpectrum(T,Sys,Mod.Basis,Mod.pairToInequiv;kwargs...)


function getChiFunction(T::FLType,Sys::Geometry,Basis::Basis_Struct,pairToInequiv::Function;kwargs...) where FLType <:AbstractFloat
    JFunc = constructJ(Sys,Basis,pairToInequiv)
    return getChiFunction(T,JFunc;kwargs...)
end

function getChiFunction(T::FLType,JFunc::Function,NCell = getNCell(JFunc),Dim = getDimfromFunc(JFunc);BZextent = 4pi,nk = 20,tol = 1e-6,verbose = true,kwargs...) where FLType <:AbstractFloat
    sinfo = analyzeSpectrum(T,JFunc,NCell,Dim;BZextent = BZextent,nk = nk,verbose = verbose)
    JFunc = sinfo.JFunc
    EV = sinfo.EV
    Lam = sinfo.LamSing
    LeftBound = sinfo.LeftBound
    RightBound = sinfo.RightBound

    cons = 1E16
    try
        # Lam = optimizeConstraint(EV,T,guess=[LeftBound,RightBound];kwargs...)
        Lam = find_zero(sinfo.constraint,[LeftBound,RightBound];kwargs...)
        # Lam = optimizeConstraint(JFunc,T,BZextent,guess=[LamSing+0.1*tol,LamSing+10+T])
        cons = sinfo.constraint(Lam)
        verbose && @info "constraint minimized to $cons for λ = $Lam"
    catch e
        verbose && @warn "constraint could not be optimized!
        (Errored with: $e). Lam set to singularity."
    end

    if abs(cons) > tol
        @warn "T= $T: constraint possibly not converged!
        constraint value: $cons \t tolerance: $tol"
    end
    
    Lam0::FLType = Lam
    T0::FLType = T
    @inline chi(q::AbstractVector)= X_q(T0,JFunc(q),Lam0)
    @inline chi(qx::Real,qy::Real)= X_q(T0,JFunc(SA[qx,qy]),Lam0)
    @inline chi(qx::Real,qy::Real,qz::Real)= X_q(T0,JFunc(SA[qx,qy,qz]),Lam0)

    return chi
end
getChiFunction(T,Sys::Geometry,Mod::Module;kwargs...) = getChiFunction(T,Sys,Mod.Basis,Mod.pairToInequiv;kwargs...)

function getZeroTChi(Sys::Geometry,Basis::Basis_Struct,pairToInequiv::Function)
    JFunc = constructJ(Sys,Basis,pairToInequiv)
    chi(x::AbstractVector) = zeroTweight(JFunc(x))
    chi(x...) = zeroTweight(JFunc(SVector(x)))
end

function getZeroTChi(JFunc::Function)
    chi(x::AbstractVector) = zeroTweight(JFunc(x))
    chi(x...) = zeroTweight(JFunc(SVector(x)))
end
getZeroTChi(Sys::Geometry,Mod::Module) = getZeroTChi(Sys,Mod.Basis,Mod.pairToInequiv)

##
