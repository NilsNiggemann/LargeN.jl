function M_inverse(T::Real,J̃::AbstractMatrix,λ::Real)
    dims = size(J̃)
    inv(J̃ ./T +λ .* Matrix(I,dims...))
end

function X_q(nb::Integer,Minv::AbstractMatrix)
    real(sum(Minv)/nb)
end

function getNCell(Jfunc::Function)
    try
        return maximum(Jfunc.alpha_vec)
    catch
        error("Function does not contain sufficient information, please provide dimension argument explicitly")
    end
end

getNCell(JMat::AbstractMatrix) = size(JMat)[1]


function X_q(T::Real,J̃::AbstractMatrix,λ::Real,nb::Integer = getNCell(J̃))
    real(sum(M_inverse(T,J̃,λ)))/nb
end

function constraint(J̃::Function,λ::Real,T::Real,BZextent::Real,d::Integer = getDimfromFunc(J̃);kwargs...)
    f(q) = real(tr(M_inverse(T,J̃(q),λ)))
    return 1/(BZextent)^d * BZIntegral(f,d,BZextent;kwargs...) -1
end

function optimizeConstraint(J̃::Function,T::Real,BZextent::Real,d::Integer= getDimfromFunc(J̃);min = 2/T+0.01, max=2/T+10,guess = (min,max),method = Roots.Bisection(),kwargs... )
    f(λ) = constraint(J̃,λ,T,BZextent,d;kwargs...)
    find_zero(f,guess,method;atol = 1e-15, rtol = 1e-15,kwargs...)
end

function BZIntegral(f::Function,d::Integer,BZextent::Real; kwargs...)
    val,err = hcubature(f, zeros(d), BZextent .*ones(d); reltol=1e-10, abstol=1e-10,maxevals =5000,kwargs...)
    return val
end

function X_q_eig(T::Real,J̃::AbstractMatrix,λ::Real,nb::Integer = getNCell(J̃))
    evals, evecs = eigen(J̃)
    sum = 0. +0im
    for alpha in eachindex(evals),beta in eachindex(evals),gamma in eachindex(evals)
        sum += (evecs[alpha,beta] *conj(evecs[gamma,beta] ))/(evals[beta]/T + λ)
    end
    return real(sum/nb)
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

function ComputeEig2D(Jfunc::Function,nk::Integer;ext = 2pi,min = -ext,max = ext)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,NCell,nk,nk)
    vec = Array{ComplexF64}(undef,NCell,NCell,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray)
        sys = eigen(Jfunc(SA[kx,ky,kz]))
        eig[:,i,j] .= sys.values
        vec[:,:,i,j] .= sys.vectors
    end
    return eig,vec
end

function ComputeEig3D(Jfunc::Function,nk::Integer;ext = 2pi,min = -ext,max = ext)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,NCell,nk,nk,nk)
    vec = Array{ComplexF64}(undef,NCell,NCell,nk,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray),(k,kz) in enumerate(karray)
        sys = eigen(Jfunc(SA[kx,ky,kz]))
        eig[:,i,j,k] .= sys.values
        vec[:,:,i,j,k] .= sys.vectors
    end
    return eig,vec
end


function ComputeEigvals2D(Jfunc::Function,nk::Integer,ext;min = -ext,max = ext)
    T = eltype(Jfunc.Jij_vec)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{T}(undef,NCell,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray)
        EV =  real.(eigvals(Jfunc(SA[kx,ky])))
        eig[:,i,j] .= EV
    end
    return eig
end

function ComputeEigvals3D(Jfunc::Function,nk::Integer,ext;min = -ext,max = ext)
    T = eltype(Jfunc.Jij_vec)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{T}(undef,NCell,nk,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray),(k,kz) in enumerate(karray)
        EV = real.(eigvals(real.(Jfunc(SA[kx,ky,kz]))))
        eig[:,i,j,k] .= real(EV)
    end
    return eig
end


@inline function constraint(eval::AbstractArray,Lam::Real,T::FloatType) where FloatType
    s = zero(FloatType)
    for e in eval
        s += 1/(e/T+Lam)
    end
    s/length(eval) - one(FloatType)
end

function optimizeConstraint(eval::AbstractArray,T::Real;guess =3/T,kwargs...)
    @inline constr(Lam) = constraint(eval,Lam,T)
    Lambda = find_zero(constr,guess,atol = 1e-15, rtol = 1e-15;kwargs...)
end

function getEvals(JFunc,BZextent = 4pi;nk = 20)
    Dim = getDimfromFunc(JFunc)
    evalFunc = (nothing,ComputeEigvals2D,ComputeEigvals3D)[Dim]
    EV= evalFunc(JFunc,nk,BZextent)
end

function analyzeSpectrum(T,Sys::Geometry,Basis::Basis_Struct,pairToInequiv::Function;BZextent = 4pi,nk = 20,verbose = true)
    JFunc = constructJ(Sys,Basis,pairToInequiv)
    EV = getEvals(JFunc,BZextent,nk = nk)

    Emin = minimum(EV)
    degeneracy = length(filter(x -> isapprox(x,Emin,atol = 1e-8),EV))#/length(EV)
    
    NCell = getNCell(JFunc)
    verbose && println("ground state degeneracy: $degeneracy. Degeneracy per unit cell: $(degeneracy*NCell/length(EV))")

    LamSing = -Emin/T
    verbose && println("-Emin/T: $LamSing")

    LeftBound = LamSing
    RightBound = LamSing+10+T
    return (JFunc = JFunc, EV = EV, LamSing = LamSing, LeftBound = LeftBound, RightBound = RightBound,constraint = x->constraint(EV,x,T))
end

analyzeSpectrum(T,Sys::Geometry,Mod::Module;kwargs...) = analyzeSpectrum(T,Sys,Mod.Basis,Mod.pairToInequiv;kwargs...)

function getChiFunction(T,Sys::Geometry,Basis::Basis_Struct,pairToInequiv::Function;BZextent = 4pi,nk = 20,tol = 1e-6,verbose = true,kwargs...)
    sinfo = analyzeSpectrum(T,Sys,Basis,pairToInequiv;BZextent = BZextent,nk = nk,verbose = verbose)

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
        cons = constraint(EV,Lam,T)
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

    @inline chi(q::AbstractVector)= X_q(T,JFunc(q),Lam)
    @inline chi(q...)= X_q(T,JFunc(SVector(q)),Lam)

    return chi
end
getChiFunction(T,Sys::Geometry,Mod::Module;kwargs...) = getChiFunction(T,Sys,Mod.Basis,Mod.pairToInequiv;kwargs...)

function getZeroTChi(Sys::Geometry,Basis::Basis_Struct,pairToInequiv::Function)
    JFunc = constructJ(Sys,Basis,pairToInequiv)
    chi(x::AbstractVector) = zeroTweight(JFunc(x))
    chi(x...) = zeroTweight(JFunc(SVector(x)))
end
getZeroTChi(Sys::Geometry,Mod::Module) = getZeroTChi(Sys,Mod.Basis,Mod.pairToInequiv)

##
