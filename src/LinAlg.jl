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
    f(q) = real(tr(M_inverse(1/T,J̃(q),λ)))
    return 1/(BZextent)^d * BZIntegral(f,d,BZextent;kwargs...) -1
end

function optimizeConstraint(J̃::Function,T::Real,BZextent::Real,d::Integer= getDimfromFunc(J̃);min = 2/T+0.01, max=2/T+10,guess = (min,max),method = Roots.Bisection(),kwargs... )
    f(λ) = constraint(J̃,λ,T,BZextent,d)
    find_zero(f,guess,method;atol = 1e-15, rtol = 1e-15,kwargs...)
end

function BZIntegral(f::Function,d::Integer,BZextent::Real; kwargs...)
    val,err = hcubature(f, zeros(d), BZextent .*ones(d); reltol=1e-10, abstol=1e-10,maxevals =5000,kwargs...)
    return val
end


function ComputeEigvals2D(Jfunc::Function,nk::Integer,ext;min = -ext,max = ext)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,NCell,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray)
        EV = eigvals(Jfunc(SA[kx,ky]))
        eig[:,i,j] .= EV
    end
    return eig
end

function ComputeEigvals3D(Jfunc::Function,nk::Integer,ext;min = -ext,max = ext)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,NCell,nk,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray),(k,kz) in enumerate(karray)
        EV = eigvals(Jfunc(SA[kx,ky,kz]))
        eig[:,i,j,k] .= EV
    end
    return eig
end


@inline function constraint(eval::AbstractArray,Lam::Real,T::Real)
    s = 0.
    for e in eval
        s += 1/(e/T+Lam)
    end
    s/length(eval) -1
end

function optimizeConstraint(eval::AbstractArray,T::Real;guess =3/T)
    @inline constr(Lam) = constraint(eval,Lam,T)
    Lambda = find_zero(constr,guess,atol = 1e-15, rtol = 1e-15)
end

function getEvals(JFunc,BZextent = 4pi;nk = 20)
    Dim = getDimfromFunc(JFunc)
    evalFunc = (nothing,ComputeEigvals2D,ComputeEigvals3D)[Dim]
    EV= evalFunc(JFunc,nk,BZextent)
end

function analyzeSpectrum(T,Sys::Geometry,Mod::Module;BZextent = 4pi,nk = 20,verbose = true)
    JFunc = constructJ(Sys,Mod)
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

function getChiFunction(T,Sys::Geometry,Mod::Module;BZextent = 4pi,nk = 20,tol = 1e-6,verbose = true)
    sinfo = analyzeSpectrum(T,Sys,Mod;BZextent = BZextent,nk = nk,verbose = verbose)

    JFunc = sinfo[:JFunc]
    EV = sinfo[:EV]
    Lam = sinfo[:LamSing]
    LeftBound = sinfo[:LeftBound]
    RightBound = sinfo[:RightBound]

    cons = 1E16
    try
        Lam = optimizeConstraint(EV,T,guess=[LeftBound,RightBound])
        # Lam = optimizeConstraint(JFunc,T,BZextent,guess=[LamSing+0.1*tol,LamSing+10+T])
        cons = constraint(EV,Lam,T)
        verbose && @info "constraint minimized to $cons for λ = $Lam"
    catch e
        verbose && @warn "constraint could not be optimized!
        (Errored with: $e). Lam set to singularity."
    end

    if abs(cons) > tol
        @warn "T= $T: constraint possibly not converged!
        constraint value: $cons \t tolerance: $tol"
    end

    @inline chi(q)= X_q(T,JFunc(q),Lam)

    return chi
end
##
