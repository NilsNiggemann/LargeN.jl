function M_inverse(β::Real,J̃::AbstractMatrix,λ::Real)
    dims = size(J̃)
    inv(β .* J̃+λ .* Matrix(I,dims...))
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


function X_q(β::Real,J̃::AbstractMatrix,λ::Real,nb::Integer = getNCell(J̃))
    real(sum(M_inverse(β,J̃,λ)))/nb
end

function constraint(J̃::Function,λ::Real,T::Real,d::Integer = getDimfromFunc(J̃);kwargs...)
    f(q) = real(tr(M_inverse(1/T,J̃(q),λ)))
    return 1/(4pi)^d * BZIntegral(f,d;kwargs...) -1
end

function optimizeConstraint(J̃::Function,T::Real,d::Integer= getDimfromFunc(J̃);min = 2/T+0.01, max=2/T+10,guess = (min,max),method = Roots.Bisection(),kwargs... )
    f(λ) = constraint(J̃,λ,T,d)
    find_zero(f,guess,method;atol = 1e-15, rtol = 1e-15,kwargs...)
end

function BZIntegral(f::Function,d::Integer;kwargs...)
    val,err = hcubature(f, zeros(d), 4pi .*ones(d); reltol=1e-15, abstol=1e-15,maxevals = 1000,kwargs...)
    return val
end


function ComputeEigvals2D(Jfunc::Function,nk::Integer,ext;min = -ext,max = ext)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,NCell,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray)
        sys = eigen(Jfunc(SA[kx,ky,kz]))
        eig[:,i,j] .= sys.values
    end
    return eig
end

function ComputeEigvals3D(Jfunc::Function,nk::Integer,ext;min = -ext,max = ext)
    karray = range(min,max,length = nk)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,NCell,nk,nk,nk)
    for (i,kx) in enumerate(karray), (j,ky) in enumerate(karray),(k,kz) in enumerate(karray)
        sys = eigen(Jfunc(SA[kx,ky,kz]))
        eig[:,i,j,k] .= sys.values
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

function getEvals(JFunc,BZextent = 4pi;kEvals = 20)
    Dim = getDimfromFunc(JFunc)
    evalFunc = (nothing,ComputeEigvals2D,ComputeEigvals3D)[Dim]
    EV= evalFunc(JFunc,kEvals,BZextent)
end

function getChiFunction(T,Sys::Geometry,Mod::Module,BZextent = 4pi;kEvals = 20,tol = 1e-7,verbose = true)
    JFunc = constructJ(Sys,Mod)
    EV = getEvals(JFunc,BZextent,kEvals = kEvals)

    Emin = minimum(EV)/T
    verbose && @info "minimum Eigenvalue/T: $Emin"
    beta = 1 /T

    LamSing = -Emin
    Lam = LamSing
    try
        Lam = optimizeConstraint(JFunc,T,guess=[LamSing,LamSing+10+T])
    catch
        verbose && @warn "constraint could not be optimized! Lambda set to singularity"
    end

    cons = constraint(JFunc,Lam,T)

    if verbose
        @info "constraint minimized to $cons for Λ = $Lam"
        if abs(cons) > tol 
            @warn "constraint possibly not converged! (constraint tolerance set to $tol)"
        end
    end

    @inline chi(q)= X_q(beta,JFunc(q),Lam)

    return chi
end
##
