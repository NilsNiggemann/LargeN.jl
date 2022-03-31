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

function X_q_eig(β::Real,J̃::AbstractMatrix,λ::Real,nb::Integer = getNCell(J̃))
    evals, evecs = eigen(J̃)
    sum = 0.
    for alpha in eachindex(evals),beta in eachindex(evals),gamma in eachindex(evals)
        sum += (evecs[alpha,beta] *conj(evecs[gamma,beta] ))/(β*evals[beta] + λ)
    end
    return real(sum/nb)
end

function constraint(J̃::Function,λ::Real,β::Real,d::Integer = getDimfromFunc(J̃);kwargs...)
    f(q) = real(tr(M_inverse(β,J̃(q),λ)))
    return 1/(4pi)^d * BZIntegral(f,d;kwargs...) -1
end

function optimizeConstraint(J̃::Function,β::Real,d::Integer= getDimfromFunc(J̃);min = 2*β+0.01, max=2*β+10,guess = (min,max),method = Roots.Bisection(),kwargs... )
    f(λ) = constraint(J̃,λ,β,d)
    find_zero(f,guess,method;kwargs...)
end

function optimizeConstraint_brute(J̃::Function,β::Real,d::Integer= getDimfromFunc(J̃);min = 2*β+0.01, max=2*β+10 ,points=800,kwargs...)
    lamRange = LinRange(min,max,points)
    # lam = LinRange(1.2*beta,1.3*beta,50000)

    f(λ) = constraint(J̃,λ,β,d;kwargs...)
    cons = f.(lamRange)
    inds = findall(isfinite,cons)
    cons = cons[inds]
    lamRange = lamRange[inds]
    Lam2 = lamRange[argmin(abs.(cons))]
    return Lam2, lamRange,cons
end


function BZIntegral(f::Function,d::Integer;kwargs...)
    val,err = hcubature(f, zeros(d), 4pi .*ones(d); reltol=1e-15, abstol=1e-15,maxevals = 1000,kwargs...)
    return val
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

"""Helper function: constraint fulfillment is often to the right of a singularity, i.e. a root of this function."""
@inline function SingularConstraint(eval::AbstractArray,Lam::Real,T::Real)
    s = 0.
    for e in eval
        s += 1/(e/T+Lam)
    end
    length(eval)/s
end

function optimizeConstraint(eval::AbstractArray,T::Real;guess =3/T)
    @inline constr(Lam) = constraint(eval,Lam,T)
    Lambda = find_zero(constr,guess,atol = 1e-15, rtol = 1e-15)
end

function optimizeSingularConstraint(eval::AbstractArray,T::Real;guess =3/T,kwargs...)
    @inline constr(Lam) = SingularConstraint(eval,Lam,T)
    Lambdas = find_zeros(constr,guess,atol = 1e-15, rtol = 1e-15;kwargs...)
end

function susc2D(Evals::AbstractArray,Evecs::AbstractArray,Lambda::Real,T::Real)
    nCell,nx,ny = size(Evals)
    Chi = Array{Float64}(undef,nx,ny)
    for i in 1:nx, j in 1:ny
        v = @view Evecs[:,:,i,j]
        e = @view Evals[:,i,j]
        s = 0.
        for a in 1:nCell, b in 1:nCell, c in 1:nCell
            s += (v[a, b]* conj(v[c, b])/(e[b]/T + Lambda))/(nx*ny)
        end
        Chi[i,j] = real(s)
    end
    Chi
end


function susc3D(Evals::AbstractArray,Evecs::AbstractArray,Lambda::Real,T::Real)
    nCell,nx,ny,nz = size(Evals)
    Chi = Array{Float64}(undef,nx,ny,nz)
    for i in 1:nx, j in 1:ny,k in 1:nz
        v = @view Evecs[:,:,i,j,k]
        e = @view Evals[:,i,j,k]
        s = 0.
        for a in 1:nCell, b in 1:nCell, c in 1:nCell
            s += (v[a, b]* conj(v[c, b])/(e[b]/T + Lambda))/(nx*ny)
        end
        Chi[i,j,k] = real(s)
    end
    Chi
end

function getEvals(JFunc,BZextent = 4pi;kEvals = 20)
    Dim = getDimfromFunc(JFunc)
    evalFunc = (nothing,ComputeEigvals2D,ComputeEigvals3D)[Dim]
    EV= evalFunc(JFunc,kEvals,BZextent)
end

function getChiFunction(T,Sys::Geometry,Mod::Module,BZextent = 4pi;kEvals = 20,singularShift = 1e-12, guess = (0.2/T,40/T),tol = 1e-7)
    JFunc = constructJ(Sys,Mod)
    EV = getEvals(JFunc,BZextent,kEvals = kEvals)
    beta = 1 /T
    LamSing = last(sort(optimizeSingularConstraint(EV,T,guess=guess)))+ singularShift
    Lam = LamSing
    @info "looking to the right of singularity at $LamSing"
    try
        Lam = optimizeConstraint(EV,T,guess=[LamSing,LamSing+5]) + singularShift
    catch
        @warn "constraint could not be optimized! Lambda set to singularity"
    end

    cons = constraint(EV,Lam,T)
    consright = constraint(EV,Lam+0.001*beta,T)
    dev = abs(cons - consright)
    @info "constraint minmized to $cons for Λ = $Lam"
    if abs(cons) > tol || dev > 10
        @warn "constraint possibly not converged! (constraint tolerance set to $tol, fluctuation is $dev)"
    end

    function chi(q)
        X_q(beta,JFunc(q),Lam)
    end
    return chi
end
##
