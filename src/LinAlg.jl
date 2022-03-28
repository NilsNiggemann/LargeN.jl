function M_inverse(β::Real,J̃::AbstractMatrix,λ::Real)
    dims = size(J̃)
    inv(β .* J̃+λ .* Matrix(I,dims...))
end

function X_q(nb::Integer,Minv::AbstractMatrix)
    real(sum(Minv)/nb)
end

function getNCell(Jfunc::Function)
    try
        return Jfunc.NCell
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

function optimizeConstraint_brute(J̃::Function,β::Real,d::Integer= getDimfromFunc(J̃);min = 2*β+0.01, max=2*β+10 ,points=800)
    lamRange = LinRange(min,max,points)
    # lam = LinRange(1.2*beta,1.3*beta,50000)

    f(λ) = constraint(J̃,λ,β,d)
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

function ComputeEig(Jfunc,n)
    k = range(-4pi,4pi,length = n)
    NCell = getNCell(Jfunc)
    eig = Array{Float64}(undef,6,n,n)
    vec = Array{ComplexF64}(undef,6,6,n,n)
    for (i,kx) in enumerate(k), (j,ky) in enumerate(k)
        sys = eigSysJ(kx,ky,coup...)
        eig[:,i,j] = sys.values
        vec[:,:,i,j] = sys.vectors
    end
    return eig,vec
end
##
