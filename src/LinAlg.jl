using LinearAlgebra,StaticArrays,Roots
using Cubature,StaticArrays
# using Brillouin,Crystalline

function M_inverse(β::Real,J̃::AbstractMatrix,λ::Real)
    dims = size(J̃)
    inv(β .* J̃+λ .* Matrix(I,dims...))
end

function X_q(nb::Integer,Minv::AbstractMatrix)
    real(sum(Minv)/nb)
end

function X_q(nb::Integer,β::Real,J̃::AbstractMatrix,λ::Real)
    real(sum(M_inverse(β,J̃,λ)))/nb
end

function X_q_eig(nb::Integer,β::Real,J̃::AbstractMatrix,λ::Real)
    evals, evecs = eigen(J̃)
    sum = 0.
    for alpha in eachindex(evals),beta in eachindex(evals),gamma in eachindex(evals)
        sum += (evecs[alpha,beta] *conj(evecs[gamma,beta] ))/(β*evals[beta] + λ)
    end
    return real(sum/nb)
end

function constraint(J̃::Function,λ::Real,β::Real,d::Integer)
    f(q) = real(tr(M_inverse(β,J̃(q),λ)))
    return 1/(2pi)^d * BZIntegral(f,d) -1
end

function optimizeConstraint(J̃::Function,β::Real,d::Integer;min = 2*β+0.01, max=2*β+10 )
    f(λ) = constraint(J̃,λ,β,d)
    find_zero(f,(min,max))
end

function optimizeConstraint_brute(J̃::Function,β::Real,d::Integer;min = 2*β+0.01, max=2*β+10 ,length=5000)
    lamRange = LinRange(min,max,length)
    # lam = LinRange(1.2*beta,1.3*beta,50000)

    f(λ) = constraint(J̃,λ,β,d)
    cons = f.(lamRange)
    filter!(isfinite,cons)
    # cons = L.constraint.(Ref(J̃),lamRange,Ref(β),d)
    Lam2 = lamRange[argmin(abs.(cons))]
end


function BZIntegral(f::Function,d::Integer)
    val,err = hcubature(f, zeros(d), 2pi .*ones(d); reltol=1e12, abstol=1e-12, maxevals=0)
    return val
end


export constraint,M_inverse,X_q
##
