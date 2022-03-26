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
        sum += (evecs[beta,alpha] *conj(evecs[beta,alpha] ))/(β*evals[beta] + λ)
    end
    return sum/nb
end

function constraint(J̃::Function,λ::Real,β::Real,d::Integer)
    f(q) = real(tr(M_inverse(β,J̃(q),λ)))
    return 1/(2pi)^d * BZIntegral(f,d) -1
end

function optimizeConstraint(J̃::Function,β::Real,d::Integer)
    f(λ) = constraint(J̃,λ,β,d)
    find_zero(f,1)
end

function optimizeConstraint_brute(J̃::Function,β::Real,d::Integer,max,length=1000)
    lamRange = LinRange(0,max,length)
    f(λ) = constraint(J̃,λ,β,d)
    constr_arr = f.(lamRange)
    const2 = abs.(constr_arr)
    ind = argmin(const2)
    # println((ind, constr_arr[ind]))
    # println(const2)
    return lamRange[ind]
end


function BZIntegral(f::Function,d::Integer)
    val,err = hcubature(f, zeros(d), 2pi .*ones(d); reltol=1e12, abstol=1e-12, maxevals=0)
    return val
end


export constraint,M_inverse,X_q
##
