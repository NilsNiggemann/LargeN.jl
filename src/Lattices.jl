
getDim(B::Basis_Struct_2D) = 2
getDim(B::Basis_Struct_3D) = 3

function precomputeJ(couplings::AbstractVector{FloatType},UnitCell::AbstractVector{T},SiteList::AbstractVector{T},PairList::AbstractVector{T},PairTypes,pairToInequiv,Basis) where {T <: Rvec, FloatType <: AbstractFloat}
    Rij_vec = SVector{getDim(Basis),FloatType}[]
    Jij_vec = FloatType[]
    alpha_vec = Int[]
    beta_vec = Int[]

    for i_site in UnitCell
        Ri = getCartesian(i_site,Basis)
        # println(Ri)
        for j_site in SiteList # site summation
            Rj = getCartesian(j_site,Basis)
            Rij = Ri - Rj
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                if abs(couplings[pair]) > eps(FloatType) 
                    α = i_site.b
                    β = j_site.b
                    push!(Rij_vec,Rij)
                    push!(Jij_vec,couplings[pair])
                    push!(alpha_vec,α)
                    push!(beta_vec,β)
                end
            end
        end
    end
    return FourierInfo(Rij_vec,Jij_vec,alpha_vec,beta_vec)
end

struct FourierInfo1{B,T}
    Rij_vec::Vector{B}
    Jij_vec::Vector{T}
    alpha_vec::Vector{Int}
    beta_vec::Vector{Int}
end
FourierInfo = FourierInfo1

makeUnitCell(Basis::Basis_Struct_2D) = [Rvec(0,0,x) for x in 1:Basis.NCell]
makeUnitCell(Basis::Basis_Struct_3D) = [Rvec(0,0,0,x) for x in 1:Basis.NCell]

precomputeJ(System::Geometry,Basis::Basis_Struct,pairToInequiv::Function) = precomputeJ(System.couplings,makeUnitCell(Basis),generateLUnitCells(System.NLen,Basis),System.PairList,System.PairTypes,pairToInequiv,Basis)

precomputeJ(System::Geometry,Mod::Module) = precomputeJ(System,Mod.Basis,Mod.pairToInequiv)

function constructJ(Infos::FourierInfo{B,T},::Val{NCell}) where {B,T,NCell}

    @unpack Rij_vec,Jij_vec,alpha_vec,beta_vec = Infos

    upperinds = findall(i -> alpha_vec[i] < beta_vec[i],collect(eachindex(alpha_vec,beta_vec)))
    diaginds = findall(i -> alpha_vec[i] == beta_vec[i],collect(eachindex(alpha_vec,beta_vec)))
    
    
    @inline function J_func(q)
        # J = zeros(ComplexF64,NCell,NCell)
        # @inbounds for (Rij,Jij,α,β) in zip(Rij_vec,Jij_vec,alpha_vec,beta_vec)
        J = MMatrix{NCell,NCell,Complex{T},NCell*NCell}(undef)
        fill!(J,0. +0im)
        @inbounds for i in upperinds
            α = alpha_vec[i]
            β = beta_vec[i]
            J[α,β] += exp(1im*dot(q,Rij_vec[i]))*Jij_vec[i]
        end
        @inbounds for i in diaginds
            α = alpha_vec[i]
            J[α,α] += cos(dot(q,Rij_vec[i]))*Jij_vec[i] 
        end
        return Hermitian(J)
    end
    
    return J_func
end

# function isInLattice(r::AbstractArray,Basis)
#     r_lattice = Basis.T *r
#     for (ib,b) in enumerate(Basis.bLatt)
#         r_curr =  r_lattice .- b
#         Intvec = round.(Int,r_curr)
#         if all(abs2.(Intvec .- r_curr) .<1E-14) # r_lattice contains only Ints
#             return true
#         end
#     end
#     return false
# end

# function isInversionSymmetric(B::Basis_Struct)
#     for b in B.b
#         !isInLattice(-b,B) && return false
#     end
#     return true
# end

constructJ(Infos::FourierInfo) = constructJ(Infos,Val(maximum(Infos.alpha_vec)))

function constructJtest(Infos::FourierInfo)

    @unpack Rij_vec,Jij_vec,alpha_vec,beta_vec = Infos
    NCell = maximum(alpha_vec)
    function J_func(q)
        J = zeros(ComplexF64,NCell,NCell)
        for (Rij,Jij,α,β) in zip(Rij_vec,Jij_vec,alpha_vec,beta_vec)
            J[α,β] += exp(1im*dot(q,Rij))*Jij
        end
        return J
    end
    dim = length(eltype(Rij_vec))[1]
    if dim == 2
        test_Couplings2D(J_func)
    elseif dim == 3
        test_Couplings3D(J_func)
    end
    return J_func
end

constructJ(System::Geometry,Mod::Module) = constructJ(precomputeJ(System,Mod))
constructJ(System::Geometry,Basis::Basis_Struct,pairToInequiv::Function) = constructJ(precomputeJ(System,Basis,pairToInequiv))
constructJtest(System::Geometry,Mod::Module) = constructJtest(precomputeJ(System,Mod))

function getDimfromFunc(Jfunc::Function)
    try
        return length(eltype(Jfunc.Rij_vec))[1]
    catch
        error("Function does not contain sufficient information, please provide dimension argument explicitly")
    end
end