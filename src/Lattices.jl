
getDim(B::Basis_Struct_2D) = 2
getDim(B::Basis_Struct_3D) = 3

function precomputeJ(couplings::AbstractVector,UnitCell::AbstractVector{T},SiteList::AbstractVector{T},PairList::AbstractVector{T},PairTypes,pairToInequiv,Basis) where T <: Rvec
    Rij_vec = SVector{getDim(Basis),Float64}[]
    Jij_vec = Float64[]
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
                if abs(couplings[pair]) > 1e-15 
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

struct FourierInfo1{B}
    Rij_vec::Vector{B}
    Jij_vec::Vector{Float64}
    alpha_vec::Vector{Int}
    beta_vec::Vector{Int}
end
FourierInfo = FourierInfo1

makeUnitCell(Basis::Basis_Struct_2D) = [Rvec(0,0,x) for x in 1:Basis.NCell]
makeUnitCell(Basis::Basis_Struct_3D) = [Rvec(0,0,0,x) for x in 1:Basis.NCell]

precomputeJ(System::Geometry,Basis::Basis_Struct,pairToInequiv::Function) = precomputeJ(System.couplings,makeUnitCell(Basis),generateLUnitCells(System.NLen,Basis),System.PairList,System.PairTypes,pairToInequiv,Basis)

precomputeJ(System::Geometry,Mod::Module) = precomputeJ(System,Mod.Basis,Mod.pairToInequiv)

function constructJ(Infos::FourierInfo)

    @unpack Rij_vec,Jij_vec,alpha_vec,beta_vec = Infos
    NCell = maximum(alpha_vec)
    function J_func(q)
        J = zeros(ComplexF64,NCell,NCell)
        @inbounds for (Rij,Jij,α,β) in zip(Rij_vec,Jij_vec,alpha_vec,beta_vec)
            if α <= β
                J[α,β] += exp(1im*dot(q,Rij))*Jij
            end
        end
        return Hermitian(J)
    end
    return J_func
end

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
constructJtest(System::Geometry,Mod::Module) = constructJtest(precomputeJ(System,Mod))

function getDimfromFunc(Jfunc::Function)
    try
        return length(eltype(Jfunc.Rij_vec))[1]
    catch
        error("Function does not contain sufficient information, please provide dimension argument explicitly")
    end
end