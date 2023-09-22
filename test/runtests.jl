using LargeN, SpinFRGLattices
using StaticArrays
using Test
##
@testset "custom Jq" begin
    J(qx,qy) = SA[
        cos(qx)+cos(qy) sin(qx-qy);
        sin(qx-qy) cos(qx)+cos(qy)
    ]
    J(q::SVector{2}) = J(q[1],q[2]) 
    @testset "JFunc" begin
        LargeN.test_Couplings2D(J)
    end
    NCell = 2
    Dim = 2
    T = 0.1
    chi = getChiFunction(T,J,NCell,Dim)
    @test chi(1,2) ≈ 0.10639739802745554 atol = 1e-13
end
##
@testset "SpinFRGLattices" begin

    Py = Pyrochlore.getPyrochlore(4,[1,0])
    
    JFunc = constructJ(Py,Pyrochlore)
    @testset "JFunc" begin
        LargeN.test_Couplings3D(JFunc)
    end
    T = 0.1
    chi = getChiFunction(T,Py,Pyrochlore)
    
    q = [1,0.,0]
    qStatic = SArr.SA[1,0.,0] # non-allocating static array
    
    @test chi(q) == chi(1,0,0) == chi(qStatic)

    @test chi(1,2,3.) ≈ 0.08045869389334914 atol = 1e-13
    
end
