isApproxHermitian(M,dig=16) = ishermitian(round.(M,digits = dig))
function test_Couplings2D(Jfunc,digits=14)
    qrange = LinRange(-1.2pi,1.2pi, 10)
    nanfailed = 0
    hermitianfailed = 0
    for qx in qrange, qy in qrange
        J = Jfunc(SA[qx,qy])
        !isApproxHermitian(J,digits) && (hermitianfailed +=1)
        any(isnan.(J)) && (nanfailed +=1)
    end
    @testset "hermiticity" begin
        @test hermitianfailed == 0
    end
    @testset "nan" begin
        @test nanfailed == 0
    end
end
function test_Couplings3D(Jfunc,digits=14)
    qrange = LinRange(-1.2pi,1.2pi, 10)
    nanfailed = 0
    hermitianfailed = 0
    for qx in qrange, qy in qrange, qz in qrange
        J = Jfunc(SA[qx,qy,qz])
        !isApproxHermitian(J,digits) && (hermitianfailed +=1)
        any(isnan.(J)) && (nanfailed +=1)
    end
    @testset "hermiticity" begin
        @test hermitianfailed == 0
    end
    @testset "nan" begin
        @test nanfailed == 0
    end
end