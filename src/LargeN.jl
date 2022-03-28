module LargeN
    using LinearAlgebra,StaticArrays,Roots,Test,SpinFRGLattices,Parameters
    using Cubature,StaticArrays
    include("FourierIntegral.jl")
    include("LinAlg.jl")
    include("Lattices.jl")
    include("UnitTests.jl")
    
    export constraint,X_q,optimizeConstraint_brute,optimizeConstraint,X_q_eig,precomputeJ,constructJ,constructJtest
end