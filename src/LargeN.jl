module LargeN
    using LinearAlgebra,StaticArrays,Roots,Test,SpinFRGLattices,Parameters,GenericSchur
    using Cubature,StaticArrays
    include("LinAlg.jl")
    include("Lattices.jl")
    include("UnitTests.jl")
    
    export constraint,X_q,optimizeConstraint,X_q_eig,precomputeJ,constructJ,constructJtest, getEvals, getChiFunction
end