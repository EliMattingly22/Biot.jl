module Biot

using LinearAlgebra, PyPlot, Interpolations




include("FieldCalc.jl")
include("Toroid_Inductance.jl")
include("MakeSolenoid.jl")
include("ToroidOptimizer.jl")
include("EddyCurrent.jl")
include("Transformer.jl")
# Write your package code here.

end
