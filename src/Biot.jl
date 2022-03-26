module Biot

using LinearAlgebra, PyPlot, Interpolations




include("FieldCalc.jl")
include("Toroid_Inductance.jl")
include("MakeSolenoid.jl")
include("ToroidOptimizer.jl")
# Write your package code here.

end
