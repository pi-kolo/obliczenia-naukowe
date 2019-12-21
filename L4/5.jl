#Piotr Kołodziejczyk
#L4Z5

include("interpolation.jl")
using Plots

f1(x) = exp(x)
interpolation.rysujNnfx(f1, 0.0, 1.0, 5, "5a1")
interpolation.rysujNnfx(f1, 0.0, 1.0, 10, "5a2")
interpolation.rysujNnfx(f1, 0.0, 1.0, 15, "5a3")

f2(x) = x^2*sin(x)
interpolation.rysujNnfx(f2, -1.0, 1.0, 5, "5b1")
interpolation.rysujNnfx(f2, -1.0, 1.0, 10, "5b2")
interpolation.rysujNnfx(f2, -1.0, 1.0, 15, "5b3")