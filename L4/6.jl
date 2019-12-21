#Piotr Kołodziejczyk
#L4Z6

include("interpolation.jl")
using Plots

f1(x) = abs(x)
interpolation.rysujNnfx(f1, -1.0, 1.0, 5, "6a1")
interpolation.rysujNnfx(f1, -1.0, 1.0, 10, "6a2")
interpolation.rysujNnfx(f1, -1.0, 1.0, 15, "6a3")

f2(x) = 1/(1+x^2)
interpolation.rysujNnfx(f2, -5.0, 5.0, 5, "6b1")
interpolation.rysujNnfx(f2, -5.0, 5.0, 10, "6b2")
interpolation.rysujNnfx(f2, -5.0, 5.0, 15, "6b3")