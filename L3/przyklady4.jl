include("functions.jl")

print(AproxFunctions.mbisekcji(x->(sin(x)-0.25x^2), 1.5, 2.0, 0.5*10^-5, 0.5*10^-5))

print(AproxFunctions.mstycznych(x->(sin(x)-0.25x^2), x->(cos(x)-0.5x), 1.5, 0.5*10^-5, 0.5^10^-5, 100))

print(AproxFunctions.msiecznych(x->(sin(x)-0.25x^2), 1.0, 2.0, 0.5*10^-5, 0.5^10^-5, 100))