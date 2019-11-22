include("functions.jl")

f(x) = 3x
g(x) = exp(x)

h(x) = 3x - exp(x)

delta = 10^-4
epsilon = 10^-4

println(AproxFunctions.mbisekcji(h, 0.0, 1.0, delta, epsilon))
println(AproxFunctions.mbisekcji(h, 1.0, 2.0, delta, epsilon))