include("interpolation.jl")

x=[-1.0, 0.0, 1.0, 2.0]
f = [-1.0, 0.0, -1.0, 2.0]

fx = interpolation.ilorazyRoznicowe(x, f)
println(interpolation.naturalna(x, fx))