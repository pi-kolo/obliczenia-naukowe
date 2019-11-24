#Piotr Kołodziejczyk
#L3Z4
include("functions.jl")
"""
Testy zaimplementowanych metod dla zadanej funkcji f(x)=sin(x)-(x/2)^2

"""

#Wymagane dane: funkcja f, pochodna pf, oraz dokładności delta i epsilon
f(x)=sin(x)-(x/2)^2
pf(x) = cos(x) - x/2
delta = 10^-5 /2
epsilon = 10^-5 /2

println(AproxFunctions.mbisekcji(f, 1.0, 2.0, delta, epsilon))

println(AproxFunctions.mstycznych(f, pf, 1.5, delta, epsilon, 100))

println(AproxFunctions.msiecznych(f, 1.0, 2.0, delta, epsilon, 100))