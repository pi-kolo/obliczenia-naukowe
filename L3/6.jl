#Piotr Kołodziejczyk
#L3Z5
include("functions.jl")
"""
Znajdowanie miejsc zerowych funkcji f1 i f2 zaimplementowanymi metodami
"""

#Dane: funkcje f1, f2 i ich pochodne f1p, f2p oraz dokładności delta i epsilon
f1(x) = exp(1-x)-1
f1p(x) = -exp(1-x)
f2(x)=x*exp(-x)
f2p(x)=-x*exp(-x)+exp(-x)

delta = 10^-5
epsilon = 10^-5

#dla f1 niech [-1,2]
println(AproxFunctions.mbisekcji(f1, -1.0, 2.0, delta, epsilon))
println(AproxFunctions.mstycznych(f1, f1p, 0.0, delta, epsilon, 100))
println(AproxFunctions.msiecznych(f1, -1.0, 2.0, delta, epsilon, 100 ))

#dla f2 [-1.5, 2]
println(AproxFunctions.mbisekcji(f2, -1.5, 2.0, delta, epsilon))
println(AproxFunctions.mstycznych(f2, f2p, 0.5, delta, epsilon, 100))
println(AproxFunctions.msiecznych(f2, -1.0, 1.0, delta, epsilon, 100))