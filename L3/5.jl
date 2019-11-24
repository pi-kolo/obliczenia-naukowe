#Piotr Kołodziejczyk
#L3Z5
include("functions.jl")
"""
Znajdowanie punktu przecięcia funkcji f i g metodą bisekcji
"""

#Dane: funkcje, dokładności delta i epsilon
f(x) = 3x
g(x) = exp(x)
h(x) = 3x - exp(x)
delta = 10^-4
epsilon = 10^-4

#rozwiązanie dla dwu przedziałów: [0,1] oraz [1,2] (otrzymanych w wyniku analizy funkcji)
println(AproxFunctions.mbisekcji(h, 0.0, 1.0, delta, epsilon))
println(AproxFunctions.mbisekcji(h, 1.0, 2.0, delta, epsilon))