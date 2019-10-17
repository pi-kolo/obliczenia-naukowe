#Piotr Kołodziejczyk
#L1Z4

#funkcja szukająca iteracyjnie po kolei takiej liczby, która spełnia warunek zadania
#start - początek przedziału poszukiwań
function notEqual(start) :: Float64
    x=nextfloat(Float64(start))
    while Float64(x*Float64(1/x))==one(Float64)
        x=nextfloat(x)
    end
    return x
end

println("x: 1<x<2 and x*(1/x) in Float64 = ", notEqual(1))

#funkcja znajdująca najmniejszą taką liczbę >0
function minNotEq() :: Float64
    x=nextfloat(Float64(0.0))
    while isinf(Float64(x*Float64(1/x))) #dla małych liczb wynik to nieskończoność
        x*=2
    end
    while Float64(x*Float64(1/x))==one(Float64)
        x=nextfloat(x)
    end
    return x
end

println("The smallest x : x*(1/x) ∉ {inf, 1} = ", minNotEq())