#Piotr Kołodziejczyk
#L1Z6

#funkcje f i g jak w zadaniu
function f(x::Float64)::Float64
    sqrt(x^2+1)-1
end

function g(x::Float64)::Float64
    x^2/(sqrt(x^2+1)+1)
end

#sprawdzenie wyników dla różnych potęg 8
for i in 1:10
    println(i, ", f=",f(8.0^(-i)),", g=", g(8.0^(-i)))
end
for i in 20:10:100
    println(i, " , f=",f(8.0^(-i)),", g=", g(8.0^(-i)))
end