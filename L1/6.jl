#Piotr Kołodziejczyk
#L1Z6

function f(x::Float64)::Float64
    sqrt(x^2+1)-1
end

function g(x::Float64)::Float64
    x^2/(sqrt(x^2+1)+1)
end

for i in 1:10
    println(i, ", ",f(8.0^(-i)),", ", g(8.0^(-i)))
end
for i in 20:10:50
    println(i,", ",f(8.0^(-i)), ", ", g(8.0^(-i)))
end