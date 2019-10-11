#Piotr Kołodziejczyk
#L1Z6

function f(x)::Float64
    sqrt(x^2+1)-1
end

function g(x)::Float64
    x^2/(sqrt(x^2+1)+1)
end

for i in 1:10
    println(f(8.0^(-i))," vs ", g(8.0^(-i)))
end