#Piotr Kołodziejczzyk
#L1Z7

function derivate(func, x, h)::Float64
    (func(x+h)-func(x))/h
end

function f(x)::Float64
    sin(x)+cos(3x)
end

truth=cos(1.0)-3sin(3.0)
for i in -1:-1:-54
    println(derivate(f,1.0,2.0^i)," ", abs(derivate(f,1.0,2.0^i)-truth))
end
