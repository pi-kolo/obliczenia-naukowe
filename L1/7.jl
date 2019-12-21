#Piotr Kołodziejczzyk
#L1Z7

# przybliżona wartość pochodnej w punkcie
# func - funkcja, której pochodna w punkcie ma zostać wyliczona
# x - punkt, w którym ma być policzona, h - odległość od x nadająca dokładność
function derivate(func::Function, x::Float64, h::Float64)::Float64
    Float64((func(x+h)-func(x))/h)
end

#badana funkcja
function f(x::Float64)::Float64
    Float64(sin(x)+cos(3x))
end

exactResult=Float64(cos(1.0)-3sin(3.0))
for i in -1:-1:-54
    println(i,", 1+h= ",1+2.0^i,", f'=", derivate(f,1.0,2.0^i),", błąd= ", abs(derivate(f,1.0,2.0^i)-exactResult))
end