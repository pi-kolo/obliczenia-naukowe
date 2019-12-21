#Piotr Kołodziejczyk
#Moduł funkcji używanych do interpolacji

using Plots

module interpolation
export ilorazyRoznicowe, warNewton, rysujNnfx, naturalna

"""
Funkcja obliczająca ilorazy różnicowe
Dane:
    x- wektor długości n+1 zawierający węzły x_0, ..., x_n
    f - wektor długości n+1 zawierający wartości interpolowanej funkcji w węzłach
Wyniki:
    fx - wektor długości n+1 zawierający obliczone ilorazy różnicowe, tak że
        fx[1] = f[x_0]
        fx[2] = f[x_0, x_1]
        ...
        fx[n] = f[x_0, x_1, ..., x_n]
"""
function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64})
    len = length(f)
    fx = Vector{Float64}(undef, len)

    for i in 1:len
        fx[i] = f[i]
    end

    for i in 2:len
        for j in len:-1:i
            fx[j] = (fx[j] - fx[j-1]) / (x[j] - x[j-i+1])
        end
    end

    return fx
end

"""
Funkcja obliczająca wartość wielomianu interpolacyjnego stopnia n w postaci Newtona N_n(x)
    w punkcie x=t za pomocą uogólnionego algorytmu Hornera
Dane:
    x -wektor długości n+1 zawierający węzły x_0, ..., x_n
    fx - wektor długości n+1 zaierający kolene ilorazy różnicowe, takie że
        fx[1] = f[x_0]
        fx[2] = f[x_0, x_1]
        ...
        fx[n] = f[x_0, x_1, ..., x_n]
    t - punkt, w którym należy obliczyć wartość wielomianu
Wyniki:
    nt - wartość wielomianu w punkcie t
"""
function warNewton(x::Vector{Float64}, fx::Vector{Float64}, t::Float64)
    leng = length(x)
    nt = fx[leng]
    for i in leng-1 : -1 : 1
        nt = fx[i] + (t-(x[i]))*nt
    end
    return nt
end

"""
Funkcja wyznaczająca postać normalną wielomiani interpolacyjnego na podstawie wektora ilorazów różnicowych i wektora węzłów
Dane:
    x - wektor długości n+1 zawierający węzły x_0, ..., x_n
    fx - wektor długości n+1 zawierający ilorazy różnicowe, takie że
        fx[1] = f[x_0]
        fx[2] = f[x_0, x_1]
        ...
        fx[n] = f[x_0, x_1, ..., x_n]
Wyniki:
    a - wetkor długości n+1 zawierający obliczone współczynniki postaci naturalnej, takie że
        a[1] = a_0
        a[2] = a_1
        ...
        a[n+1] = a_n
"""
function naturalna(x::Vector{Float64}, fx::Vector{Float64})
    
    len = length(x)
    a = Vector{Float64}(undef, len)

    a[len] = fx[len]
    for i in len-1 : -1 : 1
        a[i] =  fx[i] - a[i+1] * x[i]
        for j in i+1 : len-1
            a[j] = a[j] - a[j+1] * x[i]
        end
    end
    return a
end

"""
Funkcja interpolująca f(x) na przedziale [a,b] wielomianem Newtona stopnia n
Dane:
    f - funkcja f(x) 
    a,b - przedział interpolacji
    n - stopień wielomianu interpolacyjnego
Wyniki:
    wykres funkcji i wielomianu interpolacyjnego na [a,b]
"""
function rysujNnfx(f, a::Float64, b::Float64, n::Int, filename)
    m = n*20
    x = Vector{Float64}(undef, n+1)
    y = Vector{Float64}(undef, n+1)

    deltax = (b-a)/n
    
    for i in 1:n+1
            x[i] = a+(i-1)*deltax
            y[i] = f(x[i])
    end
    fx = ilorazyRoznicowe(x,y)
    
    args = Vector{Float64}(undef, m)
    f_values = zeros(m)
    p_values = zeros(m)
    
    deltax = (b-a)/m
    
    for i in 1:m
        args[i] = a + deltax*i
        f_values[i] = f(args[i])
        p_values[i] = warNewton(x, fx, args[i])
        
    end
    
    plot(range(a, stop=b, length=m), f_values, color="red", label="f(x)")
    plot!(range(a, stop=b, length=m), p_values, color="blue", label="p(x)")
    savefig("$filename.png")
    
end

end