#Piotr Kołodziejczyk
#Lista3

module AproxFunctions

export mbisekcji, msiecznych, mstycznych

"""
Funkcja realizująca metodę bisekcji rozwiązywania równiania f(x)=0
Parametry:
    f - badana Funkcja
    a,b - końce przedziałów
    delta,epsilon - dokładności obliczeń

Zwracana czwórka (r,v,it,err) taka że:
    r - przybliżenie pierwiastka równania f(x)=0
    v - wartość f(r)
    it - liczba wykonanych iteracji
    err - sygnalizacja błędu
        0 - brak błędu
        1 - funkcja nie zmienia znaku na przedziale [a,b]
"""

#while nie chce współpracować jeszcze xD
function mbisekcji(f, a::Float64, b::Float64, delta::Float64, epsilon::Float64)
    u=f(a)
    v=f(b)
    e=b-a
    # c, w = 0, 0 
    it=0

    if sign(u)==sign(v)
        return (0,0,0,1)
    end

    # while abs(e) > delta  
    for k in 1:200
        it++
        e=e/2
        c=a+e
        w=f(c)
        if abs(w) < epsilon || abs(e) < delta
            return (c, w, it, 0)
        end
        if sign(w) != sign(u)
            b,v = c, w
        else
            a, u=c, w
        end
    end
    return "xd"
end

"""
Funkcja rozwiązująca równanie f(x)=0 metodą Newtona (stycznych)
Parametry:
    f, pf - funkcja f(x) oraz pochodna f'(x) 
    x0  -przybliżenie początkowe
    delta, espilon - dokładności obliczeń
    maxit - maksymalna dopuszczalna liczba iteracji

Zwracana czwórka (r,v,it,err) taka że:
    r - przybliżenie pierwiastka równiania f(x)=0
    v - wartość f(r)
    it - liczba wykoanych iteracji
    err - sygnalizacja błędu
        0 - metoda zbieżna
        1 - nie osiągnięto wymaganej dokładności w maxit iteracji
        2 - pochodna bliska zeru
"""

function mstycznych(f, pf, x0::Float64, delta::Float64, epsilon::Float64, maxit::Int)
    v=f(x0)
    if abs(v) < epsilon
        return (x0, v, 0, 0)
    end
    for it in 1:maxit
        x1 = x0-v/pf(x0)
        v = f(x1)
        if abs(x1-x0) < delta || abs(v) < epsilon
            return (x1, v, it, 0)
        end
        x0=x1
    end
    return (x0, v, maxit, 1)
end

"""
Funkcja rozwiązująca rówanie f(x)=0 metodą siecznych
Parametry:
    f - funkcja f
    x0,x1 - przybliżenia początkowe
    delta, epsilon - dokładności obliczeń
    maxit - maksymalna dopuszczalna liczba iteracji

Zwracana czwórka (t, v, it, err) taka że:
    r - przybliżenie pierwiastka równania f(x)=0
    v - wartość f(r)
    it - liczba wykonanych iteracji
    err - sygnalizacja błędu
        0 - metoda zbieżna
        1 - nie osiągnięto wymaganej dokładności w maxit iteracji
"""

function msiecznych(f, x0::Float64, x1::Float64, delta::Float64, epsilon::Float64, maxit::Int)
    a, b = x0, x1
    fa, fb = f(a), f(b)
    for it in 1:maxit
        if abs(fa) > abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end
        s = (b-a)/(fb - fa);
        b, fb = a, fa
        a = a - fa*s
        fa = f(a)
        if abs(b-a) < delta || abs(fa) < epsilon
            return (a, fa, it, 0)
        end
    end
    return (a, fa, maxit, 1) 
end

end