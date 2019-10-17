#Piotr Kołodziejczyk
#L1Z1


#tablica typów liczb zmiennopozycyjnych
types=[Float16, Float32, Float64]

#funkcja wyliczająca wartość epsilonu maszynowego
#type - typ, dla którego ma być wyliczona 
function macheps(type::Type)
    x=one(type)
    while one(type)+x/2>one(type)
        x/=2
        if typeof(x)!=type
            println("Nastąpiła niejawna konwersja typu (?)")
        end
    end
    return x
end

#funkcja wyliczająca wartość eta
#type - typ, dla którego ma być wyliczona
function eta(type::Type)
    x=one(type)
    while x/2>zero(type)
        x/=2
        if typeof(x)!=type
            println("Nastąpiła niejawna konwersja typu (?)")
        end
    end
    return x
end

#funkcja wyliczająca wartość maksymalna
#type - typ, dla którego ma być wyliczona
function maxInType(type)
    x=one(type)
    while !isinf(x*2)
        x*=2
    end
    y=x/2
    while !isinf(x+y) && y>=1.0
        x+=y
        y/=2
    end
    return x
end

println("Macheps:")
foreach(type -> println(type, ": ",macheps(type)," vs. eps: ", eps(type)), types)
println("Eta:")
foreach(type -> println(type, ": ",eta(type)," vs. eps: ", nextfloat(type(0.0))), types)
println("Max number:")
foreach(type -> println(type, ": ",maxInType(type)," vs. max: ", floatmax(type)), types)


