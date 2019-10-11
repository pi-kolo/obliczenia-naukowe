#Piotr Kołodziejczyk
#L1Z1

types=[Float16, Float32, Float64]

function macheps(type)
    x=one(type)
    while one(type)+x/2>one(type)
        x/=2
        if typeof(x)!=type
            println("LOL")
        end
    end
    return x
end

function eta(type)
    x=one(type)
    while x/2>zero(type)
        x/=2
    end
    return x
end

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


