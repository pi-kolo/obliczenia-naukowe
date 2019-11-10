#Piotr Kołodziejczyk
#L1Z5 / L2Z1

#wektory ze zmienionymi wartościami
x=[2.718281828, -3.141592654, 1.414213562, 0.577215664, 0.301029995]
y= [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

# wektory w Float32
x32 = map(x->Float32(x),x)
y32 = map(x->Float32(x),y)

# wektory w Float64
x64 = map(x->Float64(x), x)
y64 = map(x->Float64(x), y)

# sposób a) - sumowanie w przód
# x,y - wektory
function addForward(x::Vector,y::Vector)
    sum=0
    for i=1:length(x)
        sum+=x[i]*y[i]
    end
    return sum
end

# sposób b) - sumowanie w tył
# x,y - wektory
function addBackward(x::Vector,y::Vector)
    sum=0
    for i=length(x):-1:1
        sum+=x[i]*y[i]
    end
    return sum
end

# sposób c) 
# x,y - wektory
function addOrderedC(x::Vector,y::Vector)
    results=[]
    for i=length(x):-1:1
        push!(results, x[i]*y[i])
    end
    sumPos=foldl(+,sort(filter(x -> x>=0, results), rev=true))
    sumNeg=foldl(+, sort(filter(x -> x<0, results)))
    return sumPos+sumNeg
end

# sposób d) 
# x,y - wektory
function addOrderedD(x::Vector,y::Vector)
    results=[]
    for i=length(x):-1:1
        push!(results, x[i]*y[i])
    end
    sumPos=foldl(+,sort(filter(x -> x>=0, results)))
    sumNeg=foldl(+, sort(filter(x -> x<0, results), rev=true))
    return sumPos+sumNeg
end



functions = [addForward, addBackward, addOrderedC, addOrderedD]
println("Float32:")
foreach(type -> println(type, ": ",type(x32,y32)), functions)
println("Float64:")
foreach(type -> println(type, ": ",type(x64,y64)), functions)

