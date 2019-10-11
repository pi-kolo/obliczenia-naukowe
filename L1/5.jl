#Piotr Kołodziejczyk
#L1Z5

x=[2.718281828, -3.141592654, 1.414213562, 0.5772156649, 0.3010299957]
y= [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

function addForward(x,y)
    sum=0
    for i=1:length(x)
        sum+=x[i]*y[i]
    end
    return sum
end

function addBackward(x,y)
    sum=0
    for i=length(x):-1:1
        sum+=x[i]*y[i]
    end
    return sum
end

function addOrdered(x,y)
    sum=0
    posRes=[]
    negRes=[]
    for i=length(x):-1:1
        res=x[i]*y[i]
        res>=0 ? push!(posRes, res) : push!(negRes, res)
    end
    #=
    reverse!(sort!(posRes))
    sort!(negRes)
    posSum=foldl(+, posRes)
    negSum=foldl(+, negRes)
    return posSum+negSum
    =#
    return posRes, negRes
end



functions = [addForward, addBackward, addOrdered]#, addReversed]
foreach(type -> println(type, ": ",type(x,y)), functions)

