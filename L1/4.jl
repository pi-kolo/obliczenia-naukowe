function notEqual(start)
    x=nextfloat(Float64(start))
    while Float64(x*Float64(1/x))==one(Float64)
        x=nextfloat(x)
    end
    return Float64(x)
end

println("x: 1<x<2 and x*(1/x) in Float64 = ", notEqual(1))

function minNotEq()
    x=nextfloat(Float64(0.0))
    while isinf(Float64(x*Float64(1/x))) #for small numbers result is infinity
        x*=2
    end
    while Float64(x*Float64(1/x))==one(Float64)
        x=nextfloat(x)
    end
    return x
end

println("The smallest x : x*(1/x) ∉ {inf, 1} = ", minNotEq())