#Piotr Kołodziejczyk
#L2Z5

# model wzrostu populacji
# p_(n+1) = p_n + r*p_n*(1-p_n)

#funkcja wyliczająca wartość w następnym okresie
function next(prev, r, type)
    return type(prev + r*prev*(1-prev))
end

#iteracja kolejnych n okresów
function consecutiveIteration(n, p0, r, type)
    prev=p0
    for i in 1:n
        curr = next(prev, r, type)
        println("&",curr, "\\\\")
        prev=curr
    end
end

#iteracja ze z obcięciem po n1 okresach
function iterationWithChange(n, n1, p0, r, type)
    prev=p0
    for i in 1:n1
        curr = next(prev, r, type)
        println(curr)
        prev=curr
    end
    prev=0.722
    for i in n1:n
        curr = next(prev,r,type)
        println(curr)
        prev=curr
    end
end

consecutiveIteration(40, 0.01, 3, Float32)
println()
iterationWithChange(40, 10, 0.01, 3, Float32)


consecutiveIteration(40, 0.01, 3, Float32)
println()
consecutiveIteration(40, 0.01, 3, Float64)
