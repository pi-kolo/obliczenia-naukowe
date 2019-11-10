#Piotr Kołodziejczyk
#L2Z6

# x_(n+1) = (x_n)^2 +c

#parametry równania
parameters = [(-2, 1), (-2,2), (-2, 1.99999999999999), (-1,1), (-1,-1), (-1.0, 0.75), (-1.0 ,0.25)]

#iteracyjna realizacja równania logistycznego
function recursion(c, x0)
    prev=x0
    for i in 1:40
        println(prev)
        curr=prev*prev+c
        prev=curr
    end
end
n=7
recursion(parameters[n][1], parameters[n][2])