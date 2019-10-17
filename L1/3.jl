#Piotr Kołodziejczyk
#L1Z3

function evenly(x)
    for i in 1:2^52-1
        if nextfloat(x)-x!=2^(-52)
            println("error")
        end
        x=nextfloat(x)
    end
end

evenly(1.0)