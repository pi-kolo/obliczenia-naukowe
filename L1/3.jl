#Piotr Kołodziejczyk
#L1Z3

#funkcja sprawdzająca, czy w przedziale (1,2) liczby rozmieszczone są równomiernie
function evenly(x)
    for i in 1:2^52-1
        if nextfloat(x)-x!=2^(-52)
            println("error")
        end
        x=nextfloat(x)
    end
end

#empiryczne przekonanie się o zapisie dwójkowym liczb
function printBits(x)
    for i in 1:5
        println(x)
        println(bitstring((x)))
        x=nextfloat(x)
    end
end

printBits(2.0)

printBits(0.5)
