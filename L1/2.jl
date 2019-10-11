#Piotr Kołodziejczyk
#L1Z2
types=[Float16, Float32, Float64]

function kahan(type)
    abs(type(3(type(4/3)-1)-1))
end

foreach(type -> println(type, ": ",kahan(type)," vs. eps: ", eps(type)), types)
