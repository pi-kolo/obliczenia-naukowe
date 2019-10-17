#Piotr Kołodziejczyk

#L1Z2

#tablica typów liczb zmiennopozycyjnych
types=[Float16, Float32, Float64]

#funkcja wyliczająca wartość 3*(4/3 - 1) -1,
#type - typ, w jakim wartość ma być wyliczona 
function kahan(type::Type)
    abs(3(type(4/3)-1)-1)
end

foreach(type -> println(type, ": ",kahan(type)," vs. eps: ", eps(type)), types)
