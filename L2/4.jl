#Piotr Kołodziejczyk
#L2Z4

using Polynomials

#współczynniki wielomianu
global p=[1, -210.0-2^(-23), 20615.0,-1256850.0,
      53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,          
      11310276995381.0, -135585182899530.0,
      1307535010540395.0,     -10142299865511450.0,
      63030812099294896.0,     -311333643161390640.0,
      1206647803780373360.0,     -3599979517947607200.0,
      8037811822645051776.0,      -12870931245150988800.0,
      13803759753640704000.0,      -8752948036761600000.0,
      2432902008176640000.0]
reverse!(p)

#dwa przedstawienia wielomianu, w postaci normalnej i iloczynowej
P = Poly(p)
p = poly(collect(1.0:20.0))

zera = roots(P)

for i in 1:20
    println("\\hline", i, "&", zera[i], "&", abs(P(zera[i])), "&", abs(p(zera[i])), "&", abs(i-zera[i]), "\\\\")
end

