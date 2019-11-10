#Piotr Kołodziejczyk
#L2Z3

using LinearAlgebra

function matcond(n::Int, c::Float64)
    # Function generates a random square matrix A of size n with
    # a given condition number c.
    # Inputs:
    #	n: size of matrix A, n>1
    #	c: condition of matrix A, c>= 1.0
    #
    # Usage: matcond(10, 100.0)
    #
    # Pawel Zielinski
            if n < 2
             error("size n should be > 1")
            end
            if c< 1.0
             error("condition number  c of a matrix  should be >= 1.0")
            end
            (U,S,V)=svd(rand(n,n))
            return U*diagm(0 =>[LinRange(1.0,c,n);])*V'
end

function hilb(n::Int)
    # Function generates the Hilbert matrix  A of size n,
    #  A (i, j) = 1 / (i + j - 1)
    # Inputs:
    #	n: size of matrix A, n>=1
    #
    #
    # Usage: hilb(10)
    #
    # Pawel Zielinski
            if n < 1
             error("size n should be >= 1")
            end
            return [1 / (i + j - 1) for i in 1:n, j in 1:n]
end

#testy dla macierzy Hilberta
for i in 1:20
    A=hilb(i)
    b=A*ones(i)
    println("\\hline", i, "& ", rank(A), "& ", cond(A), "& ", (norm(A\b-ones(i)))/norm(ones(i)), "& ", (norm(inv(A)*b - ones(i)))/norm(ones(i)), "\\\\")
end

sizes=[5,10,20]
conds = [1.0,10.0, 10.0^3, 10.0^7, 10.0^12, 10.0^16]

#testy dla losowych macierzy z zadanym wskaźnikiem uwarunkowania
for i in 1:3
    for j in 1:6
        A=matcond(sizes[i], conds[j])
        b=A*ones(sizes[i])
        println("\\hline ",sizes[i], "& ",rank(A), "& ",conds[j], "& ", (norm(A\b - ones(sizes[i])))/norm(ones(sizes[i])), "& ", (norm(inv(A)*b - ones(sizes[i])))/norm(ones(sizes[i])), "\\\\")
    end
end