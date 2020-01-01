#Piotr Kołodziejczyk

#solveWithGauss działa, w miare
#solveWithChoiceGauss działa, nawet w miarę, dla 50k nie tak długo, macierz kilka sekund ino
#solveWithLU działa, bez optymalizacji dla 50k kilkanaście sekund
#solveWithLUChoice działa (?, jakoś nagle zaczęło xD)

#so, oczyścić, uczytelnić i zakresy zrobić dokładniejsze :))

include("IOfunctions.jl")
include("matrixgen.jl")
import SparseArrays
import LinearAlgebra

function classicGauss(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64) :: SparseArrays.SparseMatrixCSC{Float64, Int64}
    B = copy(A)
    count=0
    #loop over columns
    for k in 1:n-1
        #loop over rows
        for i in k+1 : n
            z = B[i,k]/B[k,k]
            B[i,k] = 0
            #for each cell change others
            for j in k+1 : n
                B[i,j] = B[i,j] - z*B[k,j]
            end
            # println(count)
            count+=1
        end
    end
    println(count)
    return B
end 

"""
Function computing Gauss matrix (upper triangular) from matrix A
"""
function changedGauss(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
    B = copy(A)
    # v=Int64(n/l)
    count = 0
    #loop over columns
    for k in 1:n-1
        #loop over rows
        for i in k+1 : min(n,k+2*l)#min(k+l+3, n)
            z = B[i,k]/B[k,k]
            B[i,k] = z#0.0
            #for each cell change others in row
            for j in k+1 : min(n,k+2*l)
                B[i,j] -=  z*B[k,j]
            end
            b[i] -= z * b[k]
        end
    end
    return B
end

function changedGaussWithChoice(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b)
    B = copy(A)
    count = 0
    #permutations array
    p = [1:n;]
   
    #itereate over columns
    for k in 1:n-1

        #get max value in each column (from rows k to n)
        maxInCol=0
        indexMax = 0
        
        #over rows
        for i in k: min(n, k+ 2*l)
            if abs(B[p[i], k]) > maxInCol
                maxInCol = abs(B[p[i], k])
                indexMax = i
            end
        end

        p[indexMax], p[k] = p[k], p[indexMax]

        #downto in k-column
        for i in k+1 : min(n, k+ 2*l)
            z = B[p[i],k] / B[p[k],k]
            B[p[i],k] = Float64(0.0)
            for j in k+1 : min(n,k+2*l)
                B[p[i],j] -= z*B[p[k],j]
            end
            b[p[i]] -= z * b[p[k]]
        end
    end
    println("Gauss matrix made")
    return B, p, b
end

"""
Function solving A x = b
"""
function solveWithGauss(n::Int64,A::SparseArrays.SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, l::Int64) :: Vector{Float64}
    B = changedGauss(A, n, l, b)
    print("macierz wygenerowana")
    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0
        for j in i+1 : n
            sum += B[i,j]*x[j]
        end
        x[i] = (b[i] - sum)/B[i,i]
    end
    return x
end


"""
Function solving A x = b with choice of main element 
"""
function solveWithChoiceGauss(n::Int64, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, l::Int64)
    B = changedGaussWithChoice(A, n, l, b)
    matrix = B[1]
    p = B[2]
    b2 = B[3]
    x=zeros(Float64,n)
    for k in 1: n-1
        for i in k+1 : n
            b2[p[i]] -= matrix[p[i],k]*b2[p[k]]
        end
    end
    for i in n:-1:1
        sum = 0
        for j in i+1 : n
            sum+=matrix[p[i],j]*x[j]
        end
        x[i]=(b2[p[i]]-sum)/matrix[p[i],i]
    end
    return x
end
"""
Function generating matrices L and U after Gauss elimination
"""
function changedGaussLU(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    U = copy(A)
    # v=Int64(n/l)
    L = SparseArrays.spzeros(n,n)
    count = 0
    #loop over columns
    for k in 1:n-1
        L[k,k]=1.0
        #loop over rows
        for i in k+1 : min(n,k+2*l)#min(k+l+3, n)
            z = U[i,k]/U[k,k]
            L[i,k] = z
            U[i,k] = 0.0
            #for each cell change others in row
            for j in k+1 : min(n, k+2*l)
                U[i,j] -=  z*U[k,j]
            end
        end
    end
    L[n,n]=1
    return L, U
end

function changedGaussWithChoiceLU(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    U = copy(A)
    count = 0
    L = SparseArrays.spzeros(n,n)
    #permutations array
    p = [1:n;]
   
    #itereate over columns
    for k in 1:n-1
        # L[k,k]=1.0
        #get max value in each column (from rows k to n)
        maxInCol=0
        indexMax = 0
        
        #over rows
        for i in k: min(n, k+ 2*l - k%l)
            if abs(U[p[i], k]) > maxInCol
                maxInCol = abs(U[p[i], k])
                indexMax = i
            end
        end

        p[indexMax], p[k] = p[k], p[indexMax]
        # L[p[k],p[k]] = 1
        #downto in k-column
        for i in k+1 : min(n, k+ 2*l )
            z = U[p[i],k] / U[p[k],k]
            L[p[i],k] = z
            U[p[i],k] = Float64(0.0)
            for j in k+1 : min(n,k+2*l)
                U[p[i],j] -= z*U[p[k],j]
            end
        end
    end
    # L[p[1],p[1]] = 1.0
    # L[p[n],p[n]] = 1.0
    return  L, U, p
end


function solveWithLU(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b)
    x=zeros(Float64,n)
    for k in 1: n-1
        for i in k+1 : n
            b[i] -= L[i,k]*b[k]
        end
    end
    for i in n:-1:1
        sum = 0
        for j in i+1 : n
            sum+=U[i,j]*x[j]
        end
        x[i]=(b[i]-sum)/U[i,i]
    end
    return x
end

function solveWithLUChoice(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64},p, n::Int64, l::Int64, b)
    x=zeros(Float64,n)
    for k in 1: n-1
        for i in k+1 : n
            b[p[i]] -= L[p[i],k]*b[p[k]]
        end
    end
    for i in n:-1:1
        sum = 0
        for j in i+1 : n
            sum+=U[p[i],j]*x[j]
        end
        x[i]=(b[p[i]]-sum)/U[p[i],i]
    end
    return x
end
        


E = IOfunctions.readMatrix("50000/A.txt")
f = IOfunctions.readRightSideVector("50000/b.txt")

y = changedGaussWithChoiceLU(E, 50000, 4)
x = solveWithLUChoice(y[1], y[2],y[3], 50000, 4, f)
println(x)

