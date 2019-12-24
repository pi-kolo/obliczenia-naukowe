#Piotr Kołodziejczyk

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

function changedGauss(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    B = copy(A)
    # v=Int64(n/l)
    count = 0
    #loop over columns
    for k in 1:n-1
        #loop over rows
        for i in k+1 : min(k+l+3, n)
            z = B[i,k]/B[k,k]
            B[i,k] = 0.0
            #for each cell change others in row
            for j in k+1 : min(k+2*l, n)
                B[i,j] -=  z*B[k,j]
            end
        end
    end
    return B
end

function changedGaussWithChoose(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    B = copy(A)
    # v=Int64(n/l)
    count = 0
    #permutations array
    p = [1:n;]
    #max elements for each column
    s = zeros(Float64, n)

    #itereate over columns
    for k in 1:n-1
        #get max values arrays 
        maxInCol=0
        indexMax = 0
        #over rows

        for i in k:n
            if abs(B[p[i], k]) > maxInCol
                maxInCol = abs(B[p[i], k])
                indexMax = i
            end
        end

        # println("indeks:", indexMax)
        # println("na maksie: ", B[p[indexMax],k])
        # maxT = findmax(s[k:n])
        # index = maxT[2]+k-1
        
        # println("Max: " ,maxT[1], " id:", index)
        p[indexMax], p[k] = p[k], p[indexMax]
        # println(p)

        #downto in k-column
        for i in k+1 : n
            z = B[p[i],k] / B[p[k],k]
            B[p[i],k] = z#Float64(0.0)
            for j in k+1 : n
                B[p[i],j] -= z*B[p[k],j]
            end
        end
    end
    println("Gauss matrix made")
    return B, p
end

function solveWithGauss(n::Int64,A::SparseArrays.SparseMatrixCSC{Float64, Int64}, b::Vector{Float64})
    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0
        for j in i+1 : n
            sum += A[i,j]*x[j]
        end
        x[i] = (b[i] - sum)/A[i,i]
    end
    return x
end

#
function solveWithChoiceGauss(n::Int64, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, p::Vector{Int64}, b::Vector{Float64})
    x=zeros(Float64,n)
    for k in 1: n-1
        for i in k+1 : n
            b[p[i]] -= A[p[i],k]*b[p[k]]
        end
    end
    for i in n:-1:1
        sum = 0
        for j in i+1 : n
            sum+=A[p[i],j]*x[j]
        end
        x[i]=(b[p[i]]-sum)/A[p[i],i]
    end
    return x
end

# x = [6,-2,2,4, 12,-8,6,10, 3,-13,9,3,-6,4,1,-18]
# y = SparseArrays.sparse( [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4],[1:4;1:4;1:4;1:4], x)
# A = IOfunctions.readMatrix("50000/A.txt")
# # println(A)
# B=copy(A)
# matrixgen.blockmat(10, 2, 15.0, "10.txt")
A = IOfunctions.readMatrix("16/A.txt")
b = IOfunctions.readRightSideVector("16/b.txt")

matrixgen.blockmat(100, 4, 5.0, "100.txt")

C = IOfunctions.readMatrix("100.txt")
d = Array(C)*ones(Float64, 100)

# C2 = changedGaussWithChoose(C, 100, 4)
# println(solveWithChoiceGauss(100, C2[1], C2[2], d))

E = IOfunctions.readMatrix("10000/A.txt")
f = IOfunctions.readRightSideVector("10000/b.txt")

E2 = changedGauss(E, 10000, 4)
println(solveWithGauss(10000, E2, f))
# A2 = changedGaussWithChoose(A, 16, 4)
# println(solveWithChoiceGauss(16, A2[1], A2[2], b))
# println(solveWithChoiceGauss)

