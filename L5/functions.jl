#Piotr Kołodziejczyk

#funkcje działają, nad zakresami jeszcze pomyślę

module blocksys

include("IOfunctions.jl")
include("matrixgen.jl")
import SparseArrays
import LinearAlgebra

export classicGauss, modifiedGauss, modifiedGaussLU, modifiedGaussWithChoice, modifiedGaussWithChoiceLU, solveWithGauss, solveWithLU, solveWithChoiceGauss, solveWithLUChoice, LU, LUChoice

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
    Parameters:
        A - matrix (as SparseArray)
        n - size of matrix
        l - size of one block
        b - right side Vector
    Result:
        B - modified (Gauss) matrix
        b - modified right side vector
"""
function modifiedGauss(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, bv::Vector{Float64})
    B = copy(A)
    b = copy(bv)
    count = 0
    #loop over columns
    for k in 1:n-1
        #loop over rows 
        for i in k+1 : min(n, k+l+1)
            multiplier = B[i,k] / B[k,k]
            B[i,k] = Float64(0.0)
            #for each cell change others in the row
            for j in k+1 : min(n, k+l+1)
                B[i,j] -=  multiplier * B[k,j]
            end
            #when modifing matrix, right side vector must be modified too
            b[i] -= multiplier * b[k]
        end
    end
    return B, b
end

"""
Function computing Gauss matrix (upper triangular) from matrix A with choice of main element
    Parameters:
        A - matrix (as SparseArray)
        n - size of matrix
        l - size of one block
        b - right side Vector
    Result:
        B - modified (Gauss) matrix
        b - modified right side vector
"""
function modifiedGaussWithChoice(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, bv)
    B = copy(A)
    b = copy(bv)
    #permutations array
    p = [1:n;]
   
    #itereate over columns
    for k in 1:n-1

        #get max value in each column (from rows k to n)
        maxInCol=0
        indexMax = 0
        
        #over rows
        for i in k: min(n, k+l+1)
            if abs(B[p[i], k]) > maxInCol
                maxInCol = abs(B[p[i], k])
                indexMax = i
            end
        end

        #swap elements in permutation vector
        p[indexMax], p[k] = p[k], p[indexMax]

        #downto in k-column
        for i in k+1 : min(n, k+l+1)

            multiplier = B[p[i],k] / B[p[k],k]
            B[p[i],k] = Float64(0.0)
            #modify values over the row
            for j in k+1 : min(n, k+2*l)
                B[p[i],j] -= multiplier * B[p[k],j]
            end
            #modify values in right side vector
            b[p[i]] -= multiplier * b[p[k]]
        end
    end
    return B, p, b
end

"""
Function solving Ax=b with Gauss elimination
    Parameters:
        A - matrix (as SparseArray)
        n - size of matrix
        l - size of one block
        b - right side Vector
    Result:
        x - result vector
"""
function solveWithGauss(n::Int64,A::SparseArrays.SparseMatrixCSC{Float64, Int64}, bv::Vector{Float64}, l::Int64) :: Vector{Float64}
    gauss = modifiedGauss(A, n, l, bv)
    # gauss matrix
    B = gauss[1]
    #modified b vector
    b = gauss[2]
    x = zeros(Float64, n)
    #iterate over matrix to compute x-vector
    for i in n:-1:1
        sum = 0
        for j in i+1 : min(n, i+l+2)
            sum += B[i,j] * x[j]
        end
        x[i] = (b[i] - sum) / B[i,i]
    end
    return x
end


"""
Function solving Ax=b with Gauss elimination with choice of main element, using implicitly PA=LU
    Parameters:
        A - matrix (as SparseArray)
        n - size of matrix
        l - size of one block
        b - right side Vector
    Result:
        x - result vector
"""
function solveWithChoiceGauss(n::Int64, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, bv::Vector{Float64}, l::Int64)
    gauss = modifiedGaussWithChoice(A, n, l, bv)
    B = gauss[1]
    p = gauss[2]
    b = gauss[3]
    x=zeros(Float64,n)

    # solves Lz=Pb
    for k in 1: n-1
        for i in k+1 : min(n, k+2*l)
            b[p[i]] -= B[p[i],k] * b[p[k]]
        end
    end
    #solves Ux=z
    for i in n:-1:1
        sum = 0
        for j in i+1 : min(n, i+2*l)
            sum+=B[p[i],j] * x[j]
        end
        x[i] = (b[p[i]]-sum)/B[p[i],i]
    end
    return x
end
"""
Function generating matrices L and U after Gauss elimination
    Parameters:
        A - matrix as SparseArray
        n - size od matrix
        l - size of block ( l|n )
    Result:
        L (lower triangular) and U (upper triangular) matrices
"""
function modifiedGaussLU(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    U = copy(A)
    L = SparseArrays.spzeros(n,n)
    
    #loop over columns
    for k in 1:n-1
        L[k,k]=1.0
        #loop over rows
        for i in k+1 : min(n, k+l+1)
            multiplier = U[i,k] / U[k,k]
            L[i,k] = multiplier 
            U[i,k] = 0.0
            #for each cell change others in row
            for j in k+1 : min(n, k+2*l)
                U[i,j] -=  multiplier*U[k,j]
            end
        end
    end
    L[n,n]=1
    return L, U
end

"""
Function generating matrices L and U after Gauss elimination with choice of main element
    Parameters:
        A - matrix as SparseArray
        n - size od matrix
        l - size of block ( l|n )
    Result:
        L (lower triangular) and U (upper triangular) matrices
"""
function modifiedGaussWithChoiceLU(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    U = copy(A)
    L = SparseArrays.spzeros(n,n)
    #permutations array
    p = [1:n;]
   
    #itereate over columns
    for k in 1:n-1
        #get max value in each column (from rows k to n)
        maxInCol=0
        indexMax = 0
        
        #find max element in column 
        for i in k: min(n, k+l+1)
            if abs(U[p[i], k]) > maxInCol
                maxInCol = abs(U[p[i], k])
                indexMax = i
            end
        end
        #swap elements in permutation vector
        p[indexMax], p[k] = p[k], p[indexMax]
        
        #downto in k-column
        for i in k+1 : min(n, k+l+1)
            multiplier = U[p[i],k] / U[p[k],k]

            L[p[i],k] = multiplier
            U[p[i],k] = Float64(0.0)

            for j in k+1 : min(n, k + 2*l)
                U[p[i],j] -= multiplier * U[p[k],j]
            end
        end
    end
    return  L, U, p
end

"""
Function resolving Ax=b, using using L and U matrices
    Parameters:
        L - lower triangular matrix
        U - upper triangular matrix
        n - size of matrix
        l - size of one block
        b - right side Vector
    Result:
        x - result vector
"""
function solveWithLU(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, bv::Vector{Float64})
    x = zeros(Float64,n)
    b = copy(bv)
    #solves Lz = b
    for k in 1: n-1
        for i in k+1 : min(n, k+l+1)
            b[i] -= L[i,k] * b[k]
        end
    end
    #solves Ux = z
    for i in n:-1:1
        sum = 0
        for j in i+1 : min(n, i+l)
            sum +=U[i,j] * x[j]
        end
        x[i] = (b[i]-sum) / U[i,i]
    end
    return x
end

"""
Function resolving Ax=b, using using L and U matrices and p permutation vector
    Parameters:
        L - lower triangular matrix
        U - upper triangular matrix
        p - permutation vector
        n - size of matrix
        l - size of one block
        b - right side Vector
    Result:
        x - result vector
"""
function solveWithLUChoice(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, p::Vector{Int64}, n::Int64, l::Int64, bv::Vector{Float64})
    b = copy(bv)
    x = zeros(Float64,n)
    # Lz = Pb
    for k in 1: n-1
        for i in k+1 : min(n, k+2*l+5)
            b[p[i]] -= L[p[i],k] * b[p[k]]
        end
    end
    # Ux = z
    for i in n:-1:1
        sum = 0
        for j in i+1 : min(n, i+2*l)
            sum += U[p[i],j] * x[j]
        end
        x[i] = (b[p[i]]-sum) / U[p[i],i]
    end
    return x
end

function LU(n::Int64, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, bv::Vector{Float64}, l::Int64)
    LU = modifiedGaussLU(A, n, l)
    x = solveWithLU(LU[1], LU[2], n, l, bv)
    return x
end

function LUChoice(n::Int64, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, bv::Vector{Float64}, l::Int64)
    LU = modifiedGaussWithChoiceLU(A, n, l)
    x = solveWithLUChoice(LU[1], LU[2], LU[3], n, l, bv)
    return x
end

end

# E = IOfunctions.readMatrix("16/A.txt")
# f = IOfunctions.readRightSideVector("16/b.txt")

# f2 = IOfunctions.rightSideFromMatrix(E, 16, 4)

# # LU = modifiedGaussWithChoiceLU(E, 16, 4)
# # k = solveWithLUChoice(LU[1], LU[2], LU[3], 16, 4, f)
# # println(k)
# # A = modifiedGauss(E, 16, 4, f)
# # X = modifiedGauss(E,16, 4)
# x = modifiedGaussLU(E, 16, 4)
# println(solveWithLU(x[1], x[2], 16, 4, f))
# # println(solveWithGauss(16, E, f, 4))
# # Y = solveWithLU(X[1], X[2], 16, 4, f)
# # IOfunctions.saveVectorToFileDelta(Y, "wektor.txt")
# # println(Array(X[1])*Array(X[2]))
# # println(Array(Y))
# # println(f)
# # y = modifiedGaussWithChoiceLU(E, 50000, 4)
# # x = solveWithLUChoice(y[1], y[2],y[3], 50000, 4, f)
# # println(x)

