
module IOfunctions
import SparseArrays
using LinearAlgebra

export readRightSideVector, readMatrix, rightSideFromMatrix, saveVectorToFile, saveVectorToFileDelta

function readMatrix(file) :: SparseArrays.SparseMatrixCSC{Float64, Int64}
    open(file) do f
        inf = split(readline(f))
        n = parse(Int64, inf[1])
        l = parse(Int64, inf[2])
        rows = []
        columns = []
        values = []
        for ln in eachline(f)
            line = split(ln)
            push!(rows, parse(Int64, line[1]))
            push!(columns, parse(Int64, line[2]))
            push!(values, parse(Float64, line[3]))
        end
        A = SparseArrays.sparse(rows, columns, values)
        return A
    end
end

function readRightSideVector(file) :: Vector{Float64}
    open(file) do f
        n = parse(Int64, readline(f))
        values = []
        for ln in eachline(f)
            push!(values, parse(Float64, ln))
        end
        return values
    end
end

function rightSideFromMatrix(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64) :: Vector{Float64}
    x = ones(Float64, n)
    b = zeros(Float64, n)

    #iterate over rows
    for i in 1 : n
        for j in max(1, i-(2+l)): min(n, i+l)
            b[i] += A[i,j]
        end
    end
    return b
end

function saveVectorToFile(x::Vector{Float64}, filename)
    open(filename, "w") do file
        for i in 1 : length(x)
            write(file, "$(x[i])\n")
        end
    end
end

function saveVectorToFileDelta(x::Vector{Float64}, filename)
    open(filename, "w") do file
        delta = norm(ones(length(x))-x)
        write(file, "$delta\n")
        for i in 1 : length(x)
            write(file, "$(x[i])\n")
        end
    end
end

end