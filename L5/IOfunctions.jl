
module IOfunctions
import SparseArrays
export readRightSideVector, readMatrix

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

end