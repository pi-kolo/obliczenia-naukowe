include("functions.jl")
include("IOfunctions.jl")
include("matrixgen.jl")
using SparseArrays
# using blockmat


# A16 = IOfunctions.readMatrix("16/A.txt")
# b16 = IOfunctions.readRightSideVector("16/b.txt")

# matrixgen.blockmat(100, 10, 10.0, "100/A.txt")
# A100 = IOfunctions.readMatrix("100/A.txt")
# b100 = IOfunctions.rightSideFromMatrix(A100, 100, 10)

# println(blocksys.LUChoice(100, A100, b100, 10))

matrixgen.blockmat(1000, 4, 10.0, "1000/A.txt")
A1000 = IOfunctions.readMatrix("1000/A.txt")
b1000 = IOfunctions.rightSideFromMatrix(A1000, 1000, 4)

# A10000 = IOfunctions.readMatrix("10000/A.txt")
# b10000 = IOfunctions.readRightSideVector("10000/b.txt")

A50000 = IOfunctions.readMatrix("50000/A.txt")
b50000 = IOfunctions.readRightSideVector("50000/b.txt")

g  = blocksys.solveWithGauss(1000, A1000, b1000, 4)
g2 = blocksys.LU(1000, A1000, b1000, 4)


h = blocksys.solveWithChoiceGauss(50000, A50000, b50000, 4)
# b162 = IOfunctions.rightSideFromMatrix(A10000, 10000, 4)

function benchmark(n, l, repeats, func)
    totalTime = 0
    totalMemory = 0
    for i in 1:repeats
        matrixgen.blockmat(n, l, 10.0, "test/A.txt")
        A= IOfunctions.readMatrix("test/A.txt")
        b= IOfunctions.rightSideFromMatrix(A, n, l)
        (_, time, memory) = @timed func(n, A, b, l)
        totalTime += time
        totalMemory += memory
    end
    println(n, ";",totalTime/repeats, ";", totalMemory/repeats)
end

functions = [blocksys.solveWithGauss, blocksys.solveWithChoiceGauss, blocksys.LU, blocksys.LUChoice]

for i in 100:100:1000
    tt=0
    tm=0
    for j in 1:10
        matrixgen.blockmat(i, 4, 10.0, "test/A.txt")
        A= IOfunctions.readMatrix("test/A.txt")
        b= IOfunctions.rightSideFromMatrix(A, i, 4)
        (_, time, memory) = @timed Array(A)\Array(b)
        tt+=time
        tm+=memory
    end
    println("$i;$(tt/10);$(tm/10)")
end

for j in 1:4
    println("Algorytm $j")
    for k in 100:100:1000
        benchmark(k, 4, 30, functions[j])
    end
    for k in 5000:2500:40000
        benchmark(k, 4, 30, functions[j])
    end
end
