#Piotr Kołodziejczyk
#L3 testy
include("functions.jl")
using .AproxFunctions
using Test


@testset "metoda bisekcji" begin
    @test AproxFunctions.mbisekcji(x -> x^2-1, 0.3, 1.2, 0.01, 0.01)[4] == 0
    @test AproxFunctions.mbisekcji(x -> x^2-1, 0.3, 1.2, 0.01, 0.01)[1] ≈ 1.0 atol=0.01
    @test AproxFunctions.mbisekcji(x -> x^2-1, 0.3, 1.2, 0.01, 0.01)[2] ≈ 0.0 atol=0.01

    @test AproxFunctions.mbisekcji(x-> sin(2x)+exp(x)-3x^2, 0.0, 2.0, 0.001, 0.001)[4] == 0
    @test AproxFunctions.mbisekcji(x-> sin(2x)+exp(x)-3x^2, 0.0, 2.0, 0.001, 0.001)[1] ≈ 1.1373 atol=0.001

    @test AproxFunctions.mbisekcji(x-> x^2+1, 0.0, 5.0, 0.1, 0.1)[4] == 1
end

@testset "metoda Newtona" begin
    @test AproxFunctions.mstycznych(x-> x^2-1, x-> 2x, 0.5, 0.01, 0.01, 50)[1] ≈ 1.0 atol=0.01
    @test AproxFunctions.mstycznych(x-> x^2-1, x-> 2x, 0.5, 0.01, 0.01, 50)[2] ≈ 0.0 atol=0.01

    @test AproxFunctions.mstycznych(x-> sin(2x)+exp(x)-3x^2,x-> 2*cos(2x)+exp(x)-6x, 0.0, 0.001, 0.001,50)[1] ≈ -0.2774 atol=0.001
    @test AproxFunctions.mstycznych(x-> sin(2x)+exp(x)-3x^2,x-> 2*cos(2x)+exp(x)-6x, 0.0, 0.001, 0.001,50)[2] ≈ 0.0 atol=0.001

end


@testset "metoda siecznych" begin
    @test AproxFunctions.msiecznych(x -> x^2-1, 0.3, 1.2, 0.01, 0.01, 50)[4] == 0
    @test AproxFunctions.msiecznych(x -> x^2-1, 0.3, 1.2, 0.01, 0.01, 50)[1] ≈ 1.0 atol=0.01
    @test AproxFunctions.msiecznych(x -> x^2-1, 0.3, 1.2, 0.01, 0.01, 50)[2] ≈ 0.0 atol=0.01

    @test AproxFunctions.msiecznych(x-> sin(2x)+exp(x)-3x^2, 0.0, 2.0, 0.001, 0.001, 50)[4] == 0
    @test AproxFunctions.msiecznych(x-> sin(2x)+exp(x)-3x^2, 0.0, 2.0, 0.001, 0.001, 50)[1] ≈ -0.2774 atol=0.001

    @test AproxFunctions.msiecznych(x-> x^2+1, 0.0, 5.0, 0.1, 0.1, 50)[4] == 1
end
