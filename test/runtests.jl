using FiniteElements
using Test

@testset "FiniteElements.jl" begin
    a = Abc(1.2)
    @test a.hallo == 1.2
end
@testset "FiniteElements2" begin
    a = Abc(1.2)
    @test_broken a.hallo == 1.3
end
