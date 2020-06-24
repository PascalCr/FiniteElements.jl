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
@testset "PlattenConstructor" begin
    @test Platte(1.0, 1.0, 10, 10, 7874.0, 80.0, 449.0, 273.15) != nothing
end

@testset "Tempertaturiteration" begin
    p = Platte(1.0, 1.0, 10, 10, 7874.0, 80.0, 449.0, 273.15)
    pneu = tempiter(p, 1.0)
    @test pneu.temp ≈ p.temp
end
