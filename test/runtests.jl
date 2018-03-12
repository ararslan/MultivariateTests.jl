using MultivariateTests
using Compat
using Compat.Test
using Compat.DelimitedFiles
using StatsBase
using HypothesisTests

# Columns are calcium, iron, protein, vitamin A, vitamin C
nutrient = readdlm(joinpath(@__DIR__, "data", "nutrient.txt"))[:,2:end]

# Columns are information, similarities, arithmetic, picture completion
wechsler = readdlm(joinpath(@__DIR__, "data", "wechsler.txt"))[:,2:end]

# Columns are type, length, left, right, bottom, top, diag
swiss = readdlm(joinpath(@__DIR__, "data", "swiss3.txt"))
genuine = convert(Matrix{Float64}, swiss[view(swiss, :, 1) .== "real", 2:end])
counterfeit = convert(Matrix{Float64}, swiss[view(swiss, :, 1) .== "fake", 2:end])

@testset "Utility functions" begin
    MT = MultivariateTests

    @test MT.checkdims(rand(3, 2), rand(4, 2)) == (2, 3, 4)
    # Mismatched number of variables
    @test_throws DimensionMismatch MT.checkdims(rand(4, 3), rand(4, 6))
    # Empty observations
    @test_throws ArgumentError MT.checkdims(rand(0, 4), rand(2, 4))

    Sx = [ 3.0 -1.5 0.0
          -1.5  1.0 0.5
           0.0  0.5 1.0]
    Sy = fill(4.5, (3, 3))
    out = [3.5 0.5      1.5
           0.5 2.166667 1.833333
           1.5 1.833333 2.166667]
    @test MT._poolcov!(Sx, 2, Sy, 1) ≈ out atol=1e-6
    # Input is modified
    @test Sx ≈ out atol=1e-6

    # Positive semi-definite but not positive definite
    P = [ 1.0000  0.7426  0.1601 -0.7000 0.5500
          0.7426  1.0000 -0.2133 -0.5818 0.5000
          0.1601 -0.2133  1.0000 -0.1121 0.1000
         -0.7000 -0.5818 -0.1121  1.0000 0.4500
          0.5500  0.5000  0.1000  0.4500 1.0000]
    # Positive definite
    Q = [ 1.0 -0.5  0.0
         -0.5  1.0 -0.5
          0.0 -0.5  1.0]
    @test @inferred(MT.At_Binv_A(ones(5), P)) ≈ -0.8008792 atol=1e-6
    @test @inferred(MT.At_Binv_A(ones(3), Q)) ≈ 10.0
end

@testset "Generalizations of sample variance" begin
    @test @inferred(genvar(nutrient)) ≈ 2.8310418e19 atol=1e-6
    @test @inferred(totalvar(nutrient)) ≈ 2.83266877e6 atol=1e-6

    X = [1 2 5
         4 1 6
         4 0 4]
    @test @inferred(genvar(X)) ≈ 0.0
    @test @inferred(totalvar(X)) ≈ 5.0

    x = rand(Float32, 10)
    @test genvar(x) == totalvar(x) == var(x)

    it = (xᵢ for xᵢ in x)
    @test genvar(it) == totalvar(it) == var(it)
end

@testset "Partial correlation" begin
    @test @inferred(partialcor(wechsler[:,1], wechsler[:,2], wechsler[:,3:4])) ≈ 0.7118787 atol=1e-6

    w = PartialCorTest(wechsler[:,1], wechsler[:,2], wechsler[:,3:4])
    let out = sprint(show, w)
        @test contains(out, "reject h_0") && !contains(out, "fail to")
    end
    let ci = confint(w)
        @test first(ci) ≈ 0.4963917 atol=1e-6
        @test last(ci) ≈ 0.8447292 atol=1e-6
    end
    @test nobs(w) == 37
    @test dof(w) == 33
    @test pvalue(w) < 0.00001

    X = [ 2 1 0
          4 2 0
         15 3 1
         20 4 1]
    @test @inferred(partialcor(X[:,1], X[:,2], X[:,3])) ≈ 0.919145 atol=1e-6

    x = PartialCorTest(view(X,:,1), view(X,:,2), view(X,:,3))
    @test contains(sprint(show, x), "fail to reject")
    @test confint(x) == (-1.0, 1.0)
    @test nobs(x) == 4
    @test dof(x) == 1
    @test pvalue(x) ≈ 0.25776212 atol=1e-6
end

@testset "One sample Hotelling's T²" begin
    t = OneSampleHotellingT2(nutrient, [1000, 15, 60, 800, 75])
    @test nobs(t) == 737
    @test dof(t) == (5, 732)
    @test pvalue(t) ≈ 0.0 atol=eps()
    @test t.T² ≈ 1758.5413137 atol=1e-6
    @test t.F ≈ 349.7968048 atol=1e-6
    let out = sprint(show, t)
        @test contains(out, "reject h_0") && !contains(out, "fail to")
    end
end

@testset "Two sample Hotelling's T²" begin
    eq = EqualCovHotellingT2(genuine, counterfeit)
    @test nobs(eq) == (100, 100)
    @test dof(eq) == (6, 193)
    @test pvalue(eq) ≈ 0.0 atol=eps()
    @test eq.T² ≈ 2412.4506855 atol=1e-6
    @test eq.F ≈ 391.9217023 atol=1e-6
    @test contains(sprint(show, eq), "reject h_0")

    un = UnequalCovHotellingT2(genuine, counterfeit)
    @test nobs(un) == (100, 100)
    @test dof(un) == (6, 193)
    @test pvalue(un) ≈ 0.0 atol=eps()
    @test un.T² ≈ 2412.4506855 atol=1e-6
    @test un.F ≈ 391.9217023 atol=1e-6
    @test contains(sprint(show, un), "reject h_0")
end

@testset "Bartlett's test" begin
    b = BartlettsTest(genuine, counterfeit)
    @test nobs(b) == (100, 100)
    @test dof(b) == 21
    @test pvalue(b) ≈ 0.0 atol=1e-10
    @test b.L′ ≈ 121.8991235 atol=1e-6
    @test contains(sprint(show, b), "reject h_0")
end
