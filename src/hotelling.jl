# Hotelling T² test

abstract type HotellingT2Test <: HypothesisTest end

HypothesisTests.default_tail(::HotellingT2Test) = :right
HypothesisTests.pvalue(T::HotellingT2Test; tail=:right) =
    pvalue(FDist(dof(T)...), T.F, tail=tail)

function HypothesisTests.show_params(io::IO, T::HotellingT2Test, indent="")
    println(io, indent, "number of observations: ", nobs(T))
    println(io, indent, "number of variables:    ", T.p)
    println(io, indent, "T² statistic:           ", T.T²)
    println(io, indent, "transformed statistic:  ", T.F)
    println(io, indent, "degrees of freedom:     ", dof(T))
    println(io, indent, "covariance matrix:")
    Base.print_matrix(io, T.S, indent^2)
end

## One sample

struct OneSampleHotellingT2 <: HotellingT2Test
    T²::Real
    F::Real
    n::Int
    p::Int
    μ₀::Vector
    x̄::Vector
    S::Matrix
end

StatsBase.nobs(T::OneSampleHotellingT2) = T.n
StatsBase.dof(T::OneSampleHotellingT2) = (T.p, T.n - T.p)

function OneSampleHotellingT2(X::AbstractMatrix{T},
                              μ₀::AbstractVector=fill(middle(zero(T)), size(X, 2))) where T
    n, p = size(X)
    p == length(μ₀) ||
        throw(DimensionMismatch("Number of variables does not match number of means"))
    n > 0 || throw(ArgumentError("The input must be non-empty"))
    x̄ = vec(mean(X, 1))
    S = cov(X)
    T² = n * At_Binv_A(x̄ .- μ₀, S)
    F = (n - p) * T² / (p * (n - 1))
    return OneSampleHotellingT2(T², F, n, p, μ₀, x̄, S)
end

# paired
function OneSampleHotellingT2(X::AbstractMatrix{T}, Y::AbstractMatrix{S},
                              μ₀::AbstractVector=fill(middle(zero(T), zero(S)), size(X, 2))) where {T,S}
    p, nx, ny = checkdims(X, Y)
    nx == ny || throw(DimensionMismatch("Inconsistent number of observations: $nx, $ny"))
    p == length(μ₀) ||
        throw(DimensionMismatch("Number of variables does not match number of means"))
    return OneSampleHotellingT2(X .- Y, μ₀)
end

HypothesisTests.testname(::OneSampleHotellingT2) = "One sample Hotelling's T² test"
HypothesisTests.population_param_of_interest(T::OneSampleHotellingT2) =
    ("Mean vector", T.μ₀, T.x̄)

## Two sample

## Two sample: equal covariance

struct EqualCovHotellingT2 <: HotellingT2Test
    T²::Real
    F::Real
    nx::Int
    ny::Int
    p::Int
    Δ::Vector
    S::Matrix
end

function EqualCovHotellingT2(X::AbstractMatrix, Y::AbstractMatrix)
    p, nx, ny = checkdims(X, Y)
    Δ = vec(mean(X, 1) .- mean(Y, 1))
    S = _poolcov!(cov(X), nx - 1, cov(Y), ny - 1) .* (inv(nx) .+ inv(ny))
    T² = At_Binv_A(Δ, S)
    F = T² * (nx + ny - p - 1) / (p * (nx + ny - 2))
    return EqualCovHotellingT2(T², F, nx, ny, p, Δ, S)
end

StatsBase.nobs(T::EqualCovHotellingT2) = (T.nx, T.ny)
StatsBase.dof(T::EqualCovHotellingT2) = (T.p, T.nx + T.ny - T.p - 1)

HypothesisTests.testname(::EqualCovHotellingT2) =
    "Two sample Hotelling's T² test (equal covariance matrices)"
HypothesisTests.population_param_of_interest(T::EqualCovHotellingT2) =
    ("Difference in mean vectors", zeros(eltype(T.Δ), T.p), T.Δ)

## Two sample: unequal covariance

# Store the denominator degrees of freedom in the type, since the computation
# is expensive and we don't want to redo it every time the user calls dof
struct UnequalCovHotellingT2 <: HotellingT2Test
    T²::Real
    F::Real
    nx::Int
    ny::Int
    p::Int
    ν::Int
    Δ::Vector
    S::Matrix
end

function UnequalCovHotellingT2(X::AbstractMatrix, Y::AbstractMatrix)
    p, nx, ny = checkdims(X, Y)
    Δ = vec(mean(X, 1) .- mean(Y, 1))
    Sx = cov(X) ./ nx
    Sy = cov(Y) ./ ny
    ST = Sx .+ Sy
    T² = At_Binv_A(Δ, ST)
    F = (nx + ny - p - 1) * T² / (p * (nx + ny - 2))
    tmp = Symmetric(ST) \ Δ
    iν = (dot(tmp, Sx * tmp) / T²)^2 / (nx - 1) + (dot(tmp, Sy * tmp) / T²)^2 / (ny - 1)
    ν = trunc(Int, inv(iν))
    return UnequalCovHotellingT2(T², F, nx, ny, p, ν, Δ, ST)
end

StatsBase.nobs(T::UnequalCovHotellingT2) = (T.nx, T.ny)
StatsBase.dof(T::UnequalCovHotellingT2) = (T.p, T.ν)

HypothesisTests.testname(::UnequalCovHotellingT2) =
    "Two sample Hotelling's T² test (unequal covariance matrices)"
HypothesisTests.population_param_of_interest(T::UnequalCovHotellingT2) =
    ("Difference in mean vectors", zeros(eltype(T.Δ), T.p), T.Δ)
