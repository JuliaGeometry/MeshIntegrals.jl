################################################################################
#                           Integration Rules
################################################################################

function _kwargs_to_string(kwargs)
    return join([string(k) * " = " * string(v) for (k, v) in pairs(kwargs)], ", ")
end

abstract type IntegrationRule end

"""
    GaussKronrod(kwargs...)

The h-adaptive Gauss-Kronrod quadrature rule implemented by
[QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) which can be used for any
single-dimensional geometry. All standard `QuadGK.quadgk` keyword arguments are
supported.
"""
struct GaussKronrod <: IntegrationRule
    kwargs::Base.Pairs
    GaussKronrod(; kwargs...) = new(kwargs)
end

function Base.show(io::IO, rule::GaussKronrod)
    print(io, "GaussKronrod(; ", _kwargs_to_string(rule.kwargs), ")")
end

"""
    GaussLegendre(n)

An `n`'th-order Gauss-Legendre quadrature rule. Nodes and weights are
efficiently calculated using
[FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).

So long as the integrand function can be well-approximated by a polynomial of
order `2n-1`, this method should yield results with 16-digit accuracy in `O(n)`
time. If the function is know to have some periodic content, then `n` should
(at a minimum) be greater than the expected number of periods over the geometry,
e.g. `length(geometry)/λ`.
"""
struct GaussLegendre <: IntegrationRule
    n::Int64
    nodes::Vector{Float64}
    weights::Vector{Float64}

    GaussLegendre(n::Int64) = new(n, FastGaussQuadrature.gausslegendre(n)...)
end

function Base.show(io::IO, rule::GaussLegendre)
    print(io, "GaussLegendre(", rule.n, ")")
end

"""
    HAdaptiveCubature(kwargs...)

The h-adaptive cubature rule implemented by
[HCubature.jl](https://github.com/JuliaMath/HCubature.jl). All standard
`HCubature.hcubature` keyword arguments are supported.
"""
struct HAdaptiveCubature <: IntegrationRule
    kwargs::Base.Pairs
    HAdaptiveCubature(; kwargs...) = new(kwargs)
end

function Base.show(io::IO, rule::HAdaptiveCubature)
    print(io, "HAdaptiveCubature(; ", _kwargs_to_string(rule.kwargs), ")")
end
