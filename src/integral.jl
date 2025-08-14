################################################################################
#                         Master Integral Function
################################################################################

"""
    integral(f, geometry[, rule]; kwargs...)

Numerically integrate a given function `f(::Point)` over the domain defined by
a `geometry` using a particular numerical integration `rule` with floating point
precision of type `FP`.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `geometry`: a Meshes.jl `Geometry` or `Domain` that defines the integration domain
- `rule`: optionally, the `IntegrationRule` used for integration (by default
`GaussKronrod()` in 1D and `HAdaptiveCubature()` else)

# Keyword Arguments
- `diff_method::DifferentiationMethod`: manually specifies the differentiation method use to
calculate Jacobians within the integration domain.
- `FP = Float64`: manually specifies the desired output floating point precision
"""
function integral end

# Default integration rule to use if unspecified
function _default_rule(geometry)
    ifelse(Meshes.paramdim(geometry) == 1, GaussKronrod(), HAdaptiveCubature())
end

# If only f and geometry are specified, select default rule
function integral(
        f,
        geometry::Geometry,
        rule::I = _default_rule(geometry);
        kwargs...
) where {I <: IntegrationRule}
    _integral(f, geometry, rule; kwargs...)
end

function integral(
        f,
        domain::Meshes.Domain,
        rule::I = _default_rule(domain);
        kwargs...
) where {I <: IntegrationRule}
    # Discretize the Domain into primitive geometries, sum the integrals over those
    subintegral(geometry) = integral(f, geometry, rule; kwargs...)
    subgeometries = collect(Meshes.elements(Meshes.discretize(domain)))
    return sum(subintegral, subgeometries)
end

################################################################################
#                            Integral Workers
################################################################################

# GaussKronrod
function _integral(
        f,
        geometry,
        rule::GaussKronrod;
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(geometry, FP)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _check_diff_method_support(geometry, diff_method)

    # Only supported for 1D geometries
    if Meshes.paramdim(geometry) != 1
        msg = "GaussKronrod rules not supported in more than one parametric dimension."
        throw(ArgumentError(msg))
    end

    integrand(t) = f(geometry(t)) * differential(geometry, (t,), diff_method)
    return QuadGK.quadgk(integrand, zero(FP), one(FP); rule.kwargs...)[1]
end

# GaussLegendre
function _integral(
        f,
        geometry,
        rule::GaussLegendre;
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(geometry, FP)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _check_diff_method_support(geometry, diff_method)

    N = Meshes.paramdim(geometry)

    # Get Gauss-Legendre nodes and weights of type FP for a region [-1,1]ᴺ
    xs = Iterators.map(FP, rule.nodes)
    ws = Iterators.map(FP, rule.weights)
    weight_grid = Iterators.product(ntuple(Returns(ws), N)...)
    node_grid = Iterators.product(ntuple(Returns(xs), N)...)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = (1 // 2) * x + (1 // 2)

    function integrand((weights, nodes))
        # ts = t.(nodes), but non-allocating
        ts = ntuple(i -> t(nodes[i]), length(nodes))
        # Integrand function
        prod(weights) * f(geometry(ts...)) * differential(geometry, ts, diff_method)
    end

    return (1 // (2^N)) .* sum(integrand, zip(weight_grid, node_grid))
end

# HAdaptiveCubature
function _integral(
        f,
        geometry,
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(geometry, FP)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _check_diff_method_support(geometry, diff_method)

    N = Meshes.paramdim(geometry)

    integrand(ts) = f(geometry(ts...)) * differential(geometry, ts, diff_method)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = _zeros(FP, N)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(ts) = Unitful.ustrip.(integrandunits, integrand(ts))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, _zeros(FP, N), _ones(FP, N); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
