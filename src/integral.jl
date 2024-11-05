################################################################################
#                         Master Integral Function
################################################################################

"""
    integral(f, geometry[, rule]; diff_method=_default_method(geometry,FP), FP=Float64)

Numerically integrate a given function `f(::Point)` over the domain defined by
a `geometry` using a particular numerical integration `rule` with floating point
precision of type `FP`.

# Arguments
- `f`: an integrand function with a method `f(::Meshes.Point)`
- `geometry`: some `Meshes.Geometry` that defines the integration domain
- `rule`: optionally, the `IntegrationRule` used for integration (by default
`GaussKronrod()` in 1D and `HAdaptiveCubature()` else)

# Keyword Arguments
- `diff_method::DifferentiationMethod = _default_method(geometry,FP)`: the method to
use for calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired.
"""
function integral end

# If only f and geometry are specified, select default rule
function integral(
        f::Function,
        geometry::Geometry,
        rule::I = Meshes.paramdim(geometry) == 1 ? GaussKronrod() : HAdaptiveCubature();
        kwargs...
) where {I <: IntegrationRule}
    _integral(f, geometry, rule; kwargs...)
end

################################################################################
#                    Generalized (n-Dimensional) Worker Methods
################################################################################

# GaussKronrod
function _integral(
        f,
        geometry,
        rule::GaussKronrod;
        kwargs...
)
    # Pass through to dim-specific workers in next section of this file
    N = Meshes.paramdim(geometry)
    if N == 1
        return _integral_gk_1d(f, geometry, rule; kwargs...)
    elseif N == 2
        return _integral_gk_2d(f, geometry, rule; kwargs...)
    else
        _error_unsupported_combination("geometry with more than two parametric dimensions",
            "GaussKronrod")
    end
end

# GaussLegendre
function _integral(
        f,
        geometry,
        rule::GaussLegendre;
        FP::Type{T} = Float64,
        diff_method::DM = _default_method(geometry)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    N = Meshes.paramdim(geometry)

    # Get Gauss-Legendre nodes and weights for a region [-1,1]^N
    xs, ws = _gausslegendre(FP, rule.n)
    weights = Iterators.product(ntuple(Returns(ws), N)...)
    nodes = Iterators.product(ntuple(Returns(xs), N)...)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = FP(1 // 2) * x + FP(1 // 2)

    function integrand((weights, nodes))
        ts = t.(nodes)
        prod(weights) * f(geometry(ts...)) * differential(geometry, ts, diff_method)
    end

    return FP(1 // (2^N)) .* sum(integrand, zip(weights, nodes))
end

# HAdaptiveCubature
function _integral(
        f,
        geometry,
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64,
        diff_method::DM = _default_method(geometry)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    N = Meshes.paramdim(geometry)

    integrand(ts) = f(geometry(ts...)) * differential(geometry, ts, diff_method)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zeros(FP, N)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(ts) = Unitful.ustrip.(integrandunits, integrand(ts))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zeros(FP, N), ones(FP, N); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                    Specialized GaussKronrod Methods
################################################################################

function _integral_gk_1d(
        f,
        geometry,
        rule::GaussKronrod;
        FP::Type{T} = Float64,
        diff_method::DM = _default_method(geometry)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    integrand(t) = f(geometry(t)) * differential(geometry, (t,), diff_method)
    return QuadGK.quadgk(integrand, zero(FP), one(FP); rule.kwargs...)[1]
end

function _integral_gk_2d(
        f,
        geometry2d,
        rule::GaussKronrod;
        FP::Type{T} = Float64,
        diff_method::DM = _default_method(geometry2d)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    integrand(u, v) = f(geometry2d(u, v)) * differential(geometry2d, (u, v), diff_method)
    ∫₁(v) = QuadGK.quadgk(u -> integrand(u, v), zero(FP), one(FP); rule.kwargs...)[1]
    return QuadGK.quadgk(v -> ∫₁(v), zero(FP), one(FP); rule.kwargs...)[1]
end
