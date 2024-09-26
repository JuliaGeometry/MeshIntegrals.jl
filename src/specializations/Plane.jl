################################################################################
#                      Specialized Methods for Plane
#
# Why Specialized?
#   The Plane geometry is a special case, representing a planar surface with
#   infinite extent along two basis vectors. This requires a pair of domain
#   transformations mapping from the typical parametric region [0,1]² to an
#   infinite one (-∞,∞)².
################################################################################

function integral(
        f::F,
        plane::Meshes.Plane,
        rule::GaussLegendre,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, rule.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Normalize the Plane's orthogonal vectors
    uu = Meshes.unormalize(plane.u)
    uv = Meshes.unormalize(plane.v)
    plane = Meshes.Plane(plane.p, uu, uv)

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(plane(0, 0))
    function integrand(((wi, wj), (xi, xj)))
        wi * wj * f(plane(t(xi), t(xj))) * t′(xi) * t′(xj) * domainunits^2
    end
    return sum(integrand, zip(wws, xxs))
end

function integral(
        f::F,
        plane::Meshes.Plane,
        rule::GaussKronrod,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    uu = Meshes.unormalize(plane.u)
    uv = Meshes.unormalize(plane.v)
    plane = Meshes.Plane(plane.p, uu, uv)

    # Integrate f over the Plane
    domainunits = _units(plane(0, 0))^2
    integrand(u, v) = f(plane(u, v)) * domainunits
    inner∫(v) = QuadGK.quadgk(u -> integrand(u, v), FP(-Inf), FP(Inf); rule.kwargs...)[1]
    return QuadGK.quadgk(inner∫, FP(-Inf), FP(Inf); rule.kwargs...)[1]
end

function integral(
        f::F,
        plane::Meshes.Plane,
        rule::HAdaptiveCubature,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    uu = Meshes.unormalize(plane.u)
    uv = Meshes.unormalize(plane.v)
    plane = Meshes.Plane(plane.p, uu, uv)

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(plane(0, 0))
    function integrand(x::AbstractVector)
        f(plane(t(x[1]), t(x[2]))) * t′(x[1]) * t′(x[2]) * domainunits^2
    end

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zeros(FP, 2)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, -ones(FP, 2), ones(FP, 2); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
