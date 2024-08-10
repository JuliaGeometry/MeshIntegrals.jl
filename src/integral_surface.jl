################################################################################
#                          Generalized 2D Methods
################################################################################

function _integral_2d(
    f,
    geometry2d,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = FP(1/2) * x + FP(1/2)
    point(xi, xj) = geometry2d(t(xi), t(xj))

    # Integrate f over the geometry
    integrand(((wi,wj), (xi,xj))) = wi * wj * f(point(xi,xj)) * differential(geometry2d, t.([xi, xj]))
    return FP(1/4) .* sum(integrand, zip(wws,xxs))
end

function _integral_2d(
    f,
    geometry2d,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    integrand(u,v) = f(geometry2d(u,v)) * differential(geometry2d, [u,v])
    ∫₁(v) = QuadGK.quadgk(u -> integrand(u,v), FP(0), FP(1); settings.kwargs...)[1]
    return QuadGK.quadgk(v -> ∫₁(v), FP(0), FP(1); settings.kwargs...)[1]
end

function _integral_2d(
    f,
    geometry,
    settings::HAdaptiveCubature,
)
    return _integral(f, geometry, settings)
end

function _integral_2d(
    f,
    geometry,
    settings::HAdaptiveCubature,
    FP::Type{T}
) where {T<:AbstractFloat}
    return _integral(f, geometry, settings, FP)
end


################################################################################
#                  Specialized Methods for CylinderSurface
################################################################################

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a CylinderSurface with GaussLegendre not supported.")
    # TODO Planned to support in the future
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Integrate the curved sides of the CylinderSurface
    # \int ( \int f(r̄) dz ) dφ
    function sides_inner∫(φ)
        sidelength = norm(cyl(φ,FP(1)) - cyl(φ,FP(0)))
        return sidelength * QuadGK.quadgk(z -> f(cyl(φ,z)), FP(0), FP(1); settings.kwargs...)[1]
    end
    sides = (FP(2π) * cyl.radius) .* QuadGK.quadgk(φ -> sides_inner∫(φ), FP(0), FP(1); settings.kwargs...)[1]

    # Integrate the Disk at the top of the CylinderSurface
    disk_top = Meshes.Disk(cyl.top, cyl.radius)
    top = _integral_2d(f, disk_top, settings, FP)

    # Integrate the Disk at the bottom of the CylinderSurface
    disk_bottom = Meshes.Disk(cyl.bot, cyl.radius)
    bottom = _integral_2d(f, disk_bottom, settings, FP)

    return sides + top + bottom
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    Dim = Meshes.paramdim(cyl)

    integrand(t) = f(cyl(t...)) * differential(cyl, t)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = fill(FP(0.5),Dim)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zeros(FP,Dim), ones(FP,Dim); settings.kwargs...)[1]
    # Reapply units
    sides = value .* integrandunits

    # Integrate the Disk at the top of the CylinderSurface
    disk_top = Meshes.Disk(cyl.top, cyl.radius)
    top = _integral_2d(f, disk_top, settings, FP)

    # Integrate the Disk at the bottom of the CylinderSurface
    disk_bottom = Meshes.Disk(cyl.bot, cyl.radius)
    bottom = _integral_2d(f, disk_bottom, settings, FP)

    return sides + top + bottom
end

################################################################################
#                      Specialized Methods for Plane
################################################################################

function integral(
    f::F,
    plane::Meshes.Plane,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(plane(0,0))
    integrand(((wi,wj), (xi,xj))) = wi * wj * f(plane(t(xi), t(xj))) * t′(xi) * t′(xj) * domainunits^2
    return sum(integrand, zip(wws,xxs))
end

function integral(
    f::F,
    plane::Meshes.Plane,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Integrate f over the Plane
    domainunits = _units(plane(0,0))
    inner∫(v) = QuadGK.quadgk(u -> f(plane(u,v)) * domainunits, FP(-Inf), FP(Inf); settings.kwargs...)[1]
    return QuadGK.quadgk(v -> inner∫(v) * domainunits, FP(-Inf), FP(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    plane::Meshes.Plane,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(plane(0,0))
    integrand(x::AbstractVector) = f(plane(t(x[1]), t(x[2]))) * t′(x[1]) * t′(x[2]) * domainunits^2

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = FP[0.5, 0.5]
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[-1,-1], FP[1,1]; settings.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                    Specialized Methods for Triangle
################################################################################

"""
    integral(f, triangle::Meshes.Triangle, ::GaussLegendre)

Like [`integral`](@ref) but integrates over the surface of a `triangle`
by transforming the triangle into a polar-barycentric coordinate system and
using a Gauss-Legendre quadrature rule along each barycentric dimension of the
triangle.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3},
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformations:
    #   xᵢ [-1,1] ↦ R [0,1]
    #   xⱼ [-1,1] ↦ φ [0,π/2]
    uR(xᵢ) = T(1/2) * (xᵢ + 1)
    uφ(xⱼ) = T(π/4) * (xⱼ + 1)

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(((wᵢ,wⱼ), (xᵢ,xⱼ)))
        R = uR(xᵢ)
        φ = uφ(xⱼ)
        a,b = sincos(φ)
        u = R * (1 - a / (a + b))
        v = R * (1 - b / (a + b))
        return wᵢ * wⱼ * f(triangle(u, v)) * R / (a + b)^2
    end

    # Calculate 2D Gauss-Legendre integral over modified-polar-Barycentric coordinates
    # Apply a linear domain-correction factor
    return FP(π/4) * area(triangle) .* sum(integrand, zip(wws,xxs))
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` using nested
Gauss-Kronrod quadrature rules along each barycentric dimension of the triangle.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3},
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    inner∫(u) = QuadGK.quadgk(v -> f(triangle(u,v)), FP(0), FP(1-u); settings.kwargs...)[1]
    outer∫ = QuadGK.quadgk(inner∫, FP(0), FP(1); settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * area(triangle) .* outer∫
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` by
transforming the triangle into a polar-barycentric coordinate system and using
an h-adaptive cubature rule.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3},
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(Rφ)
        R,φ = Rφ
        a,b = sincos(φ)
        u = R * (1 - a / (a + b))
        v = R * (1 - b / (a + b))
        return f(triangle(u, v)) * R / (a + b)^2
    end
    intval = HCubature.hcubature(integrand, FP[0, 0], FP[1, π/2], settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * area(triangle) .* intval
end
