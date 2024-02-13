# Validate that f has a method defined for f(::Point{Dim,T})
@inline function _validate_integrand(f,Dim,T)
    if hasmethod(f, (Point{Dim,T},))
        return nothing
    else
        error("The provided Function f must have a method f(::Meshes.Point{$Dim,$T})")
    end
end

################################################################################
#                         Integration Algorithms
################################################################################

abstract type IntegrationAlgorithm end

"""
    GaussKronrod(kwargs...)

Numerically integrate using the h-adaptive Gauss-Kronrod quadrature rule implemented
by QuadGK.jl. All standard [`QuadGK.quadgk`](@ref) keyword arguments are supported.
"""
struct GaussKronrod <: IntegrationAlgorithm
    kwargs
    GaussKronrod(; kwargs...) = new(kwargs)
end

"""
    GaussLegendre(n)

Numerically integrate using an `n`'th-order Gauss-Legendre quadrature rule. nodes
and weights are efficiently calculated using FastGaussQuadrature.jl.

So long as the integrand function can be well-approximated by a polynomial of
order `2n-1`, this method should yield results with 16-digit accuracy in `O(n)`
time. If the function is know to have some periodic content, then `n` should
(at a minimum) be greater than the expected number of periods over the geometry,
e.g. `length(geometry)/lambda`.
"""
struct GaussLegendre <: IntegrationAlgorithm
    n::Int64
end

"""
    GaussKronrod(kwargs...)

Numerically integrate areas and surfaces using the h-adaptive cubature rule
implemented by HCubature.jl. All standard [`HCubature.hcubature`](@ref) keyword
arguments are supported.
"""
struct HAdaptiveCubature <: IntegrationAlgorithm
    kwargs
    HAdaptiveCubature(; kwargs...) = new(kwargs)
end

################################################################################
#                        Integral Function API
################################################################################

"""
    lineintegral(f, geometry, algorithm::IntegrationAlgorithm=GaussKronrod)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular `integration algorithm`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
"""
function lineintegral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    return lineintegral(f, geometry, GaussKronrod())
end

"""
    surfaceintegral(f, geometry, algorithm::IntegrationAlgorithm=HAdaptiveCubature)

Numerically integrate a given function `f(::Point)` over a surface `geometry`
using a particular `integration algorithm`.
"""
function surfaceintegral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    return surfaceintegral(f, geometry, HAdaptiveCubature())
end

"""
    volumeintegral(f, geometry, algorithm::IntegrationAlgorithm=HAdaptiveCubature)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular `integration algorithm`.

"""
function volumeintegral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    return volumeintegral(f, geometry, HAdaptiveCubature())
end


################################################################################
#                            surfaceintegral
################################################################################

function surfaceintegral(
    f::F,
    circle::Meshes.Circle,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Circle is parametrically 1-dimensional
    # Convert to a Disk, integrate that
    return surfaceintegral(f, Disk(circle.plane, circle.radius), settings)
end

function surfaceintegral(
    f,
    box::Meshes.Box{2,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = box(s(xi), t(xj))

    # Calculate weight-node product
    g(((wi,wj), (xi,xj))) = wi * wj * f(point(xi,xj))

    # Calculate 2D Gauss-Legendre integral of f over parametric coordinates [-1,1]^2
    # Apply a linear domain-correction factor [-1,1]^2 ↦ area(box)
    return 0.25 * area(box) .* sum(g, zip(wws,xxs))
end

function surfaceintegral(
    f,
    box::Meshes.Box{2,T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Integrate the box in parametric (u,v)-space
    innerintegral(u) = QuadGK.quadgk(v -> f(box(u,v)), 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor 1 ↦ area(box)
    return area(box) .* outerintegral
end

function surfaceintegral(
    f,
    disk::Meshes.Disk{T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = disk(s(xi), t(xj))

    # Calculate weight-node product with curvilinear correction
    g(((wi,wj), (xρ,xϕ))) = wi * wj * s(xρ) * f(point(xρ,xϕ))

    # Calculate 2D Gauss-Legendre integral of f over parametric coordinates [-1,1]²
    # Apply curvilinear domain-correction factor [-1,1]² ↦ area(disk)
    return 0.5 * area(disk) .* sum(g, zip(wws,xxs))
end

function surfaceintegral(
    f,
    disk::Meshes.Disk{T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Integrate the box in parametric (ρ,ϕ)-space
    innerintegral(ϕ) = QuadGK.quadgk(ρ -> ρ * f(disk(ρ,ϕ)), 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor π ↦ area(disk)
    return 2π .* outerintegral
end

"""
    surfaceintegral(f, triangle::Meshes.Triangle, ::GaussLegendre)

Like [`surfaceintegral`](@ref) but integrates over the surface of a `triangle`
using a Gauss-Legendre quadrature rule along each barycentric dimension of the
triangle.
"""
function surfaceintegral(
    f,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::GaussLegendre
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = triangle(s(xi), t(xj))

    # Determine output type of f at a Point inside the triangle
    # Define an applicable zero value
    fzero = zero(f(point(-0.5,-0.5)))

    # Calculate weight-node product
    function g(((wi,wj), (xi,xj)))
        if 0 < (s(xi) + t(xj)) < 1
            # Valid coordinate (inside triangle)
            return wi * wj * f(point(xi,xj))
        else
            # Invalid coordinate (outside triangle)
            return fzero
        end
    end

    # Calculate 2D Gauss-Legendre integral of f over Barycentric coordinates [-1,1]^2
    # Apply a linear domain-correction factor [-1,1]^2 ↦ area(triangle)
    return 0.5 * area(triangle) .* sum(g, zip(wws,xxs))
end

# Integrate f(::Point{Dim,T}) over a Triangle{Dim,T}
function surfaceintegral(
    f,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::GaussKronrod
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    innerintegral(u) = QuadGK.quadgk(v -> f(triangle(u,v)), 0, 1-u; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2.0 * area(triangle) .* outerintegral
end
