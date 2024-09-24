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
  FP::Type{T}=Float64
) where {F<:Function,T<:AbstractFloat}
  # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
  xs, ws = _gausslegendre(FP, settings.n)
  wws = Iterators.product(ws, ws)
  xxs = Iterators.product(xs, xs)

  # Domain transformations:
  #   xᵢ [-1,1] ↦ R [0,1]
  #   xⱼ [-1,1] ↦ φ [0,π/2]
  uR(xᵢ) = T(1 / 2) * (xᵢ + 1)
  uφ(xⱼ) = T(π / 4) * (xⱼ + 1)

  # Integrate the Barycentric triangle by transforming it into polar coordinates
  #   with a modified radius
  #     R = r ( sinφ + cosφ )
  #   s.t. integration bounds become rectangular
  #     R ∈ [0, 1] and φ ∈ [0, π/2]
  function integrand(((wᵢ, wⱼ), (xᵢ, xⱼ)))
    R = uR(xᵢ)
    φ = uφ(xⱼ)
    a, b = sincos(φ)
    u = R * (1 - a / (a + b))
    v = R * (1 - b / (a + b))
    return wᵢ * wⱼ * f(triangle(u, v)) * R / (a + b)^2
  end

  # Calculate 2D Gauss-Legendre integral over modified-polar-Barycentric coordinates
  # Apply a linear domain-correction factor
  return FP(π / 4) * area(triangle) .* sum(integrand, zip(wws, xxs))
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
  FP::Type{T}=Float64
) where {F<:Function,T<:AbstractFloat}
  # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
  #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
  inner∫(u) = QuadGK.quadgk(v -> f(triangle(u, v)), FP(0), FP(1 - u); settings.kwargs...)[1]
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
  FP::Type{T}=Float64
) where {F<:Function,T<:AbstractFloat}
  # Integrate the Barycentric triangle by transforming it into polar coordinates
  #   with a modified radius
  #     R = r ( sinφ + cosφ )
  #   s.t. integration bounds become rectangular
  #     R ∈ [0, 1] and φ ∈ [0, π/2]
  function integrand(Rφ)
    R, φ = Rφ
    a, b = sincos(φ)
    u = R * (1 - a / (a + b))
    v = R * (1 - b / (a + b))
    return f(triangle(u, v)) * R / (a + b)^2
  end
  intval = HCubature.hcubature(integrand, FP[0, 0], FP[1, π / 2], settings.kwargs...)[1]

  # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
  return 2 * area(triangle) .* intval
end
