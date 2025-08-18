function _alias_error_msg(name, N)
    "$name integrals not supported on a geometry without exactly $N parametric dimensions."
end

################################################################################
#                              Line Integral
################################################################################

"""
    lineintegral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular numerical integration `rule` with floating point precision of
type `FP`.

This is a convenience wrapper around [`integral`](@ref) that additionally enforces
a requirement that the geometry have one parametric dimension.

Rule types available:
- [`GaussKronrod`](@ref) (default)
- [`GaussLegendre`](@ref)
- [`HAdaptiveCubature`](@ref)
"""
function lineintegral(
        f,
        geometry::GeometryOrDomain,
        rule::IntegrationRule = GaussKronrod();
        kwargs...
)
    if !Meshes.iscurve(geometry)
        throw(ArgumentError(_alias_error_msg("Line", 1)))
    end

    return integral(f, geometry, rule; kwargs...)
end

################################################################################
#                              Surface Integral
################################################################################

"""
    surfaceintegral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` along a surface `geometry`
using a particular numerical integration `rule` with floating point precision of
type `FP`.

This is a convenience wrapper around [`integral`](@ref) that additionally enforces
a requirement that the geometry have two parametric dimensions.

Algorithm types available:
- [`GaussKronrod`](@ref)
- [`GaussLegendre`](@ref)
- [`HAdaptiveCubature`](@ref) (default)
"""
function surfaceintegral(
        f,
        geometry::GeometryOrDomain,
        rule::IntegrationRule = HAdaptiveCubature();
        kwargs...
)
    if !Meshes.issurface(geometry)
        throw(ArgumentError(_alias_error_msg("Surface", 2)))
    end

    return integral(f, geometry, rule; kwargs...)
end

################################################################################
#                              Volume Integral
################################################################################

"""
    volumeintegral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular numerical integration `rule` with floating point
precision of type `FP`.

This is a convenience wrapper around [`integral`](@ref) that additionally enforces
a requirement that the geometry have three parametric dimensions.

Algorithm types available:
- [`GaussKronrod`](@ref)
- [`GaussLegendre`](@ref)
- [`HAdaptiveCubature`](@ref) (default)
"""
function volumeintegral(
        f,
        geometry::GeometryOrDomain,
        rule::IntegrationRule = HAdaptiveCubature();
        kwargs...
)
    if !Meshes.issolid(geometry)
        throw(ArgumentError(_alias_error_msg("Volume", 3)))
    end

    return integral(f, geometry, rule; kwargs...)
end
