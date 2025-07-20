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
    if (N = Meshes.paramdim(geometry)) != 1
        throw(ArgumentError("Performing a line integral on a geometry \
                            with $N parametric dimensions not supported."))
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
    if (N = Meshes.paramdim(geometry)) != 2
        throw(ArgumentError("Performing a surface integral on a geometry \
                            with $N parametric dimensions not supported."))
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
    if (N = Meshes.paramdim(geometry)) != 3
        throw(ArgumentError("Performing a volume integral on a geometry \
                            with $N parametric dimensions not supported."))
    end

    return integral(f, geometry, rule; kwargs...)
end
