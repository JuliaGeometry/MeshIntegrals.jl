################################################################################
#                      Specialized Methods for PolyArea
#
# Why Specialized?
#   The PolyArea geometry defines a surface bounded by a polygon, but Meshes.jl
#   cannot define a parametric function for a polyarea because a general solution
#   for one does not exist. However, they can be partitioned into simpler elements
#   and integrated separately before finally summing the result.
################################################################################

"""
    integral(f, area::PolyArea, rule = HAdaptiveCubature();
             diff_method=FiniteDifference(), FP=Float64)

Like [`integral`](@ref) but integrates over the surface domain defined by a `PolyArea`.
The surface is first discretized into facets that are integrated independently using
the specified integration `rule`.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `area`: a `PolyArea` that defines the integration domain
- `rule = HAdaptiveCubature()`: optionally, the `IntegrationRule` used for integration

# Keyword Arguments
- `diff_method::DifferentiationMethod`: the method to use for
calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired
"""
function integral(
        f,
        area::Meshes.PolyArea,
        rule::I;
        kwargs...
) where {I <: IntegrationRule}
    # Partition the PolyArea, sum the integrals over each of those areas
    subintegral(area) = integral(f, area, rule; kwargs...)
    subgeometries = collect(Meshes.elements(Meshes.discretize(area)))
    return sum(subintegral, subgeometries)
end
