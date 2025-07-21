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
    integral(f, area::PolyArea[, rule = HAdaptiveCubature()]; kwargs...)

Like [`integral`](@ref) but integrates over the surface domain defined by a `PolyArea`.
The surface is first discretized into facets that are integrated independently using
the specified integration `rule`.
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
