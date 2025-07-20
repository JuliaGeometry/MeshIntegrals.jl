################################################################################
#                      Specialized Methods for Rope
#
# Why Specialized?
#   The Rope geometry defines a polytope whose length spans segments between
#   consecutive points. Meshes.jl does not define a parametric function for
#   Rope's, but they can be decomposed into their constituent Segment's,
#   integrated separately, and then summed.
################################################################################

"""
    integral(f, rope::Rope[, rule = GaussKronrod()]; kwargs...)

Like [`integral`](@ref) but integrates along the domain defined by `rope`. The
specified integration `rule` is applied independently to each segment formed by
consecutive points in the Rope.
"""
function integral(
        f,
        rope::Meshes.Rope,
        rule::I;
        kwargs...
) where {I <: IntegrationRule}
    # Convert the Rope into Segments, sum the integrals of those
    subintegral(segment) = integral(f, segment, rule; kwargs...)
    return sum(subintegral, Meshes.segments(rope))
end
