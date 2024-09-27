################################################################################
#                      Specialized Methods for Rope
#
# Why Specialized?
#   The Rope geometry defines a polytope whose length spans segments between
#   consecutive points. Meshes.jl does not define a parametric function for
#   Rope's, but they can be decomposed into their constituent Segment's,
#   integrated separately, and then summed.
################################################################################

function integral(
        f::F,
        rope::Meshes.Rope,
        rule::I;
        kwargs...
) where {F <: Function, I <: IntegrationRule}
    # Convert the Rope into Segments, sum the integrals of those
    return sum(segment -> integral(f, segment, rule; kwargs...), Meshes.segments(rope))
end
