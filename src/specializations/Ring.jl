################################################################################
#                      Specialized Methods for Ring
#
# Why Specialized?
#   The Ring geometry defines a polytope whose length spans segments between
#   consecutive points that form a closed path. Meshes.jl does not define a
#   parametric function for Ring's, but they can be decomposed into their
#   constituent Segment's, integrated separately, and then summed.
################################################################################

"""
    integral(f, ring::Ring, rule = GaussKronrod();
             diff_method=FiniteDifference(), FP=Float64)

Like [`integral`](@ref) but integrates along the domain defined by `ring`. The
specified integration `rule` is applied independently to each segment formed by
consecutive points in the Ring.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `ring`: a `Ring` that defines the integration domain
- `rule = GaussKronrod()`: optionally, the `IntegrationRule` used for integration

# Keyword Arguments
- `diff_method::DifferentiationMethod = FiniteDifference()`: the method to use for
calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired
"""
function integral(
        f,
        ring::Meshes.Ring,
        rule::I;
        kwargs...
) where {I <: IntegrationRule}
    # Convert the Ring into Segments, sum the integrals of those 
    subintegral(segment) = _integral(f, segment, rule; kwargs...)
    return sum(subintegral, Meshes.segments(ring))
end
