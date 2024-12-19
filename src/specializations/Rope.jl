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
    integral(f, rope::Rope, rule = GaussKronrod();
             diff_method=FiniteDifference(), FP=Float64)

Like [`integral`](@ref) but integrates along the domain defined by `rope`. The
specified integration `rule` is applied independently to each segment formed by
consecutive points in the Rope.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `rope`: a `Rope` that defines the integration domain
- `rule = GaussKronrod()`: optionally, the `IntegrationRule` used for integration

# Keyword Arguments
- `diff_method::DifferentiationMethod = FiniteDifference()`: the method to use for
calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired
"""
function integral(
        f,
        rope::Meshes.Rope,
        rule::I;
        kwargs...
) where {I <: IntegrationRule}
    # Convert the Rope into Segments, sum the integrals of those
    return sum(segment -> integral(f, segment, rule; kwargs...), Meshes.segments(rope))
end

# Use HCubature.hcubature_buffer to reduce allocations
function integral(
    f,
    rope::Meshes.Rope,
    rule::HAdaptiveCubature;
    FP::Type{T} = Float64,
    kwargs...
) where {T <: AbstractFloat}
    # Geometry information
    N = Meshes.paramdim(rope)
    segments = Meshes.segments(rope)

    # Use a sample integrand to develop and append a buffer to the given rule
    sample = first(segments)
    integrand(ts) = f(sample(ts...)) * differential(sample, ts)
    uintegrand(ts) = Unitful.ustrip.(integrand(ts))
    buffer = HCubature.hcubature_buffer(uintegrand, _zeros(FP, N), _ones(FP, N))
    rule = HAdaptiveCubature(rule.kwargs..., buffer = buffer)

    # Convert the Rope into Segments, sum the integrals of those
    _subintegral(seg) = _integral(f, seg, rule; FP = FP, kwargs...)
    return sum(_subintegral, segments)
end
