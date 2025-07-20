################################################################################
#                   Specialized Methods for BezierCurve
#
# Why Specialized?
#   The parametric function in Meshes.jl for BezierCurve accepts an argument
#   of type Meshes.BezierEvalMethod, in which the method for generating
#   parametric points along the curve is specified. These specialized methods
#   are essentially identical to the generic ones, except that they provide a
#   keyword argument and pass through the specified algorithm choice.
################################################################################

################################################################################
#                              integral
################################################################################
"""
    integral(f, curve::BezierCurve[, rule = GaussKronrod()]; kwargs...)

Like [`integral`](@ref) but integrates along the domain defined by `curve`.

# Special Keyword Arguments
- `alg = Meshes.Horner()`:  the method to use for parametrizing `curve`. Alternatively,
`alg=Meshes.DeCasteljau()` can be specified for increased accuracy, but comes with a
steep performance cost, especially for curves with a large number of control points.
"""
function integral(
        f,
        curve::Meshes.BezierCurve,
        rule::IntegrationRule;
        alg::Meshes.BezierEvalMethod = Meshes.Horner(),
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(curve, FP),
        kwargs...
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _check_diff_method_support(curve, diff_method)

    # Generate a _ParametricGeometry whose parametric function auto-applies the alg kwarg
    param_curve = _ParametricGeometry(_parametric(curve, alg), Meshes.paramdim(curve))

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_curve, rule; diff_method = diff_method, FP = FP, kwargs...)
end

################################################################################
#                              Parametric
################################################################################

# Wrap (::BezierCurve)(t, ::BezierEvalMethod) into f(t) by embedding second argument
function _parametric(curve::Meshes.BezierCurve, alg::Meshes.BezierEvalMethod)
    return t -> curve(t, alg)
end
