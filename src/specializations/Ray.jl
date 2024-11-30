################################################################################
#                   Specialized Methods for Ray
#
# Why Specialized?
#   The Ray geometry is a special case, representing a ray, originating at a point
#   and extending an infinite length in a particular direction. This requires
#   a domain transformation mapping from the typical parametric region [0,1] to
#   an infinite one (-∞,∞).
################################################################################

function integral(
        f,
        ray::Meshes.Ray,
        rule::IntegrationRule;
        kwargs...
)
    paramfunction(t) = _parametric(ray, t)
    param_ray = _ParametricGeometry(paramfunction, 1)
    return _integral(f, param_ray, rule; kwargs...)
end

############################################################################################
#                                     Parametric
############################################################################################

# Map [0, 1] ↦ [0, ∞)
function _parametric(ray::Meshes.Ray, t)
    f(t) = t / (1 - t^2)
    return ray(f(t))
end
