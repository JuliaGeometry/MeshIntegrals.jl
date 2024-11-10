################################################################################
#                  Specialized Methods for Tetrahedron
#
# Why Specialized?
#   The Tetrahedron geometry is a volumetric simplex whose parametric function
#   in Meshes.jl uses barycentric coordinates on a domain {u,v,w} with coordinates
#   that are non-negative and bound by the surface $u + v + w ≤ 1$. This requires
#   a multi-step domain transformation whose derivation is detailed in the package
#   documentation.
################################################################################

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::GaussLegendre;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _error_unsupported_combination("Tetrahedron", "GaussLegendre")
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::GaussKronrod;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Tetrahedron, diff_method)

    o = zero(FP)
    ∫uvw(u, v, w) = f(tetrahedron(u, v, w))
    ∫vw(v, w) = QuadGK.quadgk(u -> ∫uvw(u, v, w), o, FP(1 - v - w); rule.kwargs...)[1]
    ∫w(w) = QuadGK.quadgk(v -> ∫vw(v, w), o, FP(1 - w); rule.kwargs...)[1]
    ∫ = QuadGK.quadgk(∫w, o, one(FP); rule.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * Meshes.volume(tetrahedron) * ∫
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::HAdaptiveCubature;
        diff_method::DM = FiniteDifference(),  # TODO _default_method(tetrahedron),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    function parametric(t1, t2, t3)
        _constrain(t) = (t > 1e-6) ? prevfloat(t, 100) : t

        t1, t2 = _constrain.((t1, t2))

        # Take a triangular cross-section at height t3, find point in that triangle
        rem = _constrain(1 - t3)
        a = tetrahedron(0, 0, t3)
        b = tetrahedron(0, rem, t3)
        c = tetrahedron(rem, 0, t3)
        Meshes.Triangle(a, b, c)(t1, t2)
    end

    tetra = _ParametricGeometry(parametric, 3)
    return integral(f, tetra, rule; diff_method=diff_method, FP=FP)
end

################################################################################
#                               jacobian
################################################################################

_has_analytical(::Type{T}) where {T <: Meshes.Tetrahedron} = true
