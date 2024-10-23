################################################################################
#                     Jacobian and Differential Elements
################################################################################

"""
    jacobian(geometry, ts; ε=1e-6)

Calculate the Jacobian of a geometry at some parametric point `ts` using a simple
finite-difference approximation with step size `ε`.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `ε`: step size to use for the finite-difference approximation
"""
function jacobian(
        geometry::Meshes.Geometry,
        ts::V;
        ε = 1e-6
) where {V <: Union{AbstractVector, Tuple}}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end

    T = eltype(ts)
    ε = T(ε)

    # Get the partial derivative along the n'th axis via finite difference
    #   approximation, where ts is the current parametric position
    function ∂ₙr(ts, n)
        # Build left/right parametric coordinates with non-allocating iterators 
        left = Iterators.map(((i, t),) -> i == n ? t - ε : t, enumerate(ts))
        right = Iterators.map(((i, t),) -> i == n ? t + ε : t, enumerate(ts))
        # Select orientation of finite-diff
        if ts[n] < T(0.01)
            # Right
            return (geometry(right...) - geometry(ts...)) / ε
        elseif T(0.99) < ts[n]
            # Left
            return (geometry(ts...) - geometry(left...)) / ε
        else
            # Central
            return (geometry(right...) - geometry(left...)) / 2ε
        end
    end

    return ntuple(n -> ∂ₙr(ts, n), Dim)
end

function jacobian(
        bz::Meshes.BezierCurve,
        ts::V
) where {V <: Union{AbstractVector, Tuple}}
    t = only(ts)
    # Parameter t restricted to domain [0,1] by definition
    if t < 0 || t > 1
        throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
    end

    # Aliases
    P = bz.controls
    N = Meshes.degree(bz)

    # Ensure that this implementation is tractible: limited by ability to calculate
    #   binomial(N, N/2) without overflow. It's possible to extend this range by
    #   converting N to a BigInt, but this results in always returning BigFloat's.
    N <= 1028 || error("This algorithm overflows for curves with ⪆1000 control points.")

    # Generator for Bernstein polynomial functions
    B(i, n) = t -> binomial(Int128(n), i) * t^i * (1 - t)^(n - i)

    # Derivative = N Σ_{i=0}^{N-1} sigma(i)
    #   P indices adjusted for Julia 1-based array indexing
    sigma(i) = B(i, N - 1)(t) .* (P[(i + 1) + 1] - P[(i) + 1])
    derivative = N .* sum(sigma, 0:(N - 1))

    return (derivative,)
end

################################################################################
#                          Differential Elements
################################################################################

"""
    differential(geometry, ts)

Return the magnitude of the differential element (length, area, volume, etc) of
the parametric function for `geometry` at arguments `ts`.
"""
function differential(
        geometry::Meshes.Geometry,
        ts::V
) where {V <: Union{AbstractVector, Tuple}}
    J = Iterators.map(_KVector, jacobian(geometry, ts))
    return LinearAlgebra.norm(foldl(∧, J))
end
