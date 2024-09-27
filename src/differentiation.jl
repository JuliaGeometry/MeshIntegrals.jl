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
        geometry::G,
        ts::V;
        ε = 1e-6
) where {G <: Meshes.Geometry, V <: Union{AbstractVector, Tuple}}
    T = eltype(ts)

    # Get the partial derivative along the n'th axis via finite difference approximation
    #   where ts is the current parametric position (εv is a reusable buffer)
    function ∂ₙr!(εv, ts, n)
        if ts[n] < T(0.01)
            return ∂ₙr_right!(εv, ts, n)
        elseif T(0.99) < ts[n]
            return ∂ₙr_left!(εv, ts, n)
        else
            return ∂ₙr_central!(εv, ts, n)
        end
    end

    # Central finite difference
    function ∂ₙr_central!(εv, ts, n)
        εv[n] = ε
        a = ts .- εv
        b = ts .+ εv
        return (geometry(b...) - geometry(a...)) / 2ε
    end

    # Left finite difference
    function ∂ₙr_left!(εv, ts, n)
        εv[n] = ε
        a = ts .- εv
        b = ts
        return (geometry(b...) - geometry(a...)) / ε
    end

    # Right finite difference
    function ∂ₙr_right!(εv, ts, n)
        εv[n] = ε
        a = ts
        b = ts .+ εv
        return (geometry(b...) - geometry(a...)) / ε
    end

    # Allocate a re-usable ε vector
    εv = zeros(T, length(ts))

    ∂ₙr(n) = ∂ₙr!(εv, ts, n)
    return map(∂ₙr, 1:length(ts))
end

function jacobian(
        bz::Meshes.BezierCurve,
        ts::V
) where {V <: Union{AbstractVector, Tuple}}
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

    return [derivative]
end

################################################################################
#                          Differential Elements
################################################################################

"""
    differential(geometry, ts)

Calculate the differential element (length, area, volume, etc) of the parametric
function for `geometry` at arguments `ts`.
"""
function differential(
        geometry,
        ts
)
    J = jacobian(geometry, ts)

    # TODO generalize this with geometric algebra, e.g.: norm(foldl(∧, J))
    if length(J) == 1
        return norm(J[1])
    elseif length(J) == 2
        return norm(J[1] × J[2])
    elseif length(J) == 3
        return abs((J[1] × J[2]) ⋅ J[3])
    else
        error("Not implemented yet. Please report this as an Issue on GitHub.")
    end
end
