################################################################################
#                          DifferentiationMethods
################################################################################

"""
    DifferentiationMethod

A category of types used to specify the desired method for calculating derivatives.
Derivatives are used to form Jacobian matrices when calculating the differential
element size throughout the integration region.

See also [`FiniteDifference`](@ref), [`Analytical`](@ref).
"""
abstract type DifferentiationMethod end

"""
    FiniteDifference(ε=1e-6)

Use to specify use of a finite-difference approximation method with a step size
of `ε` for calculating derivatives.
"""
struct FiniteDifference{T <: AbstractFloat} <: DifferentiationMethod
    ε::T
end

# Default constructors
FiniteDifference{T}() where {T <: AbstractFloat} = FiniteDifference{T}(T(1e-6))
FiniteDifference() = FiniteDifference{Float64}()

"""
    Analytical()

Use to specify use of analytically-derived solutions for calculating derivatives.
These solutions are currently defined only for a subset of geometry types.

# Supported Geometries:
- `BezierCurve`
- `Line`
- `Plane`
- `Ray`
- `Tetrahedron`
- `Triangle`
"""
struct Analytical <: DifferentiationMethod end

"""
    AutoEnzyme()

Use to specify use of the Enzyme.jl for calculating derivatives.
"""
struct AutoEnzyme <: DifferentiationMethod
    function AutoEnzyme()
        if :Enzyme ∉ names(Main, imported = true)
            error("Please load Enzyme.jl to use AutoEnzyme().")
        end
        return new()
    end
end

# Future Support:
#   struct AutoZygote <: DifferentiationMethod end

################################################################################
#                                  Jacobian
################################################################################

"""
    jacobian(geometry, ts[, diff_method])

Calculate the Jacobian of a `geometry`'s parametric function at some point `ts`.
Optionally, direct the use of a particular differentiation method `diff_method`; by default
use analytic solutions where possible and finite difference approximations
otherwise.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `diff_method`: the desired [`DifferentiationMethod`](@ref) to use
"""
function jacobian(
        geometry::G,
        ts::Union{AbstractVector{T}, NTuple{N, T}}
) where {G <: Geometry, N, T <: AbstractFloat}
    return jacobian(geometry, ts, _default_method(G))
end

function jacobian(
        geometry::Geometry,
        ts::Union{AbstractVector{T}, NTuple{N, T}},
        diff_method::FiniteDifference
) where {N, T <: AbstractFloat}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end

    ε = T(diff_method.ε)

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

################################################################################
#                          Differential Elements
################################################################################

"""
    differential(geometry, ts[, diff_method])

Calculate the differential element (length, area, volume, etc) of the parametric
function for `geometry` at arguments `ts`. Optionally, direct the use of a
particular differentiation method `diff_method`; by default use analytic solutions where
possible and finite difference approximations otherwise.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `diff_method`: the desired [`DifferentiationMethod`](@ref) to use
"""
function differential(
        geometry::G,
        ts::V,
        diff_method::DifferentiationMethod = _default_method(G)
) where {G <: Geometry, V <: Union{AbstractVector, Tuple}}
    # Calculate the Jacobian, convert Vec -> KVector
    J = jacobian(geometry, ts, diff_method)
    J_kvecs = Iterators.map(_kvector, J)

    # Extract units from Geometry type
    Dim = Meshes.paramdim(geometry)
    units = _units(geometry)^Dim

    # Return norm of the exterior products
    element = foldl(∧, J_kvecs)
    return LinearAlgebra.norm(element) * units
end
