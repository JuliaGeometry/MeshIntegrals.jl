################################################################################
#                           Misc. Internal Tools
################################################################################

# Common error message structure
function _error_unsupported_combination(geometry, rule)
    msg = "Integrating a $geometry using a $rule rule not supported."
    throw(ArgumentError(msg))
end

################################################################################
#                           DifferentiationMethod
################################################################################

"""
    supports_autoenzyme(geometry::Geometry)
    supports_autoenzyme(type::Type{<:GeometryOrDomain})

Return whether a geometry or geometry type has a parametric function that can be
differentiated with Enzyme. See GitHub Issue #154 for more information.
"""
function supports_autoenzyme end

# Returns false for all geometries when Enzyme extension is not loaded
supports_autoenzyme(::Type{<:Any}) = false

# If provided a geometry instance, re-run with the type as argument
supports_autoenzyme(::G) where {G <: GeometryOrDomain} = supports_autoenzyme(G)

"""
    _check_diff_method_support(::Geometry, ::DifferentiationMethod) -> nothing

Throw an error if incompatible combination {geometry, diff_method} detected.
"""
function _check_diff_method_support end

# If diff_method == Enzyme, then perform check
function _check_diff_method_support(geometry::Geometry, ::AutoEnzyme)
    if !supports_autoenzyme(geometry)
        throw(ArgumentError("AutoEnzyme not supported for this geometry."))
    end
end

# If diff_method != AutoEnzyme, then do nothing
_check_diff_method_support(::Geometry, ::DifferentiationMethod) = nothing

"""
    _default_diff_method(geometry, FP)

Return an instance of the default DifferentiationMethod for a particular geometry
(or geometry type) and floating point type.
"""
function _default_diff_method(
        ::Type{G},
        ::Type{FP}
) where {G <: Geometry, FP <: AbstractFloat}
    # Enzyme only works with these FP types
    uses_Enzyme_supported_FP_type = (FP <: Union{Float32, Float64})

    if supports_autoenzyme(G) && uses_Enzyme_supported_FP_type
        return AutoEnzyme()
    else
        return FiniteDifference()
    end
end

# If provided a geometry instance, re-run with the type as argument
function _default_diff_method(
        ::G,
        ::Type{FP}
) where {G <: Geometry, FP <: AbstractFloat}
    return _default_diff_method(G, FP)
end

################################################################################
#                           Numerical Tools
################################################################################

"""
    _KVector(v::Meshes.Vec) -> Unitful.Quantity{CliffordNumbers.KVector}

Convert a `Vec` to a Unitful KVector.
"""
function _KVector(v::Meshes.Vec{Dim, T}) where {Dim, T}
    ucoords = Iterators.map(Unitful.ustrip, v.coords)
    return CliffordNumbers.KVector{1, VGA(Dim)}(ucoords...) * _units(v)
end

"""
    _units(geometry)

Return the Unitful.jl units associated with a particular `geometry`.
"""
_units(::Geometry{M, CRS}) where {M, CRS} = Unitful.unit(CoordRefSystems.lentype(CRS))
_units(::Meshes.Vec{Dim, T}) where {Dim, T} = Unitful.unit(T)

"""
    _zeros(T = Float64, N)

Return an `NTuple{N, T}` filled with zeros. This method avoids allocating an array, which
happens when using `Base.zeros`.
"""
_zeros(T::DataType, N::Int64) = ntuple(_ -> zero(T), N)
_zeros(N::Int) = _zeros(Float64, N)

"""
    _ones(T = Float64, N)

Return an `NTuple{N, T}` filled with ones. This method avoids allocating an array, which
happens when using `Base.ones`.
"""
_ones(T::DataType, N::Int64) = ntuple(_ -> one(T), N)
_ones(N::Int) = _ones(Float64, N)
