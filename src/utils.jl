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

# Return the default DifferentiationMethod instance for a particular geometry type
function _default_diff_method(
        g::Type{G}
) where {G <: Geometry}
    return FiniteDifference()
end

# Return the default DifferentiationMethod instance for a particular geometry instance
_default_diff_method(g::G) where {G <: Geometry} = _default_diff_method(G)

################################################################################
#                           Numerical Tools
################################################################################

# Calculate Gauss-Legendre nodes/weights and convert to type T
"""
    _gausslegendre(T, n)

Return FastGaussQuadrature.gausslegendre(n) in type T.
"""
function _gausslegendre(T, n)
    xs, ws = FastGaussQuadrature.gausslegendre(n)
    return T.(xs), T.(ws)
end

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

Return the Unitful units associated with a particular `geometry`.
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
