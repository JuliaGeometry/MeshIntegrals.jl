################################################################################
#                           Misc. Internal Tools
################################################################################

# Calculate Gauss-Legendre nodes/weights and convert to type T
function _gausslegendre(T, n)
    xs, ws = FastGaussQuadrature.gausslegendre(n)
    return T.(xs), T.(ws)
end

# Extract the length units used by the CRS of a Point
_units(pt::Meshes.Point{M, CRS}) where {M, CRS} = first(CoordRefSystems.units(CRS))

# Meshes.Vec -> CliffordNumber.KVector
function _clifford(v::Meshes.Vec{Dim, T}) where {Dim, T}
    units = Unitful.unit(T)
    cliffordnumber = KVector{1, VGA(Dim)}(Unitful.ustrip.(units, v.coords...))
    return (units, cliffordnumber)
end
