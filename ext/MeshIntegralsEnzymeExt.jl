module MeshIntegralsEnzymeExt

using MeshIntegrals: MeshIntegrals, AutoEnzyme, GeometryOrDomain
using Meshes: Meshes
using Enzyme: Enzyme

function MeshIntegrals.jacobian(
        geometry::Meshes.Geometry,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        ::AutoEnzyme
) where {T <: AbstractFloat}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end
    return Meshes.to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

# Supports all geometries except for those that throw errors
# See GitHub Issue #154 for more information
MeshIntegrals.supports_autoenzyme(::Type{<:GeometryOrDomain}) = true
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.BezierCurve}) = false
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.CylinderSurface}) = false
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.Cylinder}) = false
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.ParametrizedCurve}) = false

end
