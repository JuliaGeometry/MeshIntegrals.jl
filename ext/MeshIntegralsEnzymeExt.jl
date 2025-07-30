module MeshIntegralsEnzymeExt

import MeshIntegrals
import Meshes
import Enzyme

function MeshIntegrals.jacobian(
        geometry::Meshes.Geometry,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        ::MeshIntegrals.AutoEnzyme
) where {T <: AbstractFloat}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end
    return Meshes.to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

# Supports all geometries/domains by default
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.Geometry}) = true
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.Domain}) = true
# except for those that throw errors (see GitHub Issue #154)
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.BezierCurve}) = false
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.Cylinder}) = false
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.CylinderSurface}) = false
MeshIntegrals.supports_autoenzyme(::Type{<:Meshes.ParametrizedCurve}) = false

end
