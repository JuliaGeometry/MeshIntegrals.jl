module MeshIntegralsEnzymeExt

using MeshIntegrals: MeshIntegrals, AutoEnzyme
using Meshes: Geometry, to
using Enzyme: Enzyme

function MeshIntegrals.jacobian(
        geometry::Geometry, ts::V, ::AutoEnzyme) where {V <: Union{AbstractVector, Tuple}}
    to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

end
