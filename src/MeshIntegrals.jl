module MeshIntegrals
using CliffordNumbers: CliffordNumbers, VGA, ∧
using CoordRefSystems: CoordRefSystems, CRS
using Meshes: Meshes, Geometry

import FastGaussQuadrature
import HCubature
import LinearAlgebra
import QuadGK
import Unitful

include("differentiation.jl")
export DifferentiationMethod, Analytical, FiniteDifference, AutoEnzyme, jacobian

include("utils.jl")

include("integration_rules.jl")
export IntegrationRule, GaussKronrod, GaussLegendre, HAdaptiveCubature

include("integral.jl")
export integral

include("integral_aliases.jl")
export lineintegral, surfaceintegral, volumeintegral

# Integration methods specialized for particular geometries
include("specializations/BezierCurve.jl")
include("specializations/ConeSurface.jl")
include("specializations/CylinderSurface.jl")
include("specializations/FrustumSurface.jl")
include("specializations/Line.jl")
include("specializations/Plane.jl")
include("specializations/Ray.jl")
include("specializations/Ring.jl")
include("specializations/Rope.jl")
include("specializations/Tetrahedron.jl")
include("specializations/Triangle.jl")
end
