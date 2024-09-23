module MeshIntegrals
    using CoordRefSystems
    using LinearAlgebra
    using Meshes
    using Unitful

    import FastGaussQuadrature
    import HCubature
    import QuadGK

    include("utils.jl")
    export jacobian, derivative, unitdirection

    include("integration_algorithms.jl")
    export GaussKronrod, GaussLegendre, HAdaptiveCubature

    include("integral.jl")
    include("integral_aliases.jl")
    include("integral_line.jl")
    include("integral_surface.jl")
    include("integral_volume.jl")
    export integral, lineintegral, surfaceintegral, volumeintegral

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
