# Support Status

This library aims to enable users to calculate the value of integrals over all
[**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl) geometry types using
a number of numerical integration rules and techniques. However, some combinations
of geometry types and integration rules are ill-suited (and a few are simply not
yet implemented).

## General Recommendations

In general, `GaussKronrod` integration rules are recommended (and the default) for
geometries with one parametric dimension. For geometries with more than one
parametric dimension, e.g. surfaces and volumes, `HAdaptiveCubature` rules are
recommended (and the default).

While it is currently possible to apply nested `GaussKronrod` rules to numerically
integrate surfaces, this produces results that are strictly inferior to using an
equivalent `HAdaptiveCubature` rule, so support for this usage has been deprecated.
In version 16.x of MeshIntegrals.jl, using a `GaussKronrod` rule for a surface
will work but will yield a deprecation warning. Beginning with a future version
17.0, this combination will simply be unsupported and throw an error.

## The Support Matrix

The following Support Matrix captures the current state of support for all geometry/rule
combinations. Entries with a green check mark are fully supported and pass unit tests
designed to check for accuracy.

| `Meshes.Geometry` | `GaussKronrod` | `GaussLegendre` | `HAdaptiveCubature` |
|----------|----------------|---------------|---------------------|
| `Ball` in `𝔼{2}` | ⚠️ | ✅ | ✅ |
| `Ball` in `𝔼{3}` | 🛑 | ✅ | ✅ |
| `BezierCurve` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{2}` | ⚠️ | ✅ | ✅ |
| `Box` in `𝔼{≥3}` | 🛑 | ✅ | ✅ |
| `Circle` | ✅ | ✅ | ✅ |
| `Cone` | 🛑 | ✅ | ✅ |
| `ConeSurface` | ⚠️ | ✅ | ✅ |
| `Cylinder` | 🛑 | ✅ | ✅ |
| `CylinderSurface` | ⚠️ | ✅ | ✅ |
| `Disk` | ⚠️ | ✅ | ✅ |
| `Ellipsoid` | ✅ | ✅ | ✅ |
| `Frustum` | 🛑 | 🎗️ | 🎗️ |
| `FrustumSurface` | ⚠️ | ✅ | ✅ |
| `Hexahedron` | ✅ | ✅ | ✅ |
| `Line` | ✅ | ✅ | ✅ |
| `ParaboloidSurface` | ⚠️ | ✅ | ✅ |
| `ParametrizedCurve` | ✅ | ✅ | ✅ |
| `Plane` | ✅ | ✅ | ✅ |
| `Polyarea` | 🛑 | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) |
| `Pyramid` | 🛑 | 🎗️ | 🎗️ |
| `Quadrangle` | ⚠️ | ✅ | ✅ |
| `Ray` | ✅ | ✅ | ✅ |
| `Ring` | ✅ | ✅ | ✅ |
| `Rope` | ✅ | ✅ | ✅ |
| `Segment` | ✅ | ✅ | ✅ |
| `SimpleMesh` | 🛑 | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/27) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/27) |
| `Sphere` in `𝔼{2}` | ✅ | ✅ | ✅ |
| `Sphere` in `𝔼{3}` | ⚠️ | ✅ | ✅ |
| `Tetrahedron` | ⚠️ | ✅ | ✅ |
| `Triangle` | ✅ | ✅ | ✅ |
| `Torus` | ⚠️ | ✅ | ✅ |
| `Wedge` | 🛑 | 🎗️ | 🎗️ |

| Symbol | Support Level |
|--------|---------|
| ✅ | Supported |
| 🎗️ | Planned to support in the future |
| ⚠️ | Deprecated |
| 🛑 | Not supported |