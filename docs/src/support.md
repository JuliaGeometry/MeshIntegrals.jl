# Support Status

This library aims to enable users to calculate the value of integrals over all
[**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl) geometry types using
a number of numerical integration rules and techniques.

In general, `GaussKronrod` integration rules are recommended (and the default) for
geometries with one parametric dimension. For geometries with more than one
parametric dimension, e.g. surfaces and volumes, `HAdaptiveCubature` rules are
recommended (and the default).

## The Support Matrix

The following Support Matrix captures the current state of support for all geometry/rule
combinations. Entries with a green check mark are fully supported and pass unit tests
designed to check for accuracy.

| `Meshes.Geometry/Domain` | `GaussKronrod` | `GaussLegendre` | `HAdaptiveCubature` |
|----------|----------------|---------------|---------------------|
| `Ball` in `𝔼{2}` | 🛑 | ✅ | ✅ |
| `Ball` in `𝔼{3}` | 🛑 | ✅ | ✅ |
| `BezierCurve` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{2}` | 🛑 | ✅ | ✅ |
| `Box` in `𝔼{≥3}` | 🛑 | ✅ | ✅ |
| `CartesianGrid` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `CartesianGrid` in `𝔼{≥2}` | 🛑 | ✅ | ✅ |
| `Circle` | ✅ | ✅ | ✅ |
| `Cone` | 🛑 | ✅ | ✅ |
| `ConeSurface` | 🛑 | ✅ | ✅ |
| `Cylinder` | 🛑 | ✅ | ✅ |
| `CylinderSurface` | 🛑 | ✅ | ✅ |
| `Disk` | 🛑 | ✅ | ✅ |
| `Ellipsoid` | 🛑 | ✅ | ✅ |
| `Frustum` | 🛑 | ✅ | ✅ |
| `FrustumSurface` | 🛑 | ✅ | ✅ |
| `Hexahedron` | 🛑 | ✅ | ✅ |
| `Line` | ✅ | ✅ | ✅ |
| `ParaboloidSurface` | 🛑 | ✅ | ✅ |
| `ParametrizedCurve` | ✅ | ✅ | ✅ |
| `Plane` | 🛑 | ✅ | ✅ |
| `PolyArea` | 🛑 | ✅ | ✅ |
| `Pyramid` | 🛑 | ✅ | ✅ |
| `Quadrangle` | 🛑 | ✅ | ✅ |
| `Ray` | ✅ | ✅ | ✅ |
| `RegularGrid` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `RegularGrid` in `𝔼{≥2}` | 🛑 | ✅ | ✅ |
| `Ring` | ✅ | ✅ | ✅ |
| `Rope` | ✅ | ✅ | ✅ |
| `Segment` | ✅ | ✅ | ✅ |
| `SimpleMesh` | 🛑 | ✅ | ✅ |
| `Sphere` | 🛑 | ✅ | ✅ |
| `StructuredGrid` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `StructuredGrid` in `𝔼{≥2}` | 🛑 | ✅ | ✅ |
| `Tetrahedron` | 🛑 | ✅ | ✅ |
| `Triangle` | 🛑 | ✅ | ✅ |
| `Torus` | 🛑 | ✅ | ✅ |
| `Wedge` | 🛑 | ✅ | ✅ |

| Symbol | Support Level |
|--------|---------|
| ✅ | Supported |
| 🛑 | Not supported |
