# MeshIntegrals.jl

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeometry.github.io/MeshIntegrals.jl/stable/)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeometry.github.io/MeshIntegrals.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![ColPrac](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?style=flat-square)](https://github.com/SciML/ColPrac)

[![Build Status](https://github.com/JuliaGeometry/MeshIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaGeometry/MeshIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaGeometry/MeshIntegrals.jl/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/MeshIntegrals.jl)
[![Coveralls](https://coveralls.io/repos/github/JuliaGeometry/MeshIntegrals.jl/badge.svg?branch=main)](https://coveralls.io/github/JuliaGeometry/MeshIntegrals.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


**MeshIntegrals.jl** uses differential forms to enable fast and easy numerical integration of arbitrary integrand functions over domains defined via [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl) geometries. This is achieved using:
- Gauss-Legendre quadrature rules from [**FastGaussQuadrature.jl**](https://github.com/JuliaApproximation/FastGaussQuadrature.jl): `GaussLegendre(n)`
- H-adaptive Gauss-Kronrod quadrature rules from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl): `GaussKronrod(kwargs...)`
- H-adaptive cubature rules from [**HCubature.jl**](https://github.com/JuliaMath/HCubature.jl): `HAdaptiveCubature(kwargs...)`

These solvers have support for integrand functions that produce scalars, vectors, and [**Unitful.jl**](https://github.com/PainterQubits/Unitful.jl) `Quantity` types. While HCubature.jl does not natively support `Quantity` type integrands, this package provides a compatibility layer to enable this feature.

## Usage

### Basic

```julia
integral(f, geometry)
```
Performs a numerical integration of some integrand function `f` over the domain specified by `geometry`. The integrand function can be anything callable with a method `f(::Meshes.Point)`. A default integration method will be automatically selected according to the geometry: `GaussKronrod()` for 1D, and `HAdaptiveCubature()` for all others.

```julia
integral(f, geometry, rule)
```
Performs a numerical integration of some integrand function `f` over the domain specified by `geometry` using the specified integration rule, e.g. `GaussKronrod()`. The integrand function can be anything callable with a method `f(::Meshes.Point)`.

Additionally, several optional keyword arguments are defined in [the API](https://juliageometry.github.io/MeshIntegrals.jl/stable/api/) to provide additional control over the integration mechanics.

### Aliases

```julia
lineintegral(f, geometry)
surfaceintegral(f, geometry)
volumeintegral(f, geometry)
```
Alias functions are provided for convenience. These are simply wrappers for `integral` that also validate that the provided `geometry` has the expected number of parametric dimensions. Like with `integral`, a `rule` can also optionally be specified as a third argument.
- `lineintegral` is used for curve-like geometries or polytopes (e.g. `Segment`, `Ray`, `BezierCurve`, `Rope`, etc)
- `surfaceintegral` is used for surfaces (e.g. `Disk`, `Sphere`, `CylinderSurface`, etc)
- `volumeintegral` is used for (3D) volumes (e.g. `Ball`, `Cone`, `Torus`, etc)

### Example

```julia
using Meshes
using MeshIntegrals
using Unitful

# Define a Bezier curve whose path approximates a sine-wave on the xy-plane
N = 361
curve = BezierCurve([Point(t*u"m", sin(t)*u"m", 0.0u"m") for t in range(-π, π, length=N)])

# Integrand function that outputs in units of Ohms/meter
function f(p::Point)
    x, y, z = to(p)
    return (1 / sqrt(1 + cos(x / u"m")^2)) * u"Ω/m"
end

# Use recommended defaults
integral(f, curve) # -> Approximately 2π Ω

# Using aliases
lineintegral(f, curve) # -> Approximately 2π Ω
surfaceintegral(f, curve) # -> throws ArgumentError

# Specifying an integration rule and settings (loosened absolute tolerance)
integral(f, curve, GaussKronrod(atol = 1e-4u"Ω")) # -> Approximately (2π ± 1e-4) Ω
```
