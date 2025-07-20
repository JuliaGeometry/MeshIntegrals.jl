# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Adds integral methods for `PolyArea` and subtypes of `Domain` with explicit tests for `CartesianGrid`, `PolyArea`, `RegularGrid`, `SimpleMesh`, and `StructuredGrid`.

### Changed

- Generalizes the alias functions (e.g. `lineintegral`) to also accept `Domain`s.
- Updated docstrings to improve clarity and consistency.


## [0.16.3] - 2025-06-26

### Changed

- Updated compatibility entry for [CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl) to support v0.17 and v0.18 ([#175](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/175), [#180](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/180)).


## [0.16.2] - 2025-03-16

### Added

- Add official support for integrating over `Frustum`, `Pyramid`, and `Wedge` geometries ([#168](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/168)).

### Changed

- Improved unit testing with analytical solutions ([#162](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/162), [#166](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/166)).
- Made significant improvements to documentation ([#165](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/165)).
- Updated compatibility entry for [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) to support v0.53 ([#170](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/170), [#171](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/171)).
- Updates several CI dependencies ([#163](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/163), [#167](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/167), [#172](https://github.com/JuliaGeometry/MeshIntegrals.jl/pull/172)).


## [0.16.1] - 2024-12-29

### Changed

- Implemented a more efficient internal parametric transformation for `Meshes.Tetrahedron`, resulting in about an 80% integral performance improvement.

### Fixed

- Fixed a bug where `integral` would default to `diff_method=AutoEnzyme()` even when the Enzyme extension isn't loaded.


## [0.16.0] - 2024-12-14

### Added

- Added a `diff_method` keyword argument to the `integral` API, allowing the user to specify which differentiation method should be used when calculating differential element magnitudes throughout the integration domain.
- Implemented `DifferentiationMethod` types:
    - `FiniteDifference` for finite-difference approximation.
    - `AutoEnzyme` for using [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) automatic differentiation (AD) via a package extension.
- Added `diff_method` as an optional third argument to the `jacobian` and `differential` API.
- Adds standardized support for integrating over `Tetrahedron` volumes.
- Generalizes integrand functions to support any `f::Any` with a method defined for `f(::Point)`.
- Refactored specialization methods by implementing an internal `_ParametricGeometry <: Meshes.Geometry` to define geometries with custom parametric functions, standardizing support for `BezierCurve`, `Line`, `Plane`, `Ray`, `Tetrahedron`, and `Triangle`.
- Significant performance improvements:
  - Achieved an 80x improvement when integrating over `BezierCurve`.
  - Achieved an up-to-4x improvement when integrating using `HAdaptiveCubature`.

### Deprecated

- Deprecated manual specification of `GaussKronrod` rules for surfaces, i.e. geometries where `Meshes.paramdim(geometry) == 2`. A warning is now generated recommending users switch to `HAdaptiveCubature`.

### Fixed

- Refactored the unit test system.
  - Standardized `combinations.jl` tests by constructing a `TestableGeometry` package and passing it to a `@test` generation function to provide more thorough and standardized test coverage.
  - Reorganized `@testsnippet`s to exist in same source file as relevant tests.
  - Removed `:extended` tag from `Tetrahedron` now that performance is significantly improved.


## [0.15.2] - 2024-10-25

MeshIntegrals.jl is now owned by the JuliaGeometry organization!

### Added

- Added a benchmarking suite using [AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl).
- Implemented more unit tests with analytical solutions.

### Changed

- Tagged unit tests for `Meshes.Box` (4D) and `Tetrahedron` as `:extended`, removing them from automatic CI testing due to lengthy compute times.


## [0.15.1] - 2024-10-11

### Added

- Adds official support and unit testing for integrating new `Meshes.ParametrizedCurve` geometries.


## [0.15.0] - 2024-10-10

### Added

- Adds support for integrating geometries with any number of parametric dimensions by generalizing `differential` to n-dimensions via geometric algebra.
- Adds official support and unit testing for integrating `Meshes.Ellipsoid` and `Meshes.Hexahedron` geometries.

### Changed

- Refactored the unit test system.
    - Completed transition away from previous `@test` generation system which only tested unit integrands (i.e. `f(point) = 1.0`) using `Meshes.measure` as a benchmark.
    - Used [TestItems.jl](https://github.com/julia-vscode/TestItems.jl) to define independent `@testitem` packages, many with analytically-derived solutions.


## [0.14.1] - 2024-10-04

### Added

- Adds official support and unit testing for integrating `Meshes.Quadrangle` geometries.

### Changed

- Continued work transitioning away from previous `@test` generation system.

### Fixed

- Fixed a bug in `differential` finite-difference approximation that may have reduced accuracy of some integral calculations.
- Made improvements to type stability with floating point type keyword argument `FP`.


## [0.14.0] - 2024-09-28

### Changed

- Rename `IntegrationAlgorithm` to `IntegrationRule` and consolidate terminology around integration rules (versus "algorithms", "settings", etc).
- Convert floating point precision option `FP` from an optional argument to a keyword argument.
- Continued work transitioning away from previous `@test` generation system.

### Removed

- Remove `derivative` and `unitdirection` functions. All `derivative` functionality was merged into `jacobian`.
