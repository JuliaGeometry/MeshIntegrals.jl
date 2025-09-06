# Darts (Draft)

Steps
- Construct a set of geometries representing a dartboard with individual sector scores
- Develop a model of the dart trajectory with probability density distribution
- Use integration over each geometry to determine the probabilities of particular outcomes
- Calculate expected value for the throw, repeat for other distributions to compare strategies

Note: can use `@setup darts` block to hide some implementation code

```@example darts
# using Distributions
using Meshes
using MeshIntegrals
using Unitful
```

## Modeling the Dartboard

Define a dartboard coordinate system
```@example darts
dartboard_center = Point(0u"m", 0u"m", 1.5u"m")
dartboard_plane = Plane(dartboard_center, Meshes.Vec(1, 0, 0))

function point(r::Unitful.Length, ϕ)
    t = ustrip(r, u"m")
    dartboard_plane(t * sin(ϕ), t * cos(ϕ))
end
```

Model the bullseye region
```@example darts
bullseye_inner = (geometry = Circle(dartboard_plane, 6.35u"mm"), points = 50)
bullseye_outer = (geometry = Circle(dartboard_plane, 16u"mm"), points = 25)
# TODO subtract center circle region from outer circle region -- or replace with another annular geometry
```

Model the sectors
```@example darts
# Scores on the Board
ring1 = [20, 1, 18, 4, 13, 6, 10, 15, 2, 17, 3, 19, 7, 16, 8, 11, 14, 9, 12, 5]
ring2 = 3 .* ring1
ring3 = ring1
ring4 = 2 .* ring1
board_points = hcat(ring1, ring2, ring3, ring4)

# Locations
sector_width = 2π/20
phis_a = range(0, 2π, 20) .- sector_width/2
phis_b = range(0, 2π, 20) .+ sector_width/2
rs_inner = [16, 99, 107, 162]u"mm"
rs_outer = [99, 107, 162, 170]u"mm"
```

Define a struct to manage sector data
```@example darts
struct Sector{L, A}
    r_inner::L
    r_outer::L
    phi_a::A
    phi_b::A
    points::Int64
end

function to_ngon(sector::Sector; N=8)
	ϕs = range(sector.phi_a, sector.phi_b, length=N)
    arc_o = [point(sector.r_outer, ϕ) for ϕ in ϕs]
    arc_i = [point(sector.r_inner, ϕ) for ϕ in reverse(ϕs)]
    return Ngon(arc_o..., arc_i...)
end
```

## Modeling the Dart Trajectory

Define a probability distribution for where the dart will land
```
dist = MvNormal(μs, σs)
```

Integrand function is the distribution's PDF value at any particular point
```
function integrand(p::Point)
    v_error = dist_center - p
    pdf(dist, v_error)
end
```
