# Darts (Draft)

Steps
- Construct a set of geometries representing a dartboard with individual sector scores
- Develop a model of the dart trajectory with probability density distribution
- Use integration over each geometry to determine the probabilities of particular outcomes
- Calculate expected value for the throw, repeat for other distributions to compare strategies

```@example darts
using CairoMakie
using Colors
using Distributions
using Meshes
using MeshIntegrals
using Unitful
using Unitful.DefaultSymbols: mm, m

red = colorant"red"
black = colorant"black"
white = colorant"white"
green = colorant"green"

# For containing any scored landing region on the dart board
struct ScoredRegion{G, C}
    geometry::G
    points::Int64
    color::C
end

# For defining an annular region
struct Sector{L, A}
    r_inner::L
    r_outer::L
    phi_a::A
    phi_b::A
end
Sector(rs, phis) = Sector(rs..., phis...)

# Sector -> Ngon
function to_ngon(sector::Sector; N=8)
	ϕs = range(sector.phi_a, sector.phi_b, length=N)
    arc_o = [point(sector.r_outer, ϕ) for ϕ in ϕs]
    arc_i = [point(sector.r_inner, ϕ) for ϕ in reverse(ϕs)]
    return Ngon(arc_o..., arc_i...)
end

function to_makie_poly(circle::Meshes.Circle)
    return nothing # TODO
end

function to_makie_poly(ngon::Meshes.Ngon)
    return nothing # TODO
end
```

## Modeling the Dartboard

Model the dartboard
```@example darts
dartboard_center = Meshes.Point(0m, 0m, 1.5m)
dartboard_plane = Plane(dartboard_center, Meshes.Vec(1, 0, 0))

function point(r::Unitful.Length, ϕ)
    t = ustrip(m, r)
    dartboard_plane(t * sin(ϕ), t * cos(ϕ))
end

# Scores on the Board
ring1 = [20, 1, 18, 4, 13, 6, 10, 15, 2, 17, 3, 19, 7, 16, 8, 11, 14, 9, 12, 5]
ring2 = 3 .* ring1
ring3 = ring1
ring4 = 2 .* ring1
board_points = hcat(ring1, ring2, ring3, ring4)

# Colors on the board
ring1 = repeat([black, white], 10)
ring2 = repeat([red, green], 10)
ring3 = ring1
ring4 = ring2
board_colors = hcat(ring1, ring2, ring3, ring4)

# Sector geometries
sector_width = 2π/20
phis_a = range(0, 2π, 20) .- sector_width/2
phis_b = range(0, 2π, 20) .+ sector_width/2
phis = Iterators.zip(phis_a, phis_b)
rs = [ (16mm, 99mm), (99mm, 107mm), (107mm, 162mm), (162mm, 170mm) ]
board_coords = Iterators.product(phis, rs)
board_sectors = map(((phis, rs),) -> Sector(rs, phis), board_coords)
board_ngons = to_ngon.(board_sectors)

# Consolidate the Sectors
sector_data = Iterators.zip(board_ngons, board_points, board_colors)
board_regions = map(args -> ScoredRegion(args...), sector_data)

# Center region
bullseye_inner = ScoredRegion(Circle(dartboard_plane, 6.35mm), 50, red)
bullseye_outer = ScoredRegion(to_ngon(Sector((6.35mm, 16.0mm), (0.0, 2π)); N=32), 25, green)

# Get set of all regions
all_regions = vcat(vec(board_regions), bullseye_inner, bullseye_outer)

fig = Figure()
ax = LScene(fig[1, 1], scenekw=(show_axis=true,))

for region in all_regions
    poly!(ax, to_makie_poly(region.geometry), color=region.color)
end

fig
```

## Modeling the Dart Trajectory

Define a probability distribution for where the dart will land
```
# TODO
dist = MvNormal(μs, σs)
```

Integrand function is the distribution's PDF value at any particular point
```
# TODO
function integrand(p::Point)
    v_error = dist_center - p
    pdf(dist, v_error)
end
```
