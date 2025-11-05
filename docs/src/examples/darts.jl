using CairoMakie
using Colors
using Distributions
using Meshes
using MeshIntegrals
using Unitful
using Unitful.DefaultSymbols: mm, m


################################################################################
#                             Geometry Utils
################################################################################

# Define a plane (wall) on which the dartboard rests
dartboard_center = Meshes.Point(0m, 0m, 1.5m)
dartboard_plane = Plane(dartboard_center, Meshes.Vec(1, 0, 0))
point(t, ϕ) = dartboard_plane(t * sin(ϕ), t * cos(ϕ))
point(r::Unitful.Length, ϕ) = point(ustrip(u"m", r), ϕ)

# Score-earning areas where darts can land
struct ScoredRegion{G, C}
    geometry::G
    points::Int64
    color::C
end

# A sectoral geometry bounded by constant r and ϕ limits
struct Sector{L <: Unitful.Length, A}
    r_inner::L
    r_outer::L
    ϕ_a::A 
    ϕ_b::A
end
#Sector(rs, ϕs) = Sector(rs..., ϕs...)

# Sector -> Ngon
function _Ngon(sector::Sector; N=32)
	ϕs = range(sector.ϕ_a, sector.ϕ_b, length=N)
    arc_o = [point(sector.r_outer, ϕ) for ϕ in ϕs]
    arc_i = [point(sector.r_inner, ϕ) for ϕ in reverse(ϕs)]
    return Meshes.Ngon(arc_o..., arc_i...)
end


################################################################################
#                      Construct Dartboard ScoredRegions
################################################################################

# Sectorize the board
#   scores
ring_pts = [20, 1, 18, 4, 13, 6, 10, 15, 2, 17, 3, 19, 7, 16, 8, 11, 14, 9, 12, 5]
board_points = hcat(ring_pts, (3 .* ring_pts), ring_pts, (2 .* ring_pts))
#   colors
ring_c1 = repeat([colorant"black", colorant"white"], 10)
ring_c2 = repeat([colorant"red", colorant"green"], 10)
board_colors = hcat(ring_c1, ring_c2, ring_c1, ring_c2)
#   geometries
Δϕ = 2π/20
ϕas = range(0, 2π - Δϕ, length=20) .- (Δϕ / 2)
ϕbs = range(0, 2π - Δϕ, length=20) .+ (Δϕ / 2)
ϕs = Iterators.zip(ϕas, ϕbs)
rs = [ (16mm, 99mm), (99mm, 107mm), (107mm, 162mm), (162mm, 170mm) ]
board_ngons = map(((ϕs, rs),) -> _Ngon(Sector(rs..., ϕs...)), Iterators.product(ϕs, rs))

# Consolidate the Sectors
sector_data = Iterators.zip(board_ngons, board_points, board_colors)
board_regions = map(args -> ScoredRegion(args...), sector_data)

# Center region
bullseye_inner = ScoredRegion(Meshes.Circle(dartboard_plane, (6.35e-3)m), 50, colorant"red")
bullseye_outer = ScoredRegion(_Ngon(Sector(6.35mm, 16.0mm, 0.0, 2π)), 25, colorant"green")

# Get set of all regions
all_regions = vcat(vec(board_regions), bullseye_inner, bullseye_outer)


################################################################################
#                               Makie Utils
################################################################################

# To Makie-compatible GeometryBasics types
_Point2f(p::Meshes.Point) = Point2f(ustrip.(u"m", (p.coords.y, p.coords.z))...)
_Point3f(p::Meshes.Point) = Point3f(ustrip.(u"m", (p.coords.x, p.coords.y, p.coords.z))...)

# To Makie-compatible polygons
_poly(circle::Meshes.Circle; N=32) = [(_Point3f(circle(t)) for t in range(0, 1, length=N))...]
_poly(ngon::Meshes.Ngon) = [(_Point3f(pt) for pt in ngon.vertices)...]
_poly2d(circle::Meshes.Circle; N=32) = [(_Point2f(circle(t)) for t in range(0, 1, length=N))...]
_poly2d(ngon::Meshes.Ngon) = [(_Point2f(pt) for pt in ngon.vertices)...]


################################################################################
#                            Figure - Dartboard
################################################################################

# Illustrate the dartboard
fig = Figure()
ax = Axis(fig[1, 1], xlabel="y [m]", ylabel="z [m]")
ax.aspect = DataAspect()
for region in all_regions
    poly!(ax, _poly2d(region.geometry), color=region.color)
	
    # Write score label on geometry
    centerPt = centroid(region.geometry)
    center = ustrip.(u"m", [centerPt.coords.y, centerPt.coords.z])
    text!(ax, string(region.points), position=Point2f(center...), align=(:center,:center), color=:blue, fontsize=10)
end

fig
save("dartboard.png", fig)
