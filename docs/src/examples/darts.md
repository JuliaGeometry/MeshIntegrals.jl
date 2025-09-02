# Darts (Draft)

Steps
- Comstruct a set of geometries representimg a dartboard with individual sector scores
- Develop a model of the dart trajectory with probability density distribution
- Use integration over each geometry to determine the probabilities of particular outcomes
- Calculate expected value for the throw, repeat for other distributions to compare strategies

```@example darts
using Meshes
using MeshIntegrals
using Unitful
```

## Modeling the Dartboard

Model the geometries
```@example darts
center = Point(0u"m", 0u"m", 1.5u"m")
point(r, ϕ) = center + Meshes.Vec(0u"m", r*sin(ϕ)*u"m", r*cos(ϕ)*u"m")

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

Point system
```@example darts
sector_width = 2pi/20
ring1 = [20, 1, 18, 4, 13, 6, 10, 15, 2, 17, 3, 19, 7, 16, 8, 11, 14, 9, 12, 5]
ring2 = 3 .* ring1
ring3 = ring1
ring4 = 2 .* ring1

bullseye_inner = (points=50,)
bullseye_outer = (points=25,)
```

## Modeling the Dart Trajectory
